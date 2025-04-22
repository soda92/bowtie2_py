#!/usr/bin/env python3

# Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
# Copyright 2024, Translation to Python by Gemini
#
# This file is part of Bowtie 2.
#
# Bowtie 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Bowtie 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
#

# bowtie2:
#
# A wrapper script for bowtie2.  Provides various advantages over running
# bowtie2 directly, including:
#
# 1. Handling compressed inputs
# 2. Redirecting output to various files
# 3. Output directly to bam (via samtools)

import sys
import os
import argparse
import subprocess
import platform
import socket
import tempfile
import shutil
import signal
import atexit
import re
import shlex
import glob
from urllib.parse import unquote # For decoding %xx escapes in passthrough

# --- Global Variables ---
script_path = None
prog_name = None
host = socket.gethostname()
temp_dir_base = tempfile.gettempdir()
keep_temps = False
verbose_mode = False
os_is_nix = platform.system() != "Windows"
align_bin_s_base = 'bowtie2-align-s'
build_bin_base = 'bowtie2-build'
align_bin_l_base = 'bowtie2-align-l'
idx_ext_l = 'bt2l'
idx_ext_s = 'bt2'

# Will be determined based on script location and OS
align_prog_s = None
align_prog_l = None
build_bin = None

# State variables
read_fns = {} # {type: filename} e.g. {"un": "unaligned.fq"}
read_compress = {} # {type: compression_method} e.g. {"un": "gzip"}
to_delete = [] # List of temporary files/pipes to delete on exit
to_kill = [] # List of PIDs of child processes to kill on exit
old_stderr = None
log_file_handle = None

# --- Signal Handling (Mimicking Perl's %signo/@signame) ---
signal_names = {
    sig: name
    for name, sig in signal.__dict__.items()
    if name.startswith("SIG") and not name.startswith("SIG_")
}

# --- Helper Functions ---

def _info(*args):
    """Prints info messages if verbose is enabled."""
    global verbose_mode
    if verbose_mode:
        print("(INFO):", *args, file=sys.stderr)

def _error(*args):
    """Prints error messages to stderr."""
    print("(ERR):", *args, file=sys.stderr)

def _fail(*args):
    """Prints error messages and exits."""
    _error(*args)
    sys.exit(1)

def _cleanup():
    """Cleans up temporary files and kills child processes."""
    global to_delete, to_kill, keep_temps, log_file_handle, old_stderr

    # Restore stderr if redirected
    if log_file_handle:
        log_file_handle.close()
    if old_stderr:
        os.dup2(old_stderr.fileno(), sys.stderr.fileno())
        old_stderr.close()


    # Kill child processes (e.g., from mkfifo decompression)
    # Use SIGTERM first, then SIGKILL if necessary
    for pid in to_kill:
        try:
            os.kill(pid, signal.SIGTERM)
            # Optionally add a small wait and check if process is still alive
            # os.waitpid(pid, os.WNOHANG)
        except ProcessLookupError:
            pass # Process already finished
        except OSError as e:
             _error(f"Error killing process {pid}: {e}")
             # Consider trying SIGKILL here if SIGTERM failed and it's critical
             # try:
             #     os.kill(pid, signal.SIGKILL)
             # except OSError: pass


    # Delete temporary files/pipes
    if not keep_temps:
        for item in to_delete:
            try:
                if os.path.exists(item):
                    os.remove(item)
                    _info(f"Removed temporary item: {item}")
            except OSError as e:
                _error(f"Could not remove temporary item {item}: {e}")

# Register cleanup function to run on exit
atexit.register(_cleanup)

def _resolve_script_path():
    """Finds the absolute path of the script, resolving symlinks."""
    global script_path, prog_name, align_prog_s, align_prog_l, build_bin
    try:
        prog = os.path.abspath(__file__)
        while os.path.islink(prog):
            target = os.readlink(prog)
            # If the target is relative, resolve it relative to the link's directory
            if not os.path.isabs(target):
                target = os.path.join(os.path.dirname(prog), target)
            prog = os.path.abspath(target)

        script_path = os.path.join(os.path.dirname(prog), "bin")
        prog_name = os.path.basename(prog)

        # Determine platform-specific binary names
        suffix = '.exe' if not os_is_nix else ''
        align_bin_s = align_bin_s_base + suffix
        align_bin_l = align_bin_l_base + suffix
        build_bin_exe = build_bin_base + suffix

        align_prog_s = os.path.join(script_path, align_bin_s)
        align_prog_l = os.path.join(script_path, align_bin_l)
        build_bin = os.path.join(script_path, build_bin_exe)

    except Exception as e:
        _fail(f"Error resolving script path: {e}")

def _which(program):
    """Checks if a program exists in the system PATH."""
    return shutil.which(program) is not None

def _get_bt2_desc(align_prog_path):
    """Gets argument descriptions from the bowtie2-align binary."""
    desc = {}
    cmd = [align_prog_path, "--wrapper", "basic-0", "--arg-desc"]
    try:
        _info(f"Getting arg descriptions with: {' '.join(shlex.quote(c) for c in cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                opt, opt_type = line.split('\t', 1)
                desc[opt] = opt_type
            except ValueError:
                _info(f"Could not parse description line: {line}")
        return desc
    except FileNotFoundError:
        _fail(f"Failed to run command: '{' '.join(shlex.quote(c) for c in cmd)}'. Bowtie 2 aligner not found at '{align_prog_path}'.")
    except subprocess.CalledProcessError as e:
        _error(f"Command '{' '.join(shlex.quote(c) for c in cmd)}' failed with exit code {e.returncode}")
        _error(f"Stderr:\n{e.stderr}")
        _fail("Getting description of arguments failed!")
    except Exception as e:
         _fail(f"An unexpected error occurred while getting descriptions: {e}")

def _cat_file(input_filename, output_fh):
    """Copies content of potentially compressed file to output file handle."""
    cmd = None
    uncompressor_proc = None
    try:
        if input_filename.endswith(".gz"):
            if not _which("gzip"): _fail("gzip command not found, needed for .gz input.")
            cmd = ["gzip", "-dc", input_filename]
            uncompressor_proc = subprocess.Popen(cmd, stdout=output_fh, stderr=subprocess.PIPE)
        elif input_filename.endswith(".bz2"):
            if not _which("bzip2"): _fail("bzip2 command not found, needed for .bz2 input.")
            cmd = ["bzip2", "-dc", input_filename]
            uncompressor_proc = subprocess.Popen(cmd, stdout=output_fh, stderr=subprocess.PIPE)
        elif input_filename.endswith(".lz4"):
            if not _which("lz4"): _fail("lz4 command not found, needed for .lz4 input.")
            cmd = ["lz4", "-dc", input_filename]
            uncompressor_proc = subprocess.Popen(cmd, stdout=output_fh, stderr=subprocess.PIPE)
        elif input_filename.endswith(".zst"):
            if not _which("zstd"): _fail("zstd command not found, needed for .zst input.")
            cmd = ["zstd", "-dfc", input_filename]
            uncompressor_proc = subprocess.Popen(cmd, stdout=output_fh, stderr=subprocess.PIPE)
        else:
            # Plain file
            with open(input_filename, 'rb') as ifh:
                shutil.copyfileobj(ifh, output_fh)
            return # No process to wait for

        # For compressed files, wait for the process and check errors
        stdout, stderr = uncompressor_proc.communicate()
        if uncompressor_proc.returncode != 0:
            _error(f"Decompression command failed: {' '.join(cmd)}")
            _error(f"Stderr: {stderr.decode(errors='ignore')}")
            _fail(f"Could not process read file: {input_filename}")

    except FileNotFoundError:
         _fail(f"Input file not found: {input_filename}")
    except subprocess.SubprocessError as e:
         _fail(f"Subprocess error for {input_filename}: {e}")
    except IOError as e:
         _fail(f"I/O error processing {input_filename}: {e}")
    except Exception as e:
        _fail(f"Unexpected error in _cat_file for {input_filename}: {e}")


def _write_files_to_pipe_or_file(input_files, output_path, use_pipes):
    """Writes multiple files to a named pipe or regular file. Uses fork if piping."""
    global to_delete, to_kill

    if use_pipes:
        try:
            os.mkfifo(output_path, 0o700)
            _info(f"Created FIFO: {output_path}")
            to_delete.append(output_path) # Ensure FIFO is removed
        except OSError as e:
            _fail(f"mkfifo({output_path}) failed: {e}")
        except FileExistsError:
             _info(f"FIFO already exists: {output_path}. Attempting to use.")
             # Be cautious here, an old FIFO might cause issues.
             # Consider removing and recreating if necessary.
             pass # Pipe already exists, maybe from a previous run?

        # Fork a child process to write to the pipe
        try:
            pid = os.fork()
        except OSError as e:
            _fail(f"Failed to fork process for writing to pipe {output_path}: {e}")

        if pid == 0:
            # --- Child Process ---
            # Important: Child should exit cleanly or with error indication
            # No interaction with global state (to_kill, to_delete) from child needed
            child_exit_code = 1 # Assume failure
            output_fh = None
            try:
                # Opening in binary mode 'wb' is generally safer
                output_fh = open(output_path, "wb")
                for fn in input_files:
                    _cat_file(fn, output_fh)
                child_exit_code = 0 # Success
            except Exception as e:
                # Use basic print to stderr in child as logging might be complex
                print(f"(CHILD ERR): Error writing to pipe {output_path}: {e}", file=sys.stderr)
                child_exit_code = 1
            finally:
                if output_fh:
                    try:
                        output_fh.close()
                    except IOError as e:
                         print(f"(CHILD ERR): Error closing pipe {output_path}: {e}", file=sys.stderr)
                         child_exit_code = 1 # Mark as failure if close fails
                os._exit(child_exit_code) # Use os._exit in child after fork

        else:
            # --- Parent Process ---
            _info(f"Launched child process {pid} to write to {output_path}")
            to_kill.append(pid) # Register child PID for cleanup
            # Parent continues immediately, doesn't wait for child

    else: # Not using pipes, write directly to a temporary file
        _info(f"Writing concatenated files to temporary file: {output_path}")
        try:
            with open(output_path, "wb") as ofh:
                 for fn in input_files:
                     _cat_file(fn, ofh)
            to_delete.append(output_path) # Mark temp file for deletion
        except IOError as e:
            _fail(f"Can't open temporary file '{output_path}' for writing: {e}")
        except Exception as e:
            _fail(f"Error writing to temporary file {output_path}: {e}")


def _maybe_wrap_input(args, use_pipes, temp_dir):
    """Checks input lists for compressed files and replaces them with pipes/temp files."""
    updated_args = {}
    needs_wrapping = False

    # Identify file lists and check for compression types requiring wrapping
    file_lists = {
        "m1": args.get("mate1s", []) or [],
        "m2": args.get("mate2s", []) or [],
        "unpaired": args.get("unpaired", []) or [],
        "tab5": args.get("tab5_mates", []) or [],
        "tab6": args.get("tab6_mates", []) or [],
        "interleaved": args.get("interleaved_mates", []) or [],
        "bam": args.get("bam_files", []) or []
    }

    processed_files = {key: [] for key in file_lists}
    wrapped_files = {key: [] for key in file_lists}

    # Separate files needing wrapping (.bz2, .lz4, .zst) from others
    for key, files in file_lists.items():
        if not files: continue # Skip empty lists like m2 for unpaired
        for fn in files:
            # Bowtie internaly handles .gz, others need wrapping by this script
            if fn.endswith(".bz2") or fn.endswith(".lz4") or fn.endswith(".zst"):
                wrapped_files[key].append(fn)
                needs_wrapping = True
            else:
                processed_files[key].append(fn)

    if not needs_wrapping:
        _info("No input files require decompression wrapping by this script.")
        return {key: value for key, value in file_lists.items() if value} # Return original args

    _info("Some input files require decompression wrapping.")

    # Create temp files/pipes for wrapped inputs
    pid = os.getpid()
    temp_file_prefix = os.path.join(temp_dir, f"{host}_{pid}")

    # Handle unpaired reads (-U, --tab5, --tab6, --interleaved, -b)
    unpaired_keys = ["unpaired", "tab5", "tab6", "interleaved", "bam"]
    for key in unpaired_keys:
        if wrapped_files[key]:
            pipe_or_file_name = f"{temp_file_prefix}.{key}.input"
            _write_files_to_pipe_or_file(wrapped_files[key], pipe_or_file_name, use_pipes)
            processed_files[key].append(pipe_or_file_name) # Add pipe/file to list

    # Handle paired reads (-1 / -2)
    if wrapped_files["m1"] or wrapped_files["m2"]:
         # Must wrap both if either needs it, maintaining pairing
        if len(wrapped_files["m1"]) != len(wrapped_files["m2"]):
              # This case implies some pairs need wrapping, some don't.
              # The original Perl script's logic seems to lump *all* pairs needing
              # wrapping into single pipes. Let's replicate that.
              _info("Wrapping all mate 1 files needing decompression into one pipe/file.")
              pipe_or_file_name_1 = f"{temp_file_prefix}.m1.input"
              _write_files_to_pipe_or_file(wrapped_files["m1"], pipe_or_file_name_1, use_pipes)
              processed_files["m1"].append(pipe_or_file_name_1)

              _info("Wrapping all mate 2 files needing decompression into one pipe/file.")
              pipe_or_file_name_2 = f"{temp_file_prefix}.m2.input"
              _write_files_to_pipe_or_file(wrapped_files["m2"], pipe_or_file_name_2, use_pipes)
              processed_files["m2"].append(pipe_or_file_name_2)
        else:
              # Should not happen if len(m1) == len(m2) initially and we grouped correctly
              _error("Logic error in paired read wrapping.") # Should review this case


    # Return the updated lists of files (original + pipes/temps)
    return {key: value for key, value in processed_files.items() if value}


def _extract_index_name(bt2_args, ref_str_provided):
    """Finds the index base name from arguments (-x or --index)."""
    index_opt = '--index' if ref_str_provided else '-x'
    idx_basename = None
    try:
        idx = bt2_args.index(index_opt)
        if idx + 1 < len(bt2_args):
            idx_basename = bt2_args[idx+1]
        else:
             _fail(f"Option {index_opt} requires an argument.")
    except ValueError:
        _info(f"Cannot find index option ({index_opt}) in the command line.")
        return None # No index specified directly

    # Check for existence, possibly using BOWTIE2_INDEXES
    search_paths = [idx_basename]
    bowtie2_indexes_env = os.environ.get("BOWTIE2_INDEXES")
    if bowtie2_indexes_env:
        search_paths.append(os.path.join(bowtie2_indexes_env, os.path.basename(idx_basename)))
        # Maybe also consider the full path if idx_basename contained directory components?
        # search_paths.append(os.path.join(bowtie2_indexes_env, idx_basename)) # Less common?

    found_path = None
    for potential_basename in search_paths:
        # Check for either small or large index files existence
        # Use glob to find any matching index file component
        if glob.glob(potential_basename + "*.bt2") or glob.glob(potential_basename + "*.bt2l"):
             _info(f"Found Bowtie 2 index files for base: {potential_basename}")
             found_path = potential_basename
             break # Found it

    if found_path is None:
        _fail(f"Could not find Bowtie 2 index files matching base name '{idx_basename}'" +
              (f" (also checked in BOWTIE2_INDEXES='{bowtie2_indexes_env}')" if bowtie2_indexes_env else ""))

    # Important: Update the argument list if we found the index via BOWTIE2_INDEXES
    if found_path != idx_basename:
         _info(f"Using index path found in BOWTIE2_INDEXES: {found_path}")
         bt2_args[idx+1] = found_path # Modify the list in place

    return found_path


def _handle_passthrough_output(bt_proc, out_fh, read_output_config, no_unal_filter, bam_mode, cap_out_fn):
    """Processes bowtie2-align output when --passthrough is used."""
    global read_fns, read_compress

    read_fhs = {}
    fhs_to_close = []
    compressor_procs = {} # Store compressor Popen objects if needed

    try:
        # --- Setup output file handles for --al/--un etc. ---
        for read_type, base_fn in read_output_config.items():
            if not base_fn: continue # Skip if option not provided

            compression = read_compress.get(read_type)
            vol, base_spec_dir, base_fname = "", os.path.dirname(base_fn), os.path.basename(base_fn)
            if not base_spec_dir and not base_fname and os.path.isdir(base_fn):
                 base_spec_dir = base_fn
                 base_fname = None # Use default names like 'al-seqs' etc.


            if read_type.endswith("-conc") or read_type == "un-mates":
                # Paired output: needs mate 1 and mate 2 files
                fn1_raw, fn2_raw = None, None
                default_prefix = read_type + ('-mate' if not read_type.endswith('-mates') else '')

                if base_fname:
                    # Use pattern from user input
                    if '%' in base_fname:
                        fn1_raw = base_fname.replace('%', '1', 1) # Replace only first %
                        fn2_raw = base_fname.replace('%', '2', 1)
                        if fn1_raw == fn2_raw: # If only one % or none
                            base_name_part, ext_part = os.path.splitext(base_fname)
                            fn1_raw = f"{base_name_part}.1{ext_part}"
                            fn2_raw = f"{base_name_part}.2{ext_part}"
                    else:
                        # Append .1/.2 before the extension or at the end
                        base_name_part, ext_part = os.path.splitext(base_fname)
                        if ext_part:
                            fn1_raw = f"{base_name_part}.1{ext_part}"
                            fn2_raw = f"{base_name_part}.2{ext_part}"
                        else:
                            fn1_raw = f"{base_fname}.1"
                            fn2_raw = f"{base_fname}.2"
                else:
                    # Use default names
                    fn1_raw = f"{default_prefix}.1"
                    fn2_raw = f"{default_prefix}.2"

                fn1 = os.path.join(base_spec_dir, fn1_raw)
                fn2 = os.path.join(base_spec_dir, fn2_raw)

                if fn1 == fn2:
                     _fail(f"Generated mate output filenames are identical: '{fn1}'. Check pattern in --{read_type}.")

                # Open file handles (potentially through compressors)
                fh1 = _open_compressed_output(fn1, compression, read_type + "/1", compressor_procs)
                fh2 = _open_compressed_output(fn2, compression, read_type + "/2", compressor_procs)

                read_fhs[read_type] = {1: fh1, 2: fh2}
                fhs_to_close.extend([fh1, fh2])

            else: # Single-end output (--al, --un)
                if base_fname:
                     fn = os.path.join(base_spec_dir, base_fname)
                else:
                     fn = os.path.join(base_spec_dir, f"{read_type}-seqs")

                fh = _open_compressed_output(fn, compression, read_type, compressor_procs)
                read_fhs[read_type] = fh
                fhs_to_close.append(fh)

        # --- Process bowtie2-align output stream ---
        sam_header_processed = False
        while True:
            line_bytes = bt_proc.stdout.readline()
            if not line_bytes:
                break # End of stream
            line = line_bytes.decode('utf-8', errors='replace').rstrip('\n')

            if not line: continue # Skip empty lines

            # Handle SAM Header
            if line.startswith("@"):
                if not sam_header_processed and bam_mode:
                     # In BAM mode with passthrough, the header goes to samtools via out_fh
                     # No filtering applied to header lines.
                     pass # Let it fall through to the print below
                elif sam_header_processed and bam_mode:
                     # Should not happen if header is contiguous
                     _info("Warning: Encountered SAM header line after alignments.")
                # For non-BAM or if header already passed, just print if needed
                sam_header_processed = True
                print(line, file=out_fh) # Write header to main output/samtools
                continue

            # Handle Alignment lines
            sam_header_processed = True # Mark header done after first non-@ line
            fields = line.split('\t')
            if len(fields) < 11:
                _info(f"Warning: Skipping malformed SAM line: {line}")
                continue

            try:
                flag = int(fields[1])
            except ValueError:
                _info(f"Warning: Skipping SAM line with non-integer FLAG: {line}")
                continue

            is_secondary = (flag & 0x100) != 0
            is_unmapped = (flag & 0x4) != 0

            # Apply --no-unal filter
            if no_unal_filter and is_unmapped:
                # Need to consume the extra sequence line if present
                if not is_secondary and read_output_config: # Check if read outputs requested
                    _ = bt_proc.stdout.readline() # Read and discard sequence line
                continue # Skip writing this unaligned read to main output

            # Write alignment line to main output (or samtools pipe)
            print(line, file=out_fh)

            # Handle passthrough data for --al/--un etc. if not secondary alignment
            if not is_secondary and read_output_config:
                extra_line_bytes = bt_proc.stdout.readline()
                if not extra_line_bytes:
                    _error("Error: Unexpected end of stream after SAM record when expecting sequence data.")
                    break
                # Decode the URL-encoded sequence line
                try:
                     # Need to handle potential %0A (newline) etc.
                     # The Perl script uses: $l =~ s/%(..)/chr(hex($1))/eg;
                     # In Python, urllib.parse.unquote handles %xx decoding.
                     # We decode from bytes assuming UTF-8 or similar; FASTA/Q is often ASCII/UTF-8.
                     decoded_seq_line = unquote(extra_line_bytes.decode('utf-8', errors='replace').rstrip('\n'))

                except Exception as e:
                     _error(f"Error decoding passthrough sequence line: {extra_line_bytes!r} - {e}")
                     decoded_seq_line = extra_line_bytes.decode('utf-8', errors='replace').rstrip('\n') # Fallback


                # --- Distribute sequence to appropriate files based on flags ---
                is_paired = (flag & 0x1) != 0
                is_mate1 = (flag & 0x40) != 0
                is_mate2 = (flag & 0x80) != 0
                is_proper_pair = (flag & 0x2) != 0 # Concordant

                # 1. Unpaired reads (--al, --un)
                if not is_paired:
                    if is_unmapped:
                        if "un" in read_fhs: print(decoded_seq_line, file=read_fhs["un"])
                    else:
                        if "al" in read_fhs: print(decoded_seq_line, file=read_fhs["al"])

                # 2. Paired reads (--al-conc, --un-conc)
                elif is_paired and ("al-conc" in read_fhs or "un-conc" in read_fhs):
                    target_type = "al-conc" if is_proper_pair else "un-conc"
                    mate_num = 1 if is_mate1 else (2 if is_mate2 else None)
                    if mate_num and target_type in read_fhs:
                         print(decoded_seq_line, file=read_fhs[target_type][mate_num])

                # 3. Paired reads (--un-mates) - only those where *both* mates are unmapped
                # Note: Bowtie's passthrough might send reads here even if the *other* mate mapped.
                # The flag 0x8 (mate_unmapped) indicates the status of the *other* mate.
                # So, we write if *this* read is unmapped (flag 0x4) AND it's part of a pair (flag 0x1).
                # The original Perl logic `if (defined($read_fhs{"un-mates"}) && $pair && $unal && ($fl & 4) != 0)` seems to imply this.
                elif is_paired and is_unmapped and "un-mates" in read_fhs:
                     mate_num = 1 if is_mate1 else (2 if is_mate2 else None)
                     if mate_num:
                          print(decoded_seq_line, file=read_fhs["un-mates"][mate_num])


    finally:
        # Close all the extra output file handles and wait for compressors
        for fh in fhs_to_close:
            try:
                fh.close()
            except IOError as e:
                 _error(f"Error closing output file handle: {e}")

        for cmd, proc in compressor_procs.items():
             try:
                 retcode = proc.wait()
                 if retcode != 0:
                     _error(f"Compressor command '{' '.join(cmd)}' exited with code {retcode}")
                     # Stderr might have been captured elsewhere or piped from bowtie
             except Exception as e:
                  _error(f"Error waiting for compressor '{' '.join(cmd)}': {e}")


def _open_compressed_output(filename, compression, stream_name, compressor_procs):
    """Opens a file for writing, potentially through a compression pipe."""
    fh = None
    cmd = None
    mode = "wt" # Text mode by default for sequences

    if compression == "gzip":
        if not _which("gzip"): _fail("gzip command not found, needed for compressed output.")
        cmd = ["gzip", "-c"]
        mode = "wb" # Pipe expects bytes
    elif compression == "bzip2":
        if not _which("bzip2"): _fail("bzip2 command not found, needed for compressed output.")
        cmd = ["bzip2", "-c"]
        mode = "wb"
    elif compression == "lz4":
        if not _which("lz4"): _fail("lz4 command not found, needed for compressed output.")
        cmd = ["lz4", "-c"] # Note: lz4 might need levels/options
        mode = "wb"
    elif compression == "zstd":
        if not _which("zstd"): _fail("zstd command not found, needed for compressed output.")
        cmd = ["zstd", "-c"] # Note: zstd might need levels/options
        mode = "wb"

    try:
        if cmd:
            # Open the final destination file first in binary write mode
            file_fh = open(filename, "wb")
            # Start the compressor, piping its output to the file
            proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=file_fh, stderr=subprocess.PIPE)
            compressor_procs[tuple(cmd + [stream_name])] = proc # Use tuple as key
            # Return the compressor's stdin pipe (needs to be text-encoded)
            return io.TextIOWrapper(proc.stdin, encoding='utf-8', errors='replace') # Wrap stdin for text writing
            # Important: Need to close this TextIOWrapper, which flushes and closes proc.stdin.
            # Also need to wait for 'proc' later.
        else:
            # Open regular file in text mode
            return open(filename, "w", encoding='utf-8', errors='replace')

    except IOError as e:
        _fail(f"Could not open output file '{filename}' for --{stream_name}: {e}")
    except subprocess.SubprocessError as e:
        _fail(f"Could not start compressor for --{stream_name} ('{' '.join(cmd)}'): {e}")
    except Exception as e:
        _fail(f"Unexpected error opening output for --{stream_name}: {e}")


# --- Main Execution ---
def main():
    import io # Needed for TextIOWrapper

    _resolve_script_path() # Determine script location and binary paths

    # --- Argument Parser Setup ---
    parser = argparse.ArgumentParser(
        description="Wrapper for the Bowtie 2 aligner.",
        add_help=False # Handle help manually to pass it to aligner if needed
    )

    # --- Wrapper-Specific Arguments ---
    wrapper_group = parser.add_argument_group("Wrapper Options")
    wrapper_group.add_argument("--temp-directory", default=temp_dir_base,
                               help="Directory for temporary files/pipes (default: system temp)")
    wrapper_group.add_argument("--bam", action="store_true", default=False,
                               help="Output alignments in BAM format (requires samtools)")
    wrapper_group.add_argument("--no-named-pipes", action="store_true", default=False,
                               help="Use temporary files instead of named pipes for decompression")
    wrapper_group.add_argument("--ref-string", "--reference-string",
                               help="Build temporary index from this reference string")
    wrapper_group.add_argument("--keep", action="store_true", default=False,
                               help="Keep temporary index files created by --ref-string")
    wrapper_group.add_argument("--verbose", action="store_true", default=False,
                               help="Print verbose messages from the wrapper script")
    wrapper_group.add_argument("--log-file",
                               help="Redirect Bowtie 2's stderr messages to this file")
    # Hidden/internal options
    # parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS) # For bowtie2-align-debug
    # parser.add_argument("--sanitized", action="store_true", help=argparse.SUPPRESS) # For bowtie2-align-sanitized

    # --- Help Argument ---
    parser.add_argument("-h", "--help", action="store_true", help="Show this help message or Bowtie 2's help")


    # --- Parse known args first ---
    # Use parse_known_args to separate wrapper args from bowtie2 args
    # We need all args initially to handle the --un/--al style options
    all_args = sys.argv[1:]

    # --- Pre-parse for special output options and help ---
    # Manually look for help flag before full parsing
    if "-h" in all_args or "--help" in all_args:
         # Check if it's the *only* argument or only wrapper args are present
         is_wrapper_help_request = True
         temp_parser = argparse.ArgumentParser(add_help=False)
         # Add only wrapper-specific args to temp parser to check
         temp_parser.add_argument("--temp-directory")
         temp_parser.add_argument("--bam", action="store_true")
         temp_parser.add_argument("--no-named-pipes", action="store_true")
         temp_parser.add_argument("--ref-string")
         temp_parser.add_argument("--reference-string")
         temp_parser.add_argument("--keep", action="store_true")
         temp_parser.add_argument("--verbose", action="store_true")
         temp_parser.add_argument("--log-file")
         temp_parser.add_argument("-h", "--help", action="store_true")
         try:
             # Try parsing *only* known wrapper args. If anything is left, it's likely for bowtie2
             _, remaining = temp_parser.parse_known_args(all_args)
             if remaining:
                  is_wrapper_help_request = False
         except: # Any parsing error means non-wrapper args likely present
             is_wrapper_help_request = False


         if is_wrapper_help_request:
              parser.print_help()
              # Also print Bowtie 2 basic usage?
              print("\n--- Bowtie 2 Aligner Help (run with more options for details) ---")
              bt2_help_cmd = [align_prog_s, "--help"] # Use default small aligner for help
              try:
                   subprocess.run(bt2_help_cmd, check=True)
              except Exception as e:
                   print(f"Could not retrieve Bowtie 2 help: {e}", file=sys.stderr)
              sys.exit(0)
         else:
              # Pass '-h' or '--help' through to bowtie2-align
              # Let argparse handle other wrapper args later if needed
               pass # Help will be handled by bowtie2 binary


    # --- Handle --un*/--al* options manually (similar to Perl callback) ---
    passthrough_output_options = {} # Stores {type: filename}
    temp_all_args = [] # Build a new list excluding these special args
    i = 0
    while i < len(all_args):
        arg = all_args[i]
        # Regex to match pattern like --un-gz=filename or --al-conc-bz2=filename
        match = re.match(r"^--(al|un)((?:-conc|-mates)?)(?:-(gz|bz2|lz4|zst))?=(.*)", arg)
        if match:
            base_type, modifier, compression, filename = match.groups()
            read_type = base_type + (modifier if modifier else "")
            passthrough_output_options[read_type] = filename
            read_compress[read_type] = ""
            if compression == "gz": read_compress[read_type] = "gzip"
            elif compression == "bz2": read_compress[read_type] = "bzip2"
            elif compression == "lz4": read_compress[read_type] = "lz4"
            elif compression == "zst": read_compress[read_type] = "zstd"
            # Don't add this arg to temp_all_args
            i += 1 # Skip this arg
        elif arg.startswith("--un=") or arg.startswith("--al=") or \
             arg.startswith("--un-gz=") or arg.startswith("--al-gz=") or \
             arg.startswith("--un-bz2=") or arg.startswith("--al-bz2=") or \
             arg.startswith("--un-lz4=") or arg.startswith("--al-lz4=") or \
             arg.startswith("--un-zst=") or arg.startswith("--al-zst=") or \
             arg.startswith("--un-conc=") or arg.startswith("--al-conc=") or \
             arg.startswith("--un-conc-gz=") or arg.startswith("--al-conc-gz=") or \
             arg.startswith("--un-conc-bz2=") or arg.startswith("--al-conc-bz2=") or \
             arg.startswith("--un-conc-lz4=") or arg.startswith("--al-conc-lz4=") or \
             arg.startswith("--un-conc-zst=") or arg.startswith("--al-conc-zst=") or \
             arg.startswith("--un-mates=") or arg.startswith("--al-mates=") or \
             arg.startswith("--un-mates-gz=") or arg.startswith("--al-mates-gz=") or \
             arg.startswith("--un-mates-bz2=") or arg.startswith("--al-mates-bz2=") or \
             arg.startswith("--un-mates-lz4=") or arg.startswith("--al-mates-lz4=") or \
             arg.startswith("--un-mates-zst=") or arg.startswith("--al-mates-zst="):
             # Handle case without '=' if Getopt::Long allowed that (unlikely with '=s')
             # This block primarily handles the explicit =value form detected by regex above.
             # If space separation was allowed, need more complex logic here.
             # Assuming '=s' in Perl requires the equals sign or next arg.
             # The regex handles the '=' case. If separated by space:
             match_no_eq = re.match(r"^--(al|un)((?:-conc|-mates)?)(?:-(gz|bz2|lz4|zst))?$", arg)
             if match_no_eq and i + 1 < len(all_args):
                 base_type, modifier, compression = match_no_eq.groups()
                 filename = all_args[i+1]
                 read_type = base_type + (modifier if modifier else "")
                 passthrough_output_options[read_type] = filename
                 read_compress[read_type] = ""
                 if compression == "gz": read_compress[read_type] = "gzip"
                 elif compression == "bz2": read_compress[read_type] = "bzip2"
                 elif compression == "lz4": read_compress[read_type] = "lz4"
                 elif compression == "zst": read_compress[read_type] = "zstd"
                 i += 2 # Skip this arg and its value
             else:
                  # Argument doesn't match expected pattern or is missing value
                  temp_all_args.append(arg)
                  i += 1

        else:
            temp_all_args.append(arg)
            i += 1

    # --- Parse remaining arguments ---
    # Now parse the filtered list using parse_known_args
    args, bt2_args = parser.parse_known_args(temp_all_args)

    # --- Set Global Flags ---
    global verbose_mode, keep_temps, use_pipes, temp_dir
    verbose_mode = args.verbose
    keep_temps = args.keep
    use_pipes = not args.no_named_pipes
    temp_dir = args.temp_directory

    # --- Ensure temp dir exists ---
    if not os.path.isdir(temp_dir):
        try:
            os.makedirs(temp_dir, exist_ok=True)
            _info(f"Created temporary directory: {temp_dir}")
        except OSError as e:
            _fail(f"Could not create temporary directory {temp_dir}: {e}")


    # --- Handle Log File ---
    if args.log_file:
        try:
            _info(f"Redirecting Bowtie 2 stderr to: {args.log_file}")
            # Duplicate original stderr
            stderr_fileno = sys.stderr.fileno()
            old_stderr_fileno = os.dup(stderr_fileno)
            old_stderr = os.fdopen(old_stderr_fileno, 'w')

            # Open log file and redirect stderr
            log_file_handle = open(args.log_file, "w")
            os.dup2(log_file_handle.fileno(), stderr_fileno)
            # sys.stderr = log_file_handle # This also works but dup2 is more direct
        except OSError as e:
            old_stderr = None # Ensure cleanup doesn't try to use invalid handle
            log_file_handle = None
            _fail(f"Cannot redirect stderr to log file {args.log_file}: {e}")
        except Exception as e:
             _fail(f"Unexpected error setting up log file {args.log_file}: {e}")


    # --- Check if aligner binary exists ---
    # Default to small index aligner for initial checks
    if not os.path.exists(align_prog_s) or not os.access(align_prog_s, os.X_OK):
        _fail(f"Expected bowtie2 aligner ({align_bin_s_base}) to be in the same directory as the script:\n{script_path}\nAligner path checked: {align_prog_s}")

    # --- Get Bowtie 2 Argument Descriptions ---
    # Needed later? The Perl script uses this but doesn't seem essential for the logic
    # bt2_arg_desc = _get_bt2_desc(align_prog_s) # Can uncomment if needed

    # --- Process Input Files (Decompression Wrapping) ---
    # Reconstruct input args from bt2_args for _maybe_wrap_input
    # This is complex because parse_known_args doesn't give us the nice structure Getopt::Long did.
    # We need to find -1, -2, -U etc. in bt2_args.
    # A simpler approach: Assume wrap_input is needed *only* if bz2/lz4/zst files exist.
    input_args_struct = {
        "mate1s": [], "mate2s": [], "unpaired": [], "tab5_mates": [],
        "tab6_mates": [], "interleaved_mates": [], "bam_files": []
    }
    i = 0
    while i < len(bt2_args):
        arg = bt2_args[i]
        val = bt2_args[i+1] if i + 1 < len(bt2_args) else None
        arg_processed = False
        if arg == "-1" and val:
            input_args_struct["mate1s"].extend(val.split(','))
            arg_processed = True
        elif arg == "-2" and val:
            input_args_struct["mate2s"].extend(val.split(','))
            arg_processed = True
        elif (arg == "-U" or arg == "--reads") and val:
            input_args_struct["unpaired"].extend(val.split(','))
            arg_processed = True
        elif arg == "--tab5" and val:
            input_args_struct["tab5_mates"].extend(val.split(','))
            arg_processed = True
        elif arg == "--tab6" and val:
            input_args_struct["tab6_mates"].extend(val.split(','))
            arg_processed = True
        elif arg == "--interleaved" and val:
            input_args_struct["interleaved_mates"].extend(val.split(','))
            arg_processed = True
        elif arg == "-b" and val:
            input_args_struct["bam_files"].extend(val.split(','))
            arg_processed = True

        if arg_processed:
             # Remove these from bt2_args as we will re-add them later
             del bt2_args[i:i+2]
        else:
             i += 1 # Move to next arg if this one wasn't an input file list

    # Perform the wrapping/decompression if needed
    updated_input_files = _maybe_wrap_input(input_args_struct, use_pipes, temp_dir)

    # Re-add input file arguments to bt2_args with potentially updated filenames
    if updated_input_files.get("mate1s"):
         bt2_args.extend(["-1", ",".join(updated_input_files["mate1s"])])
    if updated_input_files.get("mate2s"):
         bt2_args.extend(["-2", ",".join(updated_input_files["mate2s"])])
    if updated_input_files.get("unpaired"):
         bt2_args.extend(["-U", ",".join(updated_input_files["unpaired"])])
    if updated_input_files.get("tab5_mates"):
         bt2_args.extend(["--tab5", ",".join(updated_input_files["tab5_mates"])])
    if updated_input_files.get("tab6_mates"):
         bt2_args.extend(["--tab6", ",".join(updated_input_files["tab6_mates"])])
    if updated_input_files.get("interleaved_mates"):
         bt2_args.extend(["--interleaved", ",".join(updated_input_files["interleaved_mates"])])
    if updated_input_files.get("bam_files"):
          bt2_args.extend(["-b", ",".join(updated_input_files["bam_files"])])


    # --- Handle --ref-string ---
    if args.ref_string:
        _info("Building temporary index from --ref-string")
        if not os.path.exists(build_bin) or not os.access(build_bin, os.X_OK):
             _fail(f"bowtie2-build binary not found or not executable at: {build_bin}")

        try:
            # Create a temporary FASTA file
            # Use NamedTemporaryFile for easier management
            with tempfile.NamedTemporaryFile(mode='w', suffix=".ref_str.fa", dir=temp_dir, delete=False) as tf:
                tf.write(">ref_from_string\n")
                tf.write(args.ref_string + "\n")
                temp_fasta_name = tf.name
            _info(f"Created temporary FASTA: {temp_fasta_name}")
            to_delete.append(temp_fasta_name) # Ensure FASTA file is deleted

            # Build the index
            index_basename = temp_fasta_name # Use FASTA filename as base for index
            build_cmd = [build_bin, temp_fasta_name, index_basename]
            _info(f"Running bowtie2-build: {' '.join(shlex.quote(c) for c in build_cmd)}")
            build_result = subprocess.run(build_cmd, capture_output=True, text=True)

            if build_result.returncode != 0:
                 _error(f"bowtie2-build failed (code {build_result.returncode})")
                 _error(f"Stdout:\n{build_result.stdout}")
                 _error(f"Stderr:\n{build_result.stderr}")
                 _fail("bowtie2-build returned non-0 exit status.")

            # Add --index argument for bowtie2-align
            bt2_args.extend(["--index", index_basename])

            # Register index files for deletion
            # Determine expected index extension (start with small)
            current_idx_ext = idx_ext_s
            # Check if large index files were created (unlikely for small string, but check)
            if os.path.exists(f"{index_basename}.1.{idx_ext_l}"):
                current_idx_ext = idx_ext_l

            index_files = [
                f"{index_basename}.1.{current_idx_ext}", f"{index_basename}.2.{current_idx_ext}",
                f"{index_basename}.3.{current_idx_ext}", f"{index_basename}.4.{current_idx_ext}",
                f"{index_basename}.rev.1.{current_idx_ext}", f"{index_basename}.rev.2.{current_idx_ext}"
            ]
            for idx_f in index_files:
                 if os.path.exists(idx_f): # Only add files that were actually created
                     to_delete.append(idx_f)
                 else:
                      _info(f"Warning: Expected index file not found after build: {idx_f}")


        except (IOError, OSError, subprocess.SubprocessError) as e:
             _fail(f"Error handling --ref-string: {e}")
        except Exception as e:
            _fail(f"Unexpected error during --ref-string processing: {e}")


    # --- Determine Aligner Binary and Index Extension ---
    align_prog = align_prog_s # Default to small
    idx_ext = idx_ext_s
    large_index_flag = False # Did user explicitly request --large-index? (Need to add this arg)

    # Find index name from args (needed to check for .bt2l)
    index_basename = _extract_index_name(bt2_args, bool(args.ref_string))

    # Check if large index exists or is requested
    # Need to add --large-index to argparse if we want explicit user control
    force_large = False # Replace with check for args.large_index if added
    if index_basename:
         has_large = os.path.exists(f"{index_basename}.1.{idx_ext_l}")
         has_small = os.path.exists(f"{index_basename}.1.{idx_ext_s}")

         if force_large:
             if not has_large:
                 _fail(f"User specified --large-index, but cannot find large index files ({index_basename}.1.{idx_ext_l})")
             _info(f"Using large index (forced by user): {index_basename}.*.{idx_ext_l}")
             align_prog = align_prog_l
             idx_ext = idx_ext_l
         elif has_large and not has_small:
             _info(f"Small index not found, but large index detected. Switching to large index: {index_basename}.*.{idx_ext_l}")
             align_prog = align_prog_l
             idx_ext = idx_ext_l
             # Add --large-index flag for the aligner itself? Check bowtie2 docs.
             # bt2_args.append("--large-index") # If aligner needs the flag too
         elif not has_small and not has_large:
              # _extract_index_name should have already failed if neither exists
              _fail(f"Cannot find any index files (.bt2 or .bt2l) for base: {index_basename}")
         else:
             # Small index exists (and maybe large too, but default is small)
              _info(f"Using small index: {index_basename}.*.{idx_ext_s}")
              align_prog = align_prog_s
              idx_ext = idx_ext_s
    else:
        # No index specified via -x or --ref-string, aligner might use default or fail
        _info("No index specified with -x or --ref-string. Using default aligner.")
        align_prog = align_prog_s # Stick with default small aligner

    # --- Handle Passthrough Logic ---
    passthru_mode = bool(passthrough_output_options) or args.bam # Add --no-unal check later
    cap_out_fn = "-" # Default output is stdout if passthrough needed
    final_bt2_args = list(bt2_args) # Copy args

    # Check for --no-unal (needs passthrough)
    # Need to add --no-unal flag to argparse
    no_unal_filter = False # Replace with check for args.no_unal if added
    if no_unal_filter:
        passthru_mode = True

    if passthru_mode:
        _info("Passthrough mode enabled (due to --bam, --no-unal, or --al/--un options).")
        final_bt2_args.append("--passthrough") # Tell aligner to send extra data

        # Find and remove -S/--output from bowtie2 args, store the filename
        output_arg_indices = [i for i, x in enumerate(final_bt2_args) if x == "-S" or x == "--output"]
        if output_arg_indices:
            idx = output_arg_indices[-1] # Use the last specified one
            if idx + 1 < len(final_bt2_args):
                cap_out_fn = final_bt2_args[idx+1]
                _info(f"Capturing aligner output originally intended for: {cap_out_fn}")
                # Remove -S/--output and its argument
                del final_bt2_args[idx:idx+2]
            else:
                _fail(f"Argument {final_bt2_args[idx]} requires a filename.")
        else:
             _info("Aligner output will be captured from stdout.")
             cap_out_fn = "-" # Explicitly stdout

    # --- Add --verbose to aligner args if requested ---
    if args.verbose:
        final_bt2_args.append("--verbose")

    # --- Construct Final Command ---
    # Handle debug/sanitized suffixes if needed (requires adding args)
    suffix = ""
    # if args.debug: suffix = "-debug"
    # elif args.sanitized: suffix = "-sanitized"
    align_prog_final = align_prog # Could modify with suffix here

    # Check aligner executable again now we know if it's -s or -l
    if not os.path.exists(align_prog_final) or not os.access(align_prog_final, os.X_OK):
         aligner_base = os.path.basename(align_prog_final).replace(".exe", "")
         _fail(f"Required Bowtie 2 aligner binary '{aligner_base}' not found or not executable at '{align_prog_final}'")


    cmd_list = [align_prog_final, "--wrapper", "basic-0"] + final_bt2_args
    cmd_str_for_info = ' '.join(shlex.quote(c) for c in cmd_list)
    _info(f"Bowtie 2 command: {cmd_str_for_info}")


    # --- Execute Bowtie 2 ---
    retcode = 1 # Default to error
    bowtie_proc = None

    try:
        if passthru_mode:
            # --- Passthrough Execution ---
            _info("Running Bowtie 2 aligner and capturing output...")

            # Start the Bowtie 2 process
            bowtie_proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=None, bufsize=1) # Line buffered stdout if possible

            # Determine final output destination
            out_fh = None
            samtools_proc = None
            if args.bam:
                if not _which("samtools"):
                    _fail("samtools command is needed for --bam output but not found in PATH.")
                samtools_cmd = ["samtools", "view", "-b", "-"] # Read SAM from stdin, write BAM to stdout
                if cap_out_fn == "-":
                    _info("Piping SAM output to 'samtools view -b -' -> STDOUT")
                    # Pipe bowtie stdout -> samtools stdin -> script stdout
                    samtools_proc = subprocess.Popen(samtools_cmd, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=None)
                    out_fh = io.TextIOWrapper(samtools_proc.stdin, encoding='utf-8', errors='replace') # Wrap samtools stdin
                else:
                    _info(f"Piping SAM output to 'samtools view -b - > {cap_out_fn}'")
                    try:
                        bam_file_fh = open(cap_out_fn, "wb")
                        to_delete.append(cap_out_fn) # Add to potential cleanup if script fails mid-way
                        samtools_proc = subprocess.Popen(samtools_cmd, stdin=subprocess.PIPE, stdout=bam_file_fh, stderr=None)
                        out_fh = io.TextIOWrapper(samtools_proc.stdin, encoding='utf-8', errors='replace') # Wrap samtools stdin
                    except IOError as e:
                        _fail(f"Error opening BAM output file '{cap_out_fn}': {e}")

            elif cap_out_fn == "-":
                _info("Writing filtered SAM output to STDOUT")
                out_fh = sys.stdout # Write directly to script's stdout
            else:
                _info(f"Writing filtered SAM output to file: {cap_out_fn}")
                try:
                     # Open in text mode, let print handle encoding
                     out_fh = open(cap_out_fn, "w", encoding='utf-8', errors='replace')
                     to_delete.append(cap_out_fn) # Add to potential cleanup
                except IOError as e:
                     _fail(f"Could not open output file '{cap_out_fn}' for writing: {e}")

            # Process the output stream
            _handle_passthrough_output(bowtie_proc, out_fh, passthrough_output_options, no_unal_filter, args.bam, cap_out_fn)

            # Close the final output handle (important for pipes/files)
            # sys.stdout should not be closed here.
            if out_fh is not None and out_fh is not sys.stdout:
                 try:
                     out_fh.close()
                 except Exception as e:
                      _error(f"Error closing output handle: {e}") # e.g. samtools stdin broken pipe

            # Wait for bowtie2 aligner to finish
            bowtie_proc.wait()
            retcode = bowtie_proc.returncode

            # Wait for samtools if it was used
            if samtools_proc:
                samtools_proc.wait()
                if samtools_proc.returncode != 0:
                    _error(f"samtools exited with non-zero status: {samtools_proc.returncode}")
                    # If bowtie succeeded but samtools failed, should we override retcode?
                    if retcode == 0: retcode = samtools_proc.returncode # Propagate samtools error if bowtie was ok

            # If output was written to a file successfully, remove it from deletion list
            if retcode == 0 and cap_out_fn != "-" and cap_out_fn in to_delete:
                 if args.bam or not passthrough_output_options: # Keep if BAM or non-passthrough SAM
                     to_delete.remove(cap_out_fn)


        else:
            # --- Direct Execution (no passthrough) ---
            _info("Running Bowtie 2 aligner directly...")
            # Stderr goes where Python's stderr goes (potentially log file)
            result = subprocess.run(cmd_list, check=False) # check=False to handle non-zero exit below
            retcode = result.returncode

    except KeyboardInterrupt:
         _error("Execution interrupted by user (SIGINT).")
         # Cleanup will run via atexit
         sys.exit(signal.SIGINT) # Exit with SIGINT status
    except Exception as e:
        _error(f"An unexpected error occurred during Bowtie 2 execution: {e}")
        # Perform cleanup
        sys.exit(1)
    # finally:
        # Ensure cleanup runs. atexit should handle this.
        # _cleanup()


    # --- Check Bowtie 2 Exit Status ---
    if retcode < 0:
        # Process died on a signal
        sig_num = abs(retcode)
        sig_name = signal_names.get(sig_num, f"UNKNOWN SIGNAL {sig_num}")
        _error(f"bowtie2-align died with signal {sig_num} ({sig_name})")
        sys.exit(sig_num) # Exit with signal number convention
    elif retcode > 0:
        _error(f"bowtie2-align exited with value {retcode}")
        sys.exit(retcode)
    else:
        _info("Bowtie 2 wrapper finished successfully.")
        sys.exit(0)