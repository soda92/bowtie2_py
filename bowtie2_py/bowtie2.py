#!/usr/bin/env python3

#
# Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
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
# 3. Output directly to bam (not currently supported)

import sys
import os
import getopt
import subprocess
import tempfile
import signal
import platform
from pathlib import Path

host = platform.node()
script_path = Path(__file__).resolve().parent
prog = Path(__file__).resolve()

while prog.is_file() and prog.is_symlink():
    prog = Path(os.readlink(prog)).resolve().parent / Path(prog).name
    script_path = prog.parent

vol = None  # Python doesn't really have volume in the same way
os_is_nix = os.name != 'nt'
align_bin_s = 'bowtie2-align-s' if os_is_nix else 'bowtie2-align-s.exe'
build_bin = 'bowtie2-build' if os_is_nix else 'bowtie2-build.exe'
align_bin_l = 'bowtie2-align-l' if os_is_nix else 'bowtie2-align-l.exe'
align_prog_s = script_path / align_bin_s
align_prog_l = script_path / align_bin_l
align_prog = align_prog_s
idx_ext_l = 'bt2l'
idx_ext_s = 'bt2'
idx_ext = idx_ext_s
signo = {}
signame = []

params_to_quote = {
    '-S': 1, '-U': 1, '-1': 1, '-2': 1, '-x': 1, '-b': 1,
    '--interleaved': 1, '--rg': 1, '--rg-id': 1,
    '--tab5': 1, '--tab6': 1, '--sam-acc': 1, '--bam': 1
}

def quote_params(param_list):
    quoting = 0
    for i in range(len(param_list)):
        if quoting:
            param_list[i] = f'"{param_list[i]}"'
            quoting = 0
            continue
        if param_list[i] in params_to_quote:
            quoting = 1

# Get signal info (Python's signal module is more direct)
for name in signal.valid_signals():
    try:
        signum = getattr(signal, name)
        signo[name] = signum
        if signum < len(signame):
            signame[signum] = name
        else:
            signame.extend([None] * (signum - len(signame) + 1))
            signame[signum] = name
    except AttributeError:
        pass

if not os.access(align_prog, os.X_OK):
    fail(f"Expected bowtie2 to be in same directory with bowtie2-align:\n{script_path}\n")

# Get description of arguments from Bowtie 2
def get_bt2_desc(desc):
    cmd = [str(align_prog), "--wrapper", "basic-0", "--arg-desc"]
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in process.stdout.strip().split('\n'):
            if not line.strip():
                continue
            parts = line.split('\t')
            desc[parts[0]] = parts[1]
    except subprocess.CalledProcessError as e:
        fail(f"Description of arguments failed!\n{e}")
    except FileNotFoundError:
        fail(f"Could not find executable: {cmd[0]}\n")

desc = {}
wrapped = {"1": 1, "2": 1}
get_bt2_desc(desc)

# Given an option like -1, determine whether it's wrapped
def is_wrapped(opt):
    return opt in wrapped

debug = 0
sanitized = 0
read_fns = {}
read_compress = {}
cap_out = None  # Filename for passthrough
no_unal = 0
large_idx = 0

def handle_un_or_al(option, value):
    match = re.match(r"((?:al|un)(?:-conc|-mates)?)(?:-(gz|bz2|lz4|zst))?", option)
    if match:
        name, compression_type = match.groups()
        read_fns[name] = value
        read_compress[name] = ""
        if compression_type == "gz":
            read_compress[name] = "gzip"
        elif compression_type == "bz2":
            read_compress[name] = "bzip2"
        elif compression_type == "lz4":
            read_compress[name] = "lz4"
        elif compression_type == "zst":
            read_compress[name] = "zstd"

unps = []
mate1s = []
mate2s = []
tab5_mates = []
tab6_mates = []
interleaved_mates = []
bam_files = []
to_delete = []
to_kill = []
temp_dir = tempfile.gettempdir()
bam_out = 0
ref_str = None
no_pipes = 0
keep = 0
verbose = 0
readpipe = None
log_fname = None
help_flag = 0

long_options = [
    "1=", "2=", "reads=", "U=", "tab5=", "tab6=", "interleaved=", "b=",
    "temp-directory=", "bam", "no-named-pipes", "ref-string=",
    "reference-string=", "keep", "verbose", "debug", "sanitized",
    "large-index", "no-unal", "un=", "un-gz=", "un-bz2=", "un-lz4=",
    "un-zst=", "al=", "al-gz=", "al-bz2=", "al-lz4=", "al-zst=",
    "un-conc=", "un-conc-gz=", "un-conc-bz2=", "un-conc-lz4=",
    "un-conc-zst=", "al-conc=", "al-conc-gz=", "al-conc-bz2=",
    "al-conc-lz4=", "al-conc-zst=", "un-mates=", "un-mates-gz=",
    "un-mates-bz2=", "un-mates-lz4=", "un-mates-zst=", "log-file=",
    "help"
]
try:
    opts, args = getopt.getopt(sys.argv[1:], "1:2:U:b:", long_options)
except getopt.GetoptError as err:
    print(str(err))
    sys.exit(2)

bt2_args = list(args)

for opt, arg in opts:
    if opt in ("-1"):
        mate1s.extend(arg.split(','))
    elif opt in ("-2"):
        mate2s.extend(arg.split(','))
    elif opt in ("-U", "--reads"):
        unps.extend(arg.split(','))
    elif opt == "--tab5":
        tab5_mates.extend(arg.split(','))
    elif opt == "--tab6":
        tab6_mates.extend(arg.split(','))
    elif opt == "--interleaved":
        interleaved_mates.extend(arg.split(','))
    elif opt in ("-b"):
        bam_files.extend(arg.split(','))
    elif opt == "--temp-directory":
        temp_dir = arg
    elif opt == "--bam":
        bam_out = 1
    elif opt == "--no-named-pipes":
        no_pipes = 1
    elif opt in ("--ref-string", "--reference-string"):
        ref_str = arg
    elif opt == "--keep":
        keep = 1
    elif opt == "--verbose":
        verbose = 1
    elif opt == "--debug":
        debug = 1
    elif opt == "--sanitized":
        sanitized = 1
    elif opt == "--large-index":
        large_idx = 1
    elif opt == "--no-unal":
        no_unal = 1
    elif opt == "--un":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-gz":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-bz2":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-lz4":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-zst":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-gz":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-bz2":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-lz4":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-zst":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-conc":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-conc-gz":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-conc-bz2":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-conc-lz4":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-conc-zst":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-conc":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-conc-gz":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-conc-bz2":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-conc-lz4":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--al-conc-zst":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-mates":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-mates-gz":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-mates-bz2":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-mates-lz4":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--un-mates-zst":
        handle_un_or_al(opt[2:], arg)
    elif opt == "--log-file":
        log_fname = arg
    elif opt == "--help":
        help_flag = 1
        # Placeholder for help message if needed

if verbose == 1:
    bt2_args.append("--verbose")

# If the user asked us to redirect some reads to files, or to suppress
# unaligned reads, then we need to capture the output from Bowtie 2 and pass it
# through this wrapper.
passthru = 0
if len(read_fns) > 0 or no_unal or bam_out:
    passthru = 1
    bt2_args.append("--passthrough")
    cap_out = "-"
    i = 0
    while i < len(bt2_args):
        arg = bt2_args[i]
        if arg in ("-S", "--output"):
            if i < len(bt2_args) - 1:
                cap_out = bt2_args[i + 1]
                bt2_args[i] = None
                bt2_args[i + 1] = None
                i += 1
            else:
                fail(f"{arg} takes an argument.\n")
        i += 1
    bt2_args = [arg for arg in bt2_args if arg is not None]

old_stderr = None
if log_fname:
    try:
        old_stderr_fileno = os.dup(sys.stderr.fileno())
        old_stderr = open(os.fdopen(old_stderr_fileno, 'w'))
        sys.stderr.close()
        sys.stderr = open(log_fname, "w")
    except OSError as e:
        fail(f"Cannot redirect to log file {log_fname}.\n{e}")

def which(exec_name):
    for path in os.environ["PATH"].split(os.pathsep):
        file_path = Path(path) / exec_name
        if file_path.is_file() and os.access(file_path, os.X_OK):
            return True
    return False

def cat_file(ifn, ofh):
    if ifn.endswith(".gz"):
        try:
            with gzip.open(ifn, "rt") as ifh:
                for line in ifh:
                    ofh.write(line)
        except gzip.BadGzipFile:
            fail(f"Could not open gzipped read file: {ifn}\n")
        except FileNotFoundError:
            fail(f"Could not open read file: {ifn}\n")
    elif ifn.endswith(".bz2"):
        try:
            process = subprocess.run(["bzip2", "-dc", ifn], capture_output=True, text=True, check=True)
            ofh.write(process.stdout)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            fail(f"Could not open bzip2ed read file: {ifn}\n{e}")
    elif ifn.endswith(".lz4"):
        try:
            process = subprocess.run(["lz4", "-dc", ifn], capture_output=True, text=True, check=True)
            ofh.write(process.stdout)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            fail(f"Could not open lz4ed read file: {ifn}\n{e}")
    elif ifn.endswith(".zst"):
        try:
            process = subprocess.run(["zstd", "-dfc", ifn], capture_output=True, text=True, check=True)
            ofh.write(process.stdout)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            fail(f"Could not open zstded read file: {ifn}\n{e}")
    else:
        try:
            with open(ifn, "r") as ifh:
                for line in ifh:
                    ofh.write(line)
        except FileNotFoundError:
            fail(f"Could not open read file: {ifn}\n")

def write_files(input_files, output_file, no_pipes):
    if not no_pipes:
        try:
            os.mkfifo(output_file, 0o700)
        except OSError as e:
            fail(f"mkfifo({output_file}) failed.\n{e}")

    pid = 0
    if not no_pipes:
        pid = os.fork()
        if pid:
            to_kill.append(pid)

    if pid == 0:
        try:
            with open(output_file, "w") as ofh:
                for fn in input_files:
                    cat_file(fn, ofh)
            if not no_pipes:
                os._exit(0)
        except OSError as e:
            fail(f"Can't open '{output_file}' for writing.\n{e}")
            if not no_pipes:
                os._exit(1)

def maybe_wrap_input(orig_m1s, orig_m2s, ext):
    m1s = []