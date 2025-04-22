from bowtie2_py.bowtie2 import main
import platform

def main2():
    if platform.system() != "Windows":
        print("""The current OS was not supported.
However, feel free to fork or contribute to this project:
              https://github.com/soda92/bowtie2_py
""")
        exit(-1)
    main()
