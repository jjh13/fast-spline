"""
This tool fills a very specific need for this project.

It takes a single PTX file and a dummy c file that

"""

import argparse
import shlex
import subprocess
import tempfile
import shutil
import sys
import os
import time

if len(sys.argv) < 4:
    print(f"Tool usage: "
          f"python asptx.py [input ptx] [input dummy] [output o] [optional nvcc arguments]")

input_ptx = sys.argv[1]
input_cc = sys.argv[2]
output_o = sys.argv[3]
nvcc_args = sys.argv[4:]

possible_commands = [
    'gcc', 'cicc','cicc','ptxas', 'fatbinary', 'cudafe++'
]

if not os.path.isfile(input_ptx):
    raise ValueError(f"PTX file not found {input_ptx}")

if not os.path.isfile(input_cc):
    raise ValueError(f"Dummoy file not found: {input_cc}")

from contextlib import contextmanager

@contextmanager
def cwd(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)

with tempfile.TemporaryDirectory() as tmpdirname:
    print('created temporary directory', tmpdirname)

    new_ptx = os.path.abspath(input_ptx)
    new_cc = os.path.join(tmpdirname, os.path.split(input_cc)[-1])
    new_o = os.path.abspath(output_o) #os.path.join(tmpdirname, os.path.split(output_o)[-1])
    shutil.copy(input_ptx, tmpdirname)
    shutil.copy(input_cc, tmpdirname)
    with cwd(tmpdirname):
        r = subprocess.run([
            'nvcc',
            new_cc,
            '-dc',
            '-x', 'cu',
            '-o', new_o,
            '--dryrun',
            '--keep'
        ] + nvcc_args, capture_output=True)
        commands = r.stderr.decode("utf-8").split('\n')

        for cmd in commands:
            if len(cmd) == 0:
                continue
            command_args = shlex.split(cmd[3:], posix=False)

            if command_args[0] not in possible_commands:
                ps = cmd[3:].split('=')
                variable = ps[0]
                value = '='.join(ps[1:])
                os.environ[variable] = value
            else:
                if command_args[0] == 'ptxas':
                    print("$# Injecting new ptx")
                    for idx, _ in enumerate(command_args[1:]):
                        if _[-4:] == 'ptx"':
                            os.unlink(_[1:-1])
                            shutil.copy(new_ptx, tmpdirname)
                            os.rename(os.path.split(new_ptx)[-1],
                                      _[1:-1])
                            print(f"$# Overriding {_} with \"{new_ptx}\"")
                    cmd = " ".join(command_args)
                    print(cmd)
                    os.system(cmd)
                else:
                    print(cmd[3:])
                    os.system(cmd[3:])
    # time.sleep(1000)