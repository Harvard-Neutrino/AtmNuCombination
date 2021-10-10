import subprocess
import shlex
import sys

x = sys.argv[1]
y = sys.argv[2]
z = './runf.sh'
space = ' '
lalala = z + space + x + space + y

subprocess.call(shlex.split(lalala))