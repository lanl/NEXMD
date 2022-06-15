import numpy as np
import sys
import os
import shutil
import copy
import math

import re
import glob
import subprocess
import shlex
import filecmp

from sys import path as syspath


path = os.getcwd()
print(path)

# get fds

fds = glob.glob("*")
fds.remove("convert.py")

for f in fds:
    if f[-2:] != 'py':
        fds.remove(f)

print(fds)

for f in fds:
    oldf = open(f, 'r').readlines()
    newf = open(f[:-3]+'3.py','w')

    for line in oldf:
        if 'print >>' in line:
            fname = line.split()[2]
            k = line.index(fname) + len(fname)
            k2 = line.index('print')
            newline = line[:k2]+fname[:-1]+'.wirte(' + line[k:-1] + ')\n'
            newf.write(newline)
        elif ' print ' in line:
            k = line.index('print') + 6
            newline = line[:k-6]+ 'print(' + line[k:-1] + ')\n'
            newf.write(newline)
        else:
            newf.write(line)
    newf.close()
    shutil.copy(f[:-3]+'3.py', f)
