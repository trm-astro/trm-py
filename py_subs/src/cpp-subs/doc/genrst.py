#!/usr/bin/env python

"""
This generates the rst files for sphinx to use.
Takes any lines between lines containing !!sphinx 
from 'source' and writes them out to 'dest'
""" 

import os, sys

if len(sys.argv) != 3:
    print('usage: genrst.py source dest')
    exit(1)

source, dest = sys.argv[1:]

if not os.path.exists(source):
    print('Source =', source,'does not exist')
    exit(1)

if not dest.endswith('.rst'):
    print('Destination =', dest,'does not end with .rst')
    exit(1)

out = False
with open(source) as fin, open(dest, 'w') as fout:
    for line in fin:
        if line.find('!!sphinx') > -1:
            out = not out
        elif out:
            fout.write(line)

print(source, '--->', dest)
