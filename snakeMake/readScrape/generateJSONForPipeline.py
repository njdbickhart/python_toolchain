#!/usr/bin/env python3
import sys
import glob
import re
from os.path import basename

usage = "Usage: python3 " + sys.argv[0] + " \<space delimited list of directories with BAM files\>"

fsep = re.compile('[\._]')

if length(sys.argv) == 1:
    print(usage)
    sys.exit()

print(f'\{\n\t\"repeatmasker\": \"\/path\/to\/\",\n\t\"samples\":\{')
for i in range(1, length(sys.argv)):
    dir=sys.argv[i]
    for filename in glob.iglob(dir, recursive=True):
        bname = basename(filename)
        bsegs = re.split(fsep, bname)
        print(f'\t\t\"{bsegs[0]}\" : \"filename\",')

print(f'\t\}\n\}')
