import sys
import glob
import re
from os.path import basename

usage = "Usage: python3 " + sys.argv[0] + " <space delimited list of directories with BAM files>"

reqPaths = ['fasta']

fsep = re.compile('[\._]')

if len(sys.argv) == 1:
    print(usage)
    sys.exit()

print("{")
for i in reqPaths:
    print(f' \"{i}\": \"/path/to/\",')

print(f' \"samples\":{{')
for i in range(1, len(sys.argv)):
    dir=sys.argv[i]
    for filename in glob.iglob(dir + "/**/*.bam", recursive=True):
        bname = basename(filename)
        bsegs = re.split(fsep, bname)
        print(f'  \"{bsegs[0]}\" : \"{filename}\",')

print(f' }}\n')
