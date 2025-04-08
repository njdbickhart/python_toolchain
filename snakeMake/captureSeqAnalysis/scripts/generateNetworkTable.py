# This script is designed to take consolidated mate tabulation data with a threshold copy number value
# It produces a table that is designed to be parsed by the python networkx library for visualization

from collections import defaultdict
import re
import os
import sys

usage = f'python {sys.argv[0]} <input mate table> <sample copy number file> <output table file>'

if len(sys.argv) != 4:
    print(usage)
    sys.exit(-1)

def main(args):
    print("hey")



if __name__ == "__main__":
    main(sys.argv)