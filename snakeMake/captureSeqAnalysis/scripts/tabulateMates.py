# This script is designed to count and consolidate mate maps from the vector sequence
# This is a one-shot script designed to work in a snakemake workflow
from collections import defaultdict
import re
import os
import sys

if len(sys.argv) != 4
