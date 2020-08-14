import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib.collections import BrokenBarHCollection
from itertools import cycle
from collections import defaultdict
import pandas
import numpy as np


asms = snakemake.wildcards["asm"]
print(asms)
