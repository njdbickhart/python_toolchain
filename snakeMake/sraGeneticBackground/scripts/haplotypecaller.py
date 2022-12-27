__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"
# I had to modify this script to comply with cluster-specific limitations on the temp file

import os
import sys
import logging
from multiprocessing import Pool
from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

usage = f'python3 {sys.argv[0]} <input bams> <input ref> <extra> <java_opts> <outputbase> <output> <threads>'

reference = sys.argv[2]
extra = sys.argv[3]
java_opts = sys.argv[4]

bams = sys.argv[1].split(',')
if isinstance(bams, str):
    bams = [bams]
bams = list(map("--input {}".format, bams))

#intervals = snakemake.input.get("intervals", "")
#if not intervals:
#    intervals = snakemake.params.get("intervals", "")
#if intervals:
#    intervals = "--intervals {}".format(intervals)

#known = snakemake.input.get("known", "")
#if known:
#    known = "--dbsnp " + str(known)
output = ""
vcf_output = sys.argv[6] if sys.argv[6].endswith('.vcf') else ""
if vcf_output:
    output = " --output " + str(vcf_output)

gvcf_output = sys.argv[6] if sys.argv[6].endswith('.gvcf') else ""
if gvcf_output:
    output = " --emit-ref-confidence GVCF " + " --output " + str(gvcf_output)

if (vcf_output and gvcf_output) or (not gvcf_output and not vcf_output):
    if vcf_output and gvcf_output:
        raise ValueError(
            "please set vcf or gvcf as output, not both! It's not supported by gatk"
        )
    else:
        raise ValueError("please set one of vcf or gvcf as output (not both)!")

bam_output = sys.argv[6] if sys.argv[6].endswith('.bam') else ""
if bam_output:
    bam_output = " --bam-output " + str(bam_output)

threads = int(sys.argv[7]) - 4

shell(
    "gatk --java-options '{java_opts}' HaplotypeCaller"
    " --native-pair-hmm-threads {threads}"
    " {bams}"
    " --reference {reference}"
    " {extra}"
    " {output}"
    " {bam_output}"
)
