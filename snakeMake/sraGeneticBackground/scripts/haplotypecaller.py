__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"
# I had to modify this script to comply with cluster-specific limitations on the temp file

import os
import sys
import logging
import argparse
import subprocess as sp
from multiprocessing import Pool
from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

usage = f'python3 {sys.argv[0]} <input bams> <input ref> <extra> <java_opts> <outputbase> <output> <threads>'

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Generate summary information and PDF plots of Fastq Stats"
            )
    parser.add_argument('-b', '--bam',
                        help="An indexed bam file comma separated",
                        type=str, required=True
                        )
    parser.add_argument('-F', '--fasta',
                        help="The reference fasta file",
                        type=str, required=True
                        )
    parser.add_argument('-e', '--extra',
                        help="Temp directory",
                        type=str, required=True
                        )
    parser.add_argument('-j', '--java',
                        help="java options",
                        action='append', default=[]
                        )
    parser.add_argument('-o', '--outbase',
                        help="output base name and path",
                        type=str, required=True
                        )
    parser.add_argument('-v', '--vcf',
                        help="output vcf file name",
                        type=str, required=True
                        )
    parser.add_argument('-t', '--threads',
                        help="threads",
                        type=int, required=True
                        )

    return parser.parse_args()

args = parse_user_input()

reference = args.fasta
extra = f'--tmp-dir {args.extra}'
java_opts = " ".join([f'-{x}' for x in args.java])

bams = args.bam
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
vcf_output = args.vcf if args.vcf.endswith('.vcf') else ""
if vcf_output:
    output = " --output " + str(vcf_output)

gvcf_output = args.vcf if args.vcf.endswith('.gvcf') else ""
if gvcf_output:
    output = " --emit-ref-confidence GVCF " + " --output " + str(gvcf_output)

if (vcf_output and gvcf_output) or (not gvcf_output and not vcf_output):
    if vcf_output and gvcf_output:
        raise ValueError(
            "please set vcf or gvcf as output, not both! It's not supported by gatk"
        )
    else:
        raise ValueError("please set one of vcf or gvcf as output (not both)!")

bam_output = args.vcf if args.vcf.endswith('.bam') else ""
if bam_output:
    bam_output = " --bam-output " + str(bam_output)

threads = int(args.threads) - 4

#shell(
#    "gatk --java-options '{java_opts}' HaplotypeCaller"
#    " --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP"
#    " --native-pair-hmm-threads {threads}"
#    " {bams}"
#    " --reference {reference}"
#    " --verbosity DEBUG"
#    " {extra}"
#    " {output}"
#    " {bam_output}"
#)
cmd = ['gatk', '--java-options', f'\'{java_opts}\'', 'HaplotypeCaller',
    '--pair-hmm-implementation', 'AVX_LOGLESS_CACHING_OMP',
    '--native-pair-hmm-threads', f'{threads}',
    " ".join(bams),
    '--reference', f'{reference}',
    f'{extra}',
    f'{output}',
    f'{bam_output}'
    ]
print(cmd)
cmdstr = " ".join(cmd)
print(cmdstr)
proc = sp.Popen(" ".join(cmd), shell=True, stdout=sys.stdout, stderr=sys.stderr)

retcode = proc.wait()
print(f'retcode: {retcode}')
