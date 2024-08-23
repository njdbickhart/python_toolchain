import os
import sys

usage = f'python3 {sys.argv[0]} <input chrconverter tab file> <input vcf> <output vcf> <sample name>'

if len(sys.argv) != 5:
    print(usage)
    sys.exit(-1)

convertF = sys.argv[1]
inputVCF = sys.argv[2]
outputVCF = sys.argv[3]

def parse_vcf(vcf, converter):
    """parse and analyze a vcf

    Go through a vcf line-by-line, gleaning necessary info from the
    header and then yielding a concordance score for each non-header
    line.
    """
    for line in vcf:
        if line.startswith("#"):
            pass
        else:
            segs = line.rstrip().split()
            chr = converter.get(segs[0], "None")
            lref = len(segs[3])
            lalt = len(segs[4])
            if chr == "None":
                pass # Skip missing chromosomes
            elif segs[6] != "PASS":
                pass # Skip variants that do not pass filters
            elif lref < 50 and lalt < 50:
                pass # Skip non-SV variants
            else:
                yield convert_to_sv(chr, lref, lalt, segs)

def convert_to_sv(chr, lref, lalt, segs):
    svtype = "None"
    if lref > 1:
        svtype = "<DEL>"
    elif lalt > 1:
        svtype = "<INS>"

    length = max(lref, lalt)

    return (chr, segs[1], ".", segs[3][0], svtype, segs[5], segs[6],
    f'SVLEN={length};{segs[7]}', segs[8], segs[9])

# Read in converter file and create dict
converter = dict()

with open(sys.argv[1], 'r') as input:
    for l in input:
        s = l.rstrip().split()
        converter[s[0]] = s[1]

# Parse vcf file line by line to create a headerless vcf
vcfh = open(sys.argv[2], 'r')
with open(sys.argv[3], 'w') as output:
    output.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sys.argv[4]}\n')
    for rows in parse_vcf(vcfh, converter):
        if rows[0] == "None" or rows[4] == "None" or rows[6] != "PASS":
            continue
        output.write("\t".join(rows) + "\n")


vcfh.close()
