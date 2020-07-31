from collections import defaultdict
import pysam
import re

entries = ["merQV", "merErrorRate", "merCompleteness", "baseQV",
"unmap%", "LOW_COV_PE", "LOW_NORM_COV_PE", "HIGH_SPAN_PE", "HIGH_COV_PE",
"HIGH_NORM_COV_PE", "STRECH_PE", "COMPR_PE", "HIGH_OUTIE_PE",
"HIGH_SINGLE_PE", "SVDEL", "SVDUP", "SVBND"]

descriptions = {"merQV" : "kmer-based Quality", "merErrorRate" : "kmer-based error rate",
"merCompleteness" : "Proportion of complete assembly based on kmers", "baseQV" : "SNP and INDEL Quality value",
"unmap%" : "Percentage of short-reads unmapped", "LOW_COV_PE" : "Low read COV areas",
"LOW_NORM_COV_PE" : "Low COV of normal PE reads", "HIGH_SPAN_PE" : "Regions with high numbers of inter-contig PE reads",
"HIGH_COV_PE" : "Regions with High read coverage", "HIGH_NORM_COV_PE" : "Regions with high coverage of normal PE reads",
"STRECH_PE" : "Regions with high Comp/Expansion (CE) statistics", "COMPR_PE" : "Regions with low Comp/Expansion (CE) statistics",
"HIGH_OUTIE_PE" : "Regions with high counts of improperly paired reads",
"HIGH_SINGLE_PE" : "Regions with high counts of single unmapped reads", "SVDEL" : "Number of deletion SVs",
"SVDUP" : "Number of Duplication SVs", "SVBND" : "Number of Complex SVs"}

solid = dict()
data = defaultdict(int)

# Populate merqury entries
with open(snakemake.input["merqv"], 'r') as qv:
    l = qv.readline()
    s = l.rstrip().split()
    solid["merQV"] = s[3]
    solid["merErrorRate"] = s[4]


with open(snakemake.input["complete"], 'r') as comp:
    l = comp.readline()
    s = l.rstrip().split()
    solid["merCompleteness"] = s[4]

print("loaded merqury stats")

# Populate QV and mapped reads entries
with open(snakemake.input["snpqv"], 'r') as qv:
    l = qv.readline()
    solid["baseQV"] = l.rstrip()

text = pysam.idxstats(snakemake.input["bams"])
lines = text.split(sep="\n")
mapped = 0
unmapped = 0
for i in lines:
    segs = i.split()
    if len(segs) < 4:
        continue
    mapped += int(segs[2])
    unmapped += int(segs[3])
solid["unmap%"] = "{:.2f}".format((unmapped / (mapped + unmapped)) * 100)

print("loaded QV and mapping stats")

# Populate FRC entries
with open(snakemake.input["features"], 'r') as frc:
    for l in frc:
        s = l.rstrip().split()
        data[s[1]] += 1

print("loaded FRC entries")

# Populate lumpy entries
with open(snakemake.input["lumpy"], 'r') as lump:
    for l in lump:
        if l.startswith('#'):
            continue
        s = l.rstrip().split()
        d = s[7].split(';')
        t = d[0].replace("SVTYPE=", "SV")
        data[t] += 1

print("loaded lumpy entries")

# Write out Table
ecol = 5
ccol = 5
dcol = 60
# update e and c col widths
for d in [solid, data]:
    for k, v in d.items():
        ecol = len(str(k)) if len(str(k)) > ecol else ecol
        ccol = len(str(v)) if len(str(v)) > ccol else ccol

ecol += 1
ccol += 1

# format d col string
elines = dict()
for k, v in descriptions.items():
    (s, nsubs) = re.subn(r'(.{59})', r'\1\n', v)
    descriptions[k] = s
    elines[k] = nsubs

esep = '-' * (ecol - 1)
csep = '-' * (ccol - 1)
dsep = '-' * (dcol - 1)

print(f'fixed width entry sizes: {ecol} {ccol} {dcol}')

with open(snakemake.output["table"], 'w') as out:
    # First QV Scores
    out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format("Q Scores", "Value", "Description", ecol= ecol, ccol = ccol, dcol=dcol))
    out.write('|:{}|{}:|:{}|\n'.format(esep, csep, dsep))

    for i in ["merQV", "merErrorRate", "merCompleteness", "baseQV", "unmap%"]:
        nsubs = elines[i]
        d = descriptions[i].split('\n')
        out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format(i, str(solid[i]), d[0], ecol= ecol, ccol = ccol, dcol=dcol))
        # Writing out what's left of the Description Line
        for j in range(nsubs):
            out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format("", "", d[j+1], ecol= ecol, ccol = ccol, dcol=dcol))

    print("Done with QV")
    # Next FRC
    out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format("Features", "Value", "Description", ecol= ecol, ccol = ccol, dcol=dcol))
    out.write('|:{}|{}:|:{}|\n'.format(esep, csep, dsep))

    for i in ["LOW_COV_PE", "LOW_NORM_COV_PE", "HIGH_SPAN_PE", "HIGH_COV_PE",
    "HIGH_NORM_COV_PE", "HIGH_OUTIE_PE",
    "HIGH_SINGLE_PE", "STRECH_PE", "COMPR_PE"]:
        nsubs = elines[i]
        d = descriptions[i].split('\n')
        out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format(i, str(data[i]), d[0], ecol= ecol, ccol = ccol, dcol=dcol))
        # Writing out what's left of the Description Line
        for j in range(nsubs):
            out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format("", "", d[j+1], ecol= ecol, ccol = ccol, dcol=dcol))

    print("Donw with FRC")
    # Finally SV calls
    out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format("SVs", "Value", "Description", ecol= ecol, ccol = ccol, dcol=dcol))
    out.write('|:{}|{}:|:{}|\n'.format(esep, csep, dsep))

    for i in ["SVDEL", "SVDUP", "SVBND"]:
        nsubs = elines[i]
        d = descriptions[i].split('\n')
        out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format(i, str(data[i]), d[0], ecol= ecol, ccol = ccol, dcol=dcol))
        # Writing out what's left of the Description Line
        for j in range(nsubs):
            out.write('|{0: <{ecol}}|{1: >{ccol}}|{2: <{dcol}}|\n'.format("", "", d[j+1], ecol= ecol, ccol = ccol, dcol=dcol))
