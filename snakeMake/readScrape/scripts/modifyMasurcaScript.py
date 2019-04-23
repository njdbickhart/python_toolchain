#!/usr/bin/env python3
import re

# Generate list of conversion files
files = ['gapClose.err', 'runCA2.out', 'tigStore.err', 'unitig_cov.txt', 'global_arrival_rate.txt', 'unitig_layout.txt', 'genome.uid', 'runCA1.out', \
'runCA0.out', 'superReadSequences_shr.frg', 'pe.linking.frg', 'pe.linking.fa', 'work1', 'super1.err', 'sj.cor.clean.frg', 'sj.cor.ext.reduced.fa', \
'mates_to_break.txt', 'compute_jump_coverage.txt', 'sj.cor.clean.fa', 'redundant_sj.txt', 'chimeric_sj.txt', 'super2.err', 'guillaumeKUnitigAtLeast', \
'pe.cor.fa', 'error_correct.log', 'pe_data.tmp', 'sj.renamed.fastq', 'meanAndStdevByPrefix', 'pe.renamed.fastq']

reCompiled  = dict((v, re.compile(v) for v in files)

# Now go line by line through the file to replace the files with a prefix directory for them
replace = snakemake.input[0]
replace = re.sub('/$', '', replace)

rcount = dict()
with open(snakemake.input[1], 'r') as in, open(snakemake.output[0], 'w') as out:
    for l in in:
        l = l.rstrip()
        for k, v in reCompiled:
            newS, ct = re.subn(v, replace + '/' + k, l)
            rcount[k] += ct
            l = newS
        # Special cases that are difficult to code otherwise
        l = re.sub(r'\s+CA', replace + '/CA', l)
        l = re.sub(r'"CA', '"' + replace + '/CA', l)
        out.write(l + "\n")
    print('#Replacement stats')
    for k, v in rcount:
        out.write(f'#{k}\t{v}\n')
