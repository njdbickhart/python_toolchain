#!/usr/bin/env python3

import argparse
import os
import sys
import jinja2
import markdown
import re

sepline = re.compile(r'^-+\n$')

TEMPLATE = """<!DOCTYPE html>
<html>
<head>
    <link href="http://netdna.bootstrapcdn.com/twitter-bootstrap/2.3.0/css/bootstrap-combined.min.css" rel="stylesheet">
    <style>
        body {
            font-family: sans-serif;
        }
        code, pre {
            font-family: monospace;
        }
        h1 code,
        h2 code,
        h3 code,
        h4 code,
        h5 code,
        h6 code {
            font-size: inherit;
        }
    </style>
</head>
<body>
<div class="container">
{{content}}
</div>
</body>
</html>
"""


def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Compare multiple genome assemblies with short-read alignment data" + version
            )
    parser.add_argument('-f', '--final',
                        help="The location of the final folder of the pipeline",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output',
                        help="output file basename ",
                        type=str, required=True
                        )
    parser.add_argument('-c', '--combos',
                        help="Expected assembly combinations",
                        action="append", default=[]
                        )
    parser.add_argument('-s', '--fasta',
                        help="Assembly file name",
                        action="append", default=[]
                        )
    parser.add_argument('-a', '--assembly',
                        help="Assembly label name",
                        action="append", default=[]
                        )

    return parser.parse_args(), parser

def main(args, parser):
    if len(args.fasta) < 2 or len(args.assembly) < 2:
        print("Error! Must input more than one assembly file!")
        parser.print_help()
        sys.exit()

    # Create the main index html
    mdlines = indexMd(args.fasta, args.assembly, args.combos, args.final)
    createHtml(mdlines, args.output + ".html")

    # Now create sub htmls for comparison plots
    for c in args.combos:
        asms = c.split('_')
        mdlines = compMd(c, asms, args.outbase, args.final)
        createHtml(mdlines, f'{args.final}/sum{combo}.html')

    # We should be done!
    print("Finished!")

def compMd(combo, assemblies, outbase, finalfolder):
    fname = f'{finalfolder}/sum{combo}.html'

    mdlines = list()
    mdlines.extend(f'## Assembly comparison: {assemblies[0]} vs {assemblies[1]}',
                    f'In all cases, {assemblies[1]} is the query and {assemblies[0]} is the target or reference',
                    '## Assembly alignment dotplot',
                    'This is a comparative alignment of the assemblies colored by the average percent identity of each alignment. If the assemblies are the same, you would expect a straight diagonal line with near max percent identity alignments.',
                    f'![Comparison of minimap2 alignments of each assembly to the other]({combo}/plot{combo}.png)',
                    '## All assembly alignment variants',
                    f'These are all of the discernable alignment variants identified in this assembly comparison, plotted on a log scale. Insertions and expansions indicate an increase in assembly size in {assemblies[1]} compared to {assemblies[0]}. Vice versa for deletions and contractions.',
                    f'![Identified minimap2 alignment variants found within the assembly]({combo}/vars{combo}.log_all_sizes.pdf)',
                    '## Subsets of assembly alignment variants',
                    f'These plots show distributions of variants by size. These variant sites are the same as in the larger plot above, but scaled for easier viewing. This can be of interest when identifying differences in repeat structure between assemblies.',
                    f'![Minimap2 alignment variants between 75 and 1000 bp]({combo}/vars{combo}.75-1000.pdf) ![Minimap2 alignment variants between 1000 and 500 kb]({combo}/vars{combo}.1000-500000.pdf)',
                    f'[Return to previous summary page](../{outbase}.html)'

    # Same hack to add new lines. No excuse this time
    for i, l in enumerate(mdlines):
        mdlines[i] = l + '\n'

    # Returns a list of all lines in the md file
    return mdlines

def indexMd(fastas, assemblies, combos, finalfolder):
    mdlines = list()
    mdlines.append('# Assembly Report')
    mdlines.append('---')
    mdlines.append('## Table of Contents')
    mdlines.extend(['* [Assembly quality comparison](#asmqual)',
                    '* [Assembly feature comparisons](#asmfeat)',
                    '* [Assembly error windows](#asmwin)',
                    '* [Assembly comparisons](#asmcomp)',
                    '<a name="asmqual"></a>',
                    '## Assembly quality comparisons',
                    'These are general statistics for estimating assembly error rate and quality. These tables and figures provide some information regarding an assembly\'s completeness, but may obscure smaller defects or large structural issues within the assembly. Notably, scaffold misjoins are a common error that can artificially inflate several of these statistics.',
                    '#### Assemblies',
                    'Here are the assemblies being compared, along with their short-hand labels'])

    for f, a in zip(fastas, assemblies):
        mdlines.append(f'>{a}:{f}')

    mdlines.extend(['#### NG(X) plot',
                    'This is a measure of assembly continuity. The dotted line is the 50% mark of the anticipated assembly length. The higher the assembly\'s line is at this point, the more continuous the assembly',
                    f'![Comparison of assembly continuity. Higher to the left is better]({finalfolder}/combined_ngx_plot.pdf)',
                    '#### Busco score plots',
                    '> To be added later',
                    '#### Assembly Quality Statistics',
                    'These are statistics derived from the overal continuity of the assembly and the alignment of reads/kmers to it. Better assemblies tend to have fewer contigs, higher QV values, lower error rates, and higher BUSCO scores.'])

    # Process the summary table for later insertion in the md file lines
    # The first string is the asm quality Stats
    # The second is the feature statistics
    # The last are the structural variants
    tablines = list()
    tablines.append('')
    with open(finalfolder + '/summary_table.tab', 'r') as input:
        for l in input:
            if re.match(sepline, l):
                tablines.append('')
            else:
                tablines[-1] += l

    mdlines.append(tablines[0])

    mdlines.extend(['---',
                    '<a name="asmfeat"></a>',
                    '## Alignment feature comparisons',
                    'These statistics represent smaller scale variants detected from the alignment of reads to the assembly.',
                    '#### Feature Response Curves',
                    'The following plot shows sorted lengths of the assemblies with the fewest errors. A "better" assembly "peaks" further to the left and top of the plot. These metrics do not always correlate with assembly continuity, so an assembly with a higher N50 might not perform as well in these metrics if it has more errors.',
                    f'![Feature Response Plots. Higher to the left is better]({finalfolder}/combined_frc_plot.pdf)',
                    '#### Feature Statistics',
                    'These are the errors plotted in the above Feature Response Curve image. All errors are defined in more detail by the [FRC_align](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052210) program. The fewer the number of errors detected, the better.'])

    mdlines.append(tablines[1])

    mdlines.extend(['#### Structural Variant Statistics',
                    'These structural variants represent larger (> 500 bp) potential misassemblies in the assembly. While having more of these variants is a sign of relatively poor quality, there may be a higher than expected count of these variants if the comparison read dataset is from a different individual than the reference individual used in the assembly. Alternatively, high heterozygosity in the sequenced individual can also inflate these statistics.'])

    mdlines.append(tablines[2])

    # Assembly error windows
    mdlines.extend(['---',
                    '<a name="asmwin"></a>',
                    '## Assembly error windows',
                    'These are plots of a sliding window analysis to identify regions of each assembly that have higher than normal (upper quartile) counts of errors. Only the top contigs are plotted due to space constraints. In all cases, yellow represents the upper quartile (> 25\% of all values) whereas red represents the upper 5\% of all windows'])

    for a in assemblies:
        mdlines.append(f'#### {a} ASM feature density on largest contigs')
        mdlines.append(f'![Ideogram plot of {a} assembly. Error windows are on the bottom of each Ideogram bar.]({finalfolder}/ideogram_errors.{a}.pdf)')

    mdlines.extend(['---',
                    '<a name="asmcomp"></a>',
                    '## Assembly comparisons',
                    'The following links are to pairwise comparisons of each assembly to the other. These can be informative when comparing one assembly to a reference genome of the same organism. There may also be some value in comparing assemblies between different species or breeds.')]

    # Pair wise combination links
    for c in combos:
        mdlines.append(f'#### [{c} comparison]({finalfolder}/sum{c}.html)')

    # Just because I'm lazy and realized this later, I'm going to loop through and add new lines to each md string
    for i, l in enumerate(mdlines):
        mdlines[i] = l + '\n'

    # Returns a list of all lines in the md file
    return mdlines

def createHtml(mdlines, outfile):
    extensions = ['extra', 'smarty']
    html = markdown.markdown(mdlines, extensions=extensions, output_format='html5')
    doc = jinja2.Template(TEMPLATE).render(content=html)
    with open(outfile, 'w') as out:
        out.write(doc)

if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
