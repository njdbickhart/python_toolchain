#!/usr/bin/env python3

import argparse
import os
import sys
import jinja2
import markdown
import re

sepline = re.compile(r'-+')

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
        div.container {
            border-style: outset;
        }
        div.sticky {
            position: -webkit-sticky;
            position: sticky;
            top: 0;
            padding: 5px;
            margin-right:auto;margin-left:auto;*zoom:1;
            background-color: #C4C4C4;
            border: solid;
            border-width: 2px 6px 12px 18px;
        }
        table {
            width: 100%;
            border: 1px solid #000000;
            border-collapse: collapse;
        }
        table td, th {
            border: 1px solid #000000;
        }
        table td:first-child {
            font-weight: bold;
        }
        table tbody td {
            font-size: 13 px;
        }
        table tr:nth-child(even){
            background: #D0E4F5;
        }
        table thead {
            background: #0B6FA4;
            border-bottom: 5 px solid #053047;
        }
        table thead th {
            font-size: 17px;
            font-weight: bold;
            color: #FFFFFF;
            text-align: center;
            border-left: 2px solid #053047;
        }
        img[src*="#half"]{
            display: block;
            width: 65%;
            height: 65%;
            margin-left: auto;
            margin-right: auto;
        }
        img[src*="#regular"]{
            display: block;
            margin-left: auto;
            margin-right: auto;
        }
    </style>
</head>
<body>
<div class="sticky">
<h1> {{headtitle}} </h1>
</div>
<div class="container">
{{content}}
</div>
</body>
</html>
"""

TABLE = """<!DOCTYPE html>
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
        div.sticky {
            position: -webkit-sticky;
            position: sticky;
            top: 0;
            padding: 5px;
            background-color: #C4C4C4;
            border: solid;
            border-width: 2px 6px 12px 18px;
        }
        table {
            width: 100%;
            border: 1px solid #000000;
            border-collapse: collapse;
        }
        table td, th {
            border: 1px solid #000000;
        }
        table td:first-child {
            font-weight: bold;
        }
        table tbody td {
            font-size: 13 px;
        }
        table tr:nth-child(even){
            background: #D0E4F5;
        }
        table thead {
            background: #0B6FA4;
            border-bottom: 5 px solid #053047;
        }
        table thead th {
            font-size: 17px;
            font-weight: bold;
            color: #FFFFFF;
            text-align: center;
            border-left: 2px solid #053047;
        }
        img[src*="#half"]{
            display: block;
            width: 65%;
            height: 65%;
            margin-left: auto;
            margin-right: auto;
        }
        img[src*="#regular"]{
            display: block;
            margin-left: auto;
            margin-right: auto;
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
            description = "Create a summary webpage for the pipeline"
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
    createHtml(''.join(mdlines), "Assembly Quality Control", args.output + ".html")

    # Now create sub htmls for comparison plots
    for c in args.combos:
        asms = c.split('_')
        mdlines = compMd(c, asms, args.output, args.final)
        createHtml(''.join(mdlines), f'Assembly Comparison: {asms[0]} vs {asms[1]}', f'{args.final}/sum{c}.html')

    # We should be done!
    print("Finished!")

def compMd(combo, assemblies, outbase, finalfolder):
    fname = f'{finalfolder}/sum{combo}.html'

    mdlines = list()
    mdlines.extend([
                    #f'## Assembly comparison: {assemblies[0]} vs {assemblies[1]}',
                    '<br>', f'In all cases, {assemblies[1]} is the query and {assemblies[0]} is the target or reference',
                    '---',
                    '## Assembly alignment dotplot',
                    'This is a comparative alignment of the assemblies colored by the average percent identity of each alignment. If the assemblies are the same, you would expect a straight diagonal line with near max percent identity alignments.',
                    f'![Comparison of minimap2 alignments of each assembly to the other]({combo}/plot{combo}.png)',
                    f'[Click here for an interactive version of this plot!]({combo}/int{combo}.html)',
                    '---',
                    '## All assembly alignment variants',
                    f'These are all of the discernable alignment variants identified in this assembly comparison, plotted on a log scale. Insertions and expansions indicate an increase in assembly size in {assemblies[1]} compared to {assemblies[0]}. Vice versa for deletions and contractions.',
                    #f'<embed src="{combo}/vars{combo}.log_all_sizes.pdf" type="application/pdf" width="100%" height="600px" />',
                    f'![All variant sizes]({combo}/vars{combo}.log_all_sizes.png)',
                    '---',
                    '## Subsets of assembly alignment variants',
                    f'These plots show distributions of variants by size. These variant sites are the same as in the larger plot above, but scaled for easier viewing. This can be of interest when identifying differences in repeat structure between assemblies.',
                    #f'<embed src="{combo}/vars{combo}.75-1000.pdf" type="application/pdf" width="100%" height="600px" />\n<embed src="{combo}/vars{combo}.1000-500000.pdf" type="application/pdf" width="100%" height="600px" />',
                    f'![Variants from 75 to 1000 bp]({combo}/vars{combo}.75-1000.png) ![Variants from 1000 to 500kb]({combo}/vars{combo}.1000-500000.png)',
                    '---',
                    '## Pairwise kmer spectra plots',
                    f'These plots show the same information as on the main page, but they are stacked side by side for comparison. {assemblies[1]} is the leftmost plot and {assemblies[0]} is the rightmost plot.',
                    f'![Kmer spectrum plot]({assemblies[1]}.spectra-asm.st.png) ![Kmer spectrum plot]({assemblies[0]}.spectra-asm.st.png)',
                    '---',
                    '## Pairwise ideogram error plots',
                    f'These plots show the same information as on the main page, but they are stacked side by side for comparison. {assemblies[1]} is the leftmost plot and {assemblies[0]} is the rightmost plot.',
                    f'![Ideogram error plots](ideogram_errors.{assemblies[1]}.png) ![Ideogram error plots](ideogram_errors.{assemblies[0]}.png)',
                    f'##[Return to previous summary page](../{outbase}.html)'])

    # Same hack to add new lines. No excuse this time
    for i, l in enumerate(mdlines):
        mdlines[i] = l + '\n'

    # Returns a list of all lines in the md file
    return mdlines

def indexMd(fastas, assemblies, combos, finalfolder):
    mdlines = list()
    mdlines.append('# Assembly Report')
    #mdlines.append('<p id="pex">')
    mdlines.append('---')
    mdlines.append('## Table of Contents')
    mdlines.extend(['* [Assembly quality comparison](#asmqual)',
                    '* [Assembly feature comparisons](#asmfeat)',
                    '* [Assembly kmer comparisons](#asmkmer)',
                    '* [Assembly error windows](#asmwin)',
                    '* [Assembly comparisons](#asmcomp)',
#                    '</p>',
#                    '<p id="pex">',
                    '---',
                    '<a name="asmqual"></a>',
                    '## Assembly quality comparisons',
                    '<br>',
                    'These are general statistics for estimating assembly error rate and quality. These tables and figures provide some information regarding an assembly\'s completeness, but may obscure smaller defects or large structural issues within the assembly. Notably, scaffold misjoins are a common error that can artificially inflate several of these statistics.',
                    '---',
                    '#### Assemblies',
                    'Here are the assemblies being compared, along with their short-hand labels'])

    for f, a in zip(fastas, assemblies):
        mdlines.append(f'>{a}:\t\t{f}<br><br>')

    mdlines.extend(['---', '#### NG(X) plot',
                    'This is a measure of assembly continuity. The dotted line is the 50% mark of the anticipated assembly length. The higher the assembly\'s line is at this point, the more continuous the assembly.<br>',
                    #f'<embed src="{finalfolder}/combined_ngx_plot.pdf" type="application/pdf" width="100%" height="600 px" />',
                    f'![NG(x) plot of all assemblies]({finalfolder}/combined_ngx_plot.png#regular)',
                    '---', '#### Busco score plots',
                    f'![BUSCO category plots]({finalfolder}/combined_buscos.png#regular)',
                    '---', '#### Assembly Quality Statistics',
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
                #l, _ = re.subn(r'^|', '', l)
                #l, _ = re.subn(r'|\n$', '\n', l)
                #l, _ = re.subn(':', '-', l)
                tablines[-1] += l

    mdlines.append(testHtml(tablines[0]))

    mdlines.extend(['---',
#                    '<p id="pex">',
                    '<a name="asmfeat"></a>',
                    '---', '## Alignment feature comparisons',
                    '<br>',
                    'These statistics represent smaller scale variants detected from the alignment of reads to the assembly.',
                    '---', '#### Feature Response Curves',
                    'The following plot shows sorted lengths of the assemblies with the fewest errors. A "better" assembly "peaks" further to the left and top of the plot. These metrics do not always correlate with assembly continuity, so an assembly with a higher N50 might not perform as well in these metrics if it has more errors.',
                    #f'<embed src="{finalfolder}/combined_frc_plot.pdf" type="application/pdf" width="100%" height="600px" />',
                    f'![Feature response curve]({finalfolder}/combined_frc_plot.png#regular)',
                    '---', '#### Feature Statistics',
                    'These are the errors plotted in the above Feature Response Curve image. All errors are defined in more detail by the [FRC_align](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052210) program. The fewer the number of errors detected, the better.'])

    mdlines.append(testHtml(tablines[1]))

    mdlines.extend(['---', '#### Structural Variant Statistics',
                    'These structural variants represent larger (> 500 bp) potential misassemblies in the assembly. While having more of these variants is a sign of relatively poor quality, there may be a higher than expected count of these variants if the comparison read dataset is from a different individual than the reference individual used in the assembly. Alternatively, high heterozygosity in the sequenced individual can also inflate these statistics.'])

    mdlines.append(testHtml(tablines[2]))

    # Assembly kmer comparison plots
    mdlines.extend(['---',
#                    '<p id="pex">',
                    '<a name="asmkmer"></a>',
                    '## Assembly kmer comparisons',
                    '<br>',
                    'These plots show differences in kmer composition between the input sequence read file and the assembly itself. Unique kmers to the assembly are on the far left, and Unique kmers in the reads are just right of that. If these are assemblies of diploid organisms, you should see two peaks of kmers corresponding to the heterozygous and homozygous regions of the genome, respectively. Significant deviations or presence of unique kmers in the read dataset can indicate problems in the assembly.'])

    for a in assemblies:
        mdlines.append('---')
        mdlines.append(f'#### {a} kmer spectra plot')
        #mdlines.append(f'<embed src="{finalfolder}/{a}.spectra-asm.st.pdf" type="application/pdf" width="100%" height="600px" />')
        mdlines.append(f'![Kmer spectrum plot]({finalfolder}/{a}.spectra-asm.st.png#half)')

    # Assembly error windows
    mdlines.extend(['---',
#                    '<p id="pex">',
                    '<a name="asmwin"></a>',
                    '## Assembly error windows',
                    '<br>',
                    'These are plots of a sliding window analysis to identify regions of each assembly that have higher than normal (upper quartile) counts of errors. Only the top contigs are plotted due to space constraints. In all cases, yellow represents the upper quartile (> 25\% of all values) whereas red represents the upper 5\% of all windows'])

    for a in assemblies:
        mdlines.append('---')
        mdlines.append(f'#### {a} ASM feature density on largest contigs')
        #mdlines.append(f'<embed src="{finalfolder}/ideogram_errors.{a}.pdf" type="application/pdf" width="100%" height="600px" />')
        mdlines.append(f'![Ideogram error plots]({finalfolder}/ideogram_errors.{a}.png#regular)')

    mdlines.extend(['---',
#                    '<p id="pex">',
                    '<a name="asmcomp"></a>',
                    '## Assembly comparisons',
                    'The following links are to pairwise comparisons of each assembly to the other. These can be informative when comparing one assembly to a reference genome of the same organism. There may also be some value in comparing assemblies between different species or breeds.'])

    # Pair wise combination links
    for c in combos:
        mdlines.append(f'#### [{c} comparison]({finalfolder}/sum{c}.html)<br><br>')

#    mdlines.append("</p>")
    # Just because I'm lazy and realized this later, I'm going to loop through and add new lines to each md string
    for i, l in enumerate(mdlines):
        mdlines[i] = l + '\n'

    # Returns a list of all lines in the md file
    return mdlines

def createHtml(mdlines, title, outfile):
    extensions = ['extra', 'smarty', 'tables']
    html = markdown.markdown(mdlines, extensions=extensions, output_format='html5')
    doc = jinja2.Template(TEMPLATE).render(headtitle=title, content=html)
    with open(outfile, 'w') as out:
        out.write(doc)

def testHtml(mdlines):
    extensions = ['extra', 'smarty', 'tables']
    html = markdown.markdown(mdlines, extensions=extensions, output_format='html5')
    doc = jinja2.Template(TABLE).render(content=html)
    return(doc)

if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
