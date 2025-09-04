
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template
from weasyprint import HTML
import argparse
import glob
import re
from collections import defaultdict
from itertools import cycle

import plotly.offline as pyo
import marker_loc_manhattan

ucsc = re.compile(r'(.+):(\d+)-(\d+)')

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to generate a report file from capture-seq data"
            )
    parser.add_argument('-f', '--file', 
                        help="Input final table file",
                        required=True, type=str
                        )
    parser.add_argument('-i', '--index', 
                        help="Input fasta index file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output directory",
                        required=True, type=str,
                        )
    parser.add_argument('-p', '--priority',
                        help="Priority chromosome for sorting", 
                        required=True, type=str,
                        )
    parser.add_argument('-c', '--contaminant',
                        help='Contaminants to avoid (may be entered multiple times)',
                        action='append', default=[]
                        )

    return parser.parse_args(), parser


def main(args, parser):
    # Create plot folders in case they do not exist
    os.makedirs(f'{args.output}/plots', exist_ok=True)

    # Load the CSV file
    df = pd.read_csv(args.file)

    # Modify the dataframe to create statistics
    df['CHRCount'] = df['CHROMOSOME'].astype(str).apply(lambda x: len(x.split(';')))
    df['LocCount'] = df['SITES'].astype(str).apply(lambda x: len(x.split(';')))
    df['GENE_OVERLAPS'] = df['GENE_OVERLAPS'].fillna(value="None")
    df['CONTAMINANTS'] = df['CONTAMINANTS'].fillna('None')

    # Need to make another DF for the manhattan plot
    data = defaultdict(list)
    with open(args.file, 'r') as input:
        head = input.readline()
        for l in input:
            s = l.rstrip().split(",")
            locs = s[2].split(';')
            for j in locs:
                psegs = re.split(r'[:-]', j)
                if len(psegs) < 3:
                    continue
                
                data['CHR'].append(psegs[0])
                data['BP'].append(int(psegs[1]))
                data['SNP'].append(s[0])
                data['GENE'].append(s[6])
                data['LOCS'].append(len(locs))
    mDF = pd.DataFrame(data)

    # Create temp bed file for plotting
    with open(args.file, 'r') as input, open(f'{args.output}/temp.bed', 'w') as output:
        head = input.readline()
        for l in input:
            s = l.rstrip().split(",")
            for ls in s[2].split(';'):
                bsegs = re.split(r'[:-]', ls)
                output.write("\t".join(bsegs) + "\t" + s[0] + "\n")

    # Sort DF based on priority chromosome
    df = pd.concat([df[df['CHROMOSOME'].astype(str).str.contains(args.priority)].sort_values('LocCount'), 
                    df[~df['CHROMOSOME'].astype(str).str.contains(args.priority)].sort_values(['CHROMOSOME','LocCount'])])
    


    # Summary statistics
    summary_stats = df[df['CHROMOSOME'].astype(str).str.contains(args.priority)][['CHROMOSOME', 
                                                                                            'CONTAMINANTS']].groupby(by=['CHROMOSOME', 
                                                                                                                         'CONTAMINANTS']).value_counts().reset_index(name='COUNTS').to_html(classes='table table-striped', border=0)


    # Priority Individuals
    priority_table = df[(df['CHROMOSOME'].astype(str).isin([args.priority])) & \
                        (df['LocCount'] == 1) & \
                            (df['GENE_OVERLAPS'] == 'None') & \
                                (~df['CONTAMINANTS'].isin(args.contaminant))].to_html(classes='table table-striped', border=0)



    # Distribution plot of chromosome insertion sites
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=df, x="LocCount", y= 'CHRCount', hue="CONTAMINANTS", fill=True, alpha=.5)
    plt.title('Comparison of insertion sites vs different chromosome insertions')
    plt.xlabel('Insertion Location Counts')
    plt.ylabel('Chromosome Insertion Count')
    plt.tight_layout()
    plt.savefig(f"{args.output}/plots/loc_count_vs_chrs.png")
    plt.close()

    # Pie chart for CONTAMINANTS
    contaminant_counts = df['CONTAMINANTS'].value_counts()
    plt.figure(figsize=(8, 8))
    contaminant_counts.plot.pie(autopct='%1.1f%%', startangle=90)
    plt.title('Distribution of CONTAMINANTS')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(f"{args.output}/plots/contaminants_pie_chart.png")
    plt.close()

    # Chromosome plot 
    plotSimpleChrGraph(args.index, f'{args.output}/temp.bed', f'{args.output}/plots/chromosome_layout.png')

    # Interactive Chromosome plot
    mobj = marker_loc_manhattan.ManhattanPlot(mDF)
    fig = mobj.figure(title="Genome-wide Plot of Insertion Frequency per Sample (top values have only one insertion)")


    # Start building the HTML content
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <link href="http://netdna.bootstrapcdn.com/twitter-bootstrap/2.3.0/css/bootstrap-combined.min.css" rel="stylesheet">
        
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            .sample-section { margin-bottom: 10px; border: 1px solid #ccc; border-radius: 5px; }
            .sample-header { background-color: #f2f2f2; cursor: pointer; font-weight: bold; }
            .sample-highlight { background-color: #facfc6; cursor: pointer; font-weight: bold; }
            .sample-content { display: none; padding: 10px; }
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
        </style>
        <script>
            function toggleContent(id) {
                var content = document.getElementById(id);
                if (content.style.display === "none") {
                    content.style.display = "block";
                } else {
                    content.style.display = "none";
                }
            }
        </script>
    </head>
    <body>
        <div class="sticky">
        <h1>Capture-seq Summary Report</h1>
        </div>
        <h2>Priority Individuals for Breeding Program</h2>        
        {{ priority_table | safe}}        

        <h2>Count of Individuals with Priority Chromosome Insertions</h2>        
        {{ summary_stats | safe}}        
        
        <h2>Distribution of Chromosome insertions vs multiple insertion sites</h2>
        <center>
        <img src="plots/loc_count_vs_chrs.png" alt="Location Count Distribution">
        </center>

        <h2>Genome-wide distribution of insertions</h2>
        <center>
        <img src="plots/chromosome_layout.png" alt="Chromosome distributions">
        </center>

        <h2>Interactive Genome-wide plot</h2>
        <center>
        {{ manhattan }}
        </center>

        <h2>Distribution of CONTAMINANTS</h2>
        <center>
        <img src="plots/contaminants_pie_chart.png" alt="Contaminants Pie Chart">
        </center>

        <h2>Full Sample Table</h2>
        <div class = "full-samples">
            <div class="table" onclick="toggleContent('full-table')">
                <h4>Click to open full sample table</h4>
            </div>
            <div class="table" id="full-table">
    """

    # Add collapsible sections for each sample
    for index, row in df.iterrows():
        sample_id = f"sample_{index}"
        sample_name = row['SAMPLE']
        loc_count = row['LocCount']
        chrom = str(row['CHROMOSOME'])
        divclass = 'sample-highlight' if args.priority in chrom else 'sample-header'
        row_html = row.to_frame().to_html(classes='table table-striped', border=0)
        evidence = '<h4>There are too many locations to show plots. Please check the plots folder for this sample</h4>'
        if loc_count ==1:
            files = glob.glob(f'{args.output}/plots/{sample_name}/*.png')
            print(f'{args.output}/plots/{sample_name}/*.png')
            print(files)
            if not files:
                next
            else:
                evidence = f'<img src="{files[0]}" alt="{sample_name} read based evidence">'
        html_content += f"""
        <div class="sample-section">
            <div class="{divclass}" onclick="toggleContent('{sample_id}')">
                <table>
                <tr>
                <th>Sample: {sample_name}</th> 
                <th>Chromosome(s): {chrom}</th> 
                <th>Insertions: {loc_count}</th>
                </tr>
                </table>
            </div>
            <div class="sample-content" id="{sample_id}">
                {row_html}
                {evidence}
            </div>
        </div>
        """

    # Close the HTML tags
    html_content += """
    </div>
    </div>
    </body>
    </html>
    """

    # TODO: add value_counts of chromosome additions
    template = Template(html_content)
    html_report = template.render(
        manhattan=fig.to_html(full_html=False),
        priority_table=priority_table,
        summary_stats=summary_stats
    )

    # Save the HTML content to a file
    with open(f'{args.output}/report.html', "w", encoding="utf-8") as f:
        f.write(html_report)

    print(f'Generated HTML report as {args.output}/report.html')

def plotSimpleChrGraph(fai, gene_bed, output):
    chrom_height = 1

    # Spacing between consecutive ideograms
    chrom_spacing = 1
    # Height of the gene track. Should be smaller than `chrom_spacing` in order to
    # fit correctly
    gene_height = 0.4
    # Padding between the top of a gene track and its corresponding ideogram
    gene_padding = 0.1
    # Width, height (in inches)
    figsize = (6, 8)
    # Decide which chromosomes to use
    list_chromosomes, list_length = get_chromosomes_names(fai)
    chr_assoc = {k : v for k, v in zip(list_chromosomes, list_length)}
    #list_chromosomes = [x for x, j in sorted(chr_assoc.items(), key=lambda item: item[1], reverse=True)]
    chromosome_list = list_chromosomes[:41]
    #list_length = [v for k, v in sorted(chr_assoc.items(), key=lambda item: item[1], reverse=True)]
    chromosome_size = list_length[:41]

    # Keep track of the y positions for ideograms and genes for each chromosome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}

    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing

    # Read in ideogram.txt, downloaded from UCSC Table Browser
    colcycle = cycle([ '#bd2309', '#bbb12d', '#1480fa', '#14fa2f', '#000000',
            '#faf214', '#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
            '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9' ])
    ideo = pd.DataFrame({'chrom' : chromosome_list,
                            'start' : [0 for x in range(len(chromosome_list))],
                            'width' : chromosome_size,
                            'colors': [next(colcycle) for x in range(len(chromosome_list))]})

    # Note, I am plotting exact intervals here instead of windows
    # If the plots are too washed out, I may have to return to window-based clustering
    genes = straight_plot_df(chromosome_list, chromosome_size, gene_bed, ['#EAF2F8', '#E2E73C', '#005295', '#2B0C1C'])

    #print(genes.head())

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Now all we have to do is call our function for the ideogram data...
    #print("adding ideograms...")
    for (xranges, yranges, colors) in chromosome_collections(ideo, chrom_ybase, chrom_height):
        ax.broken_barh(xranges, yranges, facecolors=colors)
        #ax.add_collection(collection)

    # ...and the gene data
    #print("adding genes...")
    for (xranges, yranges, colors) in chromosome_collections(
        genes, gene_ybase, gene_height):
        ax.broken_barh(xranges, yranges, color='black', alpha=0.5, linewidths=0)
        #ax.add_collection(collection)

    # Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    plt.xlabel('Base Pair Position')
    ax.axis('tight')

    fig.savefig(output)

def straight_plot_df(chromosome_list, chromosome_size, bed, colors):
    gtable = defaultdict(list)
    cycler = cycle(colors)
    names = set()
    with open(bed, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            if s[0] not in chromosome_list:
                continue
            gtable['chrom'].append(s[0])
            gtable['start'].append(int(s[1]))
            gtable['end'].append(int(s[2]))
            gtable['name'].append(s[3])
            names.add(s[3])
            
    ctranslate = dict()
    for i in names:
        ctranslate[i] = next(cycler)
        
    genes = pd.DataFrame(gtable)
    genes['colors'] = genes['name'].map(ctranslate)
    return genes
    
def chromosome_collections(df, y_positions, height, min_width=20000):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    df.loc[df['width'] < min_width, 'width'] = min_width
    for chrom, group in df.groupby('chrom'):
        #print(chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield (xranges, yrange, group['colors'])
    if del_width:
        del df['width']

def get_chromosomes_names(input):
    list_chromosomes = []
    list_length = []
    with open(input, 'r') as fai:
        for l in fai:
            segs = l.rstrip().split()
            list_chromosomes.append(segs[0])
            list_length.append(int(segs[1]))
    return list_chromosomes, list_length

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)