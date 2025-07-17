
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template
from weasyprint import HTML
import argparse

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to generate a report file from capture-seq data"
            )
    parser.add_argument('-f', '--file', 
                        help="Input final table file",
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

    return parser.parse_args(), parser


def main(args, parser):
    # Create plot folders in case they do not exist
    os.makedirs(f'{args.output}/plots', exist_ok=True)

    # Load the CSV file
    df = pd.read_csv(args.file)

    # Modify the dataframe to create statistics
    df['CHRCount'] = df['CHROMOSOME'].astype(str).apply(lambda x: len(x.split(';')))
    df['LocCount'] = df['SITES'].astype(str).apply(lambda x: len(x.split(';')))

    # Sort DF based on priority chromosome
    df = pd.concat([df[df['CHROMOSOME'].astype(str).str.contains(args.priority)].sort_values('LocCount'), 
                    df[~df['CHROMOSOME'].astype(str).str.contains(args.priority)].sort_values(['CHROMOSOME','LocCount'])])

    # Summary statistics
    summary_stats = df.describe(include='all').to_html(classes='table table-striped', border=0)


    # Distribution plot of CAL_DEPTH
    plt.figure(figsize=(10, 6))
    sns.histplot(df['CAL_DEPTH'].dropna(), bins=30, kde=True)
    plt.title('Distribution of CAL_DEPTH')
    plt.xlabel('CAL_DEPTH')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.savefig(f"{args.output}/plots/cal_depth_distribution.png")
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
            table {border: 2 px solid black; width: 100%; border-collapse: collapse; text-align: left}
            th, td {padding: 20 px;}
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

        <h2>Summary Statistics</h2>
        {{ summary_stats | safe}}
        
        <h2>Distribution of CAL_DEPTH</h2>
        <img src="plots/cal_depth_distribution.png" alt="CAL_DEPTH Distribution">

        <h2>Distribution of CONTAMINANTS</h2>
        <img src="plots/contaminants_pie_chart.png" alt="Contaminants Pie Chart">

        <h2>Full Sample Table</h2>
    """

    # Add collapsible sections for each sample
    for index, row in df.iterrows():
        sample_id = f"sample_{index}"
        sample_name = row['SAMPLE']
        loc_count = row['LocCount']
        chrom = str(row['CHROMOSOME'])
        divclass = 'sample-highlight' if args.priority in chrom else 'sample-header'
        row_html = row.to_frame().to_html(header=False)
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
            </div>
        </div>
        """

    # Close the HTML tags
    html_content += """
    </body>
    </html>
    """

    # TODO: add value_counts of chromosome additions
    template = Template(html_content)
    html_report = template.render(
        summary_stats=summary_stats
    )

    # Save the HTML content to a file
    with open(f'{args.output}/report.html', "w", encoding="utf-8") as f:
        f.write(html_report)

    print(f'Generated HTML report as {args.output}/report.html')

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)