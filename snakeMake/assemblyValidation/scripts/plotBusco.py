import argparse
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from collections import defaultdict
import pandas
import numpy as np
import re


def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Create a plot of the busco scores"
            )
    parser.add_argument('-o', '--output',
                        help="output file name",
                        type=str, required=True
                        )
    parser.add_argument('-a', '--asms',
                        help="Expected assembly combinations",
                        action="append", default=[]
                        )
    parser.add_argument('-b', '--buscos',
                        help="busco file",
                        action="append", default=[]
                        )

    return parser.parse_args(), parser

def main(args, parser):
    if len(args.asms) < 1 or len(args.buscos) < 1:
        print("Error! Must input at least one busco file and one assembly label!")
        parser.print_help()
        sys.exit()

    #asms = snakemake.params["asms"]
    #print(asms)

    #buscos = snakemake.input
    #print(list(buscos))

    data = defaultdict(list)

    for i, b in enumerate(args.buscos):
        with open(b, 'r') as input:
            for l in input:
                l = l.strip()
                if l.startswith('#'):
                    continue
                elif l.startswith('C'):
                    m = re.match(r'C:.+%\[S:(.+)%,D:(.+)%\],F:(.+)%,M:(.+)%,n:.+', l)
                    print(f'{args.asms[i]}\t{m.group(0)}')
                    data["Assembly"].append(args.asms[i])
                    data["CompleteSC"].append(float(m.group(1)))
                    data["CompleteDup"].append(float(m.group(2)))
                    data["Fragmented"].append(float(m.group(3)))
                    data["Missing"].append(float(m.group(4)))
                    break

    df = pandas.DataFrame(data)
    #df.set_index("Assembly")
    print(df.describe())

    fig, ax = plt.subplots()

    #ax = df[["CompleteSC", "CompleteDup", "Fragmented", "Missing"]].plot.hbar(stacked=True, edgecolor='none')
    ax = df.plot.barh(stacked=True, edgecolor='none')
    for rowNum, row in df.iterrows():
        xpos = 0
        for val in row[1:]:
            xpos += int(val)
            if val > 8:
                ax.text((xpos - int(val)) + (int(val) * 0.5) - 1, rowNum - 0.05, str(val), color = 'white')
        xpos = 0
    ax.legend(ncol=4, fontsize='small', loc='lower left', bbox_to_anchor=(0, 1))

    plt.yticks(np.arange(len(args.asms)), df["Assembly"])
    plt.savefig(args.output)

if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
