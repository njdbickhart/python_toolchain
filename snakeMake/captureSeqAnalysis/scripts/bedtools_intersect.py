import pathlib
import re
import subprocess as sp
import argparse

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to consolidate capture seq information for gene editing"
            )
    parser.add_argument('-f', '--file', 
                        help="Input refinement file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output overlap file name. Tab delimited",
                        required=True, type=str,
                        )
    parser.add_argument('-b', '--bed',
                        help="Temporary bed file",
                        required=True, type=str
                        )
    parser.add_argument('-g', '--genes',
                        help="Gene bed for overlap",
                        required=True, type=str
                        )
    return parser.parse_args(), parser


def main(args, parser):
    # Create temp bam
    with open(args.file, 'r') as input, open(args.bed, 'w') as bed:
        entries = list()
        input.readline()
        for l in input:
            s = l.rstrip().split()
            oseg = re.split(r'[:-]', s[1])
            entries.append(oseg)
        entries = sorted(entries, key=lambda x: (x[0], int(x[1])))
        for oseg in entries:
            bed.write("\t".join(oseg) + "\n")

    # Generate the intersection
    cmd = ['bedtools', '-a', args.bed, '-b', args.genes, '-wb']
    sp.run(' '.join(cmd), shell=True, stdout=open(f"{args.bed}.out", "w"))

    # Read the intersection and print out intersecting genes
    with open(f'{args.bed}.out', 'r') as input, open(args.output, 'w') as out:
        for l in input:
            s = l.rstrip().split()
            out.write("\t".join([f'{s[0]}:{s[1]}-{s[2]}', s[-1]]) + "\n")

    tempbed = pathlib.Path(f'{args.bed}.out')
    tempbed.unlink(missing_ok=True)
    print('Fini')

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)