# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 07:44:19 2020

@author: derek.bickhart-adm
"""

import argparse
from Bio import SeqIO

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A conversion tool to change from genbank to bed file formats"
            )
    parser.add_argument('-g', '--genbank', 
                        help="Input genbank file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output bed file name",
                        required=True, type=str,
                        )
    parser.add_argument('-s', '--split',
                        help="Delimiter to split ID name or feature [default = None]",
                        type=str, default=None
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    counter = 0
    shouldSplit = False
    if args.split != None:
        shouldSplit = True
    with open(args.output, 'w') as out:
        for record in SeqIO.parse(open(args.genbank, 'rU'), "genbank"):
            for feature in record.features:
                counter += 1
                if counter % 1000 == 0:
                    print(f'Processed: {counter} features...')
                if feature.type == 'gene' or feature.type == 'CDS':
                    start = feature.location.start.position
                    end = feature.location.end.position
                    chrom = ""
                    name = ""
                    orient = ""
                    val = ""
                    
                    try:
                        chrom = feature.qualifiers["gene"][0]
                    except:
                        chrom = feature.qualifiers["locus_tag"][0]
                    
                    if feature.strand < 0:
                        orient = "-"
                    else:
                        orient = "+"
                        
                    if shouldSplit:
                        nsegs = chrom.split(args.split)
                        if len(nsegs) > 1:
                            name = chrom
                            val = nsegs[1]
                            chrom = nsegs[0]
                    else:
                        name = feature.type
                    
                    out.write(f'{chrom}\t{start}\t{end}\t{name}\t{val}\t{orient}\n')

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
