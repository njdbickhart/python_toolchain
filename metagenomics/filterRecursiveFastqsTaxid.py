# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 15:45:53 2023

@author: Derek.Bickhart
"""

import argparse
import networkx as nx
import re
import sys
import gzip


def arg_parse():
    parser = argparse.ArgumentParser(
            description = "Filter host sequence from fastq files"
            )
    parser.add_argument('-f', '--fastq', 
                        help="Input fastq file to filter; can be gzipped",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name; will be gzipped",
                        required=True, type=str,
                        )
    parser.add_argument('-t', '--taxid',
                        help="Parent taxid to filter; recursively inclusive; can be selected more than once",
                        action="append", default=[],
                        )
    parser.add_argument('-c', '--centrifuge',
                        help="Centrifuge read association file",
                        required=True, type=str
                        )
    parser.add_argument('-n', '--nodes',
                        help="NCBI nodes.dmp file",
                        required=True, type=str
                        )
    return parser.parse_args(), parser

def main(args, parser):
    # Creating worker
    worker = tax_traversal(args.nodes, args.taxid)
    
    # Process the nodes file
    worker.process_nodes()
    
    # Process the taxID file
    worker.filter_taxid()
    
    # Now filter the input fastq and print out the results
    worker.filter_fastq(args.fastq, args.output)
 
def fastq_reader_fh(infile):
    name = infile.readline().rstrip()
    while True:
        seq = ""
        for s in infile:
            if s[0] == '+':
                break
            else:
                seq += s.rstrip()
        qual = ""
        for q in infile:
            if len(qual) > 0 and  q[0] == '@':
                yield name, seq, qual
                name = q.rstrip()
                break
            else:
                qual += q.rstrip()
        else:
            yield name, seq, qual
            return

def smartFile(filename : str, mode : str = 'r'):
    fh = None
    if filename.endswith('.gz') and mode == 'r':
        fh = gzip.open(filename, mode='rt')
    elif filename.endswith('.gz') and mode == 'w':
        fh = gzip.open(filename, mode='wt')
    else:
        fh = open(filename, mode)
    return fh

class tax_traversal:

    def __init__(self, nodes_file, taxlist):
        self.nodes_file = nodes_file
        self.taxlist = taxlist
        
        # Public attributes
        self.graph = None
        self.diagtree = list()
        self.set = set() 
        
        self.rids = set()
        
    def process_nodes(self):
        self.graph = nx.DiGraph()
        fields = re.compile(r'\t|\t')
        with open(self.nodes_file, 'r') as input:
            for l in input:
                s = re.split(fields, l)
                self.graph.add_edge((s[1], s[0]))
                
        # QC and validity testing
        if not nx.is_directed_acyclic_graph(self.graph):
            print('Error! Nodes.dmp Graph did not validate as a DAG!')
            sys.exit(-1)
        else:
            print('Nodes graph loaded')
      
    def filter_taxid(self):
        for i in self.taxlist:
            if self.graph.has_node(i):
                nset = self.graph.descendants(i)
                print(f'Identified {len(nset)} descendants from taxid {i}')
                self.set.update(nset)
        print(f'Total set entries comprise {len(self.set)} elements.')
        
    def identify_valid_reads(self, cfile):
        rtotal = 0
        with open(cfile, 'r') as input:
            input.readline() # Clear header
            for l in input:
                rtotal += 1
                s = l.rstrip().split()
                if s[2] not in self.set:
                    self.rids.add(s[0])
                    
        prop = len(self.rids) / rtotal * 100.0
        print(f'Retained {len(self.rids)} from {rtotal} reads \({prop})\%\)')
            
    def filter_fastq(self, fastq, outfile):
        ftotal = 0
        totreads = 0
        input = smartFile(fastq, 'r')
        output = smartFile(outfile, 'w')
        for head, seq, qual in fastq_reader_fh(input):
            totreads += 1
            if head[1:] in self.rid:
                output.write(f'{head}\n{seq}\n+\n{qual}\n')
            else:
                ftotal += 1
        input.close()
        output.close()
        print(f'Filtered {ftotal} reads out of a total of {totreads} \({ftotal / totreads}\)')
        
    
if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)