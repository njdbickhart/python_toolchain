# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:51:51 2019

@author: dbickhart
"""

import argparse
import PyPDF2
import re

# compiled, global regular expressions
author = re.compile(r'Authorship:\s*(\d{1,2})Type')
tfield = re.compile(r'Type:\s*(.{1})Official')
lognum = re.compile(r'Log Number:\s*(\d+)Authorship')
submit = re.compile(r'Date Submitted to Journal:\s*(\S*)Log Number')
accept = re.compile(r'Journal Acceptance Date:\s*(\S*)Journal or')
publish = re.compile(r'Journal Publication Date:\s*(\S*)Date Submitted')
journal = re.compile(r'Journal or Equivalent:(.+)Title of Manuscript')
title = re.compile(r'Title of Manuscript:(.+)Authors / Category')
stopfield = re.compile(r'Authors / Category')

# Dictionary to enable easy looping later
reDict = {'author' : author,
          'tfield' : tfield,
          'lognum' : lognum,
          'submit' : submit,
          'accept' : accept,
          'publish' : publish, 
          'journal' : journal,
          'title' : title
          }

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Process ARS 115 report pdf into a tabular spreadsheet"
            )
    parser.add_argument('-f', '--file', 
                        help="The input 115 pdf file",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output', 
                        help="Output tab delimited text file",
                        type=str, required=True
                        )
    return parser.parse_args()

def main(args):
    reader = PyPDF2.PdfFileReader(open(args.file, 'rb'))
    entryList = list()
    
    # Read through all pages 
    for i in range(reader.getNumPages()):
        page = reader.getPage(i)
        
        entry = None
        # Read through all lines of the page
        for l in page.extractText().splitlines():
            l = l.rstrip()
            entry = arisEntry()
                
            for f,r in reDict.items():
                if r.findall(l):
                    m = r.findall(l)
                    entry.populate(f, m[0])
                    
            entryList.append(entry)
                    
    # Now to print out
    with open(args.output, 'w') as out:
        out.write('\t'.join(['Authorship','Type', 'LogNum', 'Submission', 
                             'Acceptance', 'Published', 'Journal', 'Title']) + '\n')
        for e in entryList:
            out.write(e.printout() + '\n')
                    
    
class arisEntry:
    
    def __init__(self):
        self.author = "None"
        self.tfield = "None"
        self.lognum = "None"
        self.submit = "None"
        self.accept = "None"
        self.publish = "None"
        self.journal = "None"
        self.title = "None"
        
        self.lookup = {'author' : self.author,
                       'tfield' : self.tfield, 
                       'lognum' : self.lognum,
                       'submit' : self.submit,
                       'accept' : self.accept,
                       'publish' : self.publish,
                       'journal' : self.journal,
                       'title' : self.title
                       }
        
    def populate(self, field : str, value):
        self.lookup[field] = value
        
    def printout(self):
        return '\t'.join([self.lookup[x] for x in ['author', 'tfield', 'lognum',
                          'submit', 'accept', 'publish', 'journal', 'title']])
        """
        return '\t'.join([self.author, self.tfield, self.lognum, self.submit, 
                          self.accept, self.publish, self.journal, self.title])
        """
    
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
