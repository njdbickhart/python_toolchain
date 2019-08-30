# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:51:51 2019

@author: dbickhart
"""

import argparse
import PyPDF2
import re

# compiled, global regular expressions
author = re.compile(r'Authorship:\s*(\d{1,2})')
tfield = re.compile(r'Type:\s*(.)')
lognum = re.compile(r'Log Number:\s*(\d+)')
submit = re.compile(r'Date Submitted to Journal:\s*(\S*)$')
accept = re.compile(r'Journal Acceptance Date:\s*(\S*)$')
publish = re.compile(r'Journal Publication Date:\s*(\S*)$')
journal = re.compile(r'Journal or Equivalent:')
title = re.compile(r'Title of Manuscript')
stopfield = re.compile(r'Authors / Category')

# Dictionary to enable easy looping later
reDict = {'author' : author,
          'tfield' : tfield,
          'lognum' : lognum,
          'submit' : submit,
          'accept' : accept,
          'publish' : publish
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
        
        inEntry = False
        inJournal = False
        inTitle = False
        entry = None
        # Read through all lines of the page
        for l in page.extractText().splitlines():
            l = l.rstrip()
            if author.match(l):
                inEntry = True
                entry = arisEntry()
                
            if inEntry and inJournal:
                entry.populate('journal', l)
                inJournal = False
                
            if inEntry and inTitle:
                entry.populate('title', l)
                inTitle = False
                
            if inEntry:
                for f,r in reDict.items():
                    if r.match(l):
                        m = r.match(l)
                        entry.populate(f, m.group(1))
                
                if journal.match(l):
                    inJournal = True
                    
                if title.match(l):
                    inTitle = True
                    
                if stopfield.match(l):
                    inEntry = False
                    entryList.append(entry)
                    
    # Now to print out
    with open(args.output, 'w') as out:
        for e in entryList:
            out.write(e.printout + '\n')
                    
    
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
        return '\t'.join([self.author, self.tfield, self.lognum, self.submit, 
                          self.accept, self.publish, self.journal, self.title])
    
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
