# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:12:08 2021

@author: derek.bickhart-adm
You can submit your files to NCBI SRA in 32 easy(tm) steps!
1. Fill in metadata on the website
2. get your ftp data at the end of metadata submission
3. login using an alternative FTP client (ie. filezilla) and create a new folder in your upload directory
4. make sure that the python ftplib package is installed (/KEEP/rumen_longread_metagenome_assembly/aspera on Ceres)
5. Run this script with the ftp login information provided by NCBI, and the path to your upload copied right from their site, but with the final portion of the path terminating with the folder name you just created using another ftp client
6. Watch out for log messages throughout the script to ensure that your files are being copied
7-31. bang head against desk repeatedly
32. Breathe a sigh of relief -- your torture is over... for now :) 
"""

import argparse
import ftplib
import sys
import os

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A tool to automate painful NCBI SRA file transfers"
            )
    parser.add_argument('-f', '--file', 
                        help="Input file to store on their server",
                        action="append", default=[],
                        )
    parser.add_argument('-d', '--directory',
                        help="NCBI sra ftp upload directory",
                        required=True, type=str,
                        )
    parser.add_argument('-u', '--username',
                        help="NCBI sra username",
                        required=True, type=str,
                        )
    parser.add_argument('-p', '--password',
                        help="NCBI sra password",
                        required=True, type=str,
                        )
    parser.add_argument('-l', '--url',
                        help="NCBI sra ftp host url",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    if len(args.file) == 0:
        print("Must supply at least one file!")
        sys.exit(-1)
        
    baseLookup = {k : os.path.basename(k) for k in args.file}
    
    with ftplib.FTP(args.url, args.username, args.password) as ftp:
        ftp.cwd(args.directory)
        for file, base in baseLookup.items():
            with open(file, 'rb') as f:
                print(f'Storing {file} as {base} on server...')
                ftp.storbinary(f'STOR {base}', f)
                
    print("Finished! Hope it worked!")

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
