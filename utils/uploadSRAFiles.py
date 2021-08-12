# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:12:08 2021

@author: derek.bickhart-adm
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
