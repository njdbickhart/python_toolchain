# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:09:21 2019

@author: dbickhart
"""

import argparse;
import concurrent.futures
from queue import Queue, Empty
from threading import Thread

# Translation table for rev comp
transtable = str.maketrans('ACGTacgt', 'TGCATGCA')

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Sample \"leapfrog\" reads from long read sequence data to overcome mapping bias"
            )
    parser.add_argument('-f', '--file', 
                        help="The input longread sequence fasta file",
                        type=str, required=True
                        )
    parser.add_argument('-s', '--forstrand', 
                        help="The output forward fastq strand",
                        type=str, required=True
                        )
    parser.add_argument('-r', '--revstrand',
                        help="The output reverse fastq strand",
                        type=str, required=True
                        )
    parser.add_argument('-l', '--len',
                        help="The length of leapfrog reads and the span of inter-read space",
                        type=int, default=150
                        )
    parser.add_argument('-t', '--threads',
                        help="The number of threads to use in processing the file",
                        type=int, default=4
                        )
    return parser.parse_args()

def main(args):
    # Turning this into a single thread process before I debug
    
    qstr = 'I' * args.len
    fq1 = open(args.forstrand, 'w')
    fq2 = open(args.revstrand, 'w')
    #fq1 = SafeWriter(args.forstrand, 'w')
    #fq2 = SafeWriter(args.revstrand, 'w')
    
    rcount = list()
    fp = open(args.file, "r")
    """
    #Uncomment for parallel processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads - 2) as executor:
        for name, seq in read_fasta(fp):
            if seq < (args.len * 3):
                continue # Skip reads too short to leapfrog
            count = executor.submit(convert_to_fastq, name, seq, args.len, fq1, fq2, qstr)
            rcount.append(count.result())
    """
    for name, seq in read_fasta(fp):
        if seq < (args.len *3):
            continue
        count = convert_to_fastq(name, seq, args.len, fq1, fq2, qstr)
        rcount.append(count)
        
    fp.close()
    
    print(f'Processed {len(rcount)} reads with an average of {avg(rcount)} leapfrog pairs')
    fq1.close()
    fq2.close()
    

def avg(counts):
    return sum(counts) / len(counts)
    

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()        
        if line.startswith(">"):
            segs = line.split()
            segs[0] = segs[0].replace(">", "")
            if name: yield(name, ''.join(seq))
            name, seq = segs[0], []
        else:
            seq.append(line)
    if name: yield(name, ''.join(seq))
    
def revcomp(seq : str) -> str:
    revstr = seq[::-1]
    return revstr.translate(transtable)
    
def convert_to_fastq(name, seq, length, fq1, fq2, qualitystr) -> int:
    scount = 0 # Count of subdivisions
    for i in range(0, len(seq) - (length * 3), length):
        fr1 = seq[i:i+length]
        fr2 = revcomp(seq[i + (length*2):i + (length*3)])
        
        fq1.write(f'@{name}_{scount}\n{fr1}\n+\n{qualitystr}\n')
        fq2.write(f'@{name}_{scount}\n{fr2}\n+\n{qualitystr}\n')
        scount += 1
    return scount
    

class SafeWriter:
    """
    Props to Leonid Mednikov on Stack Overflow for this class
    """
    def __init__(self, *args):
        self.filewriter = open(*args)
        self.queue = Queue()
        self.finished = False
        Thread(name = "SafeWriter", target=self.internal_writer).start()  

    def write(self, data):
        self.queue.put(data)

    def internal_writer(self):
        while not self.finished:
            try:
                data = self.queue.get(True, 1)
            except Empty:
                continue    
            self.queue.task_done()
            self.filewriter.write(data)

    def close(self):
        self.queue.join()
        self.finished = True
        self.filewriter.close()
        
        
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
