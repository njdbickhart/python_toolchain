# This script is designed to run in a snakemake pipeline
# The goal is to filter away any "non-host tax" score that is not unclassified
from Bio import SeqIO
import sys

usage = "python3 {} [newline list of taxids] [kraken-out from Centrifuge] [input fastq file to filter] [output filtered fastq]".format(sys.argv[0])

class cententry:

    def __init__(self, scaff, taxid, score, qlen):
        self.scaff = scaff
        self.taxid = taxid
        self.score = score
        self.qlen = qlen

    def meetsCriteria(self, taxset):
        if self.taxid in taxset:
            return True
        elif self.score < self.qlen:
            return True # The score is a sum of hits weighted by length
        else:
            return False

if len(sys.argv) < 5:
    print(usage)
    sys.exit()

taxidfile = sys.argv[1]
centrifuge_file = sys.argv[2]
original_fastq = sys.argv[3]
filtered_output = sys.argv[4]

# Generate set of appropriate taxids
taxids = set()
with open(taxidfile, 'r') as input:
    for l in input:
        l = l.rstrip()
        taxids.add(l)

print("Loaded taxid file")

data = {}
print(f'Parsing Centrifuge file')
# Read in the centrifuge.out file
with open(centrifuge_file, 'r') as cent:
    # Get rid of the header
    head = cent.readline()
    for l in cent:
        l = l.rstrip()
        s = l.split()

        data[s[0]] = cententry(s[0], s[2], int(s[3]), int(s[6]))

# Now, read in the scaffold fasta file and filter it based on the criteria
passfilt = 0
failfilt = 0
with open(filtered_output, 'w') as out:
    for record in SeqIO.parse(original_fastq, "fastq"):
        if record.id in data:
            if data[record.id].meetsCriteria(taxids):
                SeqIO.write(record, out, "fastq")
                passfilt += 1
            else:
                #SeqIO.write(record, out, "fasta") # Added just to progress pipeline
                failfilt += 1
        else:
            print("Error parsing record: " + str(record.id) + '\n')

printf'Passfilter: {passfilt} Failfilter: {failfilt}\n')
