# This script is designed to run in a snakemake pipeline
# The goal is to filter away any "non-host tax" score that is not unclassified
from Bio import SeqIO

class cententry:

    def __init__(self, scaff, taxid, score, qlen):
        self.scaff = scaff
        self.taxid = taxid
        self.score = score
        self.qlen = qlen

    def meetsCriteria(self, hosttax):
        if self.taxid == hosttax:
            return True
        elif self.score < self.qlen:
            return True # The score is a sum of hits weighted by length
        else:
            return False

hosttax = snakemake.params['hosttax']
sample = snakemake.params['sample']

data = {}
# Read in the centrifuge.out file
with open(snakemake.input['cent'], 'r') as input:
    # Get rid of the header
    head = input.readline()
    for l in input:
        l = l.rsplit()
        s = l.split()

        data[s[0]] = cententry(s[0], s[2], int(s[3]), int(s[6]))

# Now, read in the scaffold fasta file and filter it based on the criteria
passfilt = 0
failfilt = 0
with open(snakemake.output['fasm'], 'w') as out:
    for record in SeqIO.parse(snakemake.input['rasm'], "fasta"):
        if record.id in data:
            if data[record.id].meetsCriteria(hosttax):
                SeqIO.write(record, out, "fasta")
                passfilt += 1
            else:
                failfilt += 1
        else:
            print("Error parsing record: " + record.id)

print(f'Passfilter: {passfilt} Failfilter: {failfilt}')
