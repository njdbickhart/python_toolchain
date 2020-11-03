# This script is designed to run in a snakemake pipeline
# The goal is to filter away any "non-host tax" score that is not unclassified
from Bio import SeqIO

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

taxidfile = snakemake.params['taxids']
sample = snakemake.params['samp']
log = open(snakemake.log[0], 'a+')
log.write("Beginning filtration")

# Generate set of appropriate taxids
taxids = set()
with open(taxidfile, 'r') as input:
    for l in input:
        l = l.rstrip()
        taxids.add(l)

log.write("Loaded taxid file")

data = {}
log.write(f'Parsing Centrifuge file for {sample}\n')
# Read in the centrifuge.out file
with open(snakemake.input['cent'], 'r') as cent:
    # Get rid of the header
    head = cent.readline()
    for l in cent:
        l = l.rstrip()
        s = l.split()

        data[s[0]] = cententry(s[0], s[2], int(s[3]), int(s[6]))

# Now, read in the scaffold fasta file and filter it based on the criteria
passfilt = 0
failfilt = 0
with open(snakemake.output['fasm'], 'w') as out:
    for record in SeqIO.parse(snakemake.input['rasm'], "fasta"):
        if record.id in data:
            if data[record.id].meetsCriteria(taxids):
                SeqIO.write(record, out, "fasta")
                passfilt += 1
            else:
                #SeqIO.write(record, out, "fasta") # Added just to progress pipeline
                failfilt += 1
        else:
            log.write("Error parsing record: " + str(record.id) + '\n')

log.write(f'Passfilter: {passfilt} Failfilter: {failfilt}\n')
log.close()
