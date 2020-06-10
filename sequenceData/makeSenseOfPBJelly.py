import json, re, sys, os
from collections import defaultdict
usage = "python {} gap_status.txt <assembly folder path> <output file>".format(sys.argv[0])

if len(sys.argv) < 3:
    print(usage)
    sys.exit(-1)

keeps = ["predictedGapSize", "fillBases", "contribSeqs", "contribBases", "spanCount", "avgSpanBases", "extendSeq1", "extendSeq2"]
data = defaultdict(dict)
gcat = dict()
with open(sys.argv[1], 'r') as gstatus:
    for l in gstatus:
        s = l.rstrip().split()
        gcat[s[0]] = s[1]
        jtemp = "{}/{}/fillingMetrics.json".format(sys.argv[2], s[0])
        if os.path.exists(jtemp):
            with open(jtemp, 'r') as jfile:
                temp = json.load(jfile)
                for k in keeps:
                    data[s[0]][k] = temp.get(k, "NA")
        else:
            for k in keeps:
                data[s[0]][k] = "Nope"
            
with open(sys.argv[3], 'w') as output:
    output.write("GapID\tGapCat\t" + "\t".join(keeps) + "\n")
    for k, v in data.items():
        order = [str(v[x]) for x in keeps]
        output.write(k + "\t" + gcat[k] + "\t" +  "\t".join(order) + "\n");
