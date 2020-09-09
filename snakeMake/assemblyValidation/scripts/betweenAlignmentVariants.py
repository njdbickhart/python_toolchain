# Originally by Mike Schatz, modified by Maria Nattestad and translated to a hipster language by Derek Bickhart!
# github.com/marianattestad/assemblytics
import os
import argparse
from collections import defaultdict

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Generate between paf alignment statistics"
            )
    parser.add_argument('-f', '--file',
                        help="A single PAF alignment file between two assemblies.",
                        type=str, required=True
                        )
    parser.add_argument('-m', '--minimum',
                        help="Minimum event size to filter",
                        type=int, default=100
                        )
    parser.add_argument('-a', '--maximum',
                        help="Maximum event size to filter",
                        type=int, default=100000
                        )
    parser.add_argument('-n', '--narrow',
                        help="How close alignments need to be to call dels",
                        type=int, default=50
                        )
    parser.add_argument('-q', '--qdist',
                        help="How far alignments need to be from each other before we discard",
                        type=int, default=100000
                        )
    parser.add_argument('-o', '--output',
                        help="Output file Name",
                        type=str, required=True
                        )

    return parser.parse_args(), parser

class pafLine:

    def __init__(self, line):
        segs = line.rstrip().split()
        if len(segs) < 12:
            self.valid = False
            # to allow filtration of malformed lines later
            return

        # Essential attributes
        self.valid = True
        self.rstart = int(segs[7])
        self.rend = int(segs[8])
        self.qstart = int(segs[2])
        self.qend = int(segs[3])
        self.rlen = int(segs[6])
        self.qlen = int(segs[1])
        self.rid = segs[5]
        self.qid = segs[0]
        self.qidx = 0
        self.qrc = True if segs[4] == '-' else False

    def getqlen(self):
        return self.qend - self.qstart

    def getrlen(self):
        return self.rend - self.rstart

class pafComp:

    def __init__(self, prev, curr):
        self.prev = prev
        self.curr = curr
        self.rdist = 0
        self.qdist = 0
        self.positions = []
        self.totdist = 0
        self.svtype = ""
        self.typeguess = ""
        self.abs_size = 0

    def compareOrient(self):
        if not self.prev.qrc and not self.curr.qrc:
            self.svtype = "FF"
            self.qdist = self.curr.qstart - self.prev.qend
            self.rdist = self.curr.rstart - self.prev.rend
            self.positions = [self.prev.rend, self.curr.rstart]
        elif self.prev.qrc and self.curr.qrc:
            self.svtype = "RR"
            self.qdist = self.prev.rstart - self.curr.rend
            self.rdist = self.curr.qend - self.prev.qstart
            self.positions = [self.prev.rstart, self.curr.rend]
        elif not self.prev.qrc and self.curr.qrc:
            self.svtype = "FR"
            self.qdist = self.curr.qend - self.prev.qend
            self.rdist = self.curr.rstart - self.prev.rend
            self.positions = [self.prev.rend, self.curr.rend]
        elif self.prev.qrc and not self.curr.qrc:
            self.svtype = "RF"
            self.qdist = self.prev.qend - self.curr.qend
            self.rdist = self.curr.rstart - self.prev.rend
            self.positions = [self.prev.rstart, self.curr.rstart]
        else:
            self.svtype = "ER"
        self.totdist = self.rdist + self.qdist

    def guessType(self, args):
        abs_event_size = abs(self.rdist - self.qdist)

        if self.prev.rid != self.curr.rid:
            self.typeguess = "Interchromosomal"
            self.rdist = 0
        else:
            if self.svtype == "FR" or self.svtype == "RF":
                self.typeguess = "Inversion"
                abs_event_size = self.rdist
            elif self.qdist > self.rdist:
                if self.rdist > -1 * args.narrow and self.rdist < args.narrow and self.qdist > -1 * args.narrow:
                    self.typeguess = "Insertion"
                else:
                    if self.rdist < 0 or self.qdist < 0:
                        self.typeguess = "Tandem_expansion"
                    else:
                        self.typeguess = "Repeat_expansion"
            elif self.qdist < self.rdist:
                if self.rdist > -1 * args.narrow and self.qdist > -1 * args.narrow and self.qdist < args.narrow:
                    self.typeguess = "Deletion"
                else:
                    if self.rdist < 0 or self.qdist < 0:
                        self.typeguess = "Tandem_contraction"
                    else:
                        self.typeguess = "Repeat_contraction"
            else:
                self.typeguess = "None"
        self.abs_size = abs_event_size

        if self.abs_size > args.maximum:
            self.typeguess = "Longrange"
            if abs(self.qdist) > args.qdist:
                self.typeguess = "None"

    def isNone(self):
        return True if self.typeguess == "None" else False

    def getOutStr(self, counter):
        ref_start = min(self.positions)
        ref_end = max(self.positions)
        ref_end = ref_start + 1 if ref_end == ref_start else ref_end
        return f'{self.prev.rid}\t{ref_start}\t{ref_end}\taqc_sv{counter}\t{self.abs_size}\t+\t{self.typeguess}\t{self.rdist}\t{self.qdist}\t-\tbetween\n'


def main(args, parser):
    # Load data
    pdata = defaultdict(lambda: defaultdict(list))
    qids = set()
    rids = set()
    with open(args.file, 'r') as paf:
        for l in paf:
            temp = pafLine(l)
            if temp.valid:
                pdata[temp.qid][temp.rid].append(temp)
                qids.add(temp.qid)
                rids.add(temp.rid)

    # start organization by query sequence
    with open(args.output, 'w') as out:
        svcounter = 0
        out.write("chrom\tstart\tstop\tname\tsize\tstrand\ttype\tref.dist\tquery.dist\tcontig_position\tmethod.found\n")
        for q in sorted(qids):
            working = []
            for r in sorted(rids):
                if r in pdata[q]:
                    working.extend(pdata[q][r])

            working.sort(key=lambda x: x.qstart)

            if len(working) > 1:
                for i in range(1,len(working)):
                    prev = working[i-1]
                    curr = working[i]

                    if prev.getrlen() >= args.minimum and curr.getrlen() >= args.minimum:
                        # Create workhorse class
                        comp = pafComp(prev, curr)

                        # Specify SV orientation
                        comp.compareOrient()

                        # Guess SV type based on context
                        comp.guessType(args)

                        # print out SV if it passes filter
                        if not comp.isNone():
                            svcounter += 1
                            out.write(comp.getOutStr(svcounter))

    print(f'Identified {svcounter} between alignment SVs')

if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
