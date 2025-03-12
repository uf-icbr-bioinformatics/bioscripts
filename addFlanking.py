#!/usr/bin/env python

import sys, csv

class Reference(object):
    fastafile = ""
    refname = ""
    ref = []
    upcase = True

    def __init__(self, fastafile):
        self.fastafile = fastafile

    def getFlanks(self, chrom, pos, size=15):
        """Return the flanking sequence around position `pos' on chromosome `chrom'. 
The flanking sequence is composed of `size' bases before and after position, so its
length is size*2+1. Assumes positions start at 1."""
        if chrom != self.refname:
            self.loadReference(chrom)
            self.refname = chrom
        frag = self.ref[pos-size-1:pos+size]
        if self.upcase:
            for i in range(size):
                frag[i] = frag[i].lower()
                frag[size+i+1] = frag[size+i+1].lower()
        return "".join(frag)

    def loadReference(self, chrom):
        chromlen = 0
        startpos = 0
        fai = self.fastafile + ".fai"
        with open(fai, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                if line[0] == chrom:
                    chromlen = int(line[1])
                    startpos = int(line[2])
        if not startpos:
            return

        self.ref = ['']*chromlen

        sys.stderr.write("Loading reference for {}...\n".format(chrom))
        with open(self.fastafile, "r") as f:
            f.seek(startpos)
            idx = 0
            while True:
                line = f.readline()
                if not line:
                    break
                if line[0] == '>':
                    break
                line = line.rstrip("\n")
                for ch in line:
                    self.ref[idx] = ch.upper()
                    idx += 1
        return self.ref

    def main(self, infile, outfile, chromCol=0, posCol=1):
        with open(infile, "r") as f, open(outfile, "w") as out:
            c = csv.reader(f, delimiter='\t')
            hdr = c.next()
            out.write("\t".join(hdr) + "\tFlanking\n")
            for line in c:
                chrom = line[chromCol]
                pos = int(line[posCol])
                f = self.getFlanks(chrom, pos)
                out.write("\t".join(line) + "\t" + f + "\n")

if __name__ == "__main__":
    args = sys.argv[1:]
    R = Reference(args[0])
    R.main(args[1], args[2])
