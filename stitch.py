#!/usr/bin/env python

import sys
import SeqUtils

class Stitcher(object):
    filename1 = None
    filename2 = None
    minoverlap = 9
    minident = 0.85
    filter1 = None
    filter2 = None
    filter1rc = None
    filter2rc = None
    filterlen = 0
    maxmismatch = 0
    _fqreader = None
    outfile = None
    
    def parseArgs(self, args):
        # sys.stderr.write("{} {}\n".format(self.filename1, self.filename2))
        if "-h" in args or "--help" in args:
            return False
        prev = ""
        for a in args:
            # sys.stderr.write("{} {}\n".format(a, prev))
            if prev == "-m":
                self.maxmismatch = int(a)
                prev = ""
            elif prev == "-p":
                self.minoverlap = int(a)
                prev = ""
            elif prev == "-i":
                self.minident = float(a)
                prev = ""
            elif prev == "-f":
                frags = a.split(":")
                if len(frags[0]) == len(frags[1]):
                    self.filter1 = frags[0]
                    self.filter2 = frags[1]
                    self.filter1rc = SeqUtils.revcomp(self.filter2)
                    self.filter2rc = SeqUtils.revcomp(self.filter1)
                    self.filterlen = len(self.filter1)
                else:
                    sys.stderr.write("Error: filters should be of the same length.\n")
                    sys.exit(1)
                prev = ""
            elif prev == "-o":
                self.outfile = a
                prev = ""
            elif a in ["-m", "-p", "-i", "-f", "-o"]:
                prev = a
            elif self.filename1 is None:
                # sys.stderr.write("1 " + a + "\n")
                self.filename1 = a
            else:
                # sys.stderr.write("2 " + a + "\n")
                self.filename2 = a
        sys.stderr.write("Fastqs:\n  {}\n  {}\n".format(self.filename1, self.filename2))
        sys.stderr.write("Stitching:\n  Min overlap = {}\n  Min ident = {}\n".format(self.minoverlap, self.minident))
        if self.filter1:
            sys.stderr.write("Filtering:\n  Left =  {}\n  Right = {}\n  Mismatch = {}\n".format(self.filter1, self.filter2, self.maxmismatch))
        self._fqreader = SeqUtils.PairedFastqReader(self.filename1, self.filename2)
        return self.filename1 and self.filename2

    def usage(self):
        sys.stdout.write("""stitch.py - stitch R1 and R2 reads from a pair of fastq files.

Usage: stitch.py [options] fastq_R1.fastq.gz fastq_R2.fastq.gz

Options:

  -o O | Write output to file O (default: standard output)
  -m M | Allow at most M mismatches in the overlap region (default: {})
  -p P | Minimum overlap length (default: {})
  -i I | Minimum % identity in the overlap region (default: {}).
    
""".format(self.maxmismatch, self.minoverlap, self.minident))
        
    def run(self):
        fq = self._fqreader
        seqid = 1
        nin = 0
        nstitched = 0
        if self.outfile:
            out = open(self.outfile, "w")
        else:
            out = sys.stdout
        try:
            while fq.isActive():
                nin += 1
                fq.nextRead()
                fq.fq2.seq = SeqUtils.revcomp(fq.fq2.seq)
                fq.fq2.qual = fq.fq2.qual[::-1]
                full = self.stitch(fq.fq1, fq.fq2)

                if full:
                    nstitched += 1
                    if self.filter1:
                        good = False
                        if self.filterMatch(full, False): # Normal
                            good = True
                        elif self.filterMatch(full, True): # RC
                            good = True
                    else:
                        good = True
                    if good:
                        out.write(">seq_{}\n{}\n".format(seqid, full))
                        seqid += 1
        finally:
            out.close()
        sys.stderr.write("Total sequences: {}\n".format(nin))
        sys.stderr.write("Stitched: {} ({}%)\n".format(nstitched, int(100.0 * nstitched / nin)))
        sys.stderr.write("Filtered: {} ({}%)\n".format(seqid, int(100.0 * seqid / nin)))

    def filterMatch(self, seq, rc):
        if rc:
            n = self.countMismatch(seq[:self.filterlen], self.filter1rc, self.filterlen) + self.countMismatch(seq[-self.filterlen:], self.filter2rc, self.filterlen)
        else:
            n = self.countMismatch(seq[:self.filterlen], self.filter1, self.filterlen) + self.countMismatch(seq[-self.filterlen:], self.filter2, self.filterlen)
        return n <= self.maxmismatch

    def countMismatch(self, s1, s2, l):
        n = 0
        for i in range(l):
            if s1[i] != s2[i]:
                n += 1
        return n

    def stitch(self, fq1, fq2):
        l1 = len(fq1.seq)
        l2 = len(fq2.seq)
        seqlen = min(l1, l2)
        bestov = 0
        bestovp = 0.0
        bestidx = 0
        for i in range(self.minoverlap, seqlen):
            ov = self.findOverlap(fq1.seq, fq2.seq, i, l1-i)
            ovperc = 1.0 * ov / i
            if ovperc > bestovp:
                bestovp = ovperc
                bestov = ov
                bestidx = i
        # sys.stdout.write("Best overlap: {} {} {} l={}\n".format(bestov, bestovp, bestidx, l1+l2-bestov))
        if bestov >= self.minoverlap and bestovp > self.minident:
            fullseq = self.makeFullSequence(fq1, fq2, bestidx, l1-bestidx)
            return fullseq
        else:
            return None

    def findOverlap(self, seq1, seq2, x, start):
        n = 0
        for i in range(x):
            if seq1[start+i] == seq2[i]:
                n += 1
        return n

    def makeFullSequence(self, fq1, fq2, p, start):
        #sys.stdout.write("p={}, start={}\n".format(p, start))
        bases = []
        for i in range(p):
            q = start+i
            if fq1.seq[q] == fq2.seq[i]:
                bases.append(fq2.seq[i])
            elif fq1.qual[q] > fq2.qual[i]:
                bases.append(fq1.seq[q])
            else:
                bases.append(fq2.seq[i])
        return fq1.seq[:start] + "".join(bases) + fq2.seq[p:]

if __name__ == "__main__":
    S = Stitcher()
    if S.parseArgs(sys.argv[1:]):
        S.run()
    else:
        S.usage()

