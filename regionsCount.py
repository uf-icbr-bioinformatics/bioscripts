#!/usr/bin/env python

import sys
import pysam

import Utils
import Script

def usage():
    sys.stderr.write("""regionsCount.py - Compute coverage in specified regions.

Usage: regionsCount.py [-z] [-q qual] [-o outfile] bamfile bedfile

This program examines a BAM file `bamfile' and computes coverage in all intervals
contained in BED file `bedfile'. For each line in the BED file the output (sent to 
standard output, or to `outfile' if specified) consists of the original contents of
the line followed by two additional columns: the total coverage in the interval and
its RPKM (number of reads in the interval in millions divided by the size of the 
interval in kB).

In vector mode (enabled with -w) the original columns of the `bedfile' are followed by
the coverage values from `bamfile' for each region, in separate columns. 

Options:
 -h         | Print this usage message.
 -o outfile | Write output to `outfile' (default: standard output).
 -q qual    | Discard reads with quality score below `qual' (default: {}).
 -z         | Discard intervals with coverage of 0.
 -w W       | Vector mode, using vector of length W.

""".format(BAMreader.qual))

P = Script.Script("regionsCount.py", version="1.0", usage=usage)

B = None
bedfile = None
outfile = None

class BAMreader(object):
    mode = "normal"
    bamfile = None
    aln = None                  # Alignment object
    zeros = True
    qual = 10
    nalignments = 0
    vectsize = 0
    rpkm_factor = 1.0

    def setBamfile(self, bamfile):
        self.bamfile = bamfile
        self.aln = pysam.AlignmentFile(self.bamfile, "rb")

    def countAlignments(self):
        self.nalignments = Utils.countReadsInBAM(self.bamfile)
        return self.nalignments

    def countAlignmentsInRegion(self, chrom, start, end):
        aiter = self.aln.fetch(chrom, start, end)
        c = 0
        for a in aiter:
            if a.mapping_quality >= self.qual:
                c += 1
        return c

    def addCountsToBED(self, bedfile, out):
        self.countAlignments()
        rpkm_factor = 1000000000.0 / self.nalignments
        with open(bedfile, "r") as f:
            for parsed in Utils.CSVreader(f):
                start = int(parsed[1])
                end = int(parsed[2])
                rl = end-start
                if rl > 0:
                    c = self.countAlignmentsInRegion(parsed[0], start, end)
                    rpkm = c * rpkm_factor / rl
                    if c > 0 or self.zeros:
                        out.write("\t".join(parsed) + "\t" + str(c) + "\t" + str(rpkm) + "\n")

    ### Vector mode

    def addVectorToBED(self, bedfile, out):
        self.countAlignments()
        self.aln = pysam.AlignmentFile(self.bamfile, "rb")
        hdr = True
        res = Utils.Resampler(1, self.vectsize)

        with open(bedfile, "r") as f:
            for parsed in Utils.CSVreader(f):
                chrom   = parsed[0]
                start   = int(parsed[1])
                end     = int(parsed[2])
                strand  = parsed[3]
                regsize = end - start
                self.rpkm_factor = 1000000000.0 / (regsize * self.nalignments)
                vec     = self.getCovVector(chrom, start, end, regsize)
                res.init(regsize, self.vectsize)
                ovec    = res.resample(vec)

                if hdr:
                    out.write("Chrom\tStart\tEnd\tStrand\tGene")
                    for i in range(self.vectsize):
                        out.write("\t{}".format(i+1))
                    out.write("\n")
                    hdr = False
                if strand == '-':
                    ovec.reverse()
                out.write("\t".join(parsed))
                for v in ovec:
                    out.write("\t" + str(v))
                out.write("\n")

    def getCovVector(self, chrom, start, end, regsize):
        vector = [0]*regsize
        for pc in self.aln.pileup(chrom, start, end):
            pos = pc.pos - start
            if pos >= 0 and pos < regsize:
                vector[pos] = pc.n * self.rpkm_factor
        return vector

def parseArgs(args):
    global B
    global bedfile
    global outfile
    nra = 0
    next = ''

    P.standardOpts(args)
    B = BAMreader()
    for a in args:
        if next == '-o':
            outfile = a
            next = ''
        elif next == '-q':
            B.qual = P.toInt(a)
            next = ''
        elif next == '-w':
            B.mode = "vector"
            B.vectsize = P.toInt(a)
            next = ''
        elif a == '-z':
            B.zeros = False
        elif a in ['-o', '-q', '-w']:
            next = a
        elif nra == 0:
            B.setBamfile(P.isFile(a))
            nra += 1
        elif nra == 1:
            bedfile = P.isFile(a)
            nra += 1
        else:
            usage()
    if nra != 2:
        usage()
        return False
    return True

if __name__ == "__main__":
    if parseArgs(sys.argv[1:]):
        if B.mode == "vector":
            if outfile:
                with open(outfile, "w") as out:
                    B.addVectorToBED(bedfile, out)
            else:
                B.addVectorToBED(bedfile, sys.stdout)
        else:
            if outfile:
                with open(outfile, "w") as out:
                    B.addCountsToBED(bedfile, out)
            else:
                B.addCountsToBED(bedfile, sys.stdout)
