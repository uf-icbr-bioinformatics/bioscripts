#!/usr/bin/env python

import sys
import csv
import numpy as np
import pysam

import Script

__doc__ = """This program receives the output of a 'samtools depth' command and 
rewrites it as a matrix. Its only argument is the number of positions in each
row of the matrix. The intervals passed to samtools should all be of the same
size."""

def readRegions(filename, strandcol=None, center=False, slop=0):
    """Read regions from tab-delimited file `filename'. Returns a list of
tuples containing chromosome, start, end. If `strandcol' is specified it
indicates the column containing the strand, which will be added as the 
fourth element of each tuple. If 'center' is True, find the midpoint of 
region from the input file and add `slop' before and after."""
    regions = []
    with open(filename, "r") as f:
        c = csv.reader(f, delimiter="\t")
        for line in c:
            if line[0][0] == '#':
                continue
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])

            if center:
                mid = start + end / 2
                start = start - slop
                end = end + slop

            if strandcol:
                regions.append((chrom, start, end, line[strandcol]))
            else:
                regions.append((chrom, start, end))

    return regions

def buildMatrix(bamfile, regions, npos):
    a = np.zeros((len(regions), npos))
    row = 0
    aln = pysam.AlignmentFile(bamfile, "rb")
    for reg in regions:
        m = 0
        for p in aln.pileup(reg[0], reg[1], reg[2]):
            col = p.pos - reg[1]
            if col < npos:
                a[row, col] = p.n
                if p.n > m:
                    m = p.n
        row += 1
    return a

def normGreenleaf(a, leftmargin=100, rightmargin=100):
    means = a.mean(axis=0)
    if leftmargin:
        leftnoise = sum(a[:leftmargin])
    if rightmargin:
        rightnoise = sum(a[-rightmargin:])
    if leftmargin and rightmargin:
        noise = 1.0 * (leftnoise + rightnoise) / (leftmargin + rightmargin)
    elif leftmargin:
        noise = 1.0 * leftnoise / leftmargin
    elif rightmargin:
        noise = 1.0 * rightnoise / rightmargin
    else:
        raise Error("At least one of leftmargin and rightmargin should be specified.")
    return a / noise

class RegionsToMatrix(Script.Command):
    _cmd = "matrix"
    regions = []
    bedfile = None
    bamfile = None
    outfile = None
    width = None
    center = False
    slop = 0

    def __init__(self):
        self.regions = []

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-c":    # take center of each region and add specified bps before and after
                self.center = True
                self.slop = int(a)
                self.width = self.slop * 2
                prev = ""
            elif prev == "-w":  # width of regions (unless specified with -c)
                self.width = int(a)
                prev = ""
            elif a in ["-c", "-w"]:
                prev = a
            elif self.bedfile is None:
                self.bedfile = a
            elif self.bamfile is None:
                self.bamfile = a
            elif self.outfile is None:
                self.outfile = a
        if self.width is None:
            P.errmsg(P.NOWIDTH)
        if self.bedfile and self.bamfile and self.outfile:
            return True
        else:
            P.errmsg(P.NOFILES)
    
    def run(self, P, args):
        self.parseArgs(args)
        self.regions = readRegions(self.bedfile, center=self.center, slop=self.slop)
        nregs = len(self.regions)
        sys.stderr.write("Read {} regions from {}\n".format(nregs, self.bedfile))
        sys.stderr.write("Building matrix ({}x{})...".format(nregs, self.width))
        a = buildMatrix(self.bamfile, self.regions, self.width)
        sys.stderr.write(" done.\n")
        sys.stderr.write("Applying Greenleaf normalization...")
        a = normGreenleaf(a)
        sys.stderr.write(" done.\n")
        sys.stderr.write("Saving matrix to {}.npy...".format(self.outfile))
        np.save(self.outfile, a)
        sys.stderr.write(" done.\n")

class PileupToMatrix(Script.Script):

    def main(self, args):
        self.standardOpts(args)
        cmd = args[0]
        c = self.findCommand(cmd)
        if c:
            c().run(self, args[1:])
        else:
            self.errmsg(self.NOCMD)
            
#        regions = readRegions(bedfile)
    #    print regions[:10]
#        sys.stderr.write("Read {} regions from {}\n".format(len(regions), bedfile))
#        a = buildMatrix(bamfile, regions, npos)
#        sys.stderr.write("Matrix created.\n")
#        np.save(matfile, a)
#        sys.stderr.write("Matrix written to {}\n".format(matfile))


P = PileupToMatrix("pileupToMatrix.py", version="1.0",
                   errors=[('NOCMD', 'Missing command', "Please specify a command."),
                           ('NOWIDTH', 'Missing width', "Please specify the width using either the -c or -w options."),
                           ('NOFILES', 'Missing input files', "Please specify a BED file and a BAM file.")])
P.addCommand(RegionsToMatrix)

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        P.main(args)
    else:
        P.errmsg(P.NOCMD)

