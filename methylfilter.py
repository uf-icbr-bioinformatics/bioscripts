#!/usr/bin/env python

## (c) 2014-2017, Alberto Riva (ariva@ufl.edu)
## DiBiG, ICBR Bioinformatics, University of Florida

import sys
from Bio import SeqIO

import Script

# Script class

def usage():
    sys.stderr.write("""methylfilter.py - Separate sequences by average methylation.

Usage: methylfilter.py input.fa [options] outdesc ...

The input file should be a fasta file in which the first sequence 
is the reference. All other sequences are compared to the reference
to determine % methylation.

Each `outdesc' argument should be a string of the form:

  min-max:filename

where min and max should be expressed as percentages, ie integer numbers
in the range 0..100. If min is omitted it defaults to 0. If max is omitted
it defaults to 100. Each outdesc specifies that sequences with methylation
values between min (inclusive) and max (exclusive) should be written to 
filename. Any number of outdesc arguments can be used.

Options:

 -h, --help             | Write this usage message.
 -v, --version          | Print version number.
 -a, --addref           | Write reference sequence at the beginning of each 
                          output file.
 -r, --report filename  | Write to `filename' a report showing, for each input
                          sequence, its % methylation and the output file it
                          was written to.
 -s, --summary filename | Write to `filename' a summary showing the number of
                          sequences written to each output file.
 -gcg                   | Do not exclude GCG sites from analysis.
 -gc                    | Output is based on GC methylation instead of CG.

""")

P = Script.Script("methylfilter.py", version="1.0", usage=usage, 
                  errors=[('BADRANGE', 'Bad range specification', "Cannot parse argument `{}'. Format should be: low-high:filename.")])

# Utility classes

class refDesc():
    """A class containing the reference sequence, its length, and a list of CG and GC positions."""
    sequence = None
    length = 0
    CGpositions = []
    GCpositions = []
    numCGs = 0
    numGCs = 0
    excludeGCG = True

    def __init__(self, ref, excludeGCG):
        length = len(ref)
        self.sequence = ref
        self.length = length
        self.excludeGCG = excludeGCG
        self.CGpositions = detectCG(ref, length, excludeGCG=self.excludeGCG)
        self.GCpositions = detectGC(ref, length, excludeGCG=self.excludeGCG)
        self.numCGs = len(self.CGpositions)
        self.numGCs = len(self.GCpositions)

class outFile():
    """A class that writes sequences whose methylation rate is in a specified range to a file."""
    mrmin = 0
    mrmax = 0
    filename = None
    stream = None
    nout = 0

    def __init__(self, filename, min, max):
        self.min = min
        self.max = max
        self.filename = filename
        self.nout = 0

    def open(self):
        self.stream = open(self.filename, "w")

    def close(self):
        self.stream.close()

    def writeSeq(self, seq, count=True):
        SeqIO.write(seq, self.stream, "fasta")
        if count:
            self.nout = self.nout + 1

class mfrun():
    """A class representing the whole run. Includes the refDesc object and the list of output file objects."""
    rd = None
    infile = None
    outfiles = []
    writeRef = False            # If true, writes the reference sequence at the beginning of each output file
    reportFile = False          # Name of report file, if desired
    reportStream = None
    summaryFile = False         # Name of summary file, if desired
    summaryStream = None
    excludeGCG = True           # 
    mode = "CG"

    def parseArgs(self, args):
        """Parse command-line arguments creating outfiles."""

        P.standardOpts(args)
        prev = False
        for arg in args:
            if prev == "r":
                self.reportFile = arg
                prev = False
            elif prev == "s":
                self.summaryFile = arg
                prev = False
            elif (arg == "--addref") or (arg == "-a"):
                self.writeRef = True
            elif (arg == "--report") or (arg == "-r"):
                prev = "r"
            elif (arg == "--summary") or (arg == "-s"):
                prev = "s"
            elif (arg == "-gcg"):
                self.excludeGCG = False
            elif arg == "-gc":
                self.mode = "GC"
            elif self.infile == None:
                self.infile = P.isFile(arg)
            else:
                pdash = arg.find('-')
                pcolon = arg.find(':', pdash)
                if (pdash == -1) or (pcolon == -1):
                    P.errmsg(P.BADRANGE, arg)
                else:
                    low = arg[0:pdash]
                    high = arg[pdash+1:pcolon]
                    fn = arg[pcolon+1:]
                    if low == "":
                        low = 0
                    else:
                        low = int(low)
                    if high == "":
                        high = 100
                    else:
                        high = int(high)
                    if high == 100: high = 101 # hack so that we can use < in comparisons and still catch 100%
                    self.outfiles.append(outFile(fn, low, high))

    def openAll(self):
        """Open all necessary output files."""
        for of in self.outfiles:
            of.open()
            if self.writeRef:
                of.writeSeq(self.rd.sequence, False) # don't count ref seq in number of written sequences
        if self.reportFile:
            print "Writing report to {}".format(self.reportFile)
            self.reportStream = open(self.reportFile, "w")
            if self.mode == "CG":
                self.reportStream.write("Sequence\t% CG Meth\tConv\tTot\t% GC Meth\tFilename\n")
            elif self.mode == "GC":
                self.reportStream.write("Sequence\t% GC Meth\tConv\tTot\t% CG Meth\tFilename\n")

    def closeAll(self):
        """Close all output files."""
        for of in self.outfiles:
            of.close()
        if self.reportStream:
            self.reportStream.close()
    
    def showOutfiles(self):
        """Show all defined output files with their ranges."""
        for of in self.outfiles:
            print "{}% - {}% -> {}".format(of.min, of.max, of.filename)

    def findOutfile(self, x):
        """Find the output file for value `x'."""
        for of in self.outfiles:
            if (x >= of.min) and (x < of.max):
                return of
        return None

    def showSummary(self):
        tot = 0
        for of in self.outfiles:
            tot = tot + of.nout
            print "{:5}  {}".format(of.nout, of.filename)
        if self.summaryFile:
            print "Writing summary to {}".format(self.summaryFile)
            with open(self.summaryFile, "w") as out:
                out.write("Filename\tMin\tMax\tNseqs\n")
                for of in self.outfiles:
                    out.write("{}\t{}\t{}\t{}\n".format(of.filename, of.min, of.max, of.nout))
        return tot
    
    def report(self, s, p, c, t, p2, o):
        if self.reportStream:
            if o:
                self.reportStream.write("{}\t{:.2f}\t{}\t{}\t{:.2f}\t{}\n".format(s.id, p, c, t, p2, o.filename))
            else:
                self.reportStream.write("{}\t{:.2f}\t{}\t{}\t{:.2f}\t{}\n".format(s.id, p, c, t, p2, "-"))

### General

def loadSequences(filename):
    return SeqIO.parse(filename, "fasta")

def detectCG(seq, length, excludeGCG=True):
    """Returns the list of C positions in CG dinucleotides in sequence `seq'.
If `excludeGCG' is True, ignores GCG positions."""
    result = []
    candidate = False

    for i in range(length):
        if (seq[i] == 'C'):
            candidate = i
        elif (seq[i] == 'G') and candidate:
            if excludeGCG:
                if (i < 2) or seq[i-2] != 'G':
                    result.append(candidate)
                candidate = False
            else:
                result.append(candidate)
                candidate = False
        else:
            candidate = False
    return result

def detectGC(seq, length, excludeGCG=True):
    """Returns the list of C positions in GC dinucleotides in sequence `seq'.
If `excludeGCG' is True, ignores GCG positions."""
    result = []
    candidate = False

    for i in range(1, length):
        if (seq[i] == 'C') and (seq[i-1] == 'G'):
            # This is a GC position. Now check position i+1 if excluding GCGs
            if excludeGCG:
                if (i == length-1) or seq[i+1] != 'G':
                    result.append(i)
            # Otherwise, simply add the position
            else:
                result.append(i)
    return result

def countCGconverted(rd, seq):
    return countConverted(seq, rd.CGpositions)

def countGCconverted(rd, seq):
    return countConverted(seq, rd.GCpositions)

def countConverted(seq, positions):
    """Check for C->T conversion in `seq' at the positions in the list `positions'. Returns a tuple containing the number of converted Cs
and the total number of Cs in the list `positions' actually present in `seq' (ie, skipping the - positions)."""
    tot = 0
    cnt = 0
    for cgpos in positions:
        if seq[cgpos] != '-':
            tot = tot + 1
            if seq[cgpos] == 'T':
                cnt = cnt + 1
    return (cnt, tot)

### Main

def main():
    run = mfrun()
    run.parseArgs(sys.argv[1:])

    if run.infile == None:
        sys.stderr.write("""No reference file specified. Please use the -h option for usage.\n""")
        sys.exit(-3)
    if len(run.outfiles) == 0:
        sys.stderr.write("""No output files specified. Please use the -h option for usage.\n""")
        sys.exit(-4)

    seqs = loadSequences(run.infile)
    rd = refDesc(seqs.next(), run.excludeGCG)   # reference sequence
    run.rd = rd

    print "Reference sequence loaded from file `{}'.".format(run.infile)
    if run.excludeGCG:
        print "Excluding GCG positions."
    else:
        print "Not excluding GCG positions."
    if run.mode == "CG":
        print "{}bp, {} CG positions.".format(rd.length, rd.numCGs)
    elif run.mode == "GC":
        print "{}bp, {} GC positions.".format(rd.length, rd.numGCs)

    try:
        run.openAll()
        print "{} output files opened:".format(len(run.outfiles))
        run.showOutfiles()
        print "Parsing sequences..."
        
        nread = 0
        for s in seqs:
            nread = nread + 1
            (CGcnt, CGtot) = countCGconverted(rd, s)
            (GCcnt, GCtot) = countGCconverted(rd, s)
            if GCtot > 0:
                CGp = 100 * (1.0 - (CGcnt * 1.0 / CGtot)) # ensure we work with floats
            else:
                CGp = 0.0
            if GCtot > 0:
                GCp = 100 * (1.0 - (GCcnt * 1.0 / GCtot))
            else:
                GCp = 0.0
            if run.mode == "CG":
                o = run.findOutfile(CGp)
                if o:
                    o.writeSeq(s)
                    run.report(s, CGp, CGcnt, CGtot, GCp, o)
                else:
                    run.report(s, CGp, CGcnt, CGtot, GCp, None)
            elif run.mode == "GC":
                o = run.findOutfile(GCp)
                if o:
                    o.writeSeq(s)
                    run.report(s, GCp, GCcnt, GCtot, CGp, o)
                else:
                    run.report(s, GCp, GCcnt, GCtot, CGp, None)

    finally:
        run.closeAll()

    print "Done. {} sequences read.".format(nread)
    print "Report:"
    nwritten = run.showSummary()
    print "{} sequences written.".format(nwritten)

if __name__ == "__main__":
    
    if (len(sys.argv) == 1):
        usage()
    else:
        main()

