#!/usr/bin/env python

import sys
import os.path
import Utils
import Script
import SeqUtils

__doc__ = """Operate on barcodes in fastq files"""

def usage(what=None):
    sys.stdout.write("""demux.py - {}

Usage: demux.py detect [options] fastq
       demux.py split  [options] fastq1 [fastq2]
       demux.py grep   [options] fastq1 [fastq2]
       demux.py distr  [options] fastq1 [fastq2]
       demux.py distr  [options] fasta

The `detect' command examines a fastq or fasta file extracting its barcodes. By default,
the input file is assumed to be in fastq format and the barcodes are extracted from 
the final part of the read header (e.g., 1:N:0:TACAGC). If the -s option is specified, 
barcodes are instead extracted from the read sequence, extracting its first S bases if
S is positive, or its last S bases if negative. You can also use a syntax similar to Python's 
slice notation to specify an arbitrary range of bases, e.g. if barcodes are in positions 
5-10 you should use -s 5:10. Note that in this case positions are 1-based, and the 
final position is included.

The number of lines examined can be controlled with the -n option. Only barcodes with 
a frequency greater than the one specified with the -p option are reported. Detected 
barcodes are written to standard output in tab-delimited format. The output is suitable 
as the barcodes file for the `split' command.

The `split' command separates the reads from a fasta or fastq file (or two fastq files in 
paired-end mode) into one file for each barcode, plus (optionally) a file for reads with 
unmatched barcodes. The barcodes file should be tab-delimited, with labels in the first 
column and barcode sequences in the second one. Use the format AAAA+BBBB for dual indexes.
Barcodes are read from the read header (assumed to be in Illumina format) by default, or 
from the sequence itself if -s is specified. If the input file is in fasta format, -s is required.

The `grep' command extracts read pairs containing at least one occurrence of a specified
pattern, in either the left or right mate of each pair. Use -lt to specify the pattern
to be searched for in the left mate, and -rt for the right mate (if not specified, defaults
to reverse-complement of -lt).

The `distr' command distributes the reads in the input file(s) into D different files, where
D is specified with the -d option. Reads are distributed in round-robin fashion, so that
at the end each output file will contain approximately the same number of reads. Input can
be in fastq or fasta format.

Options:

  -b FILE | File containing barcode sequences
  -r      | Barcodes are reverse-complemented
  -u      | Do not write sequences with unknown indexes to UND-...
  -s S    | Position of barcode in sequences, slice syntax allowed (e.g. 3:9).
  -m N    | Allow at most N mismatches in the barcode (default: {})
  -n N    | Number of reads to examine for barcode autodetection (default: {})
  -p N    | Do not show detected barcodes occurring in less than N% of reads (default: {})
  -lt P   | Search for pattern P in left reads.
  -rt P   | Search for pattern P in right reads.
  -d D    | Number of files to distribute reads to (default: {})
  -o O    | Prefix for files for distributed reads.
  -1      | In detect, if dual indexing is used (e.g. AAA+BBB), only return the first index.

""".format(__doc__, Demux.maxmismatch, Demux.ndetect, Demux.minpct, Demux.distr))

### Program object

class Demux(Script.Script):
    mode = None
    bcfile = ""
    fqleft = ""
    fqright = ""
    revcomp = False             # Reverse-complement barcodes?
    undet = True                # Write unclassified reads?
    maxmismatch = 1             # Maximum number of mismatches
    nf = 0                      # Number of input files
    ndetect = 10000             # Number of reads for barcode detection
    minpct = 1
    bcstart = 0
    bcend = None
    bcslice = None
    removePlus = False          # Only return first component of dual barcodes (e.g. AACC+GGTT)
    leftTarget = None
    rightTarget = None
    distr = 1
    distrout = None
    odir = ""                   # Output directory for split
    
    def parseArgs(self, args):
        self.nf = 0
        self.standardOpts(args)
        cmd = args[0]
        if cmd in ['split', 'detect', 'grep', 'distr']:
            self.mode = cmd
        else:
            P.errmsg(P.NOCMD)

        prev = ""
        for a in args[1:]:
            if prev == '-m':
                self.maxmismatch = self.toInt(a)
                prev = ""
            elif prev == '-b':
                self.bcfile = self.isFile(a)
                prev = ""
            elif prev == '-n':
                self.ndetect = self.toInt(a)
                prev = ""
            elif prev == '-p':
                self.minpct = self.toFloat(a)
                prev = ""
            elif prev == '-s':
                self.bcslice = Utils.parseSlice(a)
                prev = ""
            elif prev == '-lt':
                self.leftTarget = a
                prev = ""
            elif prev == '-rt':
                self.rightTarget = a
                prev = ""
            elif prev == '-d':
                self.distr = self.toInt(a)
                prev = ""
            elif prev == '-o':
                self.distrout = a
                prev = ""
            elif prev == "-D":
                self.odir = a + "/"
                prev = ""
            elif a in ['-m', '-b', '-n', '-p', '-s', '-lt', '-rt', '-o', '-d', "-D"]:
                prev = a
            elif a == '-r':
                self.revcomp = True
            elif a == '-u':
                self.undet = False
            elif a == "-1":
                self.removePlus = True
            elif self.nf == 0:
                self.fqleft = "-" if a == "-" else self.isFile(a)
                self.nf += 1
            elif self.nf == 1:
                self.fqright = self.isFile(a)
                self.nf += 1

P = Demux("demux", version="1.0", usage=usage,
          errors=[('NOCMD', 'Missing command', 'The first argument should be one of: split, detect, grep, distr.'),
                  ('BADFMT', 'Bad input file format', 'The input file should be in fasta of fastq format.'),
                  ('NOSLICE', 'Missing -s option', 'Demultiplexing FASTA files requires the -s options.') ])

### Classes

class Barcode():
    seq = ""
    name = ""
    filename = ""
    stream = None
    filename2 = ""              # For paired-end
    stream2 = None
    nhits = 0
    
    def __init__(self, name, seq):
        self.seq = seq
        self.name = name
        self.nhits = 0

    def openStream(self, filename, filename2=None, ext=".fastq.gz", odir=""):
        #print(filename, filename2, ext, odir)
        self.filename = odir + self.name + "-" + filename + ext
        self.stream = Utils.genOpen(self.filename, "wt")
        if filename2:
            self.filename2 = odir + self.name + "-" + filename2 + ext
            self.stream2 = Utils.genOpen(self.filename2, "wt")

    def closeStream(self):
        self.stream.close()
        if self.stream2:
            self.stream2.close()

    def writeRecord(self, fq1, fq2=None):
        self.nhits += 1
        if fq1.qual:
            self.stream.write(fq1.name + "\n")
            self.stream.write(fq1.seq + "\n")
            self.stream.write("+\n")
            self.stream.write(fq1.qual + "\n")
        else:                   # FASTA
            self.stream.write(">" + fq1.name + "\n")
            self.stream.write(fq1.seq + "\n")
        if fq2:
            self.stream2.write(fq2.name + "\n")
            self.stream2.write(fq2.seq + "\n")
            self.stream2.write("+\n")
            self.stream2.write(fq2.qual + "\n")

class BarcodeMgr():
    barcodeseqs = []
    barcodes = {}
    nbarcodes = 0
    nhits = 0

    def __init__(self):
        self.barcodeseqs = []
        self.barcodes = {}
        self.nbarcodes = 0
        self.nhits = 0

    def initFromFile(self, filename, rc=False, undet=False):
        with open(filename, "rt") as f:
            for line in f:
                if line[0] == '#':
                    continue
                line = line.rstrip("\r\n").split("\t")
                if len(line) < 2:
                    continue
                if line[0] == 'Other': # when reading output of `detect' command...
                    continue
                if rc:
                    bc = SeqUtils.revcomp(line[1])
                else:
                    bc = line[1]
                self.add(line[0], bc)
        if undet:
            self.add("UND", "*")

    def initForDistribute(self, n):
        for i in range(1, n+1):
            si = str(i)
            self.add(si, si)

    def add(self, name, seq):
        b = Barcode(name, seq)
        self.barcodeseqs.append(seq)
        self.barcodes[seq] = b

    def openAll(self, filename, filename2=None, ext=".fastq.gz", odir=""):
        #print(filename, filename2, ext, odir)
        for b in self.barcodes.values():
            b.openStream(filename, filename2, ext=ext, odir=odir)

    def allFilenames(self):
        result = []
        for b in self.barcodes.values():
            result.append(b.filename)
            if b.filename2:
                result.append(b.filename2)
        return result

    def closeAll(self):
        try:
            for b in self.barcodes.itervalues():
                b.closeStream()
        except AttributeError:
            for b in self.barcodes.values():
                b.closeStream()

    def findBest(self, seq, maxmismatch=1):
        maxd = len(seq)
        best = ""
        for bc in self.barcodeseqs:
            if bc != "*":
                d = SeqUtils.distance(seq, bc)
                if d < maxd:
                    maxd = d
                    best = bc
        # print(maxd)
        if maxd <= maxmismatch:
            bb = self.barcodes[best]
            self.nhits += 1
            #bb.nhits += 1
            return bb
        if "*" in self.barcodes:
            self.nhits += 1
            return self.barcodes["*"]
        # else:
        return None

    def findBestPaired(self, seq1, seq2, maxmismatch=1):
        b1 = self.findBest(seq1, maxmismatch=maxmismatch)
        b2 = self.findBest(seq2, maxmismatch=maxmismatch)
        k1 = 0
        k2 = 4
        if b1:
            if b1.seq == "*":
                k1 = 1
            else:
                k1 = 2
        if b2:
            if b2.seq == "*":
                k2 = 8
            else:
                k2 = 16
        # print(b1.seq, b2.seq)
        # print(k1, k2)
        x = k1 + k2
        if x == 18:
            if b1.seq == b2.seq:
                return b1
            else:
                return None
        elif x in [10, 18]:
            return b1
        elif x in [16, 17]:
            return b2
        else:
            return None
        
    def findOrCreate(self, seq):
        self.nhits += 1
        if seq not in self.barcodeseqs:
            self.nbarcodes += 1
            self.add("IDX{}".format(self.nbarcodes), seq)
        self.barcodes[seq].nhits += 1
        return self.barcodes[seq]

    def showCounts(self):
        hiddenb = 0             # Barcodes not shown
        hiddenh = 0             # Hits for barcodes not shown
        ranking = list(self.barcodes.values())
        ranking.sort(key=lambda b:b.nhits, reverse=True)

        try:
            sys.stdout.write("#Name\tSeq\tHits\tPct\tFile1\tFile2\n")
            for b in ranking:
                seq = b.seq
                if P.revcomp:
                    seq = SeqUtils.revcomp(seq)
                if self.nhits == 0:
                    pct = 0
                else:
                    pct = 100.0 * b.nhits / self.nhits
                if pct >= P.minpct:
                    sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(b.name, seq, b.nhits, Utils.f2dd(pct), b.filename, b.filename2 or ""))
                else:
                    hiddenb += 1
                    hiddenh += b.nhits
            if hiddenb > 0 and self.nhits > 0:
                sys.stdout.write("Other\t({})\t{}\t{}\t\t\n".format(hiddenb, hiddenh, Utils.f2dd(100.0 * hiddenh / self.nhits)))
        except IOError:
            pass

class FastqReader(SeqUtils.FastqReader):

    def demux(self, dm, filename, bcslice=None, odir=""):
        dm.openAll(filename, odir=odir)
        while True:
            self.nextRead()
            if not self.stream:
                break
            bc = self.fq.getBarcode(bcslice, removePlus=P.removePlus)
            bo = dm.findBest(bc, maxmismatch=P.maxmismatch) # Barcode object
            # print("best: " + bo.seq)
            # raw_input()
            if bo:
                self.ngood += 1
                bo.writeRecord(self.fq)
        dm.closeAll()
        sys.stderr.write("Total reads: {}\n.".format(self.nread))
        sys.stderr.write("Written: {}\n".format(self.ngood))

    def detect(self, ndetect, bcslice=None):
        """Detect the barcodes contained in the first `ndetect' reads of this fastq file."""
        dm = BarcodeMgr()
        if P.bcfile:
            dm.initFromFile(P.bcfile, rc=P.revcomp)
        while True:
            self.nextRead()
            if not self.stream:
                break
            bc = self.fq.getBarcode(bcslice, removePlus=P.removePlus)
            dm.findOrCreate(bc)
            if self.nread == ndetect:
                break
        try:
            sys.stdout.write("# Reads examined: {}\n".format(dm.nhits))
            dm.showCounts()
        except IOError:
            pass

class FastaReader(SeqUtils.FastaReader):

    def demux(self, dm, filename, bcslice=None, odir=""):
        dm.openAll(filename, ext=".fasta", odir=odir)
        while True:
            self.nextRead()
            bc = self.fq.getBarcode(bcslice, removePlus=P.removePlus)
            bo = dm.findBest(bc, maxmismatch=P.maxmismatch) # Barcode object
            # print("best: " + bo.seq)
            # raw_input()
            if bo:
                self.ngood += 1
                bo.writeRecord(self.fq)
            if not self.stream:
                break
        dm.closeAll()
        sys.stderr.write("Total reads: {}\n".format(self.nread))
        sys.stderr.write("Written: {}\n".format(self.ngood))

class PairedFastqReader(SeqUtils.PairedFastqReader):

    def demux(self, dm, filename1, filename2, bcslice=None, odir=""):
        #print(odir)
        dm.openAll(filename1, filename2=filename2, odir=odir)

        if bcslice:
            while True:
                self.nextRead()
                if not self.reader1.stream:
                    break
                bc1 = self.fq1.getBarcode(bcslice, removePlus=P.removePlus)
                bc2 = self.fq2.getBarcode(bcslice, removePlus=P.removePlus)
                # print(bc1, bc2)
                bo = dm.findBestPaired(bc1, bc2, maxmismatch=P.maxmismatch)
                # print(bo)
                # print()
                if bo:
                    self.ngood += 1
                    bo.writeRecord(self.fq1, self.fq2)

        else:
            
            while True:
                self.nextRead()
                if not self.reader1.stream:
                    break
                bc1 = self.fq1.getBarcode(None, removePlus=P.removePlus)
                bc2 = self.fq2.getBarcode(None, removePlus=P.removePlus)
                if bc1 != bc2:
                    self.nbad += 1
                    continue
                # print("barcode: " + bc1)
                # print("best: " + bo.seq)
                # raw_input()
                bo = dm.findBest(bc1, maxmismatch=P.maxmismatch)
                #print(bo)
                # print("best: " + bo.seq)
                # raw_input()
                if bo:
                    self.ngood += 1
                    bo.writeRecord(self.fq1, self.fq2)

        dm.closeAll()
        sys.stderr.write("Total reads: {}\n".format(self.nread))
        sys.stderr.write("Mismatched barcodes: {}\n".format(self.nbad))
        sys.stderr.write("Written: {}\n".format(self.ngood))

    def readMatches(self, patt1, patt2):
        f1 = self.fq1.seq.find(patt1)
        f2 = self.fq2.seq.find(patt2)
        return f1 > -1 or f2 > -1

    def grep(self, dm, filename1, filename2, patt1, patt2):
        ngood = 0
        nbad = 0
        dm.openAll(filename1, filename2)
        while True:
            self.nextRead()
            if not self.reader1.stream:
                break
            if self.readMatches(patt1, patt2):
                bc = dm.findOrCreate("MATCH")
                ngood += 1
            else:
                bc = dm.findOrCreate("MISMATCH")
                nbad += 1
            bc.writeRecord(self.fq1, self.fq2)
        dm.closeAll()
        sys.stderr.write("Total read pairs: {}\n".format(ngood+nbad))
        sys.stderr.write("Matching read pairs: {}\n".format(ngood))
        sys.stderr.write("Not matching read pairs: {}\n".format(nbad))

### Distribute

def doDistribute():
    dm = BarcodeMgr()
    dm.initForDistribute(P.distr)
    nameleft = getFastqBasename(P.fqleft)
    if P.fqright:
        nameright = getFastqBasename(P.fqright)
    else:
        nameright = None
    try:
        dm.openAll(nameleft, nameright)
        pfr = PairedFastqReader(P.fqleft, P.fqright)
        i = 0
        while True:
            pfr.nextRead()
            if not pfr.reader1.stream:
                break
            label = dm.barcodeseqs[i]
            bc = dm.barcodes[label]
            #print("writing {} to {}, {}".format(i, label, bc))
            #raw_input()
            bc.writeRecord(pfr.fq1, pfr.fq2)
            i += 1
            if i == P.distr:
                i = 0
    finally:
        dm.closeAll()

def doDistribute2():
    outfiles1 = []
    outfiles2 = []
    outstreams1 = []
    outstreams2 = []

    if not P.distrout:
        P.distrout = "out"

    sys.stderr.write("Distributing reads to:\n")
    if P.nf == 2:

        for i in range(1, P.distr + 1):
            outfiles1.append("{}.{}.R1.fastq.gz".format(P.distrout, i))
            outfiles2.append("{}.{}.R2.fastq.gz".format(P.distrout, i))

        try:
            for i in range(P.distr):
                sys.stderr.write("  {}, {}\n".format(outfiles1[i], outfiles2[i]))
                outstreams1.append(Utils.genOpen(outfiles1[i], "wt"))
                outstreams2.append(Utils.genOpen(outfiles2[i], "wt"))

            in1 = Utils.genOpen(P.fqleft, "rt")
            in2 = Utils.genOpen(P.fqright, "rt")
            i = 0
            while True:
                r = in1.readline()
                if not r:
                    break
                o1 = outstreams1[i]
                o2 = outstreams2[i]
                o1.write(r)
                o1.write(in1.readline())
                o1.write(in1.readline())
                o1.write(in1.readline())
                o2.write(in2.readline())
                o2.write(in2.readline())
                o2.write(in2.readline())
                o2.write(in2.readline())
                i += 1
                if i == P.distr:
                    i = 0
        finally:
            in1.close()
            in2.close()
            for s in outstreams1:
                s.close()
            for s in outstreams2:
                s.close()
    else:                       # single fastq file
        for i in range(1, P.distr + 1):
            outfiles1.append("{}.{}.fastq.gz".format(P.distrout, i))
        try:
            for i in range(P.distr):
                sys.stderr.write("  {}\n".format(outfiles1[i]))
                outstreams1.append(Utils.genOpen(outfiles1[i], "wt"))

            in1 = Utils.genOpen(P.fqleft, "rt")
            i = 0
            while True:
                r = in1.readline()
                if not r:
                    break
                o1 = outstreams1[i]
                o1.write(r)
                o1.write(in1.readline())
                o1.write(in1.readline())
                o1.write(in1.readline())
                i += 1
                if i == P.distr:
                    i = 0
        finally:
            in1.close()
            for s in outstreams1:
                s.close()

def doDistributeFasta():
    outfiles = []
    outstreams = []

    if not P.distrout:
        P.distrout = "out"

    sys.stderr.write("Distributing sequences to {} FASTA files.\n".format(P.distr))
    for i in range(1, P.distr + 1):
        name = "{}.{}.fasta".format(P.distrout, i)
        outfiles.append(name)
        outstreams.append(Utils.genOpen(name, "wt"))
    try:
        FR = FastaReader(P.fqleft)
        i = 0
        while True:
            r = FR.nextRead()
            if not r:
                break
            o = outstreams[i]
            o.write(">{}\n{}\n".format(FR.fq.name, FR.fq.seq))
            i += 1
            if i == P.distr:
                i = 0
    finally:
        for o in outstreams:
            o.close()
    for o in outfiles:
        sys.stdout.write("{}\n".format(o))

### Main

def getFastqBasename(filename):
    fn = os.path.split(filename)[1]
    sp = os.path.splitext(fn)
    if sp[1] == '.gz':
        sp = os.path.splitext(sp[0])
    return sp[0]

def main(args):
    P.parseArgs(args)

    if P.mode == 'split':

        bm = BarcodeMgr()
        bm.initFromFile(P.bcfile, rc=P.revcomp, undet=P.undet)

        fmt = Utils.detectFileFormat(P.fqleft)
        if fmt == "fastq":
            nameleft = getFastqBasename(P.fqleft)
            if P.fqright:
                nameright = getFastqBasename(P.fqright)
            else:
                nameright = None

            if P.nf == 2:
                fr = PairedFastqReader(P.fqleft, P.fqright)
                #print(P.odir)
                fr.demux(bm, nameleft, nameright, P.bcslice, odir=P.odir)
            else:
                fr = FastqReader(P.fqleft)
                fr.demux(bm, nameleft, P.bcslice, odir=P.odir)
            bm.showCounts()
        elif fmt == "fasta":
            if not P.bcslice:
                P.errmsg(P.NOSLICE)
            fr = FastaReader(P.fqleft)
            fr.demux(bm, getFastqBasename(P.fqleft), P.bcslice, odir=P.odir)
            bm.showCounts()
        else:
            P.errmsg(P.BADFMT)

    elif P.mode == 'detect':
        fr = FastqReader(P.fqleft)
        fr.detect(P.ndetect, P.bcslice)

    elif P.mode == 'grep':
        if P.rightTarget is None:
            P.rightTarget = SeqUtils.revcomp(P.leftTarget)
        nameleft = getFastqBasename(P.fqleft)
        if P.fqright:
            nameright = getFastqBasename(P.fqright)
        else:
            nameright = None
        bm = BarcodeMgr()
        bm.add("MATCH", "MATCH")
        bm.add("MISMATCH", "MISMATCH")

        fr = PairedFastqReader(P.fqleft, P.fqright)
        fr.grep(bm, nameleft, nameright, P.leftTarget, P.rightTarget)
    
    elif P.mode == 'distr':
        fmt = Utils.detectFileFormat(P.fqleft)
        if fmt == "fastq":
            doDistribute2()
        elif fmt == "fasta":
            doDistributeFasta()
        else:
            sys.stderr.write("File {} is in an unknown format.\n".format(P.fqleft))
        
if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        main(args)
    else:
        P.usage()

