#!/usr/bin/env python

import sys
import os.path
import Utils
import Script

def usage():
    sys.stderr.write("""demux.py - Operate on barcodes in fastq files

Usage: demux.py detect [options] fastq
       demux.py split [options] fastq1 [fastq2]
       demux.py grep [options] fastq1 [fastq2]

The `detect' command examines a fastq or fasta file extracting its barcodes. By default,
the input file is assumed to be in fastq format and the barcodes are extracted from 
the final part of the read header (e.g., 1:N:0:TACAGC). If the -s option is specified, 
barcodes are instead extracted from the read sequence, extracting its first S bases if
S is positive, or its last S bases if negative. The number of lines examined can be 
controlled with the -n option. Only barcodes with a frequency greater than the one 
specified with the -p option are reported. Detected barcodes are written to standard
output in tab-delimited format. The output is suitable as the barcodes file for the 
`split' command.

The `split' command separates the reads from a fastq file (or two fastq files in 
paired-end mode) into one file for each barcode, plus (optionally) a file for reads 
with unmatched barcodes. The barcodes file should be tab-delimited, with labels in 
the first column and barcode sequences in the second one. 

The `grep' command extracts read pairs containing at least one occurrence of a specified
pattern, in either the left or right mate of each pair. Use -lt to specify the pattern
to be searched for in the left mate, and -rt for the right mate (if not specified, defaults
to reverse-complement of -lt).

Options:

  -b FILE | File containing barcode sequences
  -r      | Barcodes are reverse-complemented
  -u      | Do not write sequences with unknown indexes to UND-...
  -s S    | Number of barcode bases at the beginning (or end if negative) of read.
  -m N    | Allow at most N mismatches in the barcode (default: {})
  -n N    | Number of reads to examine for barcode autodetection (default: {})
  -p N    | Do not show detected barcodes occurring in less than N% of reads (default: {})
  -lt P   | Search for pattern P in left reads.
  -rt P   | Search for pattern P in right reads.

""".format(Demux.maxmismatch, Demux.ndetect, Demux.minpct))

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
    bclen = None
    leftTarget = None
    rightTarget = None

    def parseArgs(self, args):
        self.nf = 0
        self.standardOpts(args)
        cmd = args[0]
        if cmd in ['split', 'detect', 'grep']:
            self.mode = cmd
        else:
            P.errmsg(P.NOCMD)

        next = ""
        for a in args[1:]:
            if next == '-m':
                self.maxmismatch = self.toInt(a)
                next = ""
            elif next == '-b':
                self.bcfile = self.isFile(a)
                next = ""
            elif next == '-n':
                self.ndetect = self.toInt(a)
                next = ""
            elif next == '-p':
                self.minpct = self.toFloat(a)
                next = ""
            elif next == '-s':
                self.bclen = self.toInt(a)
                next = ""
            elif next == '-lt':
                self.leftTarget = a
                next = ""
            elif next == '-rt':
                self.rightTarget = a
                next = ""
            elif a in ['-m', '-b', '-n', '-p', '-s', '-lt', '-rt']:
                next = a
            elif a == '-r':
                self.revcomp = True
            elif a == '-u':
                self.undet = False
            elif self.nf == 0:
                self.fqleft = "-" if a == "-" else self.isFile(a)
                self.nf += 1
            elif self.nf == 1:
                self.fqright = self.isFile(a)
                self.nf += 1

P = Demux("demux", version="1.0", usage=usage,
          errors=[('NOCMD', 'Missing command', 'The first argument should be one of: split, detect.'),
                  ('BADFMT', 'Bad input file format', 'The input file should be in fasta of fastq format.') ])

### Utils

def revcomp(seq):
    nucmap = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
              'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rc = ""
    for i in range(len(seq)-1, -1, -1):
        b = seq[i]
        if b in nucmap:
            rc += nucmap[seq[i]]
        else:
            rc += b
    return rc

def distance(s1, s2):
    d = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            d += 1
    # print "{} {}: {}".format(s1, s2, d)
    return d

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

    def openStream(self, filename, filename2=None):
        self.filename = self.name + "-" + filename + ".fastq.gz"
        self.stream = Utils.genOpen(self.filename, "w")
        if filename2:
            self.filename2 = self.name + "-" + filename2 + ".fastq.gz"
            self.stream2 = Utils.genOpen(self.filename2, "w")

    def closeStream(self):
        self.stream.close()
        if self.stream2:
            self.stream2.close()

    def writeRecord(self, fq1, fq2=None):
        self.nhits += 1
        self.stream.write(fq1.name + "\n")
        self.stream.write(fq1.seq + "\n")
        self.stream.write("+\n")
        self.stream.write(fq1.qual + "\n")
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
        with open(filename, "r") as f:
            for line in f:
                line = line.rstrip("\r\n").split("\t")
                if len(line) < 2:
                    continue
                if line[0] == 'Other': # when reading output of `detect' command...
                    continue
                if rc:
                    bc = revcomp(line[1])
                else:
                    bc = line[1]
                self.add(line[0], bc)
        if undet:
            self.add("UND", "*")

    def add(self, name, seq):
        b = Barcode(name, seq)
        self.barcodeseqs.append(seq)
        self.barcodes[seq] = b

    def openAll(self, filename, filename2=None):
        for b in self.barcodes.itervalues():
            b.openStream(filename, filename2)

    def allFilenames(self):
        result = []
        for b in self.barcodes.itervalues():
            result.append(b.filename)
            if b.filename2:
                result.append(b.filename2)
        return result

    def closeAll(self):
        for b in self.barcodes.itervalues():
            b.closeStream()

    def findBest(self, seq, maxmismatch=1):
        maxd = len(seq)
        best = ""
        for bc in self.barcodeseqs:
            if bc != "*":
                d = distance(seq, bc)
                if d < maxd:
                    maxd = d
                    best = bc
        # print maxd
        if maxd <= maxmismatch:
            bb = self.barcodes[best]
            self.nhits += 1
            #bb.nhits += 1
            return bb
        elif "*" in self.barcodes:
            self.nhits += 1
            return self.barcodes["*"]
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
        ranking = [b for b in self.barcodes.itervalues()]
        ranking.sort(key=lambda b:b.nhits, reverse=True)

        sys.stdout.write("Name\tSeq\tHits\tPct\tFile1\tFile2\n")
        for b in ranking:
            seq = b.seq
            if P.revcomp:
                seq = revcomp(seq)
            if self.nhits == 0:
                pct = 0
            else:
                pct = Utils.f2dd(100.0 * b.nhits / self.nhits)
            if pct >= P.minpct:
                sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(b.name, seq, b.nhits, pct, b.filename, b.filename2 or ""))
            else:
                hiddenb += 1
                hiddenh += b.nhits
        if hiddenb > 0 and self.nhits > 0:
            sys.stdout.write("Other\t({})\t{}\t{}\t\t\n".format(hiddenb, hiddenh, Utils.f2dd(100.0 * hiddenh / self.nhits)))

class FastqRec():
    name = ""
    seq = ""
    qual = ""

    def getBarcode(self, bclen):
        if bclen == None:
            c = self.name.rfind(":")
            if c > 0:
                return self.name[c+1:]
            else:
                return None
        elif bclen > 0:
            if len(self.seq) > bclen:
                return self.seq[:bclen]
            else:
                return None
        else:
            if len(self.seq) > bclen:
                return self.seq[bclen:]
            else:
                return None
        
class FastqReader():
    filename = ""
    fq = None
    stream = None
    nread = 0
    ngood = 0

    def __init__(self, filename):
        self.filename = filename
        self.fq = FastqRec()
        if self.filename == '-':
            self.stream = sys.stdin
        else:
            self.stream = Utils.genOpen(filename, "r")
        self.nread = 0

    def nextRead(self):
        if self.stream:
            r = self.stream.readline()
            if r == '':
                self.stream.close()
                self.stream = None
                return None
            self.fq.name = r.rstrip("\r\n")
            self.fq.seq = self.stream.readline().rstrip("\r\n")
            self.stream.readline()
            self.fq.qual = self.stream.readline().rstrip("\r\n")
            self.nread += 1

    def demux(self, dm, filename, bclen=None):
        dm.openAll(filename)
        while True:
            self.nextRead()
            if not self.stream:
                break
            bc = self.fq.getBarcode(bclen)
            bo = dm.findBest(bc, maxmismatch=P.maxmismatch) # Barcode object
            # print "best: " + bo.seq
            # raw_input()
            if bo:
                self.ngood += 1
                bo.writeRecord(self.fq)
        dm.closeAll()
        sys.stderr.write("Total reads: {}\n".format(self.nread))
        sys.stderr.write("Written: {}\n".format(self.ngood))

    def detect(self, ndetect, bclen=None):
        """Detect the barcodes contained in the first `ndetect' reads of this fastq file."""
        dm = BarcodeMgr()
        while True:
            self.nextRead()
            if not self.stream:
                break
            bc = self.fq.getBarcode(bclen)
            dm.findOrCreate(bc)
            if self.nread == ndetect:
                break
        sys.stdout.write("# Reads examined: {}\n".format(dm.nhits))
        dm.showCounts()

class FastaReader(FastqReader):

    first = True
    saved = False

    def nextRead(self):
        if self.stream:
            self.fq.seq = ""

            if self.first:
                self.saved = self.stream.readline().rstrip("\r\n")[1:]
                self.first = False

            while True:
                r = self.stream.readline()
                if r == '':
                    self.stream.close()
                    self.stream = None
                    self.fq.name = self.saved
                    return False
                r = r.rstrip("\r\n")
                if len(r) > 0 and r[0] == '>':
                    self.fq.name = self.saved
                    self.saved = r[1:]
                    return True
                else:
                    self.fq.seq += r
        else:
            return False

class PairedFastqReader():
    reader1 = None
    reader2 = None
    fq1 = None
    fq2 = None
    nread = 0
    nbad  = 0
    ngood = 0
    
    def __init__(self, filename1, filename2):
        fmt = Utils.detectFileFormat(filename1)
        if fmt == 'fasta':
            self.reader1 = FastaReader(filename1)
        elif fmt == 'fastq':
            self.reader1 = FastqReader(filename1)
        else:
            P.errmsg(P.BADFMT)
        fmt = Utils.detectFileFormat(filename2)
        if fmt == 'fasta':
            self.reader2 = FastaReader(filename2)
        elif fmt == 'fastq':
            self.reader2 = FastqReader(filename2)
        else:
            P.errmsg(P.BADFMT)

        self.fq1 = self.reader1.fq
        self.fq2 = self.reader2.fq

    def nextRead(self):
        self.reader1.nextRead()
        self.reader2.nextRead()
        self.nread += 1

    def demux(self, dm, filename1, filename2, bclen=None):
        dm.openAll(filename1, filename2)
        while True:
            self.nextRead()
            if not self.reader1.stream:
                break
            bc1 = self.fq1.getBarcode(bclen)
            bc2 = self.fq2.getBarcode(bclen)
            if bc1 != bc2:
                self.nbad += 1
                continue
            # print "barcode: " + bc1
            bo = dm.findBest(bc1, maxmismatch=P.maxmismatch)
            # print "best: " + bo.seq
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

        nameleft = getFastqBasename(P.fqleft)
        if P.fqright:
            nameright = getFastqBasename(P.fqright)
        else:
            nameright = None

        bm = BarcodeMgr()
        bm.initFromFile(P.bcfile, rc=P.revcomp, undet=P.undet)

        if P.nf == 2:
            fr = PairedFastqReader(P.fqleft, P.fqright)
            fr.demux(bm, nameleft, nameright, P.bclen)
        else:
            fr = FastqReader(P.fqleft)
            fr.demux(bm, nameleft, P.bclen)
        bm.showCounts()

    elif P.mode == 'detect':
        fr = FastqReader(P.fqleft)
        fr.detect(P.ndetect, P.bclen)

    elif P.mode == 'grep':
        if P.rightTarget == None:
            P.rightTarget = revcomp(P.leftTarget)
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
        
if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        main(args)
    else:
        P.usage()

