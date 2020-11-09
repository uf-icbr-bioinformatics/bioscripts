#!/usr/bin/env python

import sys
import Utils

class outStream():
    name = ""
    left = None
    lefts = None
    right = None
    rights = None

    def __init__(self, name):
        self.name = name

    def openStreams(self):
        self.left = self.name + "_R1.fastq.gz"
        self.right = self.name + "_R2.fastq.gz"
        self.lefts = Utils.genOpen(self.left, "w")
        self.rights = Utils.genOpen(self.right, "w")

    def closeStreams(self):
        self.lefts.close()
        self.rights.close()

    def writeRead1(self, name, seq, qual):
        self.lefts.write(name )
        self.lefts.write(seq + "+\n")
        self.lefts.write(qual)

    def writeRead2(self, name, seq, qual):
        self.rights.write(name)
        self.rights.write(seq + "+\n")
        self.rights.write(qual)

class Barcode():
    seq = ""
    length = 0
    stream = None

    def __init__(self, seq):
        self.seq = seq
        self.length = len(seq)

class BarcodeManager():
    barcodes = []
    streams = []

    def sortBarcodes(self):
        self.barcodes.sort(key=lambda b: b.length, reverse=True)

    def openAllStreams(self):
        for s in self.streams:
            s.openStreams()

    def closeAllStreams(self):
        for s in self.streams:
            s.closeStreams()

    def findBarcode(self, seq):
        offset = 1 if seq[0] == "N" else 0
        for b in self.barcodes:
            if seq[offset:b.length] == b.seq[offset:]:
                return b
        return None

def readBarcodes(filename):
    BM = BarcodeManager()
    with open(filename, "r") as f:
        while True:
            sample = f.readline()
            if not sample:
                break
            sample = sample.strip()
            os = outStream(sample)
            BM.streams.append(os)
            for i in range(4):
                seq = f.readline().strip()
                bc = Barcode(seq)
                bc.stream = os
                BM.barcodes.append(bc)
    BM.sortBarcodes()
    return BM

def demux(BM, lfq, rfq):
    nin = 0
    nout = 0
    nbad = 0
    nmism = 0
    BM.openAllStreams()
    bad = outStream("Bad")
    bad.openStreams()
    with Utils.genOpen(lfq, "r") as lf:
        with Utils.genOpen(rfq, "r") as rf:
            while True:
                hdr1 = lf.readline()
                hdr2 = rf.readline()
                if not hdr1:
                    break
                nin += 1
                seq1 = lf.readline()
                seq2 = rf.readline()
                lf.readline()
                rf.readline()
                qual1 = lf.readline()
                qual2 = rf.readline()
                bc1 = BM.findBarcode(seq1)
                bc2 = BM.findBarcode(seq2)
                if bc1 and bc2:
                    os1 = bc1.stream
                    os2 = bc2.stream
                    if os1.name == os2.name:
                        index1 = seq1[:bc1.length]
                        index2 = seq2[:bc2.length]
                        hdr1 = hdr1[:-9] + index1 + "\n"
                        hdr2 = hdr2[:-9] + index2 + "\n"
                        nout += 1
                        os1.writeRead1(hdr1, seq1[bc1.length:], qual1[bc1.length:])
                        os1.writeRead2(hdr2, seq2[bc2.length:], qual2[bc2.length:])
                    else:
                        nmism += 1
                else:
                    nbad += 1
                    bad.writeRead1(hdr1, seq1, qual1)
                    bad.writeRead2(hdr2, seq2, qual2)
#                if nin == 100000:
#                    break
    BM.closeAllStreams()
    bad.closeStreams()
    sys.stderr.write("In: {}\nOut: {}\nBad: {}\nMismatch: {}\n".format(nin, nout, nbad, nmism))

if __name__ == "__main__":
    args = sys.argv[1:]
    if args:
        bcfile = args[0]
        lfq = args[1]
        rfq = args[2]
        BM = readBarcodes(sys.argv[1])
        for b in BM.barcodes:
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(b.seq, b.length, b.stream.name, b.stream.left, b.stream.right))
        demux(BM, lfq, rfq)
    else:
        sys.stdout.write("""Usage: crisprdemux.py barcodes left right

Demultiplex two paired-end fastq files `left' and `right' based on the barcodes 
in the barcodes file. Each sample has four barcodes of varying lengths. The
barcodes file should have the following format:

  NAME1
  barcode1
  barcode2
  barcode3
  barcode4
  NAME2
  ...

""")

