import Utils

# Utils for sequences (fastq/fasta)

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
    #print("distance {} {} \n".format(s1, s2))
    d = 0
    for i in range(min(len(s1), len(s2))):
        if s1[i] != s2[i]:
            d += 1
    # print("{} {}: {}".format(s1, s2, d))
    return d

# Classes

class FastqRec():
    name = ""
    seq = ""
    qual = ""

    def getBarcode(self, bcslice, removePlus=False):
        if bcslice:
            return self.seq[bcslice]
        c = self.name.rfind(":")
        if c > 0:
            bc = self.name[c+1:]
            if removePlus:
                plus = bc.find("+")
                if plus > 0:
                    bc = bc[:plus]
            return bc
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
            r = self.stream.readline().decode()
            if r == '':
                self.stream.close()
                self.stream = None
                return None
            self.fq.name = r.rstrip("\r\n")
            self.fq.seq = self.stream.readline().decode().rstrip("\r\n")
            self.stream.readline()
            self.fq.qual = self.stream.readline().decode().rstrip("\r\n")
            self.nread += 1

    def isActive(self):
        return self.stream

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
                    self.nread += 1
                    return False
                r = r.rstrip("\r\n")
                if len(r) > 0 and r[0] == '>':
                    self.fq.name = self.saved
                    self.saved = r[1:]
                    self.nread += 1
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

    def isActive(self):
        return self.reader1.isActive() and self.reader2.isActive()
