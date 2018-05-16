# (c) 2016, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import os
import os.path
import sys
import csv
import math
import gzip
import time
import pysam
import string
import random

PYTHON_VERSION = sys.version_info[0]

def get_iterator(dict):
    if PYTHON_VERSION == 2:
        return dict.iteritems()
    else:
        return dict.items()

def genOpen(filename, mode):
    """Generalized open() function - works on both regular files and .gz files."""
    (name, ext) = os.path.splitext(filename)
    if ext == ".gz":
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def dget(key, dictionary, default=None):
    """Return the value associated with `key' in `dictionary', if present, or `default'."""
    if key in dictionary:
        return dictionary[key]
    else:
        return default

def dinc(key, dictionary, default=0):
    """Increment the value associated with `key' in dictionary. If key is not present, initialize it with `default'."""
    if key in dictionary:
        dictionary[key] += 1
    else:
        dictionary[key] = default + 1
    return dictionary[key]

def id_generator(size=6, chars=string.ascii_uppercase + string.digits, prefix=''):
    """From: http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python"""
    return prefix + ''.join(random.choice(chars) for _ in range(size))

def safeInt(v, default=None):
    try:
        return int(v)
    except ValueError:
        return default

def safeFloat(v, default=None):
    try:
        return float(v)
    except ValueError:
        return default

def convertValue(v):
    """Convert v to an int if possible, otherwise to
a float, otherwise return it unchanged."""
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v

def distance(v1, v2):
    """Euclidean distance between vectors v1 and v2."""
    def dsq(a, b):
        if a == None or b == None:
            return 0
        d = a - b
        return d*d

    return math.sqrt(sum([dsq(*p) for p in zip(v1, v2)] ))

def avgdiff(v1, v2):
    """Average difference between vectors v1 and v2 (v1 - v2)."""
    dsum = 0.0
    nd = 0
    for a, b in zip(v1, v2):
        if a == None or b == None:
            continue
        dsum += (a - b)
        nd += 1
    return dsum / nd

def colsToFloat(row, columns):
    """Returns a list containing the elements of `row' indexed by `columns' converted to floats.
Invalid entries are returned as None. If `columns' is None, use all columns."""
    if columns:
        return [ safeFloat(row[c]) for c in columns ]
    else:
        return [ safeFloat(x) for x in row ]

def _parseColspecAux(cs):
    p = cs.find("+")
    m = cs.find("-")
    if p == 0:
        v = safeInt(cs[1:])
        if v:
            return ("+", 1, v)
        else:
            return False
    elif m == 0:
        v = safeInt(cs[1:])
        if v:
            return ("-", 1, v)
        else:
            return False
    elif p > 0:
        w = "+"
        x = p
    elif m > 0:
        w = "-"
        x = m
    else:
        v = safeInt(cs)
        if v:
            return ("*", v)
        else:
            return False

    a = safeInt(cs[:x])
    b = safeInt(cs[x+1:])
    if a and b:
        return (w, a, b)
    else:
        return False

def parseColspec(s):
    """Parse a specification for a set of columns. `s' has the form C1,C2,...,Cn
where each C can be either a number (1-based column number), X-Y (from column X to
column Y inclusive), X+K (K columns starting at X)."""
    cols = []
    specs = s.split(",")
    for cs in specs:
        p = _parseColspecAux(cs)
        if p:
            if p[0] == "+":
                for i in range(p[1], p[1]+p[2]+1):
                    cols.append(i-1)
            elif p[0] == "-":
                for i in range(p[1], p[2]+1):
                    cols.append(i-1)
            elif p[0] == "*":
                cols.append(p[1]-1)
        else:
            sys.stderr.write("Incorrect field specification `{}'.\n".format(cs))
    return cols

def safeReadIntFromFile(filename, default=None, maxtries=20, delay=1):
    """Read an integer from the first line of file `filename', using `default' if the
line does not contain an integer. If the line is empty, wait `delay' seconds and
retry, up to `maxtries' times."""
    while True:
        with open(filename, "r") as f:
            c = f.readline().rstrip("\n")
        if c == '':
            maxtries -= 1
            if maxtries == 0:
                return default
            else:
                time.sleep(delay)
        else:
            return safeInt(c, default)

def safeReadLineFromFile(filename, default=None, maxtries=20, delay=1):
    """Read and return the first line of file `filename'. If the line is empty, wait `delay' seconds and
retry, up to `maxtries' times, after which return `default'."""
    while True:
        with open(filename, "r") as f:
            c = f.readline().rstrip("\n")
        if c == '':
            maxtries -= 1
            if maxtries == 0:
                return default
            else:
                time.sleep(delay)
        else:
            return c

def safeRemoveFile(filename):
    if os.path.isfile(filename):
        try:
            os.remove(filename)
        except OSError:
            pass

def filterFile(infile, outfile, column, low=None, high=None, invert=False, absolute=False, preserveHeader=True, delim='\t'):
    """Copy tab-delimited rows from `infile' to `outfile', testing the value in column `column' against `low' and `high'.
A row is copied if the value is above `low' (if specified) and below `high' (if specified). If `invert' is true, the test 
is inverted."""
    nout = 0
    colstr = (type(column).__name__ != 'int')
    with open(outfile, "w") as out:
        with open(infile, "r") as f:
            if preserveHeader or colstr:
                hdr = f.readline()
                if preserveHeader:
                    out.write(hdr)
                if colstr:
                    parsed = hdr.rstrip("\r\n").split(delim)
                    if column in parsed:
                        column = parsed.index(column)
                    else:
                        sys.stderr.write("Column `{}' not found in header of file {}.\n".format(column, infile))
                        return False
            for line in f:
                parsed = line.rstrip("\r\n").split(delim)
                x = safeFloat(parsed[column])
                if absolute:
                    x = abs(x)
                good = True
                if low and x < low:
                    good = False
                if high and x > high:
                    good = False
                if invert:
                    good = not good
                if good:
                    out.write(line)
                    nout += 1
    return nout

def filenameNoExt(s):
    return os.path.splitext(os.path.basename(s))[0]

def parseLine(line):
    return line.rstrip("\r\n").split("\t")

def fmt(n):
    return "{:,}".format(n)

def pct(x, den):
    if den == 0:
        return "0.0"
    else:
        return "{:.2f}%".format(x*100.0/den)

def f2dd(x):
    return "{:.2f}".format(x)

def f3dd(x):
    return "{:.3f}".format(x)

def f4dd(x):
    return "{:.4f}".format(x)

def up(x):
    return "<span class='upreg'>{}</span>".format(x)

def UP(x):
    return "<span class='upreg'><b>{}</b></span>".format(x)

def down(x):
    return "<span class='dnreg'>{}</span>".format(x)

def DOWN(x):
    return "<span class='dnreg'><b>{}</b></span>".format(x)

def plural(x):
    return "" if x == 1 else "s"

def fileToDict(filename, toInt=True, column=1):
    """Read `filename' and return a dictionary having as keys the strings
in the first column and as values the contents of the specified `column' 
(converted to int if `toInt' is True)."""
    result = {}
    with open(filename, "r") as f:
        for line in f:
            parsed = line.rstrip("\r\n").split("\t")
            if toInt:
                result[parsed[0]] = int(parsed[column])
            else:
                result[parsed[0]] = parsed[column]
    return result

def fileToList(filename, delimiter='\t', hdr=True):
    """Read delimited file `filename' and return its contents as a list of lists. If `hdr' is True (the
default) treat the first line as header. Returns a tuple (data, header)."""
    result = []
    header = []
    with open(filename, "r") as f:
        r = csv.reader(f, delimiter=delimiter)
        if hdr:
            header = r.next()
        for line in r:
            result.append(line)
    return (result, header)

def detectFileFormat(filename):
    """Examine the first line of file `filename' and return "fasta" if
it appears to be in FASTA format, "fastq" if it appears to be in fastq
format, "?" otherwise."""
    with genOpen(filename, "r") as f:
        line = f.readline()
        if len(line) > 0:
            if line[0] == ">":
                return "fasta"
            elif line[0] == "@":
                return "fastq"
        return "?"

### Smart CSV reader

class CSVreader():
    _reader = None
    ignorechar = '#'

    def __init__(self, stream, delimiter='\t'):
        self._reader = csv.reader(stream, delimiter=delimiter)

    def __iter__(self):
        return self

    def next(self):
        row = self._reader.next()
        while len(row) == 0 or row[0][0] == self.ignorechar:
            row = self._reader.next()
        else:
            return row

### Read a single column from a delimited file specified with the @filename:col notation.

class AtFileReader():
    filename = ""
    stream = None
    col = 0
    ignorechar = '#'

    def __init__(self, atspec, ignorechar='#'):
        self.ignorechar = ignorechar
        colonpos = atspec.rfind(":")
        if colonpos > 0:
            col = safeInt(atspec[colonpos+1:], False)
            if col:
                self.col = col - 1
                self.filename = atspec[1:colonpos]
        else:
            self.filename = atspec[1:]
        self.stream = genOpen(self.filename, "r")

    def __iter__(self):
        return self

    def close(self):
        self.stream.close()

    def next(self):
        while True:
            line = self.stream.readline()
            if line == '':
                self.close()
                raise StopIteration
            if line[0] != self.ignorechar:
                line = line.rstrip("\r\n").split("\t")
                return line[self.col]

    def readline(self):
        while True:
            line = self.stream.readline()
            if line == '':
                self.close()
                return line
            if line[0] != self.ignorechar:
                return line.rstrip("\r\n").split("\t")

### Parser for GTF files
### *** (this should be replaced with the parsers in genplot/genes)

def parseGTF(gtf):
    genes = {}
    transcripts = {}

    fn_open = gzip.open if gtf.endswith('.gz') else open
    with fn_open(gtf) as f:
        for line in f:
            if not line.startswith('#'):
                parseGTFrecord(line.rstrip("\r\n").split("\t"), genes, transcripts)
    return (genes, transcripts)

def parseGTFrecord(fields, genes, transcripts):
    what = fields[2]

    if what == 'gene':
        anndict = parseAnnotations(fields[8])
        gid = anndict['gene_id']
        genes[gid] = {'gene_id': gid,
                      'gene_biotype': anndict['gene_biotype'],
                      'gene_name': anndict['gene_name'] }

    elif what == 'transcript':
        anndict = parseAnnotations(fields[8])
        tid = anndict['transcript_id']
        transcripts[tid] = {'transcript_id': tid,
                            'gene_biotype': anndict['gene_biotype'],
                            'transcript_name': anndict['transcript_name'],
                            'gene_name': anndict['gene_name'] }

def parseAnnotations(ann):
    anndict = {}
    pieces = [ s.strip(" ") for s in ann.split(";") ]
    for p in pieces:
        pair = p.split(" ")
        if len(pair) == 2:
            anndict[pair[0]] = pair[1].strip('"')
    return anndict

### BAM utilities

def countReadsInBAM(filename):
    stats = pysam.idxstats(filename)
    nreads = 0
    for row in stats.split("\n"):
        row.rstrip("\r")
        fields = row.split("\t")
        if len(fields) > 2 and fields[0] != '*':
            nreads += int(fields[2])
    return nreads

def getChromsFromBAM(filename):
    chroms = []
    stats = pysam.idxstats(filename)
    for row in stats.split("\n"):
        fields = row.split("\t")
        if fields[0] != '*' and fields[0] != '':
            chroms.append(fields[0])
    return chroms

### resample (store a long vector into a short one, or vice-versa)

class Resampler():
    size1 = 0
    size2 = 0
    steps = 0
    p = 0.0
    q = 0.0
    vmap = []
    idxs1 = []
    idxs2 = []
    vec2  = []
    prev1 = -1
    prev2 = -1

    def __init__(self, size1, size2):
        self.init(size1, size2)

    def init(self, size1, size2):
        # print("Init, {}, {}".format(size1, size2))
        if size1 == 0 or size2 == 0:
            raise ValueError("Resampler cannot handle 0-length vectors!")
        if size1 == self.prev1 and size2 == self.prev2: # If sizes are unchanged...
            return                                      # No need to do anything.
        self.size1 = size1
        self.size2 = size2
        self.p = 1.0 / size1
        self.q = 1.0 / size2
        self.steps = (size1 * size2)
        self.vmap = []
        current = [-1, -1, 0]
        for i in range(self.steps):
            idx1 = int(i * self.q)
            idx2 = int(i * self.p)
            if idx1 == current[0] and idx2 == current[1]:
                current[2] += self.q
            else:
                current = [idx1, idx2, self.q]
                self.vmap.append(current)
        # print("Initialized vmap, length={}".format(len(self.vmap)))
        # self.idxs1 = [ int(i * self.q) for i in range(self.steps) ]
        # self.idxs2 = [ int(i * self.p) for i in range(self.steps) ]
        # print(self.idxs1)
        # print(self.idxs2)
        self.vec2 = []
        self.prev1 = size1
        self.prev2 = size2

    def resample(self, vec1):
        self.vec2 = [0.0]*self.size2

        for m in self.vmap:
            self.vec2[m[1]] = vec1[m[0]] * m[2]

        # for i in range(self.steps):
            # print("{} => {} ({})".format(wh1, wh2, i * self.q))
            # print("vec1[{}] * {} => vec2[{}]".format(int(wh1), self.q, int(wh2)))
            # print("{} {} {}".format(i, self.idxs2[i], self.idxs1[i]))
            # self.vec2[self.idxs2[i]] += vec1[self.idxs1[i]] * self.q
        return self.vec2

class CircBuf():
    """A simple circular buffer. Objects are added to the buffer with the add() method, and
retrieved with the get() method. The current() method returns the first element in the buffer
without removing it."""
    size = 0
    data = None
    aptr = 0
    bptr = 0

    def __init__(self, size=1000):
        self.size = size
        self.data = [None]*size
        self.aptr = 0
        self.bptr = 0

    def dump(self):
        sys.stderr.write("{}, a={}, b={}\n".format(self.data, self.aptr, self.bptr))

    def get(self):
        if self.aptr == self.bptr:
            return None
        v = self.data[self.aptr]
        self.data[self.aptr] = None
        self.aptr += 1
        if self.aptr == self.size:
            self.aptr = 0
        return v

    def current(self):
        if self.aptr == self.bptr:
            return None
        return self.data[self.aptr]

    def add(self, v):
        b = self.bptr
        self.data[self.bptr] = v
        self.bptr += 1
        if self.bptr == self.size:
            self.bptr = 0
        if self.bptr == self.aptr:
            raise "Error: circular buffer full!"
        return b

class Cons():
    """Let's reimplement CONS in Python..."""
    car = None
    cdr = None

    def __init__(self, car, cdr=None):
        self.car = car
        self.cdr = cdr

class LinkedList():
    """This object implements a linked list based on Cons(). Elements are added
at the front of the list with push() and removed with pop(). insert() adds an
element into a sorted list in the correct position. Specialize the smaller()
method to specify ordering."""
    top = None
    nconsed = 0

    def allocate(self, car, cdr=None):
        self.nconsed += 1
        return Cons(car, cdr)

    def deallocate(self, c):
        pass

    def preallocate(self, n):
        pass

    def smaller(self, a, b):
        return (a < b)

    def length(self):
        c = 0
        curr = self.top
        while True:
            if curr is None:
                return c
            curr = curr.cdr
            c += 1

    def current(self):
        if self.top:
            return self.top.car
        else:
            return None

    def push(self, v):
        c = self.allocate(v, self.top)
        self.top = c

    def pop(self):
        if self.top:
            c = self.top
            self.top = c.cdr
            return c
        else:
            return None

    def insert(self, v):
        """Insert value v in the appropriate place in this linked list. The 
smaller() method is used to determine ordering."""
        if self.top is None:
            self.top = self.allocate(v)
            return self.top

        p = self.top        # Current item
        if self.smaller(v, p.car):
            new = self.allocate(v, p)
            self.top = new
            return new

        return self.insertRec(p, v)

    def insertRec(self, p, v):
        while True:
            cdr = p.cdr
            if cdr is None:     # End of list?
                p.cdr = self.allocate(v)
                return
            if self.smaller(v, cdr.car): # v goes between current and next
                l = self.allocate(v, cdr)
                p.cdr = l
                return l
            # Else, move forward
            p = cdr

    def dump(self, p, out=sys.stderr):
        if p is None:
            p = self.top
        if p is None:
            out.write("()\n")
            return
        out.write("(")
        while True:
            out.write("{}".format(p.car))
            if p.cdr:
                out.write(" ")
                p = p.cdr
            else:
                out.write(")\n")
                return
    
class LinkedListPool(LinkedList):
    """This list keeps unused conses in a pool. When a new cons is needed,
it is retrieved from the pool, and when it is not needed anymore it can be
returned to the pool using the deallocate() method. The pool can be initialized
with the preallocate() method."""
    pool = None

    def preallocate(self, n):
        for i in range (n):
            self.nconsed += 1
            self.pool = Cons(None, self.pool)

    def allocate(self, car, cdr=None):
        c = self.pool
        if c:
            self.pool = c.cdr
            c.car = car
            c.cdr = cdr
            return c
        else:
            self.nconsed += 1
            return Cons(car, cdr)

    def deallocate(self, c):
        c.car = None
        c.cdr = self.pool
        self.pool = c

class GenomicWindower():
    window = 100
    chrom = None
    start = 0                   # Start of current window
    end = 0                     # End of current window
    end2 = 0                    # End of window after the current one
    data = []                   # Data contained in current window

    def __init__(self, window, out=None):
        self.window = window
        self.data = []
        self.out = out

    def open(self, echrom, epos, entry):
        self.chrom = echrom
        self.start = epos
        self.end = epos + self.window
        self.end2 = self.end + self.window
        self.data = [entry]
        
    def close(self):
        if self.data:
            self.outputWindow()

    def outputWindow(self):
        if self.out:
            self.out.write("{}\t{}\t{}\t{}\n".format(self.chrom, self.start, self.end, len(self.data)))

    def nextBlock(self, echrom, epos):
        if self.out:
            self.out.write("fixedStep chrom={} start={} step={} span={}\n".format(echrom, epos, self.window, self.window))

    def add(self, entry):
        echrom = entry[0]
        epos = entry[1]

        # Just starting? Open window.
        if self.chrom is None:
            self.open(echrom, epos, entry)
            self.nextBlock(echrom, epos)
            return entry

        # Are we starting a new chrom? Output last window of 
        # previous one and open new window.
        if self.chrom != echrom:
            self.outputWindow()
            self.open(echrom, epos, entry)
            self.nextBlock(echrom, epos)
            return entry

        # Are we still inside the window? If so, add this.
        if epos < self.end:
            self.data.append(entry)
            return entry

        # Are we in next window? Open consecutive window (without block)
        if epos < self.end2:
            self.outputWindow()
            self.open(echrom, self.end, entry)
            return entry

        # Else: new block, new window.
        self.outputWindow()
        self.open(echrom, epos, entry)
        self.nextBlock(echrom, epos)
            

def test():
    r = Resampler(9,3)
    v = range(9)
    #print(v)
    return r.resample(v)
    #return r.resample([1,1,1])
