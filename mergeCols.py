#!/usr/bin/env python

import sys
import os.path
import gzip

import Script

def usage():
    sys.stderr.write("""mergeCols.py - Merge columns from multiple files.

Usage: mergeCols.py [-o outfile] [-id 1] [-dc 2] [-cuff] [-rsem] [-na NA] [-names a,b,c...] files...

Combine gene data from multiple input files into a matrix. The program assumes that 
each input files has one column containing identifiers (e.g. gene identifiers) and 
one column containing a value associated to each identifier (e.g. FPKM). It will 
output a tab-delimited file containing one row for each identifier, with the associated
values from each input file in consecutive columns.

Options:
-o outfile | write matrix to `outfile' (default is standard output)
-id N      | column containing gene identifiers (default: 1)
-dc N      | column containing data (default: 2)
-cuff      | sets options to retrieve FPKM from cufflinks files (-id 1 -dc 10)
-rsem      | sets options to retrieve FPKM from RSEM files (-id 1 -dc 7)
-na S      | string to use for missing data (default: NA)
-names ... | comma-separated list of sample names (otherwise, use names
             of input files without extension)

""")

P = Script.Script("mergecols.py", version="1.0", usage=usage)

# Utils

def isGzip(filename):
    (name, ext) = os.path.splitext(filename)
    return (ext == ".gz")

def genOpen(filename, mode, gzipmode=False):
    """Generalized open() function - works on both regular files and .gz files."""
    if gzipmode:
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def filenameToSamplename(p):
    return os.path.splitext(os.path.basename(p))[0] # we should probably stop at the first .

class Opts():
    outfile = False
    na = "NA"
    id = 1
    dc = 2
    names = []
    infiles = []
    nfiles = 0
    valid = True

    def setNames(self):
        """Assign each missing element of the names field using the corresponding filename."""
        goodnames = [filenameToSamplename(f) for f in self.infiles]
        for i in range(min(len(self.infiles), len(self.names))):
            goodnames[i] = self.names[i]
        self.names = goodnames

    def confError(self, message):
        sys.stderr.write(message)
        self.valid = False

def parseInt(s):
    try:
        x = int(s)
        return x - 1
    except ValueError:
        return None

def parseArgs(args):
    P.standardOpts(args)
    opts = Opts()
    next = ""
    for a in args:
        if next == "-o":
            opts.outfile = a
            next = ""
        elif next == "-id":
            opts.id = P.toInt(a)
            next = ""
        elif next == "-dc":
            opts.dc = P.toInt(a)
            next = ""
        elif next == "-names":
            opts.names = a.split(",")
            next = ""
        elif next == "-na":
            opts.na = a
            next = ""
        elif a in ["-o", "-id", "-dc", "-names", "-na"]:
            next = a
        elif a == "-cuff":
            opts.id = 0
            opts.dc = 9
        elif a == "-rsem":
            opts.id = 0
            opts.dc = 6
        else:
            opts.infiles.append(P.isFile(a))
    opts.nfiles = len(opts.infiles)
    if opts.nfiles == 0:
        opts.confError("Error: no input files specified.\n")
    else:
        opts.setNames()
    return opts

def main(opts):
    table = {}
    idx = 0
    idcol = opts.id
    datacol = opts.dc
    # print "datacol={}".format(datacol)
    for infile in opts.infiles:
        sys.stderr.write("Reading {}... ".format(infile))
        nread = 0
        hdr = True              # skip header (should be replaced with cmdline option)
        with genOpen(infile, "r", gzipmode=isGzip(infile)) as f:
            for line in f:
                if hdr:
                    hdr = False
                else:
                    nread += 1
                    data = line.split("\t")
                    # print data
                    gid = data[idcol]
                    gval = data[datacol]
                    if gid in table:
                        vec = table[gid]
                    else:
                        vec = [opts.na]*opts.nfiles
                        table[gid] = vec
                    vec[idx] = gval
        sys.stderr.write("{} records.\n".format(nread))
        idx += 1
        
    if opts.outfile:
        sys.stderr.write("Writing table to {}.\n".format(opts.outfile))
        out = open(opts.outfile, "w")
    else:
        out = sys.stdout

    try:
        out.write("\t".join(["ID"] + opts.names) + "\n")
        for (gid, vec) in table.iteritems():
            out.write(gid + "\t" + "\t".join(vec) + "\n")
    finally:
        if opts.outfile:
            out.close()
        
if __name__ == "__main__":
    opts = parseArgs(sys.argv[1:])
    if opts.valid:
        main(opts)
    else:
        P.usage()
