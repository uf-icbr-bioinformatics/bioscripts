#!/usr/bin/env python

import sys
import os.path
import gzip

import Utils
import Script

def usage():
    sys.stderr.write("""mergeCols.py - Merge columns from multiple files.

Usage: mergeCols.py [-o outfile] [-id 1] [-dc 2] [-cuff] [-rsem] [-na NA] [-names a,b,c...] file1 file2 ...
       mergeCols.py -a [-o outfile] [-id 1] [-i idsfile] colspecs file

Combine gene data from multiple input files into a matrix. The program assumes that 
each input files has one column containing identifiers (e.g. gene identifiers) and 
one column containing a value associated to each identifier (e.g. FPKM). It will 
output a tab-delimited file containing one row for each identifier, with the associated
values from each input file in consecutive columns.

If -a is specified, operates in `averager' mode: the output file will consist of (possibly averaged)
columns from the input file, according to colspecs. Colspecs is a comma-separated list of entries of
the form N+M+...+Q, representing the average of the values in the specified columns. For example:
3,5+6+7,9+10 will output column 3, the average of columns 5, 6, and 7, and the average of columns 10 and 11.

Options:
-o outfile | write matrix to `outfile' (default is standard output)
-id N      | column containing gene identifiers (default: 1)
-dc N      | column containing data (default: 2)
-cuff      | sets options to retrieve FPKM from cufflinks files (-id 1 -dc 10)
-rsem      | sets options to retrieve FPKM from RSEM files (-id 1 -dc 7)
-na S      | string to use for missing data (default: NA)
-names ... | comma-separated list of sample names (otherwise, use names
             of input files without extension)
-a         | enable averager mode.
-i I       | only output rows whose IDs appear in file I.

""")

P = Script.Script("mergecols.py", version="1.0", usage=usage)

# Utils

def filenameToSamplename(p):
    return os.path.splitext(os.path.basename(p))[0] # we should probably stop at the first .

def parseColspecs(cs):
    pieces = cs.split(",")
    pieces = [ [P.toInt(x) - 1 for x in p.split("+")] for p in pieces ]
    return pieces

def getData(data, cols):
    n = 0.0
    for c in cols:
        n += float(data[c])
    return n / len(cols)

class Opts():
    outfile = False
    na = "NA"
    id = 1
    dc = 2
    names = []
    infiles = []
    nfiles = 0
    valid = True
    avgmode = False
    colspecs = None
    idsfile = None
    wanted = {}

    def setNames(self):
        """Assign each missing element of the names field using the corresponding filename."""
        goodnames = [filenameToSamplename(f) for f in self.infiles]
        for i in range(min(len(self.infiles), len(self.names))):
            goodnames[i] = self.names[i]
        self.names = goodnames

    def confError(self, message):
        sys.stderr.write(message)
        self.valid = False

    def parseArgs(self, args):
        P.standardOpts(args)
        next = ""
        for a in args:
            if next == "-o":
                self.outfile = a
                next = ""
            elif next == "-id":
                self.id = P.toInt(a) - 1
                next = ""
            elif next == "-dc":
                self.dc = P.toInt(a) - 1
                next = ""
            elif next == "-names":
                self.names = a.split(",")
                next = ""
            elif next == "-na":
                self.na = a
                next = ""
            elif next == "-i":
                self.idsfile = a
                next = ""
            elif next == "-a":
                self.avgmode = True
                self.colspecs = parseColspecs(a)
                next = ""
            elif a in ["-o", "-id", "-dc", "-names", "-na", "-i", "-a"]:
                next = a
            elif a == "-cuff":
                self.id = 0
                self.dc = 9
            elif a == "-rsem":
                self.id = 0
                self.dc = 6
            else:
                self.infiles.append(P.isFile(a))
        self.nfiles = len(self.infiles)
        if self.nfiles == 0:
            self.confError("Error: no input files specified.\n")
        else:
            self.setNames()
        return self.valid

    def run(self):
        if self.avgmode:
            if self.outfile:
                with open(self.outfile, "w") as out:
                    self.mainAvg(out)
            else:
                self.mainAvg(sys.stdout)
        else:
            self.main()

    def readWanted(self):
        self.wanted = {}        # let's make sure...
        with open(self.idsfile, "r") as f:
            for line in f:
                line = line.strip()
                self.wanted[line] = True

    def mainAvg(self, out):
        hdr = True
        if self.idsfile:
            self.readWanted()
        with Utils.genOpen(self.infiles[0], "r") as f:
            for line in f:
                if hdr:
                    hdr = False
                else:
                    data = line.split("\t")
                    gid = data[self.id].strip('"') # in case the id is surrounded by quotes
                    if self.idsfile and gid not in self.wanted:
                        continue
                    outdata = [ getData(data, cs) for cs in self.colspecs ]
                    out.write(gid)
                    for d in outdata:
                        out.write("\t{}".format(d))
                    out.write("\n")
        

    def main(self):
        table = {}
        idx = 0
        idcol = self.id
        datacol = self.dc
        # print "datacol={}".format(datacol)
        for infile in self.infiles:
            sys.stderr.write("Reading {}... ".format(infile))
            nread = 0
            hdr = True              # skip header (should be replaced with cmdline option)
            with Utils.genOpen(infile, "r") as f:
                for line in f:
                    if hdr:
                        hdr = False
                    else:
                        nread += 1
                        data = line.rstrip("\r\n").split("\t")
                        # print data
                        gid = data[idcol]
                        gval = data[datacol]
                        if gid in table:
                            vec = table[gid]
                        else:
                            vec = [self.na]*self.nfiles
                            table[gid] = vec
                        vec[idx] = gval
            sys.stderr.write("{} records.\n".format(nread))
            idx += 1

        if self.outfile:
            sys.stderr.write("Writing table to {}.\n".format(self.outfile))
            out = open(self.outfile, "w")
        else:
            out = sys.stdout

        try:
            out.write("\t".join(["ID"] + self.names) + "\n")
            for (gid, vec) in table.iteritems():
                out.write(gid + "\t" + "\t".join(vec) + "\n")
        finally:
            if self.outfile:
                out.close()
        
if __name__ == "__main__":
    opts = Opts()
    if opts.parseArgs(sys.argv[1:]):
        opts.run()
    else:
        P.usage()
