#!/usr/bin/env python

###################################################
#
# (c) 2016, Alberto Riva, ariva@ufl.edu
# DiBiG, ICBR Bioinformatics, University of Florida
#
# See the LICENSE file for license information.
###################################################

import sys
import csv
import os.path

import Script

### Program definition

def usage():
    sys.stderr.write("""cols.py - analyze columns in a delimited file.

Usage: cols.py [-h] [-c] [-0] [-m] [-d D] [-s N] [-R] [-i I] files...

With no options, print the number of columns in each of the supplied 
files (computed as the number of fields in header, by default the 
first row of the file). If no files are supplied, read standard input.

Options:

  -h, --help     | Print this help message.
  -c, --colnames | Print entries in the first line as numbered list.
  -m, --matches  | Print the number of lines matching the number of
                   fields in the header line.
  -f, --fields F | Print to standard output the columns of each input
                   file matching the items in file F (one per line).
  -r, --raw      | Print the contents of the header line, one item
                   per line.
  -n, --ncols    | Print the number of columns in each line of the input.
                   Output is tab delimited with two columns: line number
                   and number of columns.
  -s N           | Use line N as the header. If N is not a number, it
                   is interpreted as a prefix indicating lines to be
                   skipped. The header will be the first line that
                   does not start with N.
  -S N           | Print the contents of line N using the first line
                   as field names, in the format "field = value".
  -0             | Number columns starting at 0 instead of 1.
  -d D           | Use D as the delimiter (default: tab).
  -R             | Assume R-style header (number of fields in header
                   line is one less than the number of fields in the
                   rest of the file.

""")

class Cols(Script.Script):
    DELIMITER = '\t'
    MODE = 'default'
    RSTYLE = False
    ZERO = False
    FIRSTHDR = False                # First row is header regardless of ROWNUM
    ROWNUM = 1
    SKIP = None
    RAW = False
    INFILES = []
    FIELDS = []

    def parseOptions(self, args):
        C.standardOpts(args)
        next = ""
        for a in args:
            if next == "-d":
                self.DELIMITER = a
                next = ""
            elif next == "-s":
                try:
                    self.ROWNUM = int(a)
                except ValueError:
                    self.SKIP = a
                next = ""
            elif next == "-S":
                self.ROWNUM =self. toInt(a)
                self.FIRSTHDR = True
                next = ""
            elif next in ["-f", "--fields"]:
                self.MODE = 'fields'
                self.FIRSTHDR = True
                filename = self.isFile(a)
                with open(filename, "r") as f:
                    fields = f.read()
                self.FIELDS = fields.split("\n")
                next = ""
            elif a in ['-d', '-s', '-S', '-f', '--fields']:
                next = a
            elif a in ['-c', '--colnames']:
                self.MODE = 'names'
            elif a in ['-m', '--matches']:
                self.MODE = 'matches'
            elif a in ['-n', '--ncols']:
                self.MODE = 'ncols'
            elif a in ['-r', '--raw']:
                self.MODE = 'names'
                self.RAW = True
            elif a == '-0':
                self.ZERO = True
            elif a == '-R':
                self.RSTYLE = True
            else:
                self.INFILES.append(a)
        
    def processFromStream(self, filename, f):
        # print("skipping to row {}".format(self.ROWNUM))
        hdr = []
        row = []
        reader = csv.reader(f, delimiter=self.DELIMITER, quotechar='"')
        # print(MODE, ROWNUM, DELIMITER, SKIP, FIRSTHDR, RAW)

        if self.FIRSTHDR:
            hdr = reader.next()

        if self.RSTYLE:
            hdr = ["<<RowNum>>"] + hdr

        for i in range(1, self.ROWNUM):
            reader.next()

        while True:
            row = reader.next()
            if self.SKIP == None or row[0][0] != self.SKIP:
                break

        ncols = len(row)
        if self.MODE == 'default':
            sys.stderr.write("{}: {} columns.\n".format(filename, ncols))
        elif self.MODE == 'names':
            sys.stderr.write("{}: {} columns.\n".format(filename, ncols))

            if self.RAW:
                for p in range(ncols):
                    print(row[p])
            elif self.FIRSTHDR:
                for p in range(ncols):
                    print("  {} = {}".format(hdr[p], row[p]))
            else:
                idx = 0 if self.ZERO else 1
                for h in row:
                    print("  {} = {}".format(idx, h))
                    idx += 1

        elif self.MODE == 'matches':
            total = 0
            matching = 0
            for line in f:
                total += 1
                if line.count(self.DELIMITER) + 1 == ncols:
                    matching += 1
            sys.stderr.write("{}: {} columns, {}/{} matching.\n".format(filename, ncols, matching, total))
        elif self.MODE == 'fields':
            self.doFields(hdr, reader)
        elif self.MODE == 'ncols':
            ln = 1
            for line in reader:
                sys.stdout.write("{}\t{}\n".format(ln, len(line)))
                ln += 1

    def doFields(self, hdr, reader):
        columns = []
        unm = []
        ngood = 0
        nbad = 0
        for f in self.FIELDS:
            if f in hdr:
                columns.append(hdr.index(f))
                ngood += 1
            else:
                unm.append(f)
                nbad += 1
        sys.stderr.write("{} field names matched, {} not matched.\n".format(ngood, nbad))
        if nbad > 0:
            sys.stderr.write("Unmatched: {}\n".format(unm))
        if columns:
            sys.stdout.write("\t".join([hdr[c] for c in columns]) + "\n")
            for row in reader:
                outrow = [row[c] for c in columns]
                sys.stdout.write("\t".join(outrow) + "\n")

    def processOneFile(self, filename):
        ncols = 0
        if os.path.exists(filename):
            with open(filename, "r") as f:
                self.processFromStream(filename, f)
        else:
            sys.stderr.write("{}: File not found.\n".format(filename))

if __name__ == "__main__":
    C = Cols("cols.py", "1.0", usage=usage)
    C.parseOptions(sys.argv[1:])
    if C.INFILES:
        for f in C.INFILES:
            C.processOneFile(f)
    else:
        C.processFromStream("stdin", sys.stdin)
