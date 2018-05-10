#!/usr/bin/env python

import os.path
import sys
import csv
import Utils

DATA="/ufrc/data/reference/icbr/Affymetrix/HG-U133A.na22.annot.csv"

def parseFilespec(fs):
    """Parse a filespec of the form 'filename' or 'filename:c', where filename is an existing file
and c is a column number (defaulting to 1). Returns a tuple (filename, c-1) if successful, False 
otherwise."""
    c = 1
    p = fs.find(":")
    if p > -1:
        s = fs[p+1:]
        if s != '':
            c = int(s)
        fs = fs[:p]
    if os.path.isfile(fs):
        return (fs, c-1)
    else:
        return False

def parseAnnots():
    annots = {}
    sys.stderr.write("Parsing array annotations: {}\n".format(DATA))
    with open(DATA, "r") as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        hdr = reader.next()
        try:
            gpos = hdr.index("Gene Symbol")
        except:
            sys.stderr.write("Gene symbol not found in annotations??\n")
            return
        sys.stderr.write("Gene symbol found in column {}.\n".format(gpos))
        for row in reader:
            annots[row[0]] = row[gpos]
    return annots

def main(infile, outfile):
    annots = parseAnnots()
    k = annots.keys()[0]
    print(k, annots[k])
    (filename, col) = parseFilespec(infile)
    dcol = col+1
    sys.stderr.write("{} -> {}\n".format(filename, outfile))
    nin = nout = 0
    with open(outfile, "w") as out:
        with open(filename, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for line in reader:
                nin += 1
                probeset = line[col]
                if probeset in annots:
                    line[dcol:dcol] = [annots[probeset]]
                    nout += 1
                else:
                    line[dcol:dcol] = ['???']
                out.write("\t".join(line) + "\n")
    sys.stderr.write("{} probesets, {} translated, {} untranslated.\n".format(nin, nout, nin-nout))

if __name__ == "__main__":
    if len(sys.argv) == 4:
        DATA=sys.argv[3]
    main(sys.argv[1], sys.argv[2])
