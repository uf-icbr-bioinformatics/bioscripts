#!/usr/bin/env python

import sys
import pysam
import utils

def usage():
    sys.stderr.write("""bisconv.py command [args...] - Extract aligned reads from BAM file based on conversion strand.

This program examines the ZS tag of the reads in a BAM file produced by BSMAP to identify reads 
coming from conversion of the top or bottom strands. Commands:

-h                                  | print this help message.
split <bamfile> <topfile> <botfile> | write reads from top converted strand to <topfile>, and
                                    | reads from bottom converted strand to <botfile>.
zs+   <bamfile> <topfile>           | write reads from top converted strand to <topfile>. If
                                    | <topfile> is -, write to standard output.
zs-   <bamfile> <botfile>           | write reads from bottom converted strand to <botfile>. If
                                    | <botfile> is -, write to standard output.

Reads in <topfile> will only show C->T conversion, while reads in <botfile> will only show G->A
conversion.

""")

P = utils.Prog("bisconv.py", version="1.0", usage=usage)

## Code

def split(inbam, bamfile1, bamfile2):
    nin = 0
    nplus = 0
    nminus = 0
    inb = pysam.AlignmentFile(inbam, "rb")
    bam1 = pysam.AlignmentFile(bamfile1, "wb", template=inb)
    bam2 = pysam.AlignmentFile(bamfile2, "wb", template=inb)
    try:
        for read in inb.fetch():
            nin += 1
            zs = read.get_tag("ZS")
            if zs[0] == '+':
                bam1.write(read)
                nplus += 1
            else:
                bam2.write(read)
                nminus += 1
    finally:
        inb.close()
        bam1.close()
        bam2.close()
    return (nin, nplus, nminus)

def extract(inbam, outbam, flag):
    nin = 0
    nout = 0
    inb = pysam.AlignmentFile(inbam, "rb")
    bam = pysam.AlignmentFile(outbam, "wb", template=inb)
    try:
        for read in inb.fetch():
            nin += 1
            zs = read.get_tag("ZS")
            if zs[0] == flag:
                bam.write(read)
                nout += 1
    finally:
        inb.close()
        bam.close()
    return (nin, nout)

if __name__ == "__main__":
    args = sys.argv[1:]
    P.standardOpts(args)

    cmd = args[0]
    if cmd == 'split':
        (nin, nplus, nminus) = split(sys.argv[2], sys.argv[3], sys.argv[4])
        sys.stderr.write("Total: {}\nTop: {}\nBottom: {}\n".format(nin, nplus, nminus))
    elif cmd == 'zs+':
        (nin, nout) = extract(sys.argv[2], sys.argv[3], '+')
        sys.stderr.write("Total: {}\nTop: {}\n".format(nin, nout))
    else:
        extract(sys.argv[2], sys.argv[3], '-')
        sys.stderr.write("Total: {}\nBottom: {}\n".format(nin, nout))
