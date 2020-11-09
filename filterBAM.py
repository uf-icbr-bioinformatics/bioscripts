#!/usr/bin/env python

import sys
import pysam

class BAMfilter(object):
    infile = None
    outfile = None
    max_template_len = None

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-l":
                self.max_template_len = int(a)
                prev = ""
            elif a in ["-l"]:
                prev = a
            elif not self.infile:
                self.infile = a
            else:
                self.outfile = a

    def filter(self):
        nin = 0
        nout = 0
        target = 1000000
        with pysam.AlignmentFile(self.infile, "rb") as f:
            with pysam.AlignmentFile(self.outfile, "wb", template=f) as out:
                for r in f:
                    nin += 1
                    #sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(r.reference_name, r.pos, r.next_reference_name, r.mpos, r.template_length))
                    if r.reference_name == r.next_reference_name and r.template_length:
                        if abs(r.template_length) <= self.max_template_len:
                            out.write(r)
                            nout += 1
                    if nin == target:
                        sys.stderr.write("{} read, {} written ({:.2f}%)\n".format(nin, nout, 100.0 * nout / nin))
                        target += 1000000
        sys.stderr.write("{} read, {} written ({:.2f}%)\n".format(nin, nout, 100.0 * nout / nin))

if __name__ == "__main__":
    B = BAMfilter()
    B.parseArgs(sys.argv[1:])
    B.filter()
