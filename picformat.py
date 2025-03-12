#!/usr/bin/env python

import sys
import os.path

def read_metrics(filename, data=False):
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("## METRICS CLASS"):
                if data:
                    f.readline()
                return f.readline()
    return False


def main(filenames):
    first = filenames[0]
    hdr = read_metrics(first)
    if hdr:
        sys.stdout.write("Sample\t" + hdr)
        for filename in filenames:
            name = os.path.split(filename)[1].split(".")[0]
            sys.stdout.write(name + "\t" + read_metrics(filename, data=True))

if __name__ == "__main__":
    main(sys.argv[1:])
    
