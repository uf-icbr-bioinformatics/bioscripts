#!/usr/bin/env python

import sys

def maxgap(stream):
    maxg = 0
    prev = int(stream.readline().rstrip("\r\n"))
    while True:
        next = stream.readline()
        if next == '':
            break
        x = int(next.rstrip("\r\n"))
        gap = x-prev
        if gap > maxg:
            maxh = gap
        prev = x
    return maxg

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        with open(filename, "r") as f:
            maxg = maxgap(f)
    else:
        maxg = maxgap(sys.stdin)
    print "Max gap: {}".format(maxg)
