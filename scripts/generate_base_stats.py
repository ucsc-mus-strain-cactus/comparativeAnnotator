#!/usr/bin/env python2.7
import sys, os
from collections import defaultdict

def main():
    best_alignments = {}
    helper_dict = {}
    with open(sys.argv[1], 'r') as f:
        for line in f:
            line = line.split()
            qName = line[9].split("-")[0]
            matches = int(line[0])
            mismatches = int(line[1])
            rep_matches = int(line[2])
            qSize = int(line[10])
            cov = (1.0 * matches + mismatches + rep_matches) / qSize
            assert 0 <= cov <= 1.0
            if qName not in best_alignments:
                best_alignments[qName] = cov
                helper_dict[qName] = [qSize, matches + mismatches + rep_matches]
            elif best_alignments[qName] < cov:
                best_alignments[qName] = cov
                helper_dict[qName] = [qSize, matches + mismatches + rep_matches]

    with open(sys.argv[2], 'w') as outf:
        outf.write("Name\tqSize\tbasesMapped\tfracMapped\n")
        for aId, cov in best_alignments.iteritems():
            qSize, alnBases = helper_dict[aId]
            outf.write("\t".join(map(str, [aId, qSize, alnBases, cov])) + "\n")

if __name__ == "__main__":
    main()
