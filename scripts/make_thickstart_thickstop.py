#!/usr/bin/env python
import os
import sys
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
from lib.seq_lib import GenePredTranscript

for rec in sys.stdin:
    t = GenePredTranscript(rec.rstrip().split("\t"))
    if t.cds_size == 0:
        continue
    q = t.get_bed(start_offset=t.thick_start, stop_offset=t.thick_stop)
    sys.stdout.write("\t".join(map(str, q)) + "\n")
