#!/usr/bin/env python

import pysam
import sys

insam= sys.argv[1]

#samfile = pysam.AlignmentFile(insam, "rb")
samfile = pysam.AlignmentFile(insam, "r")
for aln in samfile:
    strand = '+' if not aln.is_reverse else '-'
    print(str(aln.query_name) + "\t" + str(aln.reference_name) + "\t" + str(aln.reference_start) + "\t" +  str(aln.reference_end) + "\t" + strand )
samfile.close()
sys.exit()
