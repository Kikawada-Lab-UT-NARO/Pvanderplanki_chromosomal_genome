# Determining ARId regions in the Pv5.2 genome
Input is `input/arids.fasta` obtained from our previous paper (Gusev, et al., 2014).

## Working directory
/home/yuki.yoshida/nias/analysis/reanalysis/07_ARID

## Script
- pysam_genome_rightmost_position.py

```
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
```


## Run minimap2 to determine ARId regions
```
minimap2 -ax asm5 /home/yuki.yoshida/nias/Pv_Final_assembly_2018_Olga/final_
version_from_olga/Pv11_5.0.final.fasta ARIds.fasta
python pysam_genome_rightmost_position.py ARIds.fasta.mmi.genome.sam > ARIds.bed
```
 
 The resulting `ARIds.bed` was edited manually to merge consecutive regions.
