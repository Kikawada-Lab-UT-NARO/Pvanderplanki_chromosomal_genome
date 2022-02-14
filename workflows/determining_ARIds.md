# Determining ARId regions in the Pv5.2 genome
Input is `input/arids.fasta` obtained from our previous paper (Gusev, et al., 2014) and uses `scripts/pysam_genome_rightmost_position.py` to obtain ARId regions.

## Working directory
/home/yuki.yoshida/nias/analysis/reanalysis/07_ARID

## Run minimap2 to determine ARId regions
```
minimap2 -ax asm5 /home/yuki.yoshida/nias/Pv_Final_assembly_2018_Olga/final_
version_from_olga/Pv11_5.0.final.fasta ARIds.fasta
python pysam_genome_rightmost_position.py ARIds.fasta.mmi.genome.sam > ARIds.bed
```
 
 The resulting `ARIds.bed` was edited manually to merge consecutive regions.
