# Collapsing isoforms from CTR-ONT reads with TALON
We used the TALON software to collapse isoforms from CTR-Seq ONT data.

# Input data
```

```


# Read correction
## Canu
```
 canu -correct -p 20180530_0304_20180530-CTRS-Oleg-7 genomeSize=100m useGrid=false minReadLength=250 minOverlapLength=150 corOutCoverage=all overlapper=minimap -nanopore ../20180530_0304_20180530-CTRS-Oleg-7.raw.fastq

```
 
## TranscriptClean
We first mapped the reads to the genome with minimap.

```
minimap2 -t 62 -G 10k -ax splice -uf -k14 --MD /home/yuki.yoshida/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta.mmi ${FASTQ}
```

The TranscriptClean software was run following the protocols.

```
for i in *bam.sam; 
do
  python ./TranscriptClean/TranscriptClean.py --sam $i --genome ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta --outprefix $i.TranscriptClean.out -t 64
done
```

Resulting in the following files
```
20180510_0216_20180510-CTRS-Oleg-1.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
20180510_0216_20180510-CTRS-Oleg-1.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_read_labels.tsv
20180515_0328_20180515-CTRS-Oleg-2.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
20180515_0328_20180515-CTRS-Oleg-2.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_read_labels.tsv
20180518_0759_20180518-CTRS-Oleg-3.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
20180518_0759_20180518-CTRS-Oleg-3.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_read_labels.tsv
20180522_0743_20180522-CTRS-Oleg-4.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
20180522_0743_20180522-CTRS-Oleg-4.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_read_labels.tsv
20180525_0708_20180525-CTRS-Oleg-5.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
20180525_0708_20180525-CTRS-Oleg-5.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_read_labels.tsv
20180528_0710_20180528-CTRS-Oleg-6.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
20180528_0710_20180528-CTRS-Oleg-6.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_read_labels.tsv
20180530_0304_20180530-CTRS-Oleg-7.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
20180530_0304_20180530-CTRS-Oleg-7.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_read_labels.tsv
```

# Isoform collapsing
This analysis was conduced using Pv5.0 genome assembly. The only difference with Pv5.2 was the chromosome ordering and changing Pv5.0::chr1 to Pv5.2::chr3, Pv5.0::chr3 to Pv5.2::chr1.

## Working directory
`/home/yuki.yoshida/nias/analysis/gene_structure/talonls/cov_09`

## talon_label_reads
```
talon_initialize_database --f ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.clean.v2.gtf --g Pv5.0 --o pv5.0 --a olga

for i in *clean.sam; do
  talon_label_reads --f $i --g /home/yuki.yoshida/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta --t 64 --o $i
done

talon  --f config_all.txt --db pv5.0.db --build Pv5.0 --o talon_out -t 64

talon_create_GTF  --db pv5.0.db -b Pv5.0 -a olga  --o talon_filter_transcripts

gtf2gff.pl < talon_filter_transcripts_talon.gtf --out --gff3

perl -I/home/gaou/lcl/lib/perl5 ../gene_structure/flair/fix_gtf.pl talon_filter_transcripts_talon.gtf > talon_filter_transcripts_talon.clean.gtf
```

The contents of config_all.txt is as followed
```
Oleg1,wet-1,ONT,20180510_0216_20180510-CTRS-Oleg-1.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
Oleg2,wet-2,ONT,20180515_0328_20180515-CTRS-Oleg-2.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
Oleg3,D24-1,ONT,20180518_0759_20180518-CTRS-Oleg-3.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
Oleg4,D24-2,ONT,20180522_0743_20180522-CTRS-Oleg-4.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
Oleg5,D48,ONT,20180525_0708_20180525-CTRS-Oleg-5.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
Oleg6,Pv11-T0,ONT,20180528_0710_20180528-CTRS-Oleg-6.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
Oleg7,Pv11-T48,ONT,20180530_0304_20180530-CTRS-Oleg-7.raw.fastq.mmi_MD_10kIntron.genome.sorted.sam_labeled.sam
```

# SQANTI3 validation
```
## Add polyA_peak with SAGE peak data
python ~/bin/SQANTI3/sqanti3_qc.py  --gtf talon_out_transcripts_talon.gtf ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.clean.v2.gtf ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta  --aligner_choice minimap2 --cage_peak ~/nias/analysis/single_end_homer/20200625_addedRNASEQ/mergePeaks_d50_cage.bed9 -t 64  --polyA_peak ~/nias/analysis/single_end_homer/20200625_addedRNASEQ/mergePeaks_d50_sage.bed9
```
