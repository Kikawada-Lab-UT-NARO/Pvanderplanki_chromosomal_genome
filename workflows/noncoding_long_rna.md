# Protocol for identification of non-coding lon RNA
we used RNA-Seq data from our previous paper (Gusev, 2014) and CTR-Seq ONT reads to obtain full length non-coding transcripts, not predicted in the Braker gene prediction.

## short read RNA-Seq assembly
### trinity
```
~/bin/trinityrnaseq-v2.11.0/Trinity  --seqType fq --left DRR024752.sra_1.fastq,DRR024753.sra_1.fastq,DRR024754.sra_1.fastq,DRR024755.sra_1.fastq,DRR024756.sra_1.fastq --right DRR024752.sra_2.fastq,DRR024753.sra_2.fastq,DRR024754.sra_2.fastq,DRR024755.sra_2.fastq,DRR024756.sra_2.fastq --include_supertranscripts --CPU 64 --max_memory 200G --output trinity_noLongRead
~/bin/trinityrnaseq-v2.11.0/Trinity --long_reads ONT_full_length_strandFixed.fq --seqType fq --genome_guided_bam merged.hisat2.sorted.bam --include_supertranscripts --CPU 64 --max_memory 200G --genome_guided_max_intron 1000
~/bin/trinityrnaseq-v2.11.0/Trinity --long_reads ONT_full_length_strandFixed.fq --seqType fq --genome_guided_bam merged.hisat2.sorted.bam --include_supertranscripts --CPU 64 --max_memory 200G --genome_guided_max_intron 5000
```

### rnaSPADES
```
~/bin/SPAdes-3.14.1-Linux/bin/rnaspades.py -o rnaspades --pe1-1 DRR024752.sra_1.fastq --pe2-1 DRR024753.sra_1.fastq --pe3-1 DRR024754.sra_1.fastq --pe4-1 DRR024756.sra_1.fastq --pe1-2 DRR024752.sra_2.fastq --pe2-2 DRR024753.sra_2.fastq --pe3-2 DRR024754.sra_2.fastq --pe4-2 DRR024756.sra_2.fastq --fl-rna ONT_full_length_strandFixed.fq  -t 64 -m 200
```

## long read consensus
### isONclust + medaka

```
/home/yuki.yoshida/anaconda3/envs/isoncorrect/bin/isONclust --t 64 --ont --fastq merged_full_length_unmapped.fastq --outfolder isONclust/
/home/yuki.yoshida/anaconda3/envs/isoncorrect/bin/isONclust write_fastq --N 1 --clusters isONclust/final_clusters.tsv --fastq merged_full_length_unmapped.fastq --outfolder isONclust/fastq_files

/home/yuki.yoshida/nias/analysis/non_coding/full_length/isONclust/fastq_files/run_medaka.pl

```

## hybrid assembly (mimiking LoReAn)
- map rnaSPADEs assembly to genome with minimap2
```
~/anaconda3//envs/talon/bin/minimap2 -ax splice -t 64 -a ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta transcripts.fasta.blastn.olgaCDS.1e-30.nohitlist.fna | samtools view -bS -@ 64 - | samtools sort -@ 64 -o transcripts.fasta.blastn.olgaCDS.1e-30.nohitlist.fna.mmi.genome.bam -
~/bin/bedtools2/bin/bedtools bamtobed -bed12 -split -i transcripts.fasta.blastn.olgaCDS.1e-30.nohitlist.fna.mmi.genome.bam > transcripts.fasta.blastn.olgaCDS.1e-30.nohitlist.fna.mmi.genome.bam.bed12
```
- map raw full length strand corrected ONT reads to genome with minimap2 and convert to bed12
```
~/anaconda3/envs/talon/bin/minimap2 -ax splice -t 64 -a ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta  merged_full_length_unmapped.fastq | samtools view -bS -@ 64 - | samtools sort -@ 64 -o merged_full_length_unmapped.fastq.mmi.genome.sorted.bam -
~/bin/bedtools2/bin/bedtools bamtobed -bed12 -split -i merged_full_length_unmapped.fastq.mmi.genome.sorted.bam > merged_full_length_unmapped.fastq.mmi.genome.sorted.bed12
```

- merge bam files
```
# made bed12 before bam merge
samtools merge both.bam merged_full_length_unmapped.fastq.mmi.genome.sorted.bam transcripts.fasta.blastn.olgaCDS.1e-30.nohitlist.fna.mmi.genome.bam
samtools sort -@ 64 -o both.sorted.bam both.bam
```

 - Make consensus sequences
```
# get primary alignments from bam file
samtools view  -bS -@ 64  -F 256 both.bam   | samtools sort -@ 64 -o both.No256.sorted.bam -
samtools index both.No256.sorted.bam

# convert to bed12 and sort
~/bin/bedtools2/bin/bedtools bamtobed -bed12 -split -i both.No256.sorted.bam > both.No256.sorted.bed12
~/bin/bedtools2/bin/bedtools sort -i both.No256.sorted.bed12 > both.No256.sorted.sorted.bed12

# create cluster
~/bin/bedtools2/bin/bedtools merge -c 4,4 -o count,distinct -i both.No256.sorted.sorted.bed12 > both.No256.sorted.sorted.merged.bed12
~/bin/bedtools2/bin/bedtools sort -i both.No256.sorted.sorted.merged.bed12 > both.No256.sorted.sorted.merged.sorted.bed12
perl -lane 'print if @F[3] >= 5' both.No256.sorted.sorted.merged.sorted.bed12 > both.No256.sorted.sorted.merged.sorted.over5.bed12

# make fastq files for each cluster
while read line ; do; echo $line | perl -lane 'print "samtools view -h both.No256.sorted.bam @F[0]:@F[1]-@F[2] | samtools fastq - > clusters/@F[0]_@F[1]_@F[2].fq"  '; done < both.No256.sorted.sorted.merged.sorted.bed12
cd clusters

# run racon + medaka
for i in *fq; do; echo "sh ../run_multiple_racon.sh $i" | qsub -l nodes=1:ppn=32 -N $i -d /home/yuki.yoshida/nias/analysis/non_coding/full_length/assembly/clusters; done;
cat clusters/*consensus.fa > final_consensus.fa
# remove redundant sequences
 ~/bin/cdhit-master/cd-hit-est -c 0.99  -i final_consensus.fa -o final_consensus.fa.cdhitest.fa
```

- remove rRNA and tRNA genes
```
# map to genome
~/anaconda3/envs/lorean22/bin/gmap -D ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta.gmap -d Pv11_5.0 -t 62 -O final_consensus.fa.cdhitest.fa -f 3 -B 5 -n 0 -f gff3_gene > final_consensus.fa.cdhitest.fa.gmap.genome.gff

# identify tRNA and rRNA genes
~/anaconda3/envs/lorean22/bin/gmap -D ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta.gmap -d Pv11_5.0 -t 62 -O final_consensus.fa -f 3 -B 5 -n 0 -f gff3_gene > final_consensus.fa.gmap.genome.gff

tRNAscan-SE -E -o Pv11_5.0.final.fasta.trnascan Pv11_5.0.final.fasta
barrnap --kingdom euk --threads 64 --outseq Pv11_5.0.final.fasta.barrnap.fna Pv11_5.0.final.fasta > Pv11_5.0.final.fasta.barrnap.gff3

~/bin/bedtools2/bin/bedtools intersect -a final_consensus.fa.cdhitest.fa.gmap.genome.gff -b ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta.trnascan.gff > final_consensus.fa.cdhitest.fa.gmap.genome.trnascanOverlap.bed
~/bin/bedtools2/bin/bedtools intersect -a final_consensus.fa.cdhitest.fa.gmap.genome.gff -b ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta.barrnap.gff3 > final_consensus.fa.cdhitest.fa.gmap.genome.barrnapOverlap.bed
cat final_consensus.fa.cdhitest.fa.gmap.genome.barrnapOverlap.bed final_consensus.fa.cdhitest.fa.gmap.genome.trnascanOverlap.bed  | cut -f 9 | cut -d "=" -f 3 | cut -d ";" -f 1 | sort | uniq | grep -f - -v final_consensus.fa.cdhitest.fa.gmap.genome.gff > final_consensus.cleaned.fa.gmap.gff3
```

- map CTR-Seq TSS/TTS tags to transcripts
```
# cage
ln -s ~/nias/LS3309/separated_data/CAGE/*R1.fq .
for i in *fq; do; echo "perl ~/scripts/bamqc-single.pl ../Pv11_5.1.4.genemodels.mRNA.fna $i $i.cage" | qsub -l nodes=1:ppn=32 -d /home/yuki.yoshida/nias/analysis/final_gene_set/version_5.1.4/ctrseq -N $i; done;
for i in *fq; do; grep ^@ $i | perl -ple 's/[\s\;]/\t/g' | cut -f 1,3 > $i.cage.umi &; done;
for i in *R1.fq; do; python pysam_genome_rightmost_position.py $i.cage.mem.sorted.bam $i.cage.umi > $i.cage.mem.sorted.bam.position & ; done;

# run using sage libraries 

# calculate counts and fix direction of fasta file and gff file
perl calculate_cpm.pl > calculate_cpm.pl.out
# count file = calculate_cpm.pl.out
# fasta file = calculate_cpm.pl.out.fna
# gff3  file = calculate_cpm.pl.out.gff3
# store1 calculate_cpm.pl.out.tree.store
# store2 calculate_cpm.pl.out.count.store
```

- merge with previous gene set
```
# create gene annotation/seq files
cat Pv11_5.1.3.1.genemodels.CTR.gff3 calculate_cpm.pl.out.gff3 > Pv11_5.1.4.genemodels.unsorted.gff3
cat Pv11_5.1.3.genemodels.CTR.mRNA.fna calculate_cpm.pl.out.fna > Pv11_5.1.4.genemodels.mRNA.fna

~/bin/trinityrnaseq-v2.9.1/util/align_and_estimate_abundance.pl  --transcripts ../Pv11_5.1.4.genemodels.mRNA.fna --seqType fq --prep_reference --est_method RSEM --aln_method bowtie2 --thread_count 64

# map ONT raw reads to mRNA sequences
~/anaconda3/envs/talon/bin/minimap2 -ax map-ont -t 64 Pv11_5.1.4.genemodels.mRNA.fna merged.fastq | samtools view -bS -@ 64 - | samtools sort -@ 64 -o merged.fastq.mmi.mRNA.sorted.bam -
# get ont read positions from mRNA and genome mapping

# annotate transcripts with TSS peaks
~/anaconda3/envs/talon/bin/minimap2 -ax map-ont -t 64 Pv11_5.1.4.genemodels.mRNA.fna merged.fastq | samtools view -bS -@ 64 - | samtools sort -@ 64 -o merged.fastq.mmi.mRNA.sorted.bam -

samtools flagstat merged.fastq.mmi.mRNA.sorted.bam > merged.fastq.mmi.mRNA.sorted.bam.flagstat
cat merged.fastq.mmi.mRNA.sorted.bam.flagstat
107002763 + 0 in total (QC-passed reads + QC-failed reads)
78651782 + 0 secondary
641276 + 0 supplementary
0 + 0 duplicates
102673084 + 0 mapped (95.95% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

##### Alignment to non-coding-gene-less gene set Pv_5.1.3 #####
cat merged.fastq.minimap2.mRNA.sorted.bam.flagstat
105016235 + 0 in total (QC-passed reads + QC-failed reads)
76690865 + 0 secondary
615665 + 0 supplementary
0 + 0 duplicates
98829443 + 0 mapped (94.11% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
## 1% increase !!!
###############################################################

python pysam_genome_rightmost_position.py merged.fastq.mmi.mRNA.sorted.bam > merged.fastq.mmi.mRNA.sorted.bam.pysam
perl parse_ONT_mapping_sites.pl > parse_ONT_mapping_sites.pl.out
perl annotate_peaks.pl

head /home/yuki.yoshida/nias/analysis/final_gene_set/version_5.1.4/Pv11_5.1.4.1.genemodels.gff3

```

### 121 promoter
`chr_2_segment10-chr_2:25900080.0-25900748.0`
