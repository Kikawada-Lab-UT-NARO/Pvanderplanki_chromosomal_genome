# Protocol for analyzing CTR-Seq Illumina data
This workflow was used to obtain and quantify TSS and TTS peaks for detected genes.

# hiseq data
## data merge
The input sequence data was seperated by lanes so I first merged each lane into one file
```
for i in `ls | cut -d "_" -f 1-6| sort | uniq` ;
do
  echo $i
  cat ${i}_L007_R2_001.fq ${i}_L008_R2_001.fq > ${i}_merged_R2.fq &
done;

for i in `ls | cut -d "_" -f 1-6| sort | uniq` ;
do
  echo $i
  cat ${i}_L007_R1_001.fq ${i}_L008_R1_001.fq > ${i}_merged_R1.fq &
done;
```

## hisat2 mapping

Mapped the CAGE and SAGE reads to the genome with hisat2

```
for i in `ls | grep merged | cut -d "_" -f 1-6 | sort | uniq`;
do
  echo $i
  hisat2 -p 64 --known-splicesite-infile ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.gtf.splicesites -x ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -1  ${i}_merged_R1.fq -2 ${i}_merged_R2.fq -S ${i}_merged.cage.hisat2.sam
done

```
Create bigwig files for visualization following https://github.com/suimye/cage_tutorial/blob/master/cage.counting.pipeline.b0.01.sh

```
for i in *bam; 
do
  samtools view -@ 8 -F 4 -u -q 20 $i | genomeCoverageBed -ibam /dev/stdin -5 -bg -strand + | sort -k1,1 -k2,2n > $i.gcb.fwd.txt &
done;

for i in *fwd*txt ;
do
  bedGraphToBigWig $i ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta.fai $i.ctss.fwd.bw &
 done

for i in *bam ; 
do
  samtools view -@ 8 -F 4 -u -q 20 $i | genomeCoverageBed -ibam /dev/stdin -5 -bg -strand - | sort -k1,1 -k2,2n > $i.gcb.rev.txt &
done

for i in *rev*txt ; 
do
  bedGraphToBigWig $i ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta.fai $i.ctss.rev.bw &
done
```

# Homer
We used Homer for peak detection and quantification.

```
# make tag dir for each cage-seq file
for i in `\ls | grep RN | cut -d "." -f 1-2`
do;
  echo $i
  makeTagDirectory ${i}_tagDir ${i}.hisat2.sorted.bam -genome ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -checkGC -fragLength 150 &
 done;

# find csRNATSS for each cage-seq file
for i in `\ls | grep RN | grep -v tagDir | cut -d "." -f 1-2`; 
do;
  echo $i
  findcsRNATSS.pl ${i}_tagDir/ -o ${i}_tagDir_findcsRNATSS -gtf ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.gtf -genome ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta &
done

# CAGE
## merge TSS into one file
mergePeaks *cage*.tss.txt -strand > cage_tss_peaks_merged.txt

## annotation of peaks
annotatePeaks.pl cage_tss_peaks_merged.txt ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -strand + -fragLength 1 -cpu 64 -gtf ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.gtf -d *cage_tagDir > cage_tss_peaks_merged.txt.annotated

# SAGE
## merge TTS into one file
mergePeaks *sage*.tss.txt -strand > sage_tss_peaks_merged.txt

## annotation of peaks
annotatePeaks.pl sage_tss_peaks_merged.txt ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -strand - -fragLength 1-cpu 64 -gtf ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.gtf -d *sage*tagDir > sage_tss_peaks_merged.txt.annotated

# both
annotatePeaks.pl both_csage_merged_peaks.txt ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -strand - -fragLength 1 -cpu 64 -gtf ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.gtf -d *tagDir > both_csage_merged_peaks.txt.annotated
```

Identification of DNA motifs from both peaks
```
findMotifsGenome.pl both_csage_merged_peaks.txt ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta both_csage_merged_peaks.txt.findMotifs -size -150,50 -p 64 -cache 20000

annotatePeaks.pl both_csage_merged_peaks.txt ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -strand - -fragLength 1 -cpu 64 -gtf ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.gtf -d *tagDir -m both_csage_merged_peaks.txt.findMotifs/knownResults.motif both_csage_merged_peaks.txt.findMotifs/homerMotifs.all.motifs -mbed both_csage_merged_peaks.txt.motif.annotated.bed > both_csage_merged_peaks.txt.motif.annotated.out
```

and TSS/TTS peaks seperately.

```
# find DNA motifs
## CAGE
### find all motifs
findMotifsGenome.pl cage_tss_peaks_merged.txt.annotated ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta findMotifs -size -150,50 -p 64 -cache 20000

annotatePeaks.pl cage_tss_peaks_merged.txt ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -size 200 -cpu 64 -m ./cage_tss_peaks_merged.txt.findMotifs/homerMotifs.all.motifs ./cage_tss_peaks_merged.txt.findMotifs/knownResults.motif -mbed cage_tss_peaks_merged.txt.annotated.homerALL.bed -strand + -fragLength 1 -cpu 64 -gtf ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.gtf -d *cage_tagDir > cage_tss_peaks_merged.txt.annotated.homerALL.out

#### find promoter regions with HSF motif s
cut -f 1-28,115,126,213,282,301,313,332 cage_tss_peaks_merged.txt.annotated.homerALL.out > cage_tss_peaks_merged.txt.annotated.homerALL.out.hsf
perl -lane 'print if @F[29] ne "" && @F[30] ne "" && @F[31] ne "" && @F[32] ne "" && @F[33] ne "" && @F[34] ne "" && @F[35] ne ""' cage_tss_peaks_merged.txt.annotated.homerALL.out.hsf | cut -f 11 | perl ~/scripts/get_blastp_from_stdin.pl ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.genemodels.faa.blastp.swissprot.1e-15.sorted | wc -l


### novel motifs only
annotatePeaks.pl cage_tss_peaks_merged.txt.annotated ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -size 200 -cpu 64 -m ./cage_tss_peaks_merged.txt.findMotifs/homerMotifs.all.motifs -mbed cage_tss_peaks_merged.txt.annotated.homerMotifs_all.bed > cage_tss_peaks_merged.txt.annotated.homerMotifs_all.out

### known motifs only 
annotatePeaks.pl cage_tss_peaks_merged.txt.annotated ~/nias/Pv_Final_assembly_2018_Olga/final_version_from_olga/Pv11_5.0.final.fasta -size 200 -m ./cage_tss_peaks_merged.txt.findMotifs/knownResults.motif -mbed cage_tss_peaks_merged.txt.annotated.knownMotifs.bed > cage_tss_peaks_merged.txt.annotated.knownMotifs.out
```

