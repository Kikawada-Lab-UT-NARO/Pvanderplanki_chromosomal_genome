# Analysis of CTR-Seq ONT reads

## Filtering for adaptor ligated reads and fix strands
Identify reads with Nextera adaptors on both ends
```
porechop -o 20180510_0216_20180510-CTRS-Oleg-1.raw.fastq.porechop.fastq -i 20180510_0216_20180510-CTRS-Oleg-1.raw.fastq -t 64 --check_reads 500 --end_size 200 -v 3 --discard_middle --extra_end_trim 0  > 20180510_0216_20180510-CTRS-Oleg-1.raw.fastq.porechop.log
```

Filter reads with `scripts/get_full_length_from_porechop_log.pl`
```
for i in 20180510_0216_20180510-CTRS-Oleg-1.raw.fastq 20180515_0328_20180515-CTRS-Oleg-2.raw.fastq 20180518_0759_20180518-CTRS-Oleg-3.raw.fastq 20180522_0743_20180522-CTRS-Oleg-4.raw.fastq 20180525_0708_20180525-CTRS-Oleg-5.raw.fastq 20180528_0710_20180528-CTRS-Oleg-6.raw.fastq 20180530_0304_20180530-CTRS-Oleg-7.raw.fastq;
do
  perl get_full_length_from_porechop_log.pl $i.porechop.log > $i.porechop.fastq
done;
```
