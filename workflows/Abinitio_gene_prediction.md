# Ab initio gene predictions for all dipteran genomes

The following script was names "run_gene_prediction.sh" and was used as `sh run_gene_prediction.sh ***.fasta` using the non-masked genome from ENSEMBL, NCBI, or other public genomes.

```
#!/bin/sh
set -e

# Get variables
FASTA=$1
GETWD="/home/yuki.yoshida/nias/analysis/reanalysis/18_revice/MOSGA/non-sm"

# Start conda environment
source /home/yuki.yoshida/anaconda3/etc/profile.d/conda.sh
conda activate repeat

# Set Enviromnet variables
export PATH="/home/yuki.yoshida/nias/bin/gmes_linux_64:/home/yuki.yoshida/bin/NINJA-0.95-cluster_only/NINJA:$PATH"
export PROTHINT_PATH="/home/yuki.yoshida/nias/bin/gmes_linux_64/ProtHint/bin"
export AUGUSTUS_CONFIG_PATH=/home/yuki.yoshida/anaconda3/envs/pv_genome_revise/config/


# Start Pipeline
echo "### Creating working directory "
echo "# mkdir braker_${FASTA%.dna_sm.toplevel.fa}"
mkdir -p braker_${FASTA%.dna_sm.toplevel.fa}
cd braker_${FASTA%.dna_sm.toplevel.fa}
ln -fs ${GETWD}/${FASTA} .

echo "# Current location : " `pwd`

# Run repeat masking
echo "### Run repeat masking"
echo "## Run BuildDatabase"
/home/yuki.yoshida/bin/RepeatModeler-2.0.2a/BuildDatabase -name ${FASTA%.dna_sm.toplevel.fa} -engine ncbi ${FASTA} &> log.buildDatabase
echo "## Run RepeatModeler"
/home/yuki.yoshida/bin/RepeatModeler-2.0.2a/RepeatModeler -engine ncbi -pa 16 -database ${FASTA%.dna_sm.toplevel.fa} &> log.repeatmodeler
#/home/yuki.yoshida/bin/RepeatModeler-2.0.2a/RepeatModeler -engine ncbi -pa 16 -database ${FASTA%.dna_sm.toplevel.fa} -LTRStruct &> log.repeatmodeler
echo "## Run RepeatMasker 4.1.2-p1"
RepeatMasker -e rmblast -pa 32 -lib ${FASTA%.dna_sm.toplevel.fa}-families.fa -xsmall -html -gff -xm ${FASTA} &> log.repeatmasker

MASKED=${FASTA}.masked
conda deactivate

echo "## Finished Repeat masking"
echo "## Result is ${MASKED}"


# Change to braker environment
conda activate pv_genome_revise

# Run braker
echo "### Running gene prediction"
echo "### Run Braker2"
## Braker v2.1.6
## Augustus v3.4.0
## GeneMark-ES Suite version 4.58_lic
## samtools 1.12  Using htslib 1.12
## bamtools 2.5.1
## diamond version 2.0.13

braker.pl --genome ${MASKED} --prot_seq ${GETWD}/proteins.fasta --softmasking --cores 48 --species ${FASTA%.dna_sm.toplevel.fa} --gff3 &> log.braker

AASEQ="braker/augustus.hints.aa"
conda deactivate

echo "## Finished gene prediction"
```
