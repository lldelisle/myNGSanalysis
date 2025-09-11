
# Introduction
Here you find the scripts to run a Micro-C analysis. 

To handle dependencies singularity images have been used from https://depot.galaxyproject.org/singularity/.

The pipeline follows https://micro-c.readthedocs.io/en/latest/index.html.

## Set up
### General
For github set up: generate a ssh key on the machine with ssh-keygen and put the public key on github.

```bash
cd $HOME 
git clone https://github.com/lldelisle/myNGSanalysis/tree/main/hpc/MicroC
```

Export variables for dirs in the ./bashrc.
Adapt it. 
```bash
echo '## FOR MICROC' >> $HOME/.bashrc
echo 'export SRC=/travail/ijerkovic/NGS' >> $HOME/.bashrc
echo 'export microc=$SRC/microc' >> $HOME/.bashrc 
echo 'export microcPilot=$SRC/microc/pilot' >> $HOME/.bashrc 
echo 'export microcFullData=$SRC/microc/fullData' >> $HOME/.bashrc 
echo 'export PREP=$HOME/Micro-C_pipeline/Singularity/prepare' >> $HOME/.bashrc
echo 'export RUN=$HOME/Micro-C_pipeline/Singularity/run' >> $HOME/.bashrc
source $HOME/.bashrc 
```

Create all the needed directories.
```bash
mkdir -p $SRC
mkdir -p $microc
mkdir -p $microcPilot
mkdir -p $microcFullData
mkdir -p $SRC/genomes
mkdir -p $SRC/images
mkdir -p $microcPilot/outputs
```

Make scripts executable
```bash
chmod +x $PREP/01.1_get_genome.sh
chmod +x $PREP/01.2_get_fastq.sh
chmod +x $PREP/03_fastq_table.sh
chmod +x $PREP/02_bwa_index.sh
chmod +x $RUN/04_from_fastq_to_valid_pairs_and_mcool.sh
```

### Create reference genome table
Generate a reference genome table where the first column is the genome name in the format hg38.fa.gz, the second column is the path to the fasta file, the third is http of the fastq to be downloaded.
```bash
echo -e "hg38\t$SRC/genomes/fasta/hg38.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz\nmm39\t$SRC/genomes/fasta/mm39.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz" > $SRC/genomes/genomesTable.txt
```

### Download data and genome

```bash
bash $PREP/01.1_get_genome.sh
```
Here the genome is downloaded. 
The 'genome' file is generated where the first column is chromosome name and the second column the dimension of the chrnomosome. 
Also a 'genome' file with only the numbered chromosomes is created.
```bash
bash $PREP/01.2_get_fastq.sh
```

### Index the genome
The genome is indexed with bwa.

In the script, modify the SBATCH --array=1-1 where 1-1 is the interval of rows to process in the table. Put 1-1 if you wan the hg38 genome. Put 2-2 if you want the mm39 genome
```bash
sbatch --chdir $SRC $PREP/02_bwa_index.sh
```

### Create samples fastq reference table
Generate tables for the samples sequencing data.
In the samplesFastqTable: first column is the sample name, second column is the fastq1 path, third column is the fastq2 path.
CHECK: fastq names have to end in '1.fq.gz' (for read 1), '2.fq.gz' (for read 2)

See example of the tabel in Scripts/prepare/samples_fastq_table.txt

```bash
bash $PREP/03_fastq_table.sh
```

## MicroC analysis

This script needs the 'genome' file, index file and reference fasta file.
For details check https://micro-c.readthedocs.io/en/latest/fastq_to_bam.html.

```bash
sbatch --chdir $SRC $RUN/04_from_fastq_to_valid_pairs_and_mcool.sh
```


## Compare with public data
### Get the fastq
Download the fastqs. If the public data are more deeply sequenced then your data, subset a number of reads equivalent to your depth. Then, trim to 50bp (improves alignement). Add the infos in the FastqTable.
First get the fastq:
```bash
pathToImages="$SRC/images"
wget -nc -O "$pathToImages/sra-tools:3.1.1--h4304569_2"  "https://depot.galaxyproject.org/singularity/sra-tools:3.1.1--h4304569_2"
function fasterq-dump() {
  singularity exec "$pathToImages/sra-tools:3.1.1--h4304569_2" fasterq-dump $*
}
export APPTAINER_BIND=$SRC
dirPathForFastq="${microcPilot}/fastq/"
cd $dirPathForFastq
sample=SRR29294642
fasterq-dump -o ${sample}.fastq ${sample}
gzip ${sample}_1.fastq
gzip ${sample}_2.fastq
```
Get only 50 million reads:
```bash
for r in 1 2; do
    zcat ${sample}_${r}.fastq.gz | head -n 200000000 | gzip > ${sample}_50M_${r}.fastq.gz
done
```
Trim to 50bp
```bash
toolversion="seqtk:1.4--h577a1d6_3"
wget -nc -O "$pathToImages/${toolversion}"  "https://depot.galaxyproject.org/singularity/${toolversion}"
function seqtk() {
  singularity exec "$pathToImages/${toolversion}" seqtk $*
}
for r in 1 2; do
    seqtk trimfq -l 50 ${sample}_50M_${r}.fastq.gz > ${sample}_50M_50bp_${r}.fastq.gz
done
```
Add these 2 to the table with fastqs:
```bash
echo -e "${sample}_50M_full\t${sample}_50M_1.fastq.gz\t${sample}_50M_2.fastq.gz
${sample}_50M_50bp\t${sample}_50M_50bp_1.fastq.gz\t${sample}_50M_50bp_2.fastq.gz" >> samplesFastqTable.txt
```

Given that you followed all the steps of the set-up you can run the MicroC analysis.

## Refine outputs
### Aggregate all reports in one
```bash
pathToImages="$SRC/images"
toolversion="multiqc:1.26--pyhdfd78af_0"
wget -nc -O "$pathToImages/${toolversion}"  "https://depot.galaxyproject.org/singularity/${toolversion}"
function multiqc() {
  singularity exec "$pathToImages/${toolversion}" multiqc $*
}
export APPTAINER_BIND=$SRC
cd $microcPilot/outputs
multiqc . -m pairtools --force
```
### Make a general plot
```bash
pathToImages="$SRC/images"
toolversion="pygenometracks:3.9--pyhdfd78af_0"
wget -nc -O "$pathToImages/${toolversion}"  "https://depot.galaxyproject.org/singularity/${toolversion}"
function pgt() {
  singularity exec "$pathToImages/${toolversion}" pgt $*
}
export APPTAINER_BIND=$SRC
cd $microcPilot/outputs
bin=50
testRegion="chr2:73150000-76150000"
ini_file=all_${bin}kb.ini
echo "[x-axis]" > $ini_file
for mcool in allFinalFiles/cool/*; do
    echo "[$mcool]
file = ${mcool}::/resolutions/${bin}000
depth = 2000000
min_value = 0
title = $(basename ${mcool/.mcool/})_${bin}kb
file_type = hic_matrix
[spacer]
" >> ${ini_file}
done
echo "[x-axis]" >> ${ini_file}
# Generate a basic plot on the testRegion
pgt --tracks ${ini_file} --region ${testRegion} --fontSize 6 -o ${ini_file/.ini/_testRegion.pdf}
```