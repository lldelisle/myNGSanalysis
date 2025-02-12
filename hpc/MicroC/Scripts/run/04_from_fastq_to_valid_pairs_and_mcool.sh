#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome
#SBATCH --cpus-per-task 24 # This allows to speed the indexing
#SBATCH --time 3:00:00 # This depends on the size of the fasta
#SBATCH --array 1-1 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name runMicroC # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /cecilia # This directory must exists, this is where will be the error and out files. and where it starts

###################################
#### TO SET FOR EACH ANALYSIS #####
###################################

### Specify the options for your analysis:

# number of CPU to use
# Only change if you don't want to use all CPUs allocated
nbOfThreads=${SLURM_CPUS_PER_TASK}
# Which genome to map on
genome=mm39
#bin size, in kb, for the .cool file
binSizeCoolMatrix=1
# Define a test region for a pgt plot (must be inside the captured region if it is a capture)
# chr7:155000000-158000000 SHH hg38
# chr2:174800000-177800000 HOXD hg38
# chr3:65103500-68603411 Shox2 CaptureC mm39
# chr2:73150000-76150000 HoxD mm39
# chr2:73779626-75669724 HoxD mm10
testRegion="chr2:73150000-76150000"
# bins size (in kb) for the plot:
bins="5 50 500"

### Specify the paths to the directories

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$microcPilot/outputs/"
# Where fastqs are stored:
dirPathForFastq="${microcPilot}/fastq/"
# This script will use get_qc.py
# from Micro-C:
# https://raw.githubusercontent.com/dovetail-genomics/Micro-C/refs/heads/main/get_qc.py
# If it does not exists
# The python script will be put in dirPathForScripts
dirPathForScripts="$SRC/images/"
# All samples are registered into a table where
# first column is the sample name
# second column is the path of the R1 fastq relatively to dirPathForFastq
# third column is same for R2
# Alternatively second column can be SRA number but third column must be filled by anything for example also the SRA number
filePathForTable="${dirPathForFastq}/samplesFastqTable.txt"

pathToBwaIndex=$SRC/genomes/bwaIndex/$genome
filePathForFasta="$SRC/genomes/fasta/${genome}.fa.gz"
# You need to generate a 'genome' file which is made of 2 columns. First column is the chromosom name, second column is the size of the chromosome.
# You can decide to filter chromosomes at this step:
filePathForSizesForBin="$SRC/genomes/fasta/${genome}_chrNumbered.genome"

### Specify the way to deal with dependencies:
# Here we use singularity

# Images
pathToImages="$SRC/images"

wget -nc -O $pathToImages/bwa_0.7.18.sif "http://datacache.galaxyproject.org/singularity/all/bwa:0.7.18--he4a0461_1"
function bwa() {
  singularity exec $pathToImages/bwa_0.7.18.sif bwa $*
}

wget -nc -O $pathToImages/pairtools.0.3.0 "http://datacache.galaxyproject.org/singularity/all/pairtools:0.3.0--py37h4eba2af_0"
function pairtools() {
  singularity exec $pathToImages/pairtools.0.3.0 pairtools $*
}

wget -nc -O $pathToImages/samtools.1.11.sif "http://datacache.galaxyproject.org/singularity/all/samtools:1.11--h6270b1f_0"
function samtools() {
  singularity exec $pathToImages/samtools.1.11.sif samtools $*
}
function bgzip() {
  singularity exec $pathToImages/samtools.1.11.sif bgzip $*
}

wget -nc -O $pathToImages/tabulate:0.7.5--py36_0 "https://depot.galaxyproject.org/singularity/tabulate:0.7.5--py36_0"
function python() {
  singularity exec $pathToImages/tabulate:0.7.5--py36_0 python $*
}
function tabulate() {
  singularity exec $pathToImages/tabulate:0.7.5--py36_0 tabulate $*
}

wget -nc -O $pathToImages/cooler.0.10.3 "https://depot.galaxyproject.org/singularity/cooler:0.10.3--pyhdfd78af_0"
function cooler() {
  singularity exec $pathToImages/cooler.0.10.3 cooler $*
}
function pairix() {
  singularity exec $pathToImages/cooler.0.10.3 pairix $*
}

wget -nc -O $pathToImages/pygenometracks.3.9 "https://depot.galaxyproject.org/singularity/pygenometracks:3.9--pyhdfd78af_0"
function pgt() {
  singularity exec $pathToImages/pygenometracks.3.9 pgt $*
}
function pyGenomeTracks() {
  singularity exec $pathToImages/pygenometracks.3.9 pyGenomeTracks $*
}
# Give access to path with index with fastqs etc
export APPTAINER_BIND=$SRC

# python QC script
wget -nc -O $dirPathForScripts/get_qc.py https://raw.githubusercontent.com/dovetail-genomics/Micro-C/refs/heads/main/get_qc.py

# Check installations 
v=$(bwa 2>&1)
if [[ "$v" = *"command not found" ]]
then
  echo "Bwa is not installed but required. Please install it"
  exit 1
fi
echo $v

samtools --version
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it"
  exit 1
fi

bgzip --version
if [ $? -ne 0 ]
then
  echo "bgzip is not installed but required. Please install it"
  exit 1
fi

pairtools --version
if [ $? -ne 0 ]
then
  echo "pairtools is not installed but required. Please install it"
  exit 1
fi

tabulate --version
if [ $? -ne 2 ]
then
  echo "tabulate is not installed but required. Please install it"
  exit 1
fi

python -c "import argparse;print(argparse.__version__)"
if [ $? -ne 0 ]
then
  echo "argparse is not installed but required. Please install it"
  exit 1
fi

cooler --version
if [ $? -ne 0 ]
then
  echo "cooler is not installed but required. Please install it"
  exit 1
fi

pairix --version
if [ $? -ne 0 ]
then
  echo "pairix is not installed but required. Please install it"
  exit 1
fi

pgt --version
if [ $? -ne 0 ]
then
  echo "pyGenomeTracks is not installed but required. Please install it"
  exit 1
fi


# Get the gensampleome name and fastq file from the table
sample=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
relFilePathFastqR1=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
relFilePathFastqR2=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $3}')

# The directory is created (if not existing)
pathResults=${dirPathWithResults}/${sample}/
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo ${sample}

# The analysis part takes part within the pathResults
cd ${pathResults}

# Generate SAM file
if [ ! -e ${sample}.sam ]; then
  if [ ! -e ${dirPathForFastq}/${relFilePathFastqR1} ]; then
    # If the fastq does not exists we assume it was an SRA ID
    mkdir -p ${dirPathForFastq}
    cd ${dirPathForFastq}
    # Write version to stdout:
    fasterq-dump --version
    if [ $? -ne 0 ]
    then
      echo "fasterq-dump is not installed and fastqFile not found so assumed it was a SRA ID.
Please install it for example in the conda environment (sra-tools>=2.11)."
      exit 1
    fi
    fasterq-dump -o ${sample}.fastq ${relFilePathFastqR1}
    if [ ! -s ${sample}_1.fastq ]; then
        echo "FASTQ R1 IS EMPTY"
        exit 1
    fi
    gzip ${sample}_1.fastq
    gzip ${sample}_2.fastq
    cd $pathResults
    relFilePathFastqR1=${sample}_1.fastq.gz
    relFilePathFastqR2=${sample}_2.fastq.gz
  fi
  if [ ! -s ${dirPathForFastq}/${relFilePathFastqR1} ]; then
    echo "FASTQ R1 IS EMPTY"
    exit 1
  fi
  # -5 is for split alignemnent, takes the alignemnte of the 5' read as primary. -S skips mate rescue, -P skips pairing. -T sets the minimus mapping quality (we want all reads to compute the stats)
  bwa mem -5SP -T0 -t${nbOfThreads} $pathToBwaIndex ${dirPathForFastq}/${relFilePathFastqR1} \
    ${dirPathForFastq}/${relFilePathFastqR2} > ${sample}.sam
else
  echo "${sample}.sam already exists"
fi

# Record valid ligation events.
# --min-mapq is the Mapping quality threshold for defining an alignment as a multi-mapping alignment. 
# --walks-policy is to handle multi mapping alignements.
# 5unique is used to report the 5’-most unique alignment on each side, if present (one or both sides may map to different locations on the genome, producing more than two alignments per DNA molecule)
if [ ! -e ${sample}.parsed.pairsam ]; then
  pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${nbOfThreads} --nproc-out ${nbOfThreads} --chroms-path "$filePathForSizesForBin" "${sample}.sam" > "${sample}.parsed.pairsam"
else
  echo "${sample}.parsed.pairsam already exists"
fi

# Sort the ${sample}.parsed.pairsam.
if [ ! -e ${sample}.sorted.pairsam ]; then
  pairtools sort --nproc ${nbOfThreads} --tmpdir=$(mktemp -d) "${sample}.parsed.pairsam" > "${sample}.sorted.pairsam"
else
  echo "${sample}.sorted.pairsam already exists"
fi

# Remove PCR duplicates
# duplicate pairs are marked as DD in “pair_type” and as a duplicate in the sam entries.
if [ ! -e ${sample}.dedup.pairsam ]; then
  pairtools dedup --nproc-in ${nbOfThreads} --nproc-out ${nbOfThreads} --mark-dups --output-stats "stats.txt" --output "${sample}.dedup.pairsam" "${sample}.sorted.pairsam"
else
  echo "dedup.pairsam already exists"
fi

# Generate .pair and final BAM file
if [ ! -e ${sample}.mapped.pairs ]; then
  pairtools split --nproc-in ${nbOfThreads} --nproc-out ${nbOfThreads} --output-pairs "${sample}.mapped.pairs" --output-sam "${sample}.unsorted.bam" "${sample}.dedup.pairsam"
else
  echo "${sample}.mapped.pairs already exists"
fi

## Sort and index the final BAM file
# Can be used to generate a coverage
# samtools sort -@${nbOfThreads} -T temp/temp.bam -o "${sample}.mapped.PT.bam" "${sample}unsorted.bam" 
# samtools index "$sample_output_dir/${sample}.mapped.PT.bam"

# Run QC script
python $dirPathForScripts/get_qc.py -p "${sample}.stats.txt"  > ${sample}.pretty.stats.txt

# # Save the stats in a common file for all samples
# mkdir $microc/output/stats_all/
# touch $microc/output/stats_all/stats_all.txt
# echo -e "$sample" >> $microc/output/stats_all/stats_all.txt
# cat "${sample}.stats.txt" >> "$microc/output/stats_all/"

# Generate contact matrix for .pairs with cooler
# Bgzip the pairs:
if [ ! -e ${sample}.pairs.gz ]; then
  bgzip -c "${sample}.mapped.pairs" > "${sample}.pairs.gz"
else
  echo "${sample}.pairs.gz already exists"
fi
# Index them
if [ ! -e ${sample}.pairs.gz.px2 ]; then
  pairix "${sample}.pairs.gz"
else
  echo "${sample}.pairs.gz already indexed"
fi
# Create the cooler at smaller resolution
if [ ! -e ${sample}_raw.${binSizeCoolMatrix}kb.cool ]; then
  cooler cload pairix -p 16 ${filePathForSizesForBin}:${binSizeCoolMatrix}000 "${sample}.pairs.gz" ${sample}_raw.${binSizeCoolMatrix}kb.cool
else
  echo "${sample}_raw.${binSizeCoolMatrix}kb.cool already exists"
fi
# Zoomify
if [ ! -e ${sample}.mcool ]; then
   singularity exec $pathToImages/cooler.0.10.3 cooler zoomify --balance -p ${nbOfThreads} --resolutions ${binSizeCoolMatrix}000N -o ${sample}.mcool --balance-args '--nproc 24 --cis-only' ${nb}${sample}_raw.${binSizeCoolMatrix}kb.cool
else
  echo "${sample}.mcool already exists"
fi

# Plot generation
ini_file=${sample}.ini
echo "[x-axis]" > $ini_file
for bin in $bins; do
  echo "[${sample}_${bin}kb]
file = ${sample}.mcool::/resolutions/${bin}000
depth = 2000000
min_value = 0
title = ${sample}_${bin}kb
file_type = hic_matrix
[spacer]
" >> ${ini_file}
done
echo "[x-axis]" >> ${ini_file}
# Generate a basic plot on the testRegion
pgt --tracks ${ini_file} --region ${testRegion} --fontSize 6 -o ${ini_file/.ini/_testRegion.pdf}

# Copy all final files to specific directories
mkdir -p ${dirPathWithResults}/allFinalFiles/cool
cp *.mcool ${dirPathWithResults}/allFinalFiles/cool/
mkdir -p ${dirPathWithResults}/allFinalFiles/reports
cp *.stats.txt ${dirPathWithResults}/allFinalFiles/reports/
mkdir -p ${dirPathWithResults}/allFinalFiles/pairs
cp ${sample}.pairs.gz ${dirPathWithResults}/allFinalFiles/pairs/
mkdir -p ${dirPathWithResults}/allFinalFiles/visualisation
cp *.pdf ${dirPathWithResults}/allFinalFiles/visualisation
