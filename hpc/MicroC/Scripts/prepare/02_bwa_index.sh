#!/bin/bash
set -e

#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu 50G # The memory needed depends on the size of the genome
#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --array=2-2 # Put here the row/rows from the table that need to be processed

#################
#### SET UP #####
#################

# Set paths and dirs
# The table genomesTable.txt has to be already generated. (See README.md)
# first column is the genome name
# second column is the absolute path for fasta
pathToGenomesTable="$SRC/genomes/genomesTable.txt"
mkdir -p $SRC/genomes/bwaIndex/
pathToBwaIndex="$SRC/genomes/bwaIndex/__genome__"
pathToImages="$SRC/images"

#################
#### SCRIPT #####
#################

# load singularity
module load singularity

# Pull images
# bwa
wget -nc -O $pathToImages/bwa_0.7.18.sif "http://datacache.galaxyproject.org/singularity/all/bwa:0.7.18--he4a0461_1"
function bwa() {
  singularity exec $pathToImages/bwa_0.7.18.sif bwa $*
}
# samtools
wget -nc -O $pathToImages/samtools.1.11.sif "http://datacache.galaxyproject.org/singularity/all/samtools:1.11--h6270b1f_0"
function samtools() {
  singularity exec $pathToImages/samtools.1.11.sif samtools $*
}

# To index the genome are needed bwa and samtools
# Check they are properly installed
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

# Get the genome name and fasta file from the genomesTable
genome=$(cat ${pathToGenomesTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
filePathForFasta=$(cat ${pathToGenomesTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')

# Adapt pathToBwaIndex to the name of the genome:
pathToBwaIndex=${pathToBwaIndex/__genome__/${genome}}

# Index genome
if [ ! -e ${filePathForFasta}.fai ]; then
    samtools faidx "$filePathForFasta"
fi
  bwa index -p "$pathToBwaIndex" "$filePathForFasta"
else
    echo "bwa index seems to already exists. If you want to regenerate it. Please remove it before running the job."
fi

# samtools faidx $filePathForFasta
# cut -f1,2 "${pathToBwaIndex}.fai" > "${pathToBwaIndex}.genome"
# grep -E '^chr([1-9]|1[0-9]|2[0-2])\t' "${pathToBwaIndex}.genome" > "${pathToBwaIndex}.genome"
# bwa index -p "$pathToBwaIndex" "$filePathForFasta"