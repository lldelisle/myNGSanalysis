#!/bin/bash

#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 4G # The memory needed depends on the size of the genome
#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --array=2-2 # Put here the row/rows from the table that need to be processed
#SBATCH --job-name getFasta # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /home/users/d/delislel/scratch/alphaTC1/ # This directory must exist, this is where will be the error and out files
## Specific to UNIGE:
#SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days
#SBATCH --account herrerap

# This script get fasta and bgzip it

# Adapted from a script from Cecilia Carmignoto

##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the paths
# All genomes are registered into a table where
# first column is the genome name
# second column is the path where the fasta should end
# third column is the URL
filePathForTable=$HOME/scripts/alphatc1-clone-9/prepare/genomes_table.txt

### Specify the way to deal with dependencies:

# You can choose to use singularity, then you need to define each function:
# pathToImages=/cvmfs/singularity.galaxyproject.org/all/
# I don't know why on bamboo on clusters it does not work...
pathToImages=/home/users/d/delislel/scratch/images/
mkdir -p $pathToImages
if [ ! -e "$pathToImages/samtools:1.11--h6270b1f_0" ]; then
    wget -nc -O "$pathToImages/samtools:1.11--h6270b1f_0" "http://datacache.galaxyproject.org/singularity/all/samtools:1.11--h6270b1f_0"
fi
function samtools() {
    singularity exec "$pathToImages/samtools:1.11--h6270b1f_0" samtools $*
}
# bgzip from samtools
function bgzip() {
    singularity exec "$pathToImages/samtools:1.11--h6270b1f_0" bgzip $*
}

# And additional binds:
export APPTAINER_BIND=$HOME/scratch/,$(realpath $HOME/scratch)


##################################
####### BEGINING OF SCRIPT #######
##################################

# Check everything is set correctly:
if [ ! -z ${condaEnvName} ]; then
    # Conda environment:
    # This line is to adapt the conda to the shell
    source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
    # We check if the conda environment exists
    exists=$(conda info --envs | awk -v ce=${condaEnvName} '$1==ce{print}' | wc -l)
    # It if does not exists an error is raised
    if [ $exists -ne 1 ]; then
    echo "conda environment ${condaEnvName} does not exists. Create it before."
    exit 1
    fi
    # Activate the conda environment
    conda activate ${condaEnvName}
fi

# Check all softwares are present and write version to stdout:
samtools --version
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it"
  exit 1
fi

# Get the genome name, the genome path and fasta file from the genomesTable.txt
genomeName=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
filePathForFasta=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
genomeURL=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $3}')

# Download genome
mkdir -p $(dirname $filePathForFasta)
wget -nc -O $filePathForFasta $genomeURL
gunzip -c $filePathForFasta > genomeUnzipped.fa
bgzip genomeUnzipped.fa 
mv genomeUnzipped.fa.gz $filePathForFasta

# Generate the genome file only with numbered chromosomes:
if [ ! -e ${filePathForFasta}.fai ]; then
  samtools faidx "$filePathForFasta"
fi
