#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome
#SBATCH --cpus-per-task 24 # This allows to speed the indexing
#SBATCH --time 3:00:00 # This depends on the size of the fasta
#SBATCH --array=2-2 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name star_index # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /home/users/d/delislel/scratch/alphaTC1/ # This directory must exists, this is where will be the error and out files
## Specific to baobab:
#SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days
#SBATCH --account herrerap

##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
# number of CPU to use
# Only change if you don't want to use all CPUs allocated
nbOfThreads=${SLURM_CPUS_PER_TASK}

### Specify the paths
# All genomes are registered into a table where
# first column is the genome name
# second column is the absolute path for fasta
filePathForTable="$HOME/scripts/alphatc1-clone-9/prepare/genomes_table.txt"

# STAR index may be incompatible between versions
# Therefore it may be interested to keep a record of the STAR version used
# Substitute the name of your genome by __genome__
basenamePathForB2Index="$PWD/genomes/STAR_2.7.11b/__genome__"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load GCC/8.3.0
# module load STAR/2.7.4a
# or
# module purge
# module load Miniconda3/4.9.2
# or use from home:
# export PATH=$PATH:/home/ldelisle/softwares/STAR-2.7.9a/bin/Linux_x86_64/
# You can choose to use a conda environment to solve star dependencies
# Comment it if you will use module load
# condaEnvName=rna202210

# You can choose to use singularity, then you need to define each function:
pathToImages=/cvmfs/singularity.galaxyproject.org/all/
# I don't know why on bamboo on clusters it does not work...
# pathToImages=/home/users/d/delislel/scratch/images/
if [ ! -e "$pathToImages/star:2.7.11b--h5ca1c30_5" ]; then
    wget -nc -O "$pathToImages/star:2.7.11b--h5ca1c30_5" "http://datacache.galaxyproject.org/singularity/all/star:2.7.11b--h5ca1c30_5"
fi
function STAR() {
    singularity exec "$pathToImages/star:2.7.11b--h5ca1c30_5" STAR $*
}

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

v=$(STAR --version)
if [ $? -ne 0 ]
then
  echo "STAR is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
echo "STAR version $v"

# Get the genome name and fasta file from the table
genome=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
filePathForFasta=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')

# Adapt dirPathForSTARIndex to the name of the genome:
dirPathForSTARIndex=${dirPathForSTARIndex/__genome__/${genome}}

if [ ! -e ${dirPathForSTARIndex}/chrLength.txt ]; then
    STAR --runMode genomeGenerate --genomeFastaFiles ${filePathForFasta} \
        --genomeDir ${dirPathForSTARIndex} --runThreadN $nbOfThreads
else
    echo "STAR index seems to already exists. If you want to regenerate it. Please remove it before running the job."
fi
