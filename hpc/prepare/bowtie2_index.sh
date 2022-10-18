#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome
#SBATCH --cpus-per-task 24 # This allows to speed the indexing
#SBATCH --time 3:00:00 # This depends on the size of the fasta
#SBATCH --array=1-2 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name bowtie2_index # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/ # This directory must exists, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

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
filePathForTable="/home/ldelisle/scripts/scitas_sbatchhistory/2022/20220921_test_index/table_genomes.txt"

# Substitute the name of your genome by __genome__
basenamePathForB2Index="/home/ldelisle/genomes/bowtie2/__genome__"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load GCC/8.3.0
# module load Bowtie2/2.3.5.1
# or
# module purge
# module load GCC/10.3.0
# module load Bowtie2/2.4.4
# or
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve bowtie2 dependencies
# Comment it if you will use module load
condaEnvName=hic202209

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

bowtie2 --version
if [ $? -ne 0 ]
then
  echo "Bowtie2 is not installed but required. Please install it from github (just untar)"
  exit 1
fi

# Get the genome name and fasta file from the table
genome=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
filePathForFasta=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')

# Adapt basenamePathForB2Index to the name of the genome:
basenamePathForB2Index=${basenamePathForB2Index/__genome__/${genome}}

if [ ! -e ${basenamePathForB2Index}.rev.2.bt2 ]; then
    bowtie2-build --thread ${nbOfThreads} ${filePathForFasta} ${basenamePathForB2Index}
    ln -s ${filePathForFasta} ${basenamePathForB2Index}.fa
else
    echo "bowtie2 index seems to already exists. If you want to regenerate it. Please remove it before running the job."
fi
