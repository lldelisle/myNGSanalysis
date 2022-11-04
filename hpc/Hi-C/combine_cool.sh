#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1 # This only uses a single CPU
#SBATCH --time 12:00:00 # This depends on the number of cool files
#SBATCH --array=1-14 # Put here the rows from the table that need to be processed
#SBATCH --job-name cHiC_GW_merge # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/cHi-C_GW/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

# This script run hicup with Bowtie2
# convert hicup bam to juicebox format valid pairs
# If needed will subset to pairs both in capture region
# Will sort the pairs to build matrices in cool format
# Balance with cooler
# Do some plots with pyGenomeTracks


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the paths

# Put in dirPathWithCool the directory
# where the output cool file will be
dirPathWithCool="$PWD/allFinalFiles/cool/"
# All samples to merge are registered into a table where
# first column is the merged cool name
# second column is the path of the RAW cool files to sum
# this is relative to the dirPathWithCool
filePathForTable="/home/ldelisle/scripts/scitas_sbatchhistory/2022/20221104_gastruloids_cHiC_GW/cHi-C_merged_table.txt"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module load Miniconda3/4.9.2
#

# We will use a conda environment to solve cooler and HiCExplorer dependencies
# We create a conda environment with hicexplorer as cooler
# Is a  dependency of hicexplorer
# You need to install it before running the script
# To create it simply load all the modules required to get conda on your instance
# Then choose a name for your conda environment (here pgt3.7)
# and then:
# conda create -y -n pgt3.7 -c bioconda -c conda-forge pygenometracks=3.7 hicexplorer=3.7.2
# condaEnvName=pgt3.7
# Alternatively you can use conda to solve all dependencies:
# conda create -y -n hic202209 -c bioconda -c conda-forge pygenometracks hicup hicexplorer=3.7.2
condaEnvName=hic202209


##################################
####### BEGINING OF SCRIPT #######
##################################

# Check everything is set correctly:
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
# Check all softwares are present and write version to stdout:
cooler --version
# Check cooler is installed:
if [ $? -ne 0 ]
then
  echo "Cooler is not installed but required. Please install it for example in the conda environment"
  exit 1
fi
hicSumMatrices --version
# Check hicSumMatrices is installed:
if [ $? -ne 0 ]
then
  echo "HiCExplorer is not installed but required. Please install it for example in the conda environment"
  exit 1
fi

# Get the output name and input files from the table
output=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
inputs_comma=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')

cd ${dirPathWithCool}
hicSumMatrices --matrices $(echo $inputs_comma | tr "," " ") -o ${output}
cp ${output} ${output/.cool/_raw.cool}
cooler balance --mad-max 5 --min-nnz 10 --min-count 0 --ignore-diags 2 --tol 1e-05 --max-iters 200 --cis-only -f ${output}
