#!/bin/bash

#SBATCH -o slurm-%x-%A.out # Template for the std output of the job uses the job name and the job id
#SBATCH -e slurm-%x-%A.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 4G # The memory needed depends on the size of the gtf
#SBATCH --cpus-per-task 1 # This script uses a single CPU
#SBATCH --time 30:00 # This depends on the size of the gtf
#SBATCH --job-name prepare_gtf # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/ # This directory must exists, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
# If the filePathForInputGTF does not exists,
# Ensembl gtf will be downloaded using this info:
ensemblVersion=104
genome=mm39

### Specify the paths
filePathForInputGTF="/scratch/ldelisle/RNAseq/Mus_musculus.GRCm39.104.gtf.gz"
dirPathForResults="/scratch/ldelisle/RNAseq/"
# This script will use filterGtfLikeAnouk.R and mergeGenes_overlap.R
# The last version of the scripts are available at
# https://raw.githubusercontent.com/lldelisle/toolBoxForMutantAndWTGenomes/main/scripts/filterGtfLikeAnouk.R
# and
# https://raw.githubusercontent.com/lldelisle/toolBoxForMutantAndWTGenomes/main/scripts/mergeGenes_overlap.R
# If they do not exist
# The R scripts will be put in dirPathForScripts
dirPathForScripts="/home/ldelisle/scripts/Scitas_template_bashScript/"

### Specify the way to deal with dependencies:
# (R package rtracklayer)
# The module load depends on the installation
# # For baobab:
# module purge
# module load GCC/11.3.0  OpenMPI/4.1.4
# module load R/4.2.1
# or
# module purge
# module load Miniconda3/4.9.2
# or
module purge
module load gcc/8.4.0
module load r/4.0.2

# You can choose to use a conda environment to solve r-rtracklayer dependencies
# Comment it if you will use module load
# condaEnvName=rna202210

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
R --version
if [ $? -ne 0 ]
then
  echo "R is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
# Check plotly and rmarkdown and tidyverse are installed:
Rscript -e "library(rtracklayer)"
if [ $? -ne 0 ]
then
  echo "rtracklayer R package is missing."
  exit 1
fi

mkdir -p ${dirPathForResults}
cd ${dirPathForResults}

if [ ! -e ${filePathForInputGTF} ]; then
    # It will be downloaded:
    if [ $genome = "mm10" ]; then
        versionOfGtf="Mus_musculus.GRCm38.$ensemblVersion"
        species="mus_musculus"
    elif [ $genome = "galGal6" ]; then
        # ftp://ftp.ensembl.org/pub/release-95/gtf/gallus_gallus/Gallus_gallus.GRCg6a.95.gtf.gz
        versionOfGtf="Gallus_gallus.GRCg6a.$ensemblVersion"
        species="gallus_gallus"
    elif [ $genome = "mm39" ]; then
        versionOfGtf="Mus_musculus.GRCm39.$ensemblVersion"
        species="mus_musculus"
    else
        echo "unknown genome"
        exit 1
    fi
    wget "ftp://ftp.ensembl.org/pub/release-${ensemblVersion}/gtf/${species}/${versionOfGtf}.gtf.gz" -O $filePathForInputGTF
    if [ ! -e ${filePathForInputGTF} ]; then
        echo "Download failed"
        exit 1
    fi
fi

wget "https://raw.githubusercontent.com/lldelisle/toolBoxForMutantAndWTGenomes/main/scripts/filterGtfLikeAnouk.R" \
  -O ${dirPathForScripts}/filterGtfLikeAnouk.R -nc
wget "https://raw.githubusercontent.com/lldelisle/toolBoxForMutantAndWTGenomes/main/scripts/mergeGenes_overlap.R" \
  -O ${dirPathForScripts}/mergeGenes_overlap.R -nc

cp ${filePathForInputGTF} .
# Filter for read-through transcripts and non-coding transcripts in coding genes.
Rscript ${dirPathForScripts}/filterGtfLikeAnouk.R $(basename ${filePathForInputGTF}) 
# Merge genes with same name which overlap
Rscript ${dirPathForScripts}/mergeGenes_overlap.R FilteredTranscriptsOf$(basename $(basename ${filePathForInputGTF} .gz) .gtf)_ExonsOnly_UCSC.gtf
