#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1 # This allows to speed the mapping part of the pipeline
#SBATCH --time 04:00:00 # This depends on the size of the bedgraphs
#SBATCH --array=1-3 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name RNAseq_mergeRep # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/RNAseq/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

# This script merge coverage of replicates.
# It uses deeptools bigwigAverage


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

binsize=50 # 50 is default value. Decreasing will increase computation time.

### Specify the paths (same as for script 01)

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$PWD/"
# All merges are registered into a table where
# first column is the sample name
# second column is the name of all replicates separated by comma
filePathForTable="/home/ldelisle/scripts/scitas_sbatchhistory/2022/20221018_testRNA/table_RNA.txt"

### Specify the way to deal with dependencies:

# The only dependency is deeptools version >=3.5.2

# Or if you are using conda
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve deeptools>=3.5.2 dependency
# You can create it with: conda create -n deeptools3.5.2 python=3.9 deeptools=3.5.2
# Comment it if you will use module load
condaEnvName=deeptools3.5.2
######


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
python --version
if [ $? -ne 0 ]
then
  echo "python is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
bigwigAverage --version
if [ $? -ne 0 ]
then
  echo "deeptools is not installed or version is below 3.5.2. Please install it for example in the conda environment."
  exit 1
fi

sample=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
samples=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{split($2,a,",");for(j in a){print a[j]}}')
n=$(echo $samples | wc -w)

# Each sample is processed into an independent directory:
pathResults=${dirPathWithResults}/${sample}/

# The directory is created (if not existing)
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo $sample

# The analysis part takes part within the pathResults
cd $pathResults

# Three variable are initialized:
doStrand=1
noGlobal=0
allBWfwd=""
allBWrev=""
allBW=""

for s in $samples; do
    if [ -e ${dirPathWithResults}/${s}/${s}_Uniq_norm.bw ]; then
        allBW="${allBW} ${dirPathWithResults}/${s}/${s}_Uniq_norm.bw"
    else
        echo "Uniq norm bw not found for $s"
        noGlobal=1
    fi
    if [ -e ${dirPathWithResults}/${s}/${s}_Uniq_fwd_norm.bw ]; then
        allBWfwd="${allBWfwd} ${dirPathWithResults}/${s}/${s}_Uniq_fwd_norm.bw"
    else
        doStrand=0
    fi
    if [ -e ${dirPathWithResults}/${s}/${s}_Uniq_rev_norm.bw ]; then
        allBWrev="${allBWrev} ${dirPathWithResults}/${s}/${s}_Uniq_rev_norm.bw"
    else
        doStrand=0
    fi
done

if [ $noGlobal -eq 1 ]; then
    if [ -z "$allBW" ]; then
        echo "No combined strand coverage"
    else
        echo "Some combined strand coverage are missing will not compute combined strand coverage average."
    fi
else
    bigwigAverage --bigwigs $allBW --skipNAs --binSize $binsize -o ${sample}_Uniq_norm.bw
fi
if [ $doStrand = 1 ]; then
    bigwigAverage --bigwigs $allBWfwd --skipNAs --binSize $binsize -o ${sample}_Uniq_fwd_norm.bw
    bigwigAverage --bigwigs $allBWrev --skipNAs --binSize $binsize -o ${sample}_Uniq_rev_norm.bw
else
    if [ -z "$allBWfwd" ]; then
        echo "No specific strand coverage"
    else
        echo "Some specific strand coverage are missing will not compute strand specific coverage."
    fi
fi

mkdir -p ${dirPathWithResults}/allFinalFiles/bw
cp *bw ${dirPathWithResults}/allFinalFiles/bw/
