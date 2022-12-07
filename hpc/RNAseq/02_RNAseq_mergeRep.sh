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
# It uses unionbedgraph from bedtools + awk


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
# Which genome has been used for mapping
genome=mm39

### Specify the paths (same as for script 01)

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$PWD/"
# All merges are registered into a table where
# first column is the sample name
# second column is the name of all replicates separated by comma
filePathForTable="/home/ldelisle/scripts/scitas_sbatchhistory/2022/20221018_testRNA/table_RNA.txt"
filePathForFasta="/home/ldelisle/genomes/fasta/${genome}.fa"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load GCC/8.3.0 #required for samtools, bowtie2 and bedtools
# module load OpenMPI/3.1.4
# module load SAMtools/1.10
# module load BEDTools/2.28.0
# export PATH=$PATH:/home/users/d/darbellf/live/softwares/ucsc_utilities/ #to load bedGraphToBigWig

# Or if you are using conda
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve bedtools/bedgraphtobigwig dependencies
# You can create it with: conda create -n rna202209 cutadapt samtools star cufflinks bedtools ucsc-bedgraphtobigwig "sra-tools>=2.11"
# Comment it if you will use module load
condaEnvName=rna202209
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
bedtools --version
if [ $? -ne 0 ]
then
  echo "bedtools is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
bedGraphToBigWig 2>&1
if [ $? -eq 127 ]
then
  echo "bedGraphToBigWig is not installed but required. Please install it for example in the conda environment."
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
allBDGfwd=""
allBDGrev=""
allBDG=""

for s in $samples; do
    if [ -e ${dirPathWithResults}/${s}/${s}_Uniq_norm_sorted.bedgraph ]; then
        allBDG="${allBDG} ${dirPathWithResults}/${s}/${s}_Uniq_norm_sorted.bedgraph"
    elif  [ -e ${dirPathWithResults}/${s}/${s}_Uniq_norm.bedgraph.gz ]; then
        zcat ${dirPathWithResults}/${s}/${s}_Uniq_norm.bedGraph | grep -v track | LC_ALL=C sort -k1,1 -k2,2n > ${dirPathWithResults}/${s}/${s}_Uniq_norm_sorted.bedGraph
        allBDG="${allBDG} ${dirPathWithResults}/${s}/${s}_Uniq_norm_sorted.bedgraph"
    else
        echo "Uniq norm bedgraph not found for $s"
        exit 1
    fi
    if [ -e ${dirPathWithResults}/${s}/${s}_Uniq_fwd_norm_sorted.bedgraph ]; then
        allBDGfwd="${allBDGfwd} ${dirPathWithResults}/${s}/${s}_Uniq_fwd_norm_sorted.bedgraph"
    elif  [ -e ${dirPathWithResults}/${s}/${s}_Uniq_fwd_norm.bedgraph.gz ]; then
        zcat ${dirPathWithResults}/${s}/${s}_Uniq_fwd_norm.bedGraph | grep -v track | LC_ALL=C sort -k1,1 -k2,2n > ${dirPathWithResults}/${s}/${s}_Uniq_fwd_norm_sorted.bedGraph
        allBDGfwd="${allBDGfwd} ${dirPathWithResults}/${s}/${s}_Uniq_fwd_norm_sorted.bedgraph"
    else
        doStrand=0
    fi
    if [ -e ${dirPathWithResults}/${s}/${s}_Uniq_rev_norm_sorted.bedgraph ]; then
        allBDGrev="${allBDGrev} ${dirPathWithResults}/${s}/${s}_Uniq_rev_norm_sorted.bedgraph"
    elif  [ -e ${dirPathWithResults}/${s}/${s}_Uniq_rev_norm.bedgraph.gz ]; then
        zcat ${dirPathWithResults}/${s}/${s}_Uniq_rev_norm.bedGraph | grep -v track | LC_ALL=C sort -k1,1 -k2,2n > ${dirPathWithResults}/${s}/${s}_Uniq_rev_norm_sorted.bedGraph
        allBDGrev="${allBDGrev} ${dirPathWithResults}/${s}/${s}_Uniq_rev_norm_sorted.bedgraph"
    else
        doStrand=0
    fi
done
if [ $n -ne 1 ]; then
    # When there are more than 1 bedgraph
    # unionbedg is used to cut the genome in blocks
    # where all bedgraph have constant values
    # Then awk is used to perform the mean
    bedtools unionbedg -i $allBDG | awk -v sn=$sample -v n=$n -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"average of "n" normalized Uniq reads of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > ${sample}_Uniq_norm.bedGraph
    bedGraphToBigWig ${sample}_Uniq_norm.bedGraph
    gzip ${sample}_Uniq_norm.bedGraph
    if [ $doStrand = 1 ]; then
        bedtools unionbedg -i $allBDGfwd | awk -v sn=$sample -v n=$n -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"average of "n" normalized Uniq reads forward of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > ${sample}_Uniq_fwd_norm.bedGraph
        bedGraphToBigWig ${sample}_Uniq_fwd_norm.bedGraph
        gzip ${sample}_Uniq_fwd_norm.bedGraph
        bedtools unionbedg -i $allBDGrev | awk -v sn=$sample -v n=$n -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"average of "n" normalized Uniq reads reverse of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > ${sample}_Uniq_rev_norm.bedGraph
        bedGraphToBigWig ${sample}_Uniq_rev_norm.bedGraph
        gzip ${sample}_Uniq_rev_norm.bedGraph
    fi
else
    cat $allBDG | awk -v sn=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"normalized Uniq reads of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}$4!=0{print}' > ${sample}_Uniq_nor.bedGraph
    bedGraphToBigWig ${sample}_Uniq_norm.bedGraph
    gzip ${sample}_Uniq_norm.bedGraph
    if [ $doStrand = 1 ]; then
        cat $allBDGfwd | awk -v sn=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"normalized Uniq reads forward of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > ${sample}_Uniq_fwd_norm.bedGraph
        bedGraphToBigWig ${sample}_Uniq_fwd_norm.bedGraph
        gzip ${sample}_Uniq_fwd_norm.bedGraph
        cat $allBDGrev | awk -v sn=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"normalized Uniq reads reverse of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > ${sample}_Uniq_rev_norm.bedGraph
        bedGraphToBigWig ${sample}_Uniq_rev_norm.bedGraph
        gzip ${sample}_Uniq_rev_norm.bedGraph
    fi
fi

mkdir -p ${dirPathWithResults}/allFinalFiles/bw
cp *bw ${dirPathWithResults}/allFinalFiles/bw/
