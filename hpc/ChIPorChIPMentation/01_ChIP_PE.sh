#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 24 # This allows to speed the mapping part of the pipeline
#SBATCH --time 04:00:00 # This depends on the size of the fastqs
#SBATCH --array=1-2 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name ChIP_PE_test # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/ChIP/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

# This script remove adapters 
# make alignment with Bowtie2
# select MAPQ30 alignments
# Compute coverage joining two pairs
# Normalize by nb of tags


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
# number of CPU to use
# Only change if you don't want to use all CPUs allocated
nbOfThreads=${SLURM_CPUS_PER_TASK}
# Which genome to map on
# /!\ This script will use bowtie2 and bowtie2 is not 'alt-aware' so do not use a genome with alt contigs
genome=mm10
adapterSeq="TruSeq"  # comment the non-relevant adapter type
# adapterSeq="NextSeq"

### Specify the paths

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$PWD/"
# Where fastqs are stored:
dirPathForFastq="${dirPathWithResults}/fastq/"
# This script will use fromMacs2BdgToSimplifiedBdgAndBw.sh
# The last version of the script is available at
# https://raw.githubusercontent.com/lldelisle/scriptsForAmandioEtAl2021/main/scripts/fromMacs2BdgToSimplifiedBdgAndBw.sh
# If it does not exists
# The python script will be put in dirPathForScripts
dirPathForScripts="/home/ldelisle/scripts/Scitas_template_bashScript/"
# All samples are registered into a table where
# first column is the sample name
# second column is the path of the R1 fastq relatively to dirPathForFastq
# third column is same for R2
# Alternatively second column can be SRA number but third column must be filled by anything for example also the SRA number
filePathForTable="/home/ldelisle/scripts/scitas_sbatchhistory/2022/TODO/table_ChIP_PE.txt"

basenamePathForB2Index="/home/ldelisle/genomes/bowtie2/${genome}"
filePathForFasta="/home/ldelisle/genomes/fasta/${genome}.fa"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve cutadapt/bowtie2/samtools/macs2/bedgraphtobigwig dependencies
# You can create it with: conda create -n chip202209 cutadapt samtools bowtie2 macs2 ucsc-bedgraphtobigwig
# If you want to use sra you also need sra-tools>=2.11
# Comment it if you will use module load
condaEnvName=atac202209
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
v=$(cutadapt --version)
if [ $? -ne 0 ]
then
  echo "cutadapt is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
echo "cutadapt $v"
bowtie2 --version
if [ $? -ne 0 ]
then
  echo "Bowtie2 is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
if [ ! -e ${basenamePathForB2Index}.rev.2.bt2 ]; then
  echo "${basenamePathForB2Index}.rev.2.bt2 does not exists. Either your variable basenamePathForB2Index is wrong or your bowtie2 index did not full run."
  exit 1
fi
samtools --version
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
macs2 --version
if [ $? -ne 0 ]
then
  echo "macs2 is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
bedGraphToBigWig 2>&1
if [ $? -eq 127 ]
then
  echo "bedGraphToBigWig is not installed but required. Please install it for example in the conda environment."
  exit 1
fi

sample=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
relFilePathFastqR1=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
relFilePathFastqR2=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')

# Each sample is processed into an independent directory:
pathResults=${dirPathWithResults}/${sample}/

# The directory is created (if not existing)
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo $sample

# The analysis part takes part within the pathResults
cd $pathResults

# Only run cutadapt if the report does not exists
if [ ! -e ${sample}_report-cutadapt_PE.txt ]; then
  fastqR1=${dirPathForFastq}/$relFilePathFastqR1
  fastqR2=${dirPathForFastq}/$relFilePathFastqR2
  if [ ! -e $fastqR1 ]; then
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
    fastqR1=${dirPathForFastq}/${sample}_1.fastq.gz
    fastqR2=${dirPathForFastq}/${sample}_2.fastq.gz
  fi
  if [ ! -s $fastqR1 ]; then
    echo "FASTQ R1 IS EMPTY"
    exit 1
  fi
  if [ $adapterSeq = "TruSeq" ]; then
    cutadapt -j $nbOfThreads -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -q 30 -m 15 -o ${sample}-cutadapt_R1.fastq.gz -p ${sample}-cutadapt_R2.fastq.gz \
        $fastqR1 $fastqR2 > ${sample}_report-cutadapt_PE.txt
  else
    if [ $adapterSeq = "NextSeq" ]; then
      cutadapt -j $nbOfThreads -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        -q 30 -m 15 -o ${sample}-cutadapt_R1.fastq.gz -p ${sample}-cutadapt_R2.fastq.gz \
        $fastqR1 $fastqR2 > ${sample}_report-cutadapt_PE.txt
    else
      echo "YOU NEED TO WRITE THE CODE"
      exit 1
    fi
  fi
fi

# Only run mapping if report does not exists
if [ ! -e ${sample}_mapping_stats.txt ];then
  bowtie2 -p $nbOfThreads -x $basenamePathForB2Index -1 ${sample}-cutadapt_R1.fastq.gz -2 ${sample}-cutadapt_R2.fastq.gz 2> ${sample}_mapping_stats.txt  | samtools view --threads $nbOfThreads -Su - | samtools sort --threads $nbOfThreads -o ${sample}_mapped_sorted.bam
fi

# Filtering
if [ ! -e ${sample}_mapped_sorted_q30.bam ]; then
  samtools view --threads $nbOfThreads -b ${sample}_mapped_sorted.bam -q 30 > ${sample}_mapped_sorted_q30.bam
fi

# Compute chromosome sizes:
if [ ! -e ${filePathForFasta}.fai ]; then
  samtools faidx ${filePathForFasta}
fi

# Get the bash script:
mkdir -p $dirPathForScripts
wget "https://raw.githubusercontent.com/lldelisle/scriptsForAmandioEtAl2021/main/scripts/fromMacs2BdgToSimplifiedBdgAndBw.sh" \
  -O ${dirPathForScripts}/fromMacs2BdgToSimplifiedBdgAndBw.sh -nc

# Use macs2 with BAMPE
if [ ! -e ${sample}_macs.bedGraph.gz ]; then
  macs2 callpeak -t ${sample}_mapped_sorted_q30.bam -n ${sample}_macs --call-summits --format BAMPE -B 2> ${sample}_macs.log
  # The output of macs2 can go over the size of chromosomes so a special care needs to be taken before coverting to bigwig
  bash ${dirPathForScripts}/fromMacs2BdgToSimplifiedBdgAndBw.sh ${sample}_macs_treat_pileup.bdg ${sample}_macs "macs2 (BAMPE) of ${sample}" ${filePathForFasta}.fai &
fi

# I wait for gzip to finish
wait

if [ ! -e ${sample}_macs_norm.bedGraph.gz ]; then
  nbfrags=$(grep "fragments after filtering in treatment" ${sample}_macs.log | awk '{print $NF}')
  zcat ${sample}_macs.bedGraph.gz | awk -v s=$sample -v n=$nbfrags -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\""s" BEMPE normalized by million tags\" visibility=full autoScale=on windowingFunction=mean"}NR>1{$4=$4/n*1e6; print}' > ${sample}_macs_norm.bedGraph
  bedGraphToBigWig ${sample}_macs_norm.bedGraph ${filePathForFasta}.fai ${sample}_macs_norm.bw
  gzip ${sample}_macs_norm.bedGraph &
fi

# I wait for gzip to finish
wait

echo "Everything is done"
find . -size 0 -delete

mkdir -p ${dirPathWithResults}/allFinalFiles/reports/
cp ${sample}_report-cutadapt.txt ${dirPathWithResults}/allFinalFiles/reports/
cp ${sample}_mapping_stats.txt ${dirPathWithResults}/allFinalFiles/reports/
cp ${sample}_macs.log ${dirPathWithResults}/allFinalFiles/reports/

mkdir -p ${dirPathWithResults}/allFinalFiles/bw_peaks
cp ${sample}_macs_peaks.narrowPeak ${dirPathWithResults}/allFinalFiles/bw_peaks/
cp ${sample}*bw ${dirPathWithResults}/allFinalFiles/bw_peaks/
