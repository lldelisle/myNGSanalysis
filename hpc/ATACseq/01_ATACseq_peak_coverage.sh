#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 24 # This allows to speed the mapping part of the pipeline
#SBATCH --time 04:00:00 # This depends on the size of the fastqs
#SBATCH --array=1-8 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name ATACseq # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /home/users/d/delislel/scratch/alphaTC1/ATACseq/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
#SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days
#SBATCH --account herrerap

# This script run cutadapt with Nextera PE adapters
# alignment with Bowtie2 --very-sensitive (end-to-end) --no-discordant --dovetail -X 1000
# remove duplicates with picard
# convert bam to bed then macs2 for coverage and peak calling with --keep-dup all --shift -100 --extsize 200
# Normalize the coverage by the million of reads in peaks (summit +-500bp)


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
# number of CPU to use
# Only change if you don't want to use all CPUs allocated
nbOfThreads=${SLURM_CPUS_PER_TASK}
# Which genome to map on
# /!\ This script will use bowtie2 and bowtie2 is not 'alt-aware' so do not use a genome with alt contigs
genome=mm39

### Specify the paths

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$PWD/"
# Where fastqs are stored:
dirPathForFastq="${dirPathWithResults}/fastq/"
# This script will use fromMacs2BdgToSimplifiedBdgAndBw.sh
# The last version of the script is available at
# https://github.com/lldelisle/myNGSanalysis/raw/refs/heads/main/hpc/helper_scripts/fromMacs2BdgToSimplifiedBdgAndBw.sh
# If it does not exists
# The python script will be put in dirPathForScripts
dirPathForScripts="$HOME/scripts/Scitas_template_bashScript/"
# All samples are registered into a table where
# first column is the sample name
# second column is the path of the R1 fastq relatively to dirPathForFastq
# third column is same for R2
# Alternatively second column can be SRA number but third column must be filled by anything for example also the SRA number
filePathForTable="$HOME/scripts/alphatc1-clone-9/ATACseq/table_ATAC.txt"

basenamePathForB2Index="/cvmfs/data.galaxyproject.org/byhand/${genome}/bowtie2_index/${genome}/${genome}"
filePathForFasta="/cvmfs/data.galaxyproject.org/byhand/${genome}/sam_indexes/${genome}/${genome}.fa"

picardCommand="picard"
# For baobab:
# picardCommand="java -jar $EBROOTPICARD/picard.jar"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve cutadapt/bowtie2/samtools/bedtools/picard/macs2/bedgraphtobigwig dependencies
# I created mine with: conda create -n atac202209 cutadapt samtools bedtools bowtie2 picard macs2 ucsc-bedgraphtobigwig
# Comment it if you will use module load
# condaEnvName=atac202209

# You can choose to use singularity, then you need to define each function:
pathToImages=/cvmfs/singularity.galaxyproject.org/all/
# I don't know why on bamboo on clusters it does not work...
# pathToImages=/home/users/d/delislel/scratch/images/

# In order to run the download only once
# This block is only executed for the first sample
images_needed="bowtie2:2.5.4--he96a11b_5 cutadapt:5.0--py312h0fa9677_0 samtools:1.11--h6270b1f_0 bedtools:2.31.1--hf5e1c6e_2 ucsc-bedgraphtobigwig:472--h9b8f530_1 picard:3.3.0--hdfd78af_0 macs2:2.2.9.1--py38h0020b31_1"
if [ ${SLURM_ARRAY_TASK_ID} = ${SLURM_ARRAY_TASK_MIN} ]; then
  for image in $images_needed; do
    if [ ! -e "$pathToImages/$image" ]; then
      wget -nc -O "$pathToImages/$image" "http://datacache.galaxyproject.org/singularity/all/$image"
    fi
  done
else
  echo "Waiting for the first sample to get all images."
  i=0
  while [[ "$i" -lt 20 ]]; do
    sleep 1m
    all_exists=1
    for image in $images_needed; do
      if [ ! -e "$pathToImages/$image" ]; then
        all_exists=0
      fi
    done
    if [ $all_exists = "1" ]; then
      break
    fi
    ((i++))
  done
  if [ $all_exists = "0" ]; then
    echo "After 20 minutes some images are not downloaded"
    echo "Checkout what happened or download them before running"
    exit 1
  fi
fi
function bowtie2() {
    singularity exec "$pathToImages/bowtie2:2.5.4--he96a11b_5" bowtie2 $*
}
function cutadapt() {
  singularity exec "$pathToImages/cutadapt:5.0--py312h0fa9677_0" cutadapt $*
}
function samtools() {
  singularity exec "$pathToImages/samtools:1.11--h6270b1f_0" samtools $*
}
function bedtools() {
  singularity exec "$pathToImages/bedtools:2.31.1--hf5e1c6e_2" bedtools $*
}
function bedGraphToBigWig() {
  singularity exec "$pathToImages/ucsc-bedgraphtobigwig:472--h9b8f530_1" bedGraphToBigWig $*
}
# bedGraphToBigWig is used in another script:
export -f bedGraphToBigWig
export pathToImages=$pathToImages
function picard() {
  singularity exec "$pathToImages/picard:3.3.0--hdfd78af_0" picard $*
}
function macs2() {
  singularity exec "$pathToImages/macs2:2.2.9.1--py38h0020b31_1" macs2 $*
}
# And additional binds:
export APPTAINER_BIND=$HOME/scratch/,$(realpath $HOME/scratch),/cvmfs/data.galaxyproject.org/

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
v=$($picardCommand MarkDuplicates --version 2>&1)
if [ $? -eq 127 ]
then
  echo "picard is not installed but required. Or picardCommand is wrongly set. Please install it for example in the conda environment."
  exit 1
fi
echo "picard MarkDuplicates version $v"
bedtools --version
if [ $? -ne 0 ]
then
  echo "bedtools is not installed but required. Please install it for example in the conda environment."
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
  cutadapt -j $nbOfThreads -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
    -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
    -q 30 -m 15 -o ${sample}-cutadapt_R1.fastq.gz -p ${sample}-cutadapt_R2.fastq.gz \
    $fastqR1 $fastqR2 > ${sample}_report-cutadapt_PE.txt
fi

# Only run mapping if report does not exists
if [ ! -e ${sample}_mapping_stats.txt ];then
  bowtie2 -p $nbOfThreads -x $basenamePathForB2Index --very-sensitive --no-unal --no-mixed --no-discordant --dovetail -X 1000 -1 ${sample}-cutadapt_R1.fastq.gz -2 ${sample}-cutadapt_R2.fastq.gz 2> ${sample}_mapping_stats.txt  | samtools view --threads $nbOfThreads -Su - | samtools sort --threads $nbOfThreads -o ${sample}_mapped_sorted.bam
fi
# Compute chromosome sizes:
if [ ! -e ${filePathForFasta}.fai ]; then
  samtools faidx ${filePathForFasta}
fi

if [ ! -e ${sample}_mapped_sorted_cp_q30_noM.bam ]; then
  allChrsExceptMito=$(samtools view -H ${sample}_mapped_sorted.bam | grep "@SQ" | awk '{gsub("SN:", "", $2); print $2}' | grep -v chrM)
  samtools index ${sample}_mapped_sorted.bam
  samtools view --threads $nbOfThreads -b ${sample}_mapped_sorted.bam -q 30 -f 0x2 $allChrsExceptMito > ${sample}_mapped_sorted_cp_q30_noM.bam
fi

if [ ! -e ${sample}_mapped_sorted_cp_q30_noM_rmdup.bam ]; then
  # Add a readgroup
  samtools addreplacerg -r "@RG\tID:ReadGroup1\tSM:${sample}" -o ${sample}_mapped_sorted_cp_q30_noM_with_RG.bam ${sample}_mapped_sorted_cp_q30_noM.bam
  ${picardCommand} MarkDuplicates SORTING_COLLECTION_SIZE_RATIO=0.15 I=${sample}_mapped_sorted_cp_q30_noM_with_RG.bam O=${sample}_mapped_sorted_cp_q30_noM_rmdup.bam M=${sample}_mapped_sorted_cp_q30_noM_rmdup.log REMOVE_DUPLICATES=true AS=true
fi

## To get tags for homer analysis
##export PATH=$PATH:/home/mayran/software/Homer/bin/
##mkdir -p ${dirPathWithResults}tag
##if [ ! -e ${dirPathWithResults}tag/${sample}.mapped_sorted_cp_q30_noM_rmdup.bam ]; then
##  makeTagDirectory ${padirPathWithResultsth}tag/${sample}.mapped_sorted_cp_q30_noM_rmdup.bam ${sample}_mapped_sorted_cp_q30_noM_rmdup.bam &
##fi
## end

# Get the bash script:
mkdir -p $dirPathForScripts
wget "https://github.com/lldelisle/myNGSanalysis/raw/refs/heads/main/hpc/helper_scripts/fromMacs2BdgToSimplifiedBdgAndBw.sh" \
  -O ${dirPathForScripts}/fromMacs2BdgToSimplifiedBdgAndBw.sh -nc

if [ ! -e ${sample}_macs_likeATAC.bedGraph.gz ]; then
  bedtools bamtobed -i ${sample}_mapped_sorted_cp_q30_noM_rmdup.bam > ${sample}_mapped_sorted_cp_q30_noM_rmdup.bed
  macs2 callpeak -t ${sample}_mapped_sorted_cp_q30_noM_rmdup.bed --nomodel --keep-dup all --shift -100 --extsize 200 -n ${sample} --call-summits -B 2> ${sample}_macs_likeATAC.log
  # The output of macs2 can go over the size of chromosomes so a special care needs to be taken before coverting to bigwig
  bash ${dirPathForScripts}/fromMacs2BdgToSimplifiedBdgAndBw.sh ${sample}_treat_pileup.bdg ${sample}_macs_likeATAC "macs2 like ATAC of ${sample}" ${filePathForFasta}.fai &
fi

# I wait for gzip to finish
wait

if [ ! -e ${sample}_macs_likeATAC_norm1.bedGraph.gz ]; then
  likeATAC=$(cat ${sample}_macs_likeATAC.log | awk '$0~/total tags/{print $NF}')
  zcat ${sample}_macs_likeATAC.bedGraph.gz | awk -v s=$sample -v n=$likeATAC -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\""s" like ATAC normalized by million reads\" visibility=full autoScale=on windowingFunction=mean"}NR>1{$4=$4/n*1e6; print}' > ${sample}_macs_likeATAC_norm1.bedGraph
  bedGraphToBigWig ${sample}_macs_likeATAC_norm1.bedGraph ${filePathForFasta}.fai ${sample}_macs_likeATAC_norm1.bw
  gzip ${sample}_macs_likeATAC_norm1.bedGraph &
fi

if [ ! -e ${sample}_macs_likeATAC_norm2.bedGraph.gz ]; then
  if [ ! -e ${sample}_readsInPeaks.txt ]; then
    bedtools slop -i ${sample}_summits.bed -g ${filePathForFasta}.fai -b 500 > ${sample}_macs_likeATAC_summits_1kb.bed
    bedtools merge -i ${sample}_macs_likeATAC_summits_1kb.bed > ${sample}_macs_likeATAC_summits_1kb_merged.bed
    bedtools coverage -a ${sample}_macs_likeATAC_summits_1kb_merged.bed -b ${sample}_mapped_sorted_cp_q30_noM_rmdup.bed > ${sample}_macs_likeATAC_summits_1kb_merged_coverage.txt
    cat ${sample}_macs_likeATAC_summits_1kb_merged_coverage.txt | awk '{S+=$4}END{print S}' > ${sample}_readsInPeaks.txt
  fi
  readsInPeaks=$(cat ${sample}_readsInPeaks.txt)
  zcat ${sample}_macs_likeATAC.bedGraph.gz | awk -v s=$sample -v n=$readsInPeaks -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\""s" like ATAC normalized by million reads in peaks\" visibility=full autoScale=on windowingFunction=mean"}NR>1{$4=$4/n*1e6; print}' > ${sample}_macs_likeATAC_norm2.bedGraph
  bedGraphToBigWig ${sample}_macs_likeATAC_norm2.bedGraph ${filePathForFasta}.fai ${sample}_macs_likeATAC_norm2.bw
  gzip ${sample}_macs_likeATAC_norm2.bedGraph &
fi

# I wait for all gzip to finish
wait

echo "Everything is done"
find . -size 0 -delete

mkdir -p ${dirPathWithResults}/allFinalFiles/reports/
cp ${sample}_report-cutadapt_PE.txt ${dirPathWithResults}/allFinalFiles/reports/
cp ${sample}_mapping_stats.txt ${dirPathWithResults}/allFinalFiles/reports/
cp ${sample}_mapped_sorted_cp_q30_noM_rmdup.log ${dirPathWithResults}/allFinalFiles/reports/
cp ${sample}_macs_likeATAC.log ${dirPathWithResults}/allFinalFiles/reports/

mkdir -p ${dirPathWithResults}/allFinalFiles/bw_peaks
cp ${sample}_peaks.narrowPeak ${dirPathWithResults}/allFinalFiles/bw_peaks/
cp ${sample}*bw ${dirPathWithResults}/allFinalFiles/bw_peaks/

mkdir -p ${dirPathWithResults}/allFinalFiles/bam/
cp ${sample}_mapped_sorted_cp_q30_noM_rmdup.bam ${dirPathWithResults}/allFinalFiles/bam/
