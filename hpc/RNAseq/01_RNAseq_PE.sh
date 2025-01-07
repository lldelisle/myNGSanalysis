#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 16 # This allows to speed the mapping part of the pipeline
#SBATCH --time 12:00:00 # This depends on the size of the fastqs
#SBATCH --array=2-19 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name RNAseq # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /home/users/d/delislel/scratch/NET_project_RNAseq/ # This directory must exist, this is where will be the error and out files
## Specific to UNIGE:
#SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

# This script run cutadapt to remove adapters and bad quality bases
# make alignment, counting and coverage with STAR ENCODE parameters
# Evaluate FPKM with cufflinks


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
# number of CPU to use
# Only change if you don't want to use all CPUs allocated
# I tryed 36 and it failed:
# BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file ./_STARtmp//BAMsort/19/48
# SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.
if [ ${SLURM_CPUS_PER_TASK} -gt 16 ]; then
  nbOfThreads=16
else
  nbOfThreads=${SLURM_CPUS_PER_TASK}
fi
# Which genome to map on
# /!\ This script will use STAR and STAR is not 'alt-aware' so do not use a genome with alt contigs
# genome=mm39
genome=hg38.analysisSet
adapterSeq="TruSeq"  # comment the non-relevant adapter type
# adapterSeq="NextSeq"

### Specify the paths

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$PWD/"
# Where fastqs are stored:
dirPathForFastq="${dirPathWithResults}/fastq/"
# All samples are registered into a table where
# first column is the sample name
# second column is the path of the R1 fastq relatively to dirPathForFastq
# third column is same for R2
# Alternatively second column can be SRA number but third column must be filled by anything for example also the SRA number
# fourth column is the strandness of the library: forward or reverse or unstranded
filePathForTable="/home/users/d/delislel/scripts/net_project/RNAseq/table.txt"
# filePathForGTF="${dirPathWithResults}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.104_ExonsOnly_UCSC.gtf"
filePathForGTF="${dirPathWithResults}/gencode.v47.primary_assembly.annotation.gtf"
dirPathForSTARIndex="/home/users/d/delislel/genomes/STARIndex_2.7.11a/${genome}/"
filePathForFasta="/home/users/d/delislel/genomes/fasta/${genome}.fa"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve cutadapt/star/samtools/cufflinks/bedtools/bedgraphtobigwig dependencies
# You can create it with: conda create -n rna202209 cutadapt samtools star cufflinks bedtools ucsc-bedgraphtobigwig "sra-tools>=2.11"
# If you want to use sra you also need sra-tools>=2.11
# Comment it if you will use module load
# condaEnvName=rna202209
######

# You can choose to use singularity, then you need to define each function:
function cutadapt() {
  singularity exec "/cvmfs/singularity.galaxyproject.org/all/cutadapt:5.0--py312h0fa9677_0" cutadapt $*
}
function samtools() {
  singularity exec "/cvmfs/singularity.galaxyproject.org/all/samtools:1.9--h91753b0_8" samtools $*
}
function STAR() {
  singularity exec "/cvmfs/singularity.galaxyproject.org/all/star:2.7.11b--h5ca1c30_4" STAR $*
}
function cufflinks() {
  singularity exec "/cvmfs/singularity.galaxyproject.org/all/cufflinks:2.2.1--py36_2" cufflinks $*
}
# Cufflinks is used in another script:
export -f cufflinks
function bedtools() {
  singularity exec "/cvmfs/singularity.galaxyproject.org/all/bedtools:2.31.1--hf5e1c6e_2" bedtools $*
}
# bedtools is used in another script:
export -f bedtools
function bedGraphToBigWig() {
  singularity exec "/cvmfs/singularity.galaxyproject.org/all/ucsc-bedgraphtobigwig:472--h9b8f530_1" bedGraphToBigWig $*
}
# bedtools is used in another script:
export -f bedGraphToBigWig

# And additional binds:
export APPTAINER_BIND=$dirPathWithResults

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

v=$(STAR --version)
if [ $? -ne 0 ]
then
  echo "STAR is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
echo "STAR version $v"

if [ ! -e ${dirPathForSTARIndex}/chrLength.txt ]; then
  echo "The dirPathForSTARIndex is not valid."
  exit 1
fi

samtools --version
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
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
relFilePathFastqR1=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
relFilePathFastqR2=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')
stranded=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $4}')

# In order to run the fai and chrM.gtf only once
# This block is only executed for the first sample
if [ ${SLURM_ARRAY_TASK_ID} = ${SLURM_ARRAY_TASK_MIN} ]; then
  # Compute chromosome sizes:
  if [ ! -e ${filePathForFasta}.fai ]; then
    samtools faidx ${filePathForFasta}
  fi
  # Generate the chrM.gtf:
  if [ ! -e ${dirPathWithResults}/${genome}_chrM.gtf ]; then
    mkdir -p ${dirPathWithResults}
    chrM_size=$(cat ${filePathForFasta}.fai  | grep -w chrM | cut -f 2)
    echo -e "chrM\tchrM_gene\texon\t1\t${chrM_size}\t.\t+\t.\tgene_id \"chrM_gene_plus\"; transcript_id \"chrM_tx_plus\"; exon_id \"chrM_ex_plus\";" > ${dirPathWithResults}/${genome}_chrM.gtf
    echo -e "chrM\tchrM_gene\texon\t1\t${chrM_size}\t.\t-\t.\tgene_id \"chrM_gene_minus\"; transcript_id \"chrM_tx_minus\"; exon_id \"chrM_ex_minus\";" >> ${dirPathWithResults}/${genome}_chrM.gtf
  fi
fi

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

# Map with STAR
if [ ! -e Aligned.sortedByCoord.out.bam ];then
  if [ "$stranded" = "unstranded" ]; then 
    #I need to add --outSAMstrandField intronMotif because it is not stranded library
    options="--outSAMstrandField intronMotif --outWigStrand Unstranded"
  else
    options="--outWigStrand Stranded"
  fi
  STAR --runThreadN $nbOfThreads --genomeDir ${dirPathForSTARIndex} \
    --readFilesIn ${sample}-cutadapt_R1.fastq.gz ${sample}-cutadapt_R2.fastq.gz \
    $options --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbOverhang '99' --sjdbGTFfile $filePathForGTF \
    --quantMode GeneCounts \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outWigType bedGraph 
  samtools index Aligned.sortedByCoord.out.bam
fi

# Cufflinks is run in parallel
if [ -e ${filePathForFasta} ] && [ -e ${dirPathWithResults}/${genome}_chrM.gtf ] && [ -e $filePathForGTF ] && [ -s Aligned.sortedByCoord.out.bam ]; then
  if [ ! -e cufflinks/isoforms.fpkm_tracking ];then
    echo "mkdir -p cufflinks" >>cufflinks_${sample}.sh
    echo "cufflinks -p ${nbOfThreads} -o cufflinks --max-bundle-length 10000000 --multi-read-correct --library-type \"fr-firststrand\" -b ${filePathForFasta} --no-effective-length-correction -M ${dirPathWithResults}/${genome}_chrM.gtf -G ${filePathForGTF} Aligned.sortedByCoord.out.bam" >> cufflinks_${sample}.sh
    echo "" >> cufflinks_${sample}.sh
    if [ "$stranded" = "unstranded" ]; then 
      sed -i 's/fr-firststrand/fr-unstranded/g' cufflinks_${sample}.sh
    elif [ "$stranded" = "forward" ]; then
      sed -i 's/fr-secondstrand/fr-unstranded/g' cufflinks_${sample}.sh
    fi
    echo "Launching cufflinks"
    bash cufflinks_${sample}.sh &
  fi
else
  echo "cufflinks not launch because some files are missing."
fi

if [ ! -e htseqCount_${sample}.txt ] && [ -e ReadsPerGene.out.tab ];then 
  echo "write htseqCount" # compile htseq counts-like from STAR counts
  if [ "$stranded" = "unstranded" ]; then 
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$2}' > htseqCount_${sample}.txt
  elif [ "$stranded" = "forward" ]; then
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$3}' > htseqCount_${sample}.txt
  else
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$4}' > htseqCount_${sample}.txt
  fi
fi

# Coverage
# Sort files
for f in *.out.bg; do
  output=${f/.bg/.sorted.bg}
  if [ ! -e $output ]; then
    bedtools sort -i $f > $output
  fi
done
# Convert files
if [ "$stranded" = "unstranded" ]; then 
  if [ ! -e ${sample}_Uniq_norm.bw ]; then
    bedGraphToBigWig Signal.Unique.str1.out.sorted.bg ${filePathForFasta}.fai ${sample}_Uniq_norm.bw
  fi
elif [ "$stranded" = "forward" ]; then
  if [ ! -e ${sample}_Uniq_fwd_norm.bw ]; then
    bedGraphToBigWig Signal.Unique.str1.out.sorted.bg ${filePathForFasta}.fai ${sample}_Uniq_fwd_norm.bw
  fi
  if [ ! -e ${sample}_Uniq_rev_norm.bw ]; then
  bedGraphToBigWig Signal.Unique.str2.out.sorted.bg ${filePathForFasta}.fai ${sample}_Uniq_rev_norm.bw
  fi
else
  if [ ! -e ${sample}_Uniq_fwd_norm.bw ]; then
    bedGraphToBigWig Signal.Unique.str2.out.sorted.bg ${filePathForFasta}.fai ${sample}_Uniq_fwd_norm.bw
  fi
  if [ ! -e ${sample}_Uniq_rev_norm.bw ]; then
  bedGraphToBigWig Signal.Unique.str1.out.sorted.bg ${filePathForFasta}.fai ${sample}_Uniq_rev_norm.bw
  fi
fi

wait
echo "Everything is done"
find . -size 0 -delete


mkdir -p ${dirPathWithResults}/allFinalFiles/reports/
cp ${sample}_report-cutadapt_PE.txt ${dirPathWithResults}/allFinalFiles/reports/
cp Log.final.out ${dirPathWithResults}/allFinalFiles/reports/${sample}_STAR_logFinal.txt

mkdir -p ${dirPathWithResults}/allFinalFiles/bam/
cp Aligned.sortedByCoord.out.bam ${dirPathWithResults}/allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
cp Aligned.sortedByCoord.out.bam.bai ${dirPathWithResults}/allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam.bai

mkdir -p ${dirPathWithResults}/allFinalFiles/bw
cp *bw ${dirPathWithResults}/allFinalFiles/bw/

mkdir -p ${dirPathWithResults}/allFinalFiles/counts_FPKM
cp cufflinks/genes.fpkm_tracking ${dirPathWithResults}/allFinalFiles/counts_FPKM/FPKM_${sample}.txt
cp cufflinks/isoforms.fpkm_tracking ${dirPathWithResults}/allFinalFiles/counts_FPKM/FPKM_${sample}_isoforms.txt
cp htseqCount_${sample}.txt ${dirPathWithResults}/allFinalFiles/counts_FPKM/
