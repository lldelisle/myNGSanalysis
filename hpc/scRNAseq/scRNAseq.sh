#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 16 # This allows to speed the mapping part of the pipeline
#SBATCH --time 12:00:00 # This depends on the size of the fastqs
#SBATCH --array=1-1 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name scRNAseq # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /home/users/d/delislel/scratch/scRNAseq # This directory must exist, this is where will be the error and out files
## Specific to UNIGE:
#SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days
#SBATCH --account andreygu

# This script is designed for GEM v3 or v4 (GEM-X) or Multiome
# This script runs STARsolo without filtering
# Filter with dropletutils
# Compute Loom with velocytopy
# And normalized coverage with STAR

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
genome=mm39

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
gitDir=$HOME/scripts/myproject/
filePathForTable="$gitDir/fastq_to_matrices/table_scRNA.txt"
filePathForGTF="${dirPathWithResults}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.113_ExonsCDSOnly_UCSC.gtf"
# filePathForGTF="${dirPathWithResults}/gencode.v47.primary_assembly.annotation.gtf"
dirPathForSTARIndex="/cvmfs/data.galaxyproject.org/byhand/rnastar/2.7.4a/mm39/mm39/dataset_d85576e1-d6db-411e-898c-79af8214faa4_files/"
filePathForFasta="/cvmfs/data.galaxyproject.org/byhand/${genome}/sam_indexes/${genome}/${genome}.fa"

dirPathForSTARIndex="$PWD/../genomes/${genome}/STAR_2.7.11b/${genome}"
filePathForFasta="$PWD/../genomes/${genome}/${genome}.fa"

# Comment what is not your case
# mode=GEMv4
# mode=GEMv3
mode=Multiome

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve star/samtools/bedtools/bioconductor-dropletutils dependencies
# You can create it with: conda create -n scrna bedtools samtools star "sra-tools>=2.11" bioconductor-dropletutils
# If you want to use sra you also need sra-tools>=2.11
# Comment it if you will use module load
# condaEnvName=scrna

# You can choose to use singularity, then you need to define each function:
pathToImages=/cvmfs/singularity.galaxyproject.org/all/

# In order to run the download only once
# This block is only executed for the first sample
images_needed="samtools:1.11--h6270b1f_0 star:2.7.11b--h5ca1c30_5 bedtools:2.31.1--hf5e1c6e_2 ucsc-bedgraphtobigwig:472--h9b8f530_1 bioconductor-dropletutils:1.8.0--r40h5f743cb_0 velocyto.py:0.17.17--py312h1f1cfbb_7"
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
function STAR() {
  singularity exec "$pathToImages/star:2.7.11b--h5ca1c30_5" STAR $*
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
function Rscript() {
  singularity exec "$pathToImages/bioconductor-dropletutils:1.8.0--r40h5f743cb_0" Rscript $*
}
function R() {
  singularity exec "$pathToImages/bioconductor-dropletutils:1.8.0--r40h5f743cb_0" R $*
}
function velocyto() {
  singularity exec "$pathToImages/velocyto.py:0.17.17--py312h1f1cfbb_7" velocyto $*
}
# And additional binds:
export APPTAINER_BIND=$HOME/scratch/,$(realpath $HOME/scratch),/cvmfs/data.galaxyproject.org/,/srv/beegfs


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
R --version
if [ $? -ne 0 ]
then
  echo "R is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
# Check DropletUtils installed:
Rscript -e "library(DropletUtils)"
if [ $? -ne 0 ]
then
  echo "DropletUtils R package is missing."
  exit 1
fi
velocyto --version
if [ $? -ne 0 ]
then
  echo "velocyto is not installed but required. Please install it for example in the conda environment."
  exit 1
fi

sample=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
relFilePathFastqR1=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
relFilePathFastqR2=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')

dropletutilsRscript="${dirPathWithResults}/dropletutils.R"

if [ "$mode" = "GEMv4" ]; then
  barcodeWhiteList="${dirPathWithResults}/3M-3pgex-may-2023.txt"
elif [ "$mode" = "GEMv3" ]; then
  barcodeWhiteList="${dirPathWithResults}/3M-february-2018.txt"
elif [ "$mode" = "Multiome" ]; then
  barcodeWhiteList="${dirPathWithResults}/737K-arc-v1_GEX.txt"
else
  echo "Unknown mode"
  exit 1
fi

if [ ! -e $barcodeWhiteList ]; then
  # In order to write the R script and 
  # download the barcode WL only once
  # This block is only executed for the first sample
  if [ ${SLURM_ARRAY_TASK_ID} = ${SLURM_ARRAY_TASK_MIN} ]; then
    # This script is inspired from https://github.com/galaxyproject/tools-iuc/blob/b1797a2dee3977cdf40d3cf413ab9ec1e0cb3f26/tools/dropletutils/scripts/dropletutils.Rscript
    echo "filein <- commandArgs(TRUE)[1]
fileout <- commandArgs(TRUE)[2]
fdr_threshold <- 0.01
lower <- 100
suppressWarnings(suppressPackageStartupMessages(require(DropletUtils)))
set.seed(100)
sce <- read10xCounts(filein, col.names = T, type = \"sparse\")
e_out <- emptyDrops(m = counts(sce), lower = lower)
e_out\$is_cell <- e_out\$FDR <= fdr_threshold
e_out\$is_cellandlimited <- e_out\$is_cell & e_out\$Limited
e_out\$is_cellandlimited[is.na(e_out\$is_cellandlimited)] <- FALSE
kept_barcodes <- rownames(e_out)[e_out\$is_cellandlimited]
sce_filtered <- sce[, kept_barcodes]
write10xCounts(fileout, counts(sce_filtered),
               type = \"sparse\",
               gene.symbol = rowData(sce_filtered)\$Symbol,
               overwrite = TRUE)
message(paste(\"Cells and Limited:\", sum(e_out\$is_cellandlimited)))" > ${dropletutilsRscript}
    if [ "$mode" = "GEMv4" ]; then
      wget "https://github.com/Teichlab/scg_lib_structs/raw/dc91d3940718d5e73566de3aa73ebaef33d5cfd8/data/10X-Genomics/3M-3pgex-may-2023.txt.gz" -O tmp.txt.gz
    elif [ "$mode" = "GEMv3" ]; then
      wget "https://zenodo.org/record/3457880/files/3M-february-2018.txt.gz" -O tmp.txt.gz
    elif [ "$mode" = "Multiome" ]; then
      wget "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/gex_737K-arc-v1.txt.gz" -O tmp.txt.gz
    fi
    gunzip tmp.txt.gz
    mv tmp.txt $barcodeWhiteList
  else
    echo "Waiting for the first sample to download the barcode WL file."
    i=0
    while [[ "$i" -lt 20 ]]; do
      sleep 1m
      if [ -e $barcodeWhiteList ]; then
        break
      fi
      ((i++))
    done
    if [ ! -e $barcodeWhiteList ]; then
      echo "After 20 minutes the barcode whitelist file was not created"
      echo "Checkout what happened or generate it before running"
      exit 1
    fi
  fi
fi

# In order to run the fai only once
# This block is only executed for the first sample
if [ ${SLURM_ARRAY_TASK_ID} = ${SLURM_ARRAY_TASK_MIN} ]; then
  # Compute chromosome sizes:
  if [ ! -e ${filePathForFasta}.fai ]; then
    samtools faidx ${filePathForFasta}
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

# Check if you need to download the fastqs:
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
# Map with STARsolo without filtering
if [ ! -e Aligned.sortedByCoord.out.bam ];then
  STAR --runThreadN $nbOfThreads --genomeDir ${dirPathForSTARIndex} \
    --readFilesIn ${fastqR2} ${fastqR1} \
    $options --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbOverhang '99' --sjdbGTFfile $filePathForGTF \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist $barcodeWhiteList \
    --soloUMIlen 12 \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter None \
    --quantMode GeneCounts \
    --outSAMattributes NH HI AS nM GX GN CB UB \
    --outSAMmapqUnique 60 \
    --outWigType bedGraph \
    --outWigStrand Stranded
  samtools index Aligned.sortedByCoord.out.bam
fi
if [ ! -e filtered/matrix.mtx ]; then
  Rscript ${dropletutilsRscript} Solo.out/Gene/raw/ filtered
fi

# Coverage
# Sort files
for f in *.out.bg; do
  output=${f/.bg/.sorted.bg}
  if [ ! -e $output ]; then
    bedtools sort -i $f > $output
  fi
done
# Convert to bigwig
if [ ! -e ${sample}_Uniq_fwd_norm.bw ]; then
  bedGraphToBigWig Signal.Unique.str1.out.sorted.bg ${filePathForFasta}.fai ${sample}_Uniq_fwd_norm.bw
fi
if [ ! -e ${sample}_Uniq_rev_norm.bw ]; then
  bedGraphToBigWig Signal.Unique.str2.out.sorted.bg ${filePathForFasta}.fai ${sample}_Uniq_rev_norm.bw
fi

# Run velocyto
if [ ! -e ${sample}.loom ]; then
  # Velocyto expect a specific architecture
  mkdir -p ${sample}/outs/filtered_gene_bc_matrices/whatever/
  ln -s $PWD/Aligned.sortedByCoord.out.bam ${sample}/outs/possorted_genome_bam.bam
  ln -s $PWD/filtered/barcodes.tsv ${sample}/outs/filtered_gene_bc_matrices/whatever/barcodes.tsv
  velocyto  run10x  -t 'uint16'  --samtools-threads ${nbOfThreads} '-vv' ${sample} ${filePathForGTF} 
  mv $sample/velocyto/*.loom ${sample}.loom
fi

mkdir -p ${dirPathWithResults}/allFinalFiles/reports/
cp Log.final.out ${dirPathWithResults}/allFinalFiles/reports/${sample}_STAR_logFinal.txt

mkdir -p ${dirPathWithResults}/allFinalFiles/bam/
cp Aligned.sortedByCoord.out.bam ${dirPathWithResults}/allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
cp Aligned.sortedByCoord.out.bam.bai ${dirPathWithResults}/allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam.bai

mkdir -p ${dirPathWithResults}/allFinalFiles/matrices/
if [ -e ${dirPathWithResults}/allFinalFiles/matrices/${sample} ]; then
  rm -r ${dirPathWithResults}/allFinalFiles/matrices/${sample}
fi
cp -r filtered ${dirPathWithResults}/allFinalFiles/matrices/${sample}

mkdir -p ${dirPathWithResults}/allFinalFiles/loom/
cp ${sample}.loom ${dirPathWithResults}/allFinalFiles/loom/

mkdir -p ${dirPathWithResults}/allFinalFiles/bw
cp *bw ${dirPathWithResults}/allFinalFiles/bw/
