#!/bin/bash

#SBATCH -o slurm-%x-%A.out # Template for the std output of the job uses the job name and the job id
#SBATCH -e slurm-%x-%A.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 10G # The memory needed depends on the number of samples
#SBATCH --cpus-per-task 1 # This does not use multiple CPUs
#SBATCH --time 02:00:00 # This depends on the number of samples
#SBATCH --job-name RNAseq_R_pc # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /home/users/d/delislel/scripts/whatever/RNAseq/ # This directory must exist, this is where will be the error and out files
## Specific to UNIGE:
#SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days
#SBATCH --account=herrerap

# This script run Rscripts to merge counts and FPKM
# Generate first plots
# And DESeq2 analysis


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the paths
# Put in dirPathWithConfigFiles the directory
# where all config files will be written
# this can be on your home
dirPathWithConfigFiles="$PWD/"
# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample (should match the one in 01_RNAseq_XX.sh)
dirPathWithResults="$HOME/scratch/whatever/RNAseq/"
# filePathForGTF="${dirPathWithResults}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC.gtf"
filePathForGTF="${dirPathWithResults}/gencode.v47.primary_assembly.annotation.gtf"
# This is where should be clone the 2 github repositories:
# https://github.com/lldelisle/toolBoxForMutantAndWTGenomes
# https://github.com/lldelisle/rnaseq_rscripts
dirPathWithDependencies="$HOME/scripts/"
# A samples plan is a tabular table with at least one column named 'sample'
# This should match the table of 01_RNAseq_XX.sh
# The sample names should not start with letter and not contain parenthesis or -
filePathForSamplesPlan="${dirPathWithConfigFiles}/samples_plan_both.txt"
# Sometimes some chromosomes are excluded from analysis:
# chrsToRemove="chrX,chrY,chrM"
chrsToRemove="chrM"
# For the plots you can use columns in your samplesplan
column1=Tissue
column2=Library
# If you don't have
# column1=sample
# column2=sample
# If you don't want to do deseq2 just comment:
deseqcolumn=Tissue

# If you are running 2 analysis in the same directory use suffix for config files and output directories:
suffix="_both"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load GCC/8.2.0-2.31.1  
# module load OpenMPI/3.1.3
# module load R/3.6.0

# You can choose to use a conda environment to solve R and associated packages
# You can create it with: conda create -n RNAseq_R_202302 r-base r-colorspace bioconductor-deseq2 r-ggplot2 r-pheatmap r-rcolorbrewer  bioconductor-rtracklayer
# Comment it if you will use module load
# condaEnvName=RNAseq_R_202302
######

# You can choose to use singularity, then you need to define each function:
pathToImages=/home/share/andrey_lab/singularity/
# Or use a personnal space
# pathToImages=/home/users/d/delislel/scratch/images/

verse_with_more_packages_version=4.4.3_0

if [ ! -e ${pathToImages}/verse_with_more_packages_${verse_with_more_packages_version}.sif ]; then
  # Get the docker image
  cwd=$PWD
  cd ${pathToImages}
  export APPTAINER_CACHEDIR=$PWD/.cache
  singularity pull docker://lldelisle/verse_with_more_packages:${verse_with_more_packages_version}
  cd $cwd
fi
function R() {
  singularity exec "${pathToImages}/verse_with_more_packages_${verse_with_more_packages_version}.sif" R $*
}
function Rscript() {
  singularity exec "${pathToImages}/verse_with_more_packages_${verse_with_more_packages_version}.sif" Rscript $*
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
R --version
if [ $? -ne 0 ]
then
  echo "R is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
# Check  are installed:
Rscript -e "library(colorspace);library(DESeq2);library(ggplot2);library(pheatmap);library(RColorBrewer);library(rtracklayer)"
if [ $? -ne 0 ]
then
  echo "Some R packages are missing check rtracklayer, colorspace, DESeq2, ggplot2, pheatmap and RColorBrewer are installed."
  exit 1
fi

# Check the 2 github:
if [ ! -e ${dirPathWithDependencies}/rnaseq_rscripts/ ]; then
  echo "${dirPathWithDependencies}rnaseq_rscripts/ does not exists please clone https://github.com/lldelisle/rnaseq_rscripts"
  exit 1
fi
cd ${dirPathWithDependencies}/rnaseq_rscripts/
echo "Version of rnaseq_rscripts"
git rev-parse HEAD

if [ ! -e ${dirPathWithDependencies}/toolBoxForMutantAndWTGenomes/ ]; then
  echo "${dirPathWithDependencies}/toolBoxForMutantAndWTGenomes/ does not exists please clone https://github.com/lldelisle/toolBoxForMutantAndWTGenomes"
  exit 1
fi
cd ${dirPathWithDependencies}/toolBoxForMutantAndWTGenomes/
echo "Version of toolBoxForMutantAndWTGenomes"
git rev-parse HEAD


## START
cd "${dirPathWithConfigFiles}"

# This script will merge the tables

# First generate tables
if [ ! -e ${dirPathWithResults}/mergedTables${suffix}/AllCufflinks_Simplified_subset.txt ]; then
  configFile="configFileRNAseq_step1${suffix}.R"

  echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"

#### STEP 1 - MERGE TABLES ####
# If the merged tables are not already generated:
outputFolderForStep1 <- \"${dirPathWithResults}/mergedTables${suffix}/\"
# Needed for DESeq2: Do you want to merge counts? T=yes F or commented=no
mergeCounts <- T
# Optional: subset the count table Do you want to remove some genes from the
# count table
subsetCounts <- T
# If the table with counts have already been generated and you just want to
# remove some genes.
# initialTableWithCount<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTables${suffix}/AllHTSeqCounts.txt'
# If you provide the initialTableWithCount you need to provide the name of the
# column with the ensembl id.
# geneIDColInInitialTable<-'Ens_ID'
# List of genes id to remove (one per line with no header).
genesToRmFromCounts <- \"$PWD/genes${chrsToRemove//,/_}.txt\"
# Optional:
mergeFPKM <- T
# By default cufflinks split the transcripts which do not overlap in different
# locus and so different lines, put T if you want to sum the FPKM for non
# overlapping transcripts (put F if not).
oneLinePerEnsemblID <- T
# Optional: subset the FPKM table Do you want to remove some genes from the
# FPKM table
subsetFPKM <- T
chrToRemove <- c(\"${chrsToRemove//,/\",\"}\")
# Anouk method: Genes that have the less variable rank should have the same
# expression.
normFPKMWithAnoukMethod <- T
# If the table with FPKM have already been generated and you just want to
# normalize it.
# initialTableWithFPKM<-'${dirPathWithDependencies}/rnaseq_rscripts/example/mergedTables${suffix}/AllCufflinks_Simplified.txt'
# Usually, it is recommanded to remove mitochondrial genes before doing the
# normalization. In some cases, it can also be useful to remove the sex
# chromosomes (put c('chrX','chrY','chrM')).  If you do not want to remove any
# gene put NA or comment the line.
chrToRemoveBeforeNormWithAnoukMethod <- c(\"${chrsToRemove//,/\",\"}\")
# Default is 1000, you can change here.
nbOfGenesWithAnoukMethod <- 1000
# If you want to keep the genes used in the normalization from Anouk, they will
# be written in a file.
keepGenesUsedForNorm <- F
" > ${configFile}

  if [ ! -e genes${chrsToRemove//,/_}.txt ]; then
    Rscript ${dirPathWithDependencies}/toolBoxForMutantAndWTGenomes/scripts/getGeneListFromChrAndGTF.R $filePathForGTF ${chrsToRemove} ./
    mv genesIn${chrsToRemove}from* genes${chrsToRemove//,/_}.txt
  fi
  # Adjust the samplesplan
  if [ ! $(grep "htseq_count_file" $filePathForSamplesPlan) ]; then
    cat $filePathForSamplesPlan | awk -v pa=$dirPathWithResults 'BEGIN{print "htseq_count_file"}NR==1{for (i=1;i<=NF;i++){if($i=="sample"){col=i}}}NR>1{print pa"/allFinalFiles/counts_FPKM/htseqCount_"$col".txt"}' > htseqCol.txt
    paste -d "\t" $filePathForSamplesPlan htseqCol.txt > ${filePathForSamplesPlan}_withPaths
    rm htseqCol.txt
  fi
  if [ ! $(grep "cufflinks_file" $filePathForSamplesPlan) ]; then
    cat $filePathForSamplesPlan | awk -v pa=$dirPathWithResults 'BEGIN{print "cufflinks_file"}NR==1{for (i=1;i<=NF;i++){if($i=="sample"){col=i}}}NR>1{print pa"/allFinalFiles/counts_FPKM/FPKM_"$col".txt"}' > CuffCol.txt
    if [ -e ${filePathForSamplesPlan}_withPaths ];then
      mv ${filePathForSamplesPlan}_withPaths tmp
      paste -d "\t" tmp CuffCol.txt > ${filePathForSamplesPlan}_withPaths
      rm tmp
    else
      paste -d "\t" $filePathForSamplesPlan CuffCol.txt > ${filePathForSamplesPlan}_withPaths
    fi
    rm CuffCol.txt
  fi
  if [ -e ${filePathForSamplesPlan}_withPaths ]; then
    cat $configFile | sed "s#$filePathForSamplesPlan#${filePathForSamplesPlan}_withPaths#" > ${configFile}_withPaths
    configFile=${configFile}_withPaths
  fi

  Rscript ${dirPathWithDependencies}/rnaseq_rscripts/step1-generateTables.R $configFile
  # copy to GEO:
  mkdir -p ${dirPathWithResults}/toGEO
  cp ${dirPathWithResults}/mergedTables${suffix}/AllCufflinks_Simplified.txt ${dirPathWithResults}/toGEO/AllCufflinks_Simplified${suffix}.txt
  cp ${dirPathWithResults}/mergedTables${suffix}/AllCufflinks_Simplified_subset.txt ${dirPathWithResults}/toGEO/AllCufflinks_Simplified_subset${suffix}.txt
  cp ${dirPathWithResults}/mergedTables${suffix}/AllHTSeqCounts_subset.txt ${dirPathWithResults}/toGEO/AllHTSeqCounts_subset${suffix}.txt
fi


# Then basic plots

if [ ! -e ${dirPathWithResults}/plots_default${suffix}/PC1.pdf ]; then
  configFile="configFileRNAseq_step3${suffix}.R"

  echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"



#### STEP 3 - PLOTS #### Required You can put here either the FPKM norm values (subsetted or not)
#### or the count norm values obtained after DESeq2
tableWithNormalizedExpression <- \"${dirPathWithResults}/mergedTables${suffix}/AllCufflinks_Simplified_subset.txt\"
# In case you are using a file with both raw counts and FPKM you need to choose
# which values you want to plot.  If set to T, only columns called FPKM_sample
# will be used.
useFPKM <- T
# Optional
outputFolder <- \"${dirPathWithResults}/plots_default${suffix}/\"
# By default pdf is used as output. If set to T, png will be used.
usePng <- F
# You can provide color for each value of each factor in you samples plan to
# have constistant graphs for PCA, correlation and genes.


### Common to PCA and clustering ### In DESeq2 they restrict to the 500 most
### variant genes. If you want to keep all genes, comment the line or put
### 1000000.
restrictToNMoreVariantGenes <- 500

### PCA ### Put here the number of PC you want to see (0=do not perform PCA,
### 1=Only look at first component, 2=look at the 2 first etc...)
nbOfPC <- 3
# If the nbOfPC is greater than 1, you will have a barchart of each PC and you
# may want to use different parameters to identify your samples using the
# column names of the samples plan.
PCA1D <- list(fill = \"${column1}\", color = \"${column2}\")
# Possible personalizations are : fill is for the color of the bar, alpha is
# for the transparency, color is for the color of the border of the bar
# linetype is for the type of border of the bar If you do not want to use one
# of the parameter, just remove it from the list.

# If the nbOfPC is greater than 2, you will have projection in 2 dimension and
# to identify your sample you may want to use the column names of the samples
# plan.
PCA2D <- list(color = \"${column1}\", shape = \"${column2}\")
# Possible personalizations are : color is for the color of the symbol, alpha
# is for the transparency, shape is for the shape of the symbol. You can also
# choose 2 colors fill and color. If so, the fill will be inside and the color
# will be the border. If you do not want to use one of the parameter, just
# remove it from the list.

# Do you want to have the contribution of each gene to each PC (T=yes, F=no).
getGeneContributionToPCA <- F

### Clustering ###

# Do you want to perform a correlation matrix and clustering (T=yes, F=no)
plotMatrixAndClustering <- T

### Genes ###
# One gene per line. The first line of the gene file should correspond to a
# column in the expression file.
# fileWithGenes <- \"~/rnaseq_rscripts/example/genesHoxDandAround.txt\"
# By default, the title of the plot is the id provided in the fileWithGenes but
# you can add a meaning full name like gene_short_name if it is provided in the
# tableWithNormalizedExpression.
# geneIDToAdd<-'gene_short_name'
# By default, the values of expression plotted are log2(1+expression) (when T,
# if F the raw expression will be plotted.)
useLogExpression <- T
# By default, each gene is plotted on an adjusted scale. If
# useSameYmaxForAllGenes is T, all genes will be plotted with the same y axis.
useSameYmaxForAllGenes <- T
# A factor which will be used as x axis.
xaxisForGenes <- \"Tissue\"
plotGenesPara <- list(color = \"Line\", shape = \"Replicate\")
# Possible personalizations are :
# color is for the color of the symbol, 
# alpha is for the transparency,
# shape is for the shape of the symbol.
# You can also choose 2 colors fill and color. If so, the fill will be inside and the color will be the border.
# If you do not want to use one of the parameter, just remove it from the list.

# If you only want a heatmap and not one gene per one gene. Put it to T.
doNotPlotGeneByGene <- F
# If you want to have a heatmap with the genes provided in the list. All values
# of plotGenesPara will be used to annotate the samples.
addGlobalHeatmap <- T
# By default, genes are clustered by euclidean distance and complete
# clustering. If you want to keep the original order. Put keepGeneOrder to T.
keepGeneOrder <- F
# By default, samples are not clustered.
clusterSamples <- F
" > ${configFile}
  Rscript ${dirPathWithDependencies}/rnaseq_rscripts/step3-graphClusteringPCAGenes.R $configFile
fi


if [ ! -e ${dirPathWithResults}/plots_default${suffix}/PC1.png ]; then
  configFile="configFileRNAseq_step3${suffix}_png.R"
  sed "s/usePng <- F/usePng <- T/" configFileRNAseq_step3${suffix}.R > ${configFile}
  Rscript ${dirPathWithDependencies}/rnaseq_rscripts/step3-graphClusteringPCAGenes.R $configFile
fi

# Same graphs but only on protein coding genes:

# Extract protein coding:
dir=${dirPathWithResults}/mergedTables${suffix}/
if [ ! -e ${dir}/AllHTSeqCounts_subset_onlypc.txt ]; then
  Rscript ${dirPathWithDependencies}/toolBoxForMutantAndWTGenomes/scripts/subsetForProteinCoding.R "${filePathForGTF}" ${dir}/AllHTSeqCounts_subset.txt Ens_ID
fi
if [ ! -e ${dir}/AllCufflinks_Simplified_subset_onlypc.txt ]; then
  Rscript ${dirPathWithDependencies}/toolBoxForMutantAndWTGenomes/scripts/subsetForProteinCoding.R "${filePathForGTF}" ${dir}/AllCufflinks_Simplified_subset.txt gene_id
fi

# Then basic plots using protein_coding

if [ ! -e ${dirPathWithResults}/plots_default_onlypc${suffix}/PC1.pdf ]; then
  configFile="configFileRNAseq_step3_pc${suffix}.R"

  echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"



#### STEP 3 - PLOTS #### Required You can put here either the FPKM norm values (subsetted or not)
#### or the count norm values obtained after DESeq2
tableWithNormalizedExpression <- \"${dirPathWithResults}/mergedTables${suffix}/AllCufflinks_Simplified_subset_onlypc.txt\"
# In case you are using a file with both raw counts and FPKM you need to choose
# which values you want to plot.  If set to T, only columns called FPKM_sample
# will be used.
useFPKM <- T
# Optional
outputFolder <- \"${dirPathWithResults}/plots_default_onlypc${suffix}/\"
# By default pdf is used as output. If set to T, png will be used.
usePng <- F
# You can provide color for each value of each factor in you samples plan to
# have constistant graphs for PCA, correlation and genes.


### Common to PCA and clustering ### In DESeq2 they restrict to the 500 most
### variant genes. If you want to keep all genes, comment the line or put
### 1000000.
restrictToNMoreVariantGenes <- 500

### PCA ### Put here the number of PC you want to see (0=do not perform PCA,
### 1=Only look at first component, 2=look at the 2 first etc...)
nbOfPC <- 3
# If the nbOfPC is greater than 1, you will have a barchart of each PC and you
# may want to use different parameters to identify your samples using the
# column names of the samples plan.
PCA1D <- list(fill = \"${column1}\", color = \"${column2}\")
# Possible personalizations are : fill is for the color of the bar, alpha is
# for the transparency, color is for the color of the border of the bar
# linetype is for the type of border of the bar If you do not want to use one
# of the parameter, just remove it from the list.

# If the nbOfPC is greater than 2, you will have projection in 2 dimension and
# to identify your sample you may want to use the column names of the samples
# plan.
PCA2D <- list(color = \"${column1}\", shape = \"${column2}\")
# Possible personalizations are : color is for the color of the symbol, alpha
# is for the transparency, shape is for the shape of the symbol. You can also
# choose 2 colors fill and color. If so, the fill will be inside and the color
# will be the border. If you do not want to use one of the parameter, just
# remove it from the list.

# Do you want to have the contribution of each gene to each PC (T=yes, F=no).
getGeneContributionToPCA <- F

### Clustering ###

# Do you want to perform a correlation matrix and clustering (T=yes, F=no)
plotMatrixAndClustering <- T

### Genes ###
# One gene per line. The first line of the gene file should correspond to a
# column in the expression file.
# fileWithGenes <- \"~/rnaseq_rscripts/example/genesHoxDandAround.txt\"
# By default, the title of the plot is the id provided in the fileWithGenes but
# you can add a meaning full name like gene_short_name if it is provided in the
# tableWithNormalizedExpression.
# geneIDToAdd<-'gene_short_name'
# By default, the values of expression plotted are log2(1+expression) (when T,
# if F the raw expression will be plotted.)
useLogExpression <- T
# By default, each gene is plotted on an adjusted scale. If
# useSameYmaxForAllGenes is T, all genes will be plotted with the same y axis.
useSameYmaxForAllGenes <- T
# A factor which will be used as x axis.
xaxisForGenes <- \"Tissue\"
plotGenesPara <- list(color = \"Line\", shape = \"Replicate\")
# Possible personalizations are :
# color is for the color of the symbol, 
# alpha is for the transparency,
# shape is for the shape of the symbol.
# You can also choose 2 colors fill and color. If so, the fill will be inside and the color will be the border.
# If you do not want to use one of the parameter, just remove it from the list.

# If you only want a heatmap and not one gene per one gene. Put it to T.
doNotPlotGeneByGene <- F
# If you want to have a heatmap with the genes provided in the list. All values
# of plotGenesPara will be used to annotate the samples.
addGlobalHeatmap <- T
# By default, genes are clustered by euclidean distance and complete
# clustering. If you want to keep the original order. Put keepGeneOrder to T.
keepGeneOrder <- F
# By default, samples are not clustered.
clusterSamples <- F
" > ${configFile}
  Rscript ${dirPathWithDependencies}/rnaseq_rscripts/step3-graphClusteringPCAGenes.R $configFile
fi

if [ ! -e ${dirPathWithResults}/plots_default_onlypc${suffix}/PC1.png ]; then
  configFile="configFileRNAseq_step3_pc${suffix}_png.R"
  sed "s/usePng <- F/usePng <- T/" configFileRNAseq_step3_pc${suffix}.R > ${configFile}
  Rscript ${dirPathWithDependencies}/rnaseq_rscripts/step3-graphClusteringPCAGenes.R $configFile
fi

if [ ! -z ${deseqcolumn} ]; then
  if [ ! -e ${dirPathWithResults}/DESeq2${suffix}/DESeq2Analysis.txt ]; then
    configFile="configFileRNAseq_step2${suffix}.R"

    echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"


#### STEP 2 - DESEQ 2 ANALYSIS ####
# Required
tableWithCounts <- \"${dirPathWithResults}/mergedTables${suffix}/AllHTSeqCounts_subset.txt\"
# Specify here the name of the column which contains the gene IDs (they need to
# be unique).
geneIDColCounts <- \"Ens_ID\"
# For the DESeq2 analysis you need to specify a factor on which you want to do
# the analysis: This needs to be a name of a column of the samplesPlan file.
factor <- \"${deseqcolumn}\"

# Optional
# This can be table from cufflinks or cuffdiff or Biomart to annotate genes.
# You will need to choose a file with at least one column with the Ensembl Gene IDs.
tableWithAnnotations <- \"${dirPathWithResults}/mergedTables${suffix}/AllCufflinks_Simplified_norm.txt\"
# Specify here the name of the column which contains the gene IDs (it must
# match with the content of the geneID from the table with counts).
geneIDColInAnnotations <- \"gene_id\"
# You can also provide a gtf:
gtfFile <- \"${filePathForGTF}\"
# Default test is Wald but you can change to likelihood ratio test (LRT) with
# reduced formula ~1. Put F to keep Wald and put T to use LRT.
changeTest <- F
outputDESeqTable <- \"${dirPathWithResults}/DESeq2${suffix}/DESeq2Analysis.txt\"
# If you want to have another table with only significant genes abs(l2FC) > 1.5
# and corrected p-value < 0.05
outputSignificantTable <- T

" > ${configFile}
    Rscript ${dirPathWithDependencies}/rnaseq_rscripts/step2-DESeq2.R $configFile
  fi

  if [ ! -e ${dirPathWithResults}/DESeq2_pc${suffix}/DESeq2Analysis.txt ]; then
    configFile="configFileRNAseq_step2_pc${suffix}.R"

    echo "
### Required for all steps ###
RNAseqFunctionPath <- \"${dirPathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- \"${filePathForSamplesPlan}\"


#### STEP 2 - DESEQ 2 ANALYSIS ####
# Required
tableWithCounts <- \"${dirPathWithResults}/mergedTables${suffix}/AllHTSeqCounts_subset_onlypc.txt\"
# Specify here the name of the column which contains the gene IDs (they need to
# be unique).
geneIDColCounts <- \"Ens_ID\"
# For the DESeq2 analysis you need to specify a factor on which you want to do
# the analysis: This needs to be a name of a column of the samplesPlan file.
factor <- \"${deseqcolumn}\"

# Optional
# This can be table from cufflinks or cuffdiff or Biomart to annotate genes.
# You will need to choose a file with at least one column with the Ensembl Gene IDs.
tableWithAnnotations <- \"${dirPathWithResults}/mergedTables${suffix}/AllCufflinks_Simplified_norm.txt\"
# Specify here the name of the column which contains the gene IDs (it must
# match with the content of the geneID from the table with counts).
geneIDColInAnnotations <- \"gene_id\"
# You can also provide a gtf:
gtfFile <- \"${filePathForGTF}\"
# Default test is Wald but you can change to likelihood ratio test (LRT) with
# reduced formula ~1. Put F to keep Wald and put T to use LRT.
changeTest <- F
outputDESeqTable <- \"${dirPathWithResults}/DESeq2_pc${suffix}/DESeq2Analysis.txt\"
# If you want to have another table with only significant genes abs(l2FC) > 1.5
# and corrected p-value < 0.05
outputSignificantTable <- T

" > ${configFile}
    Rscript ${dirPathWithDependencies}/rnaseq_rscripts/step2-DESeq2.R $configFile
  fi
fi
