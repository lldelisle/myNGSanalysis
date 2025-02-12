# Get fastq data from s3 with AWS
# The data are stored in a s3 bucket. We will use the AWS CLI to download them.

# Set directories
pathToFastq="$microcPilot/fastq/"
pathToImages="$SRC/images"

# Pull image
# aws
if [ ! -e $pathToImages/aws-cli_2.23.6.sif ]; then
    singularity pull docker://amazon/aws-cli:2.23.6
    mv aws-cli_2.23.6.sif $pathToImages/
fi

function aws() {
  singularity exec $pathToImages/aws-cli_2.23.6.sif aws $*
}

# Give access to pathToFastq
export APPTAINER_BIND=$SRC

# check it is properly installed
aws --version 

# Configure with the data provided form the company
aws configure
# Aws_access_key_id：
# Aws_secret_access_key：
# Region
# Defaul output format:

# Test if the data are found and summarize them
aws s3 ls --summarize --human-readable --recursive s3://musgzjor-598731762349/F24A430002451_MUSgzjoR_24DEC2024

# Download
aws s3 sync s3://musgzjor-598731762349/F24A430002451_MUSgzjoR_24DEC2024/ $pathToFastq > permanent_transfert.log 2> permanent_transfert.err 