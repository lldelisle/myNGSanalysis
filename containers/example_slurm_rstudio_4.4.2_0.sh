#!/bin/sh
#SBATCH --time=12:00:00 # This is the maximum length of your R session
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G # Here change to the memory you think you will use
#SBATCH --output=rstudio-server.job.%j.out
#SBATCH --error=rstudio-server.job.%j.err

### CUSTOMIZE THIS PART ###

# Put here the path of your sif file:
SIF=$HOME/verse_with_more_packages_4.4.2_0.sif

# Put here the path where you want your packages to be installed
# to be able to have them the next time you start a new container
# But best practices is to change the image!
# Choose a path that is specific to this sif file
LIBDIR=${HOME}/R/rocker-rstudio/4.4.2_0

# Apptainer always mounts the HOME directory
# Put here the other path you want your container accesses
# For example the scratch if you are on baobab/SCITAS
# It must starts with comma and then add the paths separated by comma
# For SCITAS use:
# othermounts=",/scratch/$(id -un)/"

# Optional but useful for copy-paste:
name_host_machine=jed.epfl.ch

### END ###

mkdir -p $LIBDIR
# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
workdir=$(mktemp -d)

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
# libraries used by R) from spawning more threads than the number of processors
# allocated to the job.
#
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment

cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export OMP_NUM_THREADS=${SLURM_JOB_CPUS_PER_NODE}
export R_LIBS_USER=${LIBDIR}
exec /usr/lib/rstudio-server/bin/rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export APPTAINER_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server${othermounts}"

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0

export APPTAINERENV_USER=$(id -un)
export APPTAINERENV_PASSWORD=$(openssl rand -base64 15)

readonly PORT=$(shuf -i8040-9999 -n1)
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:

   ssh -N -f -L 8787:$(hostname -i):${PORT} ${APPTAINERENV_USER}@${name_host_machine}

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: ${APPTAINERENV_USER}
   password: ${APPTAINERENV_PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f ${SLURM_JOB_ID}

In case of unexpected behaviour check $workdir
END

singularity exec --cleanenv $SIF \
    /usr/lib/rstudio-server/bin/rserver --www-port ${PORT} \
            --www-address=$(hostname -i) \
            --server-user=$APPTAINERENV_USER \
            --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --auth-stay-signed-in-days=30 \
            --auth-timeout-minutes=0 \
            --rsession-path=/etc/rstudio/rsession.sh
printf 'rserver exited' 1>&2
