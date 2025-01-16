username_on_galaxyduboule=ldelisle
sbatch_scripts_on_DubouleServer=mountDuboule/Lucille/rstudio_Lucille.sh
# port on galaxyduboule.epfl.ch
# Hocine 8000 to 8009
# Anaïs 8010 to 8019
# Alex 8020 to 8029
# Lucille 8030 to 8040
PORT=8030
# port on your computer
my_local_port=8787
# Make the ssh tunnel
echo "Openning tunnel"
ssh -N -f -L ${my_local_port}:127.0.1.1:${PORT} ${username_on_galaxyduboule}@galaxyduboule.epfl.ch
# Ask for the container version:
echo "Which version of verse_with_more_packages you want to use?"
read sif_version
echo "Check if $sif_version exists"
exists=$(ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch "if [ -e \"mountDuboule/UPDUB COMMON/Sif_Images/verse_with_more_packages_${sif_version}.sif\" ]; then echo 'Y'; else echo 'N'; fi" 2> /dev/null)
echo $exists
if [ "$exists" = "N" ]; then
    echo "The version you asked does not exist in UPDUB COMMON"
    echo "You need to pull it or ask Lucille"
    exit 1
fi
# Ask for the project:
echo "Which is the name of your project (no space, leave empty if do not want to create a specific directory)?"
read project_name
# Run the sbatch
run_batch=$(ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch "sbatch -J rstudio_docker_${sif_version}_${project_name} $sbatch_scripts_on_DubouleServer $PORT $sif_version $project_name")
echo $run_batch

jobid="${run_batch##* }"

echo "Your job has id $jobid"

# Thanks to ChatGPT edition:

# Check job status in a loop with timeout
check_squeue="squeue -o '%.8i %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S'"
timeout=300 # 5 minutes timeout
elapsed_time=0

echo "Checking job status..."

while [ $elapsed_time -le $timeout ]
do
    ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch $check_squeue 2> /dev/null
    status=$(ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch "squeue -j $jobid -h -o '%T'")
    if [ -z $status ]; then
        echo "Something went wrong, your job has not been submitted or failed."
        exit 1
    fi
    echo "Your job is $status"
    if [ $status = "RUNNING" ]; then
        echo "Your job is running you can access RStudio server in your browser."
        exit 0
    fi
    sleep 10s
    elapsed_time=$((elapsed_time+10))
done

if [ $elapsed_time -ge $timeout ]; then
    echo "Your job is still PENDING"
    echo "You can check its status by"
    echo "ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch "squeue -j $jobid -h -o '%T'""
    echo "You can check all running jobs by"
    echo "ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch $check_squeue 2> /dev/null"
    echo "Once your job is running you can access RStudio server in your browser."
fi
