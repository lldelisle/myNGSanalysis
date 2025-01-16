username_on_scitas=ldelisle
name_host_machine=jed.epfl.ch
# port on your computer
my_local_port=8787
sbatch_scripts_on_host_machine=rstudio_version_project_SCITAS.sh
# Ask for the container version:
read -p "Enter the version of verse_with_more_packages you want to use [4.4.2_0]:" sif_version
sif_version=${sif_version:-4.4.2_0}
# Ask for the project:
read -p "Enter the name of your project (no space, leave empty if do not want to create a specific directory):" project_name
# Ask for the memory:
read -p "Enter the memory you need in GB [30]: " mem
mem=${mem:-30}
# Ask for the time:
read -p "Enter the time you need in hh:mm:ss [12:00:00]: " maxtime
maxtime=${maxtime:-12:00:00}
# Ask for the cpus:
read -p "Enter the number of CPU you need [1]: " cpus
cpus=${cpus:-1}
# Run the sbatch
run_batch=$(ssh ${username_on_scitas}@${name_host_machine} "sbatch -J rstudio_docker_${sif_version}_${project_name} -t $maxtime --mem-per-cpu=${mem}G -c $cpus $sbatch_scripts_on_host_machine $sif_version $project_name")
echo $run_batch

jobid="${run_batch##* }"

echo "Your job has been submitted and has job id $jobid"

get_ssh_tunnel_command="cat rstudio-server.job.${jobid}.err | grep ssh"

# Thanks to ChatGPT edition:

# Check job status in a loop with timeout
check_squeue="squeue -o '%.8i %.6u %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S'"
timeout=300 # 5 minutes timeout
elapsed_time=0

echo "Checking job status..."

while [ $elapsed_time -le $timeout ]
do
    ssh ${username_on_scitas}@${name_host_machine} $check_squeue 2> /dev/null
    status=$(ssh ${username_on_scitas}@${name_host_machine} "squeue -j $jobid -h -o '%T'")
    if [ -z $status ]; then
        echo "Something went wrong, your job has not been submitted or failed."
        # exit 1
    fi
    echo "Your job is $status"
    if [ $status = "RUNNING" ]; then
        tunnel_command=$(ssh ${username_on_scitas}@${name_host_machine} $get_ssh_tunnel_command 2> /dev/null)
        while [ -z "$tunnel_command" ]; do
            echo "RStudio has not started yet, probably because the singularity image is not here"
            echo "Here is the end of the error file"
            ssh ${username_on_scitas}@${name_host_machine} tail rstudio-server.job.${jobid}.err
            echo ""
            echo ""
            sleep 2s
            tunnel_command=$(ssh ${username_on_scitas}@${name_host_machine} $get_ssh_tunnel_command 2> /dev/null)
        done
        echo "Openning the tunnel"
        ${tunnel_command/8787:/${my_local_port}:}
        echo "Your job is running you can access RStudio server in your browser."
        echo "At http://localhost:${my_local_port}"
        echo "Your username and passwords are"
        ssh ${username_on_scitas}@${name_host_machine} cat rstudio-server.job.${jobid}.err | grep -B1 "password"
        exit 0
    fi
    sleep 10s
    elapsed_time=$((elapsed_time+10))
done

if [ $elapsed_time -ge $timeout ]; then
    echo "Your job is still PENDING"
    echo "You can check its status by"
    echo "ssh ${username_on_scitas}@${name_host_machine} \"squeue -j $jobid -h -o '%T'\""
    echo ""
    echo "You can check all running/queuing jobs by"
    echo "ssh ${username_on_scitas}@${name_host_machine} \"$check_squeue\" 2> /dev/null"
    echo ""
    echo "Once your job is running you need to open a tunnel"
    echo "Launch in your terminal the following command:"
    echo "ssh ${username_on_scitas}@${name_host_machine} $get_ssh_tunnel_command | sed \"s/8787:/${my_local_port}:/\" 2> /dev/null"
    echo ""
    echo "Launch in the terminal the returned command to open the tunnel and" 
    echo "you will be able to access RStudio server in your browser."
    echo "At http://localhost:${my_local_port}"
    echo "To get the user and password you need to run"
    echo "ssh ${username_on_scitas}@${name_host_machine} cat rstudio-server.job.${jobid}.err | grep -B1 \"password\""
fi
