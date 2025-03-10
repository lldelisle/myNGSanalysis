$username_on_galaxyduboule = "yourlogin"
$sbatch_scripts_on_DubouleServer = "nas/lab.data/common.scripts/rstudio_version_project_fixed_port.sh"
# port on 192.168.202.69
# Hocine 8000 to 8009
# Anaïs 8010 to 8019
# Yoann 8020 to 8029
# Lucille 8030 to 8039
# Ivana 8040 to 8049
# Cecilia 8050 to 8059
# ...
$PORT = "8030"
# port on your computer
$my_local_port = "8787"

# Ask for the container version:
$default = "4.4.2_0"
if (!($sif_version = Read-Host "Enter the version of verse_with_more_packages you want to use [$default]:")) { $sif_version = $default }

# Ask for the project:
$project_name = Read-Host "Enter the name of your project (no space, leave empty if do not want to create a specific directory):"

# Ask for the memory:
$default = "30"
if (!($mem = Read-Host "Enter the memory you need in GB [$default]:")) { $mem = $default }

# Ask for the time:
$default = "12:00:00"
if (!($maxtime = Read-Host "Enter the time you need in hh:mm:ss [$default]:")) { $maxtime = $default }

# Ask for the cpus:
$default = "1"
if (!($cpus = Read-Host "Enter the number of CPU you need [$default]:")) { $cpus = $default }

# Run the sbatch
$run_batch = ssh ${username_on_galaxyduboule}@192.168.202.69 "sbatch -J rstudio_docker_${sif_version}_${project_name} -t $maxtime --mem-per-cpu=${mem}G -c $cpus $sbatch_scripts_on_DubouleServer $PORT $sif_version $project_name"

Write-Host $run_batch

$jobid=$run_batch.Split()[-1]

Write-Host "Your job has id $jobid"

$get_ssh_tunnel_command = "cat rstudio-server.job.${jobid}.err | grep ssh"

# Thanks to ChatGPT edition:

# Check job status in a loop with timeout
$check_squeue = "squeue -o '%.8i %.6u %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S'"
$timeout = 300 # 5 minutes timeout
$elapsed_time = 0

Write-Host "Checking job status..."

while ($elapsed_time -lt $timeout) {
    ssh ${username_on_galaxyduboule}@192.168.202.69 $check_squeue 2> $null
    $status = ssh ${username_on_galaxyduboule}@192.168.202.69 "squeue -j $jobid -h -o '%T'"
    if ($status -eq "") {
        Write-Host "Something went wrong, your job has not been submitted or failed."
        return
    }
    Write-Host "Your job is $status"
    if ($status -eq "RUNNING") {
        $tunnel_command = ssh ${username_on_galaxyduboule}@192.168.202.69 $get_ssh_tunnel_command 2> /dev/null
        while ($tunnel_command -eq "") {
            Write-Host "RStudio has not started yet, probably because the singularity image is not here"
            Write-Host "Here is the end of the error file"
            ssh ${username_on_galaxyduboule}@192.168.202.69 tail rstudio-server.job.${jobid}.err
            Write-Host ""
            Write-Host ""
            Start-Sleep -Seconds 2s
            $tunnel_command = ssh ${username_on_galaxyduboule}@192.168.202.69 $get_ssh_tunnel_command 2> /dev/null
        }
        Write-Host "Your job is running you can access RStudio server in your browser at http://localhost:${my_local_port}."
        Write-Host "Your username and passwords are"
        ssh ${username_on_galaxyduboule}@192.168.202.69 cat rstudio-server.job.${jobid}.err | grep -B1 "password"

        Write-Host "Do not forget to open the tunnel in a separate powershell with the following command"
        Write-Host "ssh -N -f -L ${my_local_port}:127.0.1.1:${PORT} ${username_on_galaxyduboule}@192.168.202.69"
        return
    }
    Start-Sleep -Seconds 10
    $elapsed_time += 10
}

if ($elapsed_time -ge $timeout) {
    Write-Host "Your job is still PENDING"
    Write-Host "You can check its status by"
    Write-Host "ssh ${username_on_galaxyduboule}@192.168.202.69 "squeue -j $jobid -h -o '%T'""
    Write-Host ""
    Write-Host "You can check all running/queuing jobs by"
    Write-Host "ssh ${username_on_galaxyduboule}@192.168.202.69 $check_squeue 2> $null"
    Write-Host ""
    Write-Host "Once your job is running you can access RStudio server in your browser"
    Write-Host "At http://localhost:${my_local_port}."
    Write-Host "To get the user and password you need to run"
    Write-Host "ssh ${username_on_galaxyduboule}@192.168.202.69 cat rstudio-server.job.${jobid}.err | grep -B1 password"
    Write-Host "Do not forget to open the tunnel in a separate powershell with the following command"
    Write-Host "ssh -N -f -L ${my_local_port}:127.0.1.1:${PORT} ${username_on_galaxyduboule}@192.168.202.69"
}
