$username_on_galaxyduboule = "mayran"
$sbatch_scripts_on_DubouleServer = "mountDuboule/Alex/docker/rstudio_version_project.sh"
# port on galaxyduboule.epfl.ch
# Hocine 8000 to 8009
# Anais 8010 to 8019
# Alex 8020 to 8029
# Lucille 8030 to 8039
$PORT = "8020"
# port on your computer
$my_local_port = "8787"

# Ask for the container version:
$sif_version = Read-Host "Which version of verse_with_more_packages you want to use?"

Write-Host "Check if $sif_version exists"
$exists = ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch "if [ -e 'mountDuboule/UPDUB COMMON/Sif_Images/verse_with_more_packages_${sif_version}.sif' ]; then echo Y; else echo N; fi" 2> $null
Write-Host "$exists"
if ( "$exists" -eq "N" ) {
    Write-Host "The version you asked does not exist in UPDUB COMMON"
    Write-Host "You need to pull it or ask Lucille"
    Exit
}
# Ask for the project:
$project_name = Read-Host "Which is the name of your project (no space, leave empty if do not want to create a specific directory)?"
# Run the sbatch
$run_batch = ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch "sbatch -J rstudio_docker_${sif_version}_${project_name} $sbatch_scripts_on_DubouleServer $PORT $sif_version $project_name"

Write-Host $run_batch

$jobid=$run_batch.Split()[-1]

Write-Host "Your job has id $jobid"

# Thanks to ChatGPT edition:

# Check job status in a loop with timeout
$check_squeue = "squeue -o '%.8i %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S'"
$timeout = 300 # 5 minutes timeout
$elapsed_time = 0

Write-Host "Checking job status..."
while ($elapsed_time -lt $timeout) {
    ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch $check_squeue 2> $null
    $status = ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch "squeue -j $jobid -h -o '%T'"
    if ($status -eq "") {
        Write-Host "Something went wrong, your job has not been submitted or failed."
        return
    }
    Write-Host "Your job is $status"
    if ($status -eq "RUNNING") {
        Write-Host "Your job is running you can access RStudio server in your browser at localhost:${my_local_port}."
        Write-Host "Do not forget to open the tunnel in a separate powershell with the following command"
        Write-Host "ssh -N -f -L ${my_local_port}:127.0.1.1:${PORT} ${username_on_galaxyduboule}@galaxyduboule.epfl.ch"
        return
    }
    Start-Sleep -Seconds 10
    $elapsed_time += 10
}

if ($elapsed_time -ge $timeout) {
    Write-Host "Your job is still PENDING"
    Write-Host "You can check its status by"
    Write-Host "ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch "squeue -j $jobid -h -o '%T'""
    Write-Host "You can check all running jobs by"
    Write-Host "ssh ${username_on_galaxyduboule}@galaxyduboule.epfl.ch $check_squeue 2> $null"
    Write-Host "Once your job is running you can access RStudio server in your browser at localhost:${my_local_port}."
    Write-Host "Do not forget to open the tunnel in a separate powershell with the following command"
    Write-Host "ssh -N -f -L ${my_local_port}:127.0.1.1:${PORT} ${username_on_galaxyduboule}@galaxyduboule.epfl.ch"
}
