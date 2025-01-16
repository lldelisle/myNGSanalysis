# Make the ssh tunnel
echo "Openning tunnel"
ssh -N -f -L 8787:127.0.1.1:8000 rekaik@galaxyduboule.epfl.ch
# Ask for the container version:
echo "Which version of verse_with_more_packages you want to use?"
read sif_version
echo "Check if $sif_version exists"
exists=$(ssh rekaik@galaxyduboule.epfl.ch "if [ -e \"Duboule-server/UPDUB COMMON/Sif_Images/verse_with_more_packages_${sif_version}.sif\" ]; then echo 'Y'; else echo 'N'; fi" 2> /dev/null)
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
ssh rekaik@galaxyduboule.epfl.ch "sbatch -J rstudio_docker_${sif_version}_${project_name} Duboule-server/Hocine/rstudio_Hocine.sh $sif_version $project_name"
ssh rekaik@galaxyduboule.epfl.ch 'squeue -o "%.8i %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S"' 2> /dev/null
sleep 1s
ssh rekaik@galaxyduboule.epfl.ch 'squeue -o "%.8i %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S"' 2> /dev/null
sleep 3s
ssh rekaik@galaxyduboule.epfl.ch 'squeue -o "%.8i %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S"' 2> /dev/null
sleep 10s
ssh rekaik@galaxyduboule.epfl.ch 'squeue -o "%.8i %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S"' 2> /dev/null
sleep 60s
ssh rekaik@galaxyduboule.epfl.ch 'squeue -o "%.8i %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S"' 2> /dev/null
sleep 120s
ssh rekaik@galaxyduboule.epfl.ch 'squeue -o "%.8i %.38j %.4t %.10m %.8M %.10l %.6D %.3C %.12r %.20S"' 2> /dev/null
sleep 240s
echo "DONE"
