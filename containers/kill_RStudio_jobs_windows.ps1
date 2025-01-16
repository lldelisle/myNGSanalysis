$username_on_server = "mayran"
$address_of_server = "galaxyduboule.epfl.ch"

$check_squeue = "squeue -u $username_on_server -o '%.8i %.38j %.8T %.10m %.8M %.10l %.20S'"
Write-Host "Here are your jobs in the queue"
ssh ${username_on_server}@${address_of_server} $check_squeue 2> $null
Write-Host "Don't forget to Exit the RStudio Session ('power' button in the top right) before cancelling jobs"
$jobs_to_cancel = Read-Host "Enter the ids of the jobs you want to cancel (separated by comma):"
ssh ${username_on_server}@${address_of_server} scancel -f $jobs_to_cancel
Write-Host "DONE"
