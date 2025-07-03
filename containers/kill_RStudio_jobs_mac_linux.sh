username_on_server=your_login
name_host_machine=jed.epfl.ch
name_host_machine=192.168.202.69

check_squeue="squeue -u $username_on_server -o '%.8i %.38j %.8T %.10m %.8M %.10l %.20S'"
echo "Here are your jobs in the queue"
ssh ${username_on_server}@${name_host_machine} $check_squeue 2> /dev/null
echo "Don't forget to Exit the RStudio Session ('power' button in the top right) before cancelling the corresponding jobs"
echo ""
read -p "Enter the ids of the jobs you want to cancel (separated by comma):" jobs_to_cancel
ssh ${username_on_server}@${name_host_machine} scancel -f $jobs_to_cancel
echo "DONE"
