# Useful files when working on a cluster

Here is a collection of files I put on my home of a cluster with slurm scheduler.

- The [bash_profile](./.bash_profile) indicates to add the bin directory in HOME to the PATH.

- The [bashrc](./.bashrc) contains simply an alias `squeueMe` with an output I like.

- The [tmux.conf](./.tmux.conf) contains an option to be able to use the mouse to change pane/windows and scroll down and also an option so that when you split the window you keep the same directory.

In the [bin](./bin) directory there are multiple useful 'functions':

- `lastMonth` allows to get the status and the memory usage of the jobs ran the last Month
- `lastWeek` is the same for the last week
- `multicol` allows to print on multiple columns and is used in most of the below functions:
- `wLastMin` allows to track the files that haved changed in the last minute. It highlights the changes in files. To stop it it Ctl + C.
- `wLast10Min` is the same for files that have changed in the last 10 minutes.
- `lsLastMin` is a helper for:
- `wLastMinLs` is like `wLastMin` but in addition it gives the size of the files.
- `myqstat` gives an overview of the state of the cluster with a summary on the nodes, a summary of your jobs and a summary of the pending jobs sorted by decreasing priority. It is mainly a helper for:
- `wqstat` allows to follow myqstat and highlight the differences.
