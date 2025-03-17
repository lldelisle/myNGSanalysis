---
marp: true
author: Lucille Delisle
title: Reproducible analysis with R
style: |
  .columns {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 1rem;
  }
  .columns32 {
    display: grid;
    grid-template-columns: repeat(5, minmax(0, 1fr));
    gap: 1rem;
  }
  .item:nth-child(-n + 1) { grid-column: span 3; }
  .item:nth-child(n + 2) { grid-column: span 2; }
  .item {
    padding: 10px;
  }
---

# Reproducible analysis with R

Lucille Delisle

2025/03/18

---

# Who am I?

Lucille Delisle

PhD in biology

Bioinformatician since 2015 (nearly 10 years at EPFL in Duboule lab)

50% in Andrey lab and 50% in Herrera lab since March 2025

Highly motivated by open data open software open trainings

---

# Who are you?

- Which lab are you?
<small>Andrey, Braun, Herrera, Neerman-Arbez, Nef, Ruiz Altaba, Wu, Zdobnov, Frayling, Bertrand</small> 

- Are you PhD? Post-doc/Collab? Master? Other?

- What is your bioinformatic level?
    - None (never used R never used command lines)
    - Beginner (used a bit of R)
    - Intermediate (At ease with R)
    - Advanced (At ease with R and CLI but want to know about reproducibility)
    - Expert (At ease with R and command lines and containers)

---

# Schedule

- Theorical part

- Hands-on

  - Containers
  - git?

- Free part if time permits

---

# Reproducibility

Start from the same input data and be able to rerun the analysis to the figure/tables.

What is needed:

- input data
- scripts (should be stored in git)
- and???

---

# Reproducibility

Start from the same input data and be able to rerun the analysis to the figure/tables.

What is needed:

- input data
- scripts (should be stored in git)
- tool versions (R + all packages used)

---

# How to manage tool versions

Issues:
    - you may have multiple projects on going
    - you may want to collaborate on a project
    - you cannot reproduce results
    - you want to install a new packages

I think the solution is "Containers"

---

### What are containers?

I've been inspired by: [SIB course material](https://sib-swiss.github.io/containers-introduction-training/latest/) and [Carpentries course](https://carpentries-incubator.github.io/docker-introduction/).

Containers are:
- Virtualized environment: an isolated file system accessible from a host computer
- Other than virtual machines (VMs), containers have specific purposes, and carry only essential information to perform their task
- In IT terms:
  - containers share the kernel with the host OS (limits the reproducibility)
  - VMs bring the entire operating system

---

### Important vocabulary

<div class="columns32">
<div class="item item--1">
<!-- ![bg left:33%](../containers-cookie-cutter.png) -->
<!-- <style>
p { columns: 2; }
</style> -->
<!-- ![width:600px](../containers-cookie-cutter.png) -->
<img src="../containers-cookie-cutter.png">
<!-- ![width:600px](../containers-cookie-cutter.png) -->
</div>
<div class="item item--2">

Developers write a **Dockerfile** (how to make a cookie cutter). You can see examples on [my github](https://github.com/lldelisle/lldelisle-docker)

They run on their computer or on the cloud a command that generates a **Container Image** (you build the cookie cutter). This image is then published for example on dockerhub.

</div>
</div>

---

### Important vocabulary

<div class="columns32">
<div class="item item--1">
<!-- ![bg left:33%](../containers-cookie-cutter.png) -->
<!-- <style>
p { columns: 2; }
</style> -->
<!-- ![width:600px](../containers-cookie-cutter.png) -->
<img src="../containers-cookie-cutter.png">
<!-- ![width:600px](../containers-cookie-cutter.png) -->
</div>
<div class="item item--2">

Everyone can then download this image on its computer or on any server and then can create multiple **containers** (cookies).

</div>
</div>


You can personalize each cookie but then you loose the reproducibility. 
It is better to modify the Dockerfile -> image -> container

Container has a short life so don't spend time on customizing something that will last so shortly...

---

### What means 'isolated' in 'isolated file system'?

By default there in no communication between what is on your computer and what is in your container.

If you want to have access in your container to files in your computer, you need to specify before which are the directories which will be 'mounted'.

---

## Application to RStudio server

---

### R packages

#### Origin

R packages can have 3 origins:
- CRAN (The Comprehensive R Archive Network): 21599 available packages. Any package of good quality (documentation + tests) can go to CRAN. Each developer deals with the releases of its own package.
- Bioconductor: 2300 packages which are manually reviewed and related to biology. There are release cycles every 6 months that ensure that at a fixed release (for example Bioconductor release 3.19) all packages are compatible and work together. The releases are attached to R version: Bioconductor releases 3.17, 3.18 require R version 4.3, Bioconductor releases 3.19, 3.20 requires R version 4.4 (R version 4.5 is ).
- Github/Gitlab/others

---

#### Dependencies

Some R packages have no dependencies (they only use basic R functions). Most of them have dependencies from other packages and requires a minimal version of these packages. Some of them depends not only R packages but requires some system libraries to be able to compile a code (written in C or in Fortran for example) to speedup the calculations.

---

### Current situation for people already using R

#### Local RStudio

Everyone has RStudio on its laptop with fixed R version, fixed version of softwares. Sometimes you find your R too old or you need to upgrade a package because you want to use another package but which is not compatible with the old one or you want to use a new functionality... And you have some RDS files from previous version so you use them and then Lucille arrives and say you need to use the same version from start to end...

#### R on HPC?

Not necessarily the same version as what you have on your computer. Modules may not be updated as fast as necessary...

---

### Docker vs Singularity/Apptainer

Docker is the most popular container software. It has a GUI on Linux, MacOS and Windows. It is great for container development. It offers a large repository (docker hub) with a lot of base images including ubuntu, rstudio server... It requires administration access which is not possible for shared servers. So you can use docker on your laptop (see [this section](#use-on-your-computer-the-same-image-as-you-are-using-on-the-server)).

Singularity was the name of the original project but disagreement between the founder and the developers lead to two different projects: one adopted to the Linux Foundation changed its name to Apptainer while SingularityCE is the second one still developped at Sylabs. Apptainer/Singularity work without administration access. Image files are '.sif' (a file format develpped by SingularityCE). You can easily create a '.sif' from a docker image (and this is what most of developpers do). So on servers/linux/WSL2 we can use apptainer.

---

# Hands-on

---

## HPC

We are going to use the HPC of the university.

### Pricing

The [cost model](https://doc.eresearch.unige.ch/hpc/hpc_clusters?s[]=pricing#cost_model) has just changed. It was free. It is now free up to 100k CPU hours per PI, then 1.5cts per CPU hours.

---

## Intro

For UNIGE people everything is simplified thanks to the Open-on-demand.

Useful links to have more info:

- [Open-on-demand documentation from UNIGE](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand).
- [Open-on-demand architecture documentation](https://osc.github.io/ood-documentation/latest/architecture.html#container-context)

Open-on-demand is a service provided by the HPC in order to simplify access to HPC resources.

At UNIGE HPC has 3 clusters: baobab (old), bamboo (new), yggdrasil (astro). For the moment Open-on-demand is deployed on baobab (accessible without VPN) and bamboo (only accessible with VPN).

---

## Access the open-on-demand dashboard

Simply click on the [baobab](https://openondemand.baobab.hpc.unige.ch/) or [bamboo](https://openondemand.baobab.hpc.unige.ch/) link.

You will be redirected to Switch edu-ID that is used to login and link your HPC account.

For this you need to have a switch edu ID and have linked it to your UNIGE account. See [here](https://plone.unige.ch/distic/pub/compte-switch-edu-id/comment-creer-compte-switch-edu-id#EN).

---

## Use the default RStudio / default image

You directly arrives on your dashboard and you can click on RStudio.

For the training we will use an image available on HPC disk (the first one from the example values), partition shared-cpu, running time 02:00:00 memory 4 GB and 1 core.

(For a real analysis see [here](#available-images) to understand which are the possible images available but the first time you download and convert an image from docker it takes 5-10 minutes)

---

## Your first (?) HPC job / Rstudio on HPC

Your job is blue = 'Queued' and will become green 'Running'.

If your job is blue for long. You can check when it will be scheduled by clicking on Clusters > Shell Access and then

```bash
squeue -l -u $(id -u)
```

You can also check the occupancy of the cluster by

```bash
cat <(echo 'queue CPU_used CPU_free Node_used Node_free') \
<(sinfo -h -o "%R %C %A" | sed 's#\([^/]*\)/\([^/]*\)/.* \([^/]*\)/\([^/]*\)#\1 \2 \3 \4#g' | sort) \
| column -t
```

When your job is green you can click on "Connect to RStudio Server".

---

## Try it

```r
library(ggpubr)

set.seed(1234)
wdata = data.frame(
   sex = factor(rep(c("F", "M"), each=200)),
   weight = c(rnorm(200, 55), rnorm(200, 58)))
head(wdata, 4)

# Density plot with mean lines and marginal rug
ggdensity(wdata, x = "weight",
   add = "mean", rug = TRUE,
   color = "sex", fill = "sex",
   palette = c("#00AFBB", "#E7B800"))
```

---

## Check what is in 'Files'

The 'HOME' directory from the cluster is automatically mounted on the container you are currently using.

Your 'scratch' is also mounted.

---

## What is hidden


Small explanation on how it works and what is hidden behind the 'RStudio app' of open-on-demand that you ran.
- The code of the 'RStudio app' is not available to us currently (controlled by the HPC). The current version (2025-03-17) has been released 1 year ago and match what is [here](https://github.com/BioinfoSupport/ood-rstudio-baobab/tree/92e1182efc2cbd17a7095e166d20ef3bb0e616ec).
- You have created a container with RStudio server and with the packages that are in the docker image you chose. These packages are installed in the container. You can list them by:

```r
rownames(installed.packages())
```

---

## What is hidden


- The RStudio server from your container is running on the compute node and is available on a specific port. To avoid conflict between different users or different containers from the same user, this port is chosen randomly among not used ports.
- Clicking on the button 'Connect to RStudio Server' open a new window that enables you to 'see and interact' with the webpage on the port of the compute node.

---

## Can I Install new packages?

### Where the packages are installed?

```r
.libPaths()
```

The first directory is a directory that is on your home available by any R running with this version on this cluster.

Other directories are on the `/usr/` so on the container.

### Where would go a new packages if I install it?

You cannot write on `/usr/` as you are not a super user in this container.

Your new packages installation would go to the first path and you loose the reproducibility... but I have a solution... but before...

---

### When you're done with your work on RStudio

Do not forget to terminate your session (top right) when you are done.

I suggest to NEVER save your workspace if you don't want to wait for minutes next time you start a RStudio + Reproducibility.

To avoid paying for nothing and most importantly, allow other users to use the CPU. You should cancel your job(s) when you are done.

To do this, come back to your dashboard and click on 'Delete'.

---

## Use the last version of the OOD RStudio for more customization

I hope it will be available soon as the new default RStudio but today we need to set it up.

The instructions are [here](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand#sandbox).

First we need to enable custom Open-On-Demand applications.

From the dashboard, go to Clusters > Shell Access

In the new window, copy paste

```bash
mkdir -p $HOME/ondemand/dev
```

Restart your environnment "Help (?) > Restart web server"

---

## Use the last version of the OOD RStudio for more customization

Go to "develop (`</>`) > My Sandbox app"

Click on 'New App'

Click on 'Clone Existing App'

Directory name: `RStudio2025-03`

Git remote: `https://github.com/BioinfoSupport/ood-rstudio-baobab.git`

---

## Use the last version of the OOD RStudio for more customization

From your dashboard you should now have a new 'RStudio' but with a different logo, if you click on it it should be named "RStudio version: 3487caa"

The first fields are the same:
- Singularity image: `/acanas/m-BioinfoSupport/singularity/ngs_v1.1.sif` (available on HPC disk)
- Partition: `shared-cpu`
- Running time: `02:00:00`
- Memory: `4`
- Number of core: ` 1`

---

## Use the last version of the OOD RStudio for more customization

With this new application, you can:

- choose which account to charge for this job (useful if you are shared between labs).
- receive an email when your job start
- mount a smb (NASAC)
- get a more isolated session

---

### Mount your NASAC

To be able to automatically mount your nasac in the compute node and in your container you need to store your credentials in a file.

From the dashboard, go to Clusters > Shell Access

```bash
nano .credentials
```

Write
```txt
<username>
ISIS
<password>
```

Save with Ctl + O
Below it is written "File Name to Write: .credentials" hit enter

Then Crl + X to exit

---

### Get a more isolated RStudio

- In order to allow you to install new packages in an isolated place and be able to use them the day after, the bash script has created a new directory in your HOME of the server with the name of your project, and has told RStudio server to use packages in this directory and that new packages should be installed there.

- In order to be able to run RStudio on 2 independent projects in parallel


---

### Create a new container RStudio with all these options

The NASAC is mounted in `/run/user`, check it with 'Files'.

Check where packages are installed with

```r
.libPaths()
```

---

# Other RStudio images


The base images with rstudio server are described [here](https://rocker-project.org/images/), they are published by 'rocker' and based on ubuntu.

I used as template the `verse` image which has already Rstudio server + some tidyverse packages (dplyr, devtools, ggplot2, knitr, rmarkdown, stringr, tidyr...).

All my images are described [here](https://github.com/lldelisle/lldelisle-docker/blob/main/verse_with_more_packages/CHANGELOG.md).

A 4.4.2 version has been released in November 2024 and I am currently collecting the packages that the next image (4.4.3) should have, feel free to contribute.

Of course, the best would be that either we manage to do an image that fit everyone needs or everyone is able to push its own image to dockerhub (one image per project for example).

---

# Use on your computer the same image as you are using on the server

I wrote a documentation [here](https://github.com/lldelisle/myNGSanalysis/tree/main/containers#use-on-your-computer-the-same-image-as-you-are-using-on-the-server)

---

# Gitlab of UNIGE

Best place to store your code (and images) associated to it.

This is to avoid having:

- Analysisv1.R
- Analysisv2.R
- Analysisv3.R

or

- Analysis_202412.R
- Analysis_20250102.R
...

---

# Git

Git is a distributed version control system for managing
source code.

- Git is a software
- Git is useful when you have multiple versions of a file (works
better with text files).
  - You can go backward to a previous version.
  - You can check what changed between 2 versions.
- Git is distributed, that means that each user have all versions.

---

## Create a project

Go to [gitlab UNIGE](https://gitlab.unige.ch/).

Log in.

Create a project (New project > Create a blank project)

Project name: `demo-202503`
Project URL - Pick a group or namespace choose your user.
Visibility Level: keep `Private`

Click on `Create project`

---

## Set a SSH key on baobab/bamboo

From the dashboard of OpenOnDemand, Clusters > Shell Access

We will generate a SSH key = it ressembles password. There is a public key you can give to gitlab/github and a private key that is stored on a computer and must never be given.

```bash
ssh-keygen
```

Put enter to all questions.

Once done print the public key to the screen:

```bash
cat .ssh/id_rsa.pub
```

Or `.ssh/id_ed25519.pub`

---

## Give the public ssh key to gitlab UNIGE

On [gitlab.unige.ch](https://gitlab.unige.ch/), click on your "picture" then to "Preferences"

Click to "SSH Keys"

Click on "Add a new key", paste in the box the whole key starting from `ssh` to `login2.cluster`.

Title: `HPC baobab` or `HPC bamboo`

Expiration date X (no expiration)

---

## Connect it in RStudio

In the upper right corner in "Project: (None)"

New project

Choose Version Control > Git

Repository URL: `git@gitlab.unige.ch:<first name>.<last name>/demo-202503.git`

(You can find it from your project page in `Code \/`, clone with SSH)

Project directory name is automatically filled.

Create project as subdirectory of: I personally prefer to create directory `scripts` in my HOME where I clone all my git repositories.

---

## Create your first R script

File > New File > R Script

File > Save As > first.R

Write anything

```r
# Load ggplot package
library(ggplot2)
# Define variables
output.directory <- "plots"
# Create the output.directory
dir.create(output.directory, showWarnings = FALSE, recursive = TRUE)
# Do a dot plot
ggplot(mpg, aes(displ, hwy, colour = class)) + 
  geom_point()
ggsave(file.path(output.directory, "my_first_plot.png"), width = 7, height = 7)

```

Run line by line or 'Source'

---

## See your changes locally

Go to 'Files'

## Send your changes in gitlab

In order to send your changes to gitlab, go to the 'Git' tab.

Click on Diff

Click on the checkbox of each of your files. Write a commit message: THIS IS SUPER IMPORTANT `add a R script that plot a dotplot`

Commit and Close (= create a snapshot)

Push and Close (= send all my snapshots to gitlab)