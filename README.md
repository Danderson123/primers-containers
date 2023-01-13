# Containers Tutorial

Task for the Primers for Predocs containers session on 18th January 2023.

## Singularity

Singularity is a tool that allows us to create and run software containers and greatly improves the portability and reproducibility of our analyses.

## Setup

1) Log into the virtual machine as described by the course instructions.
2) We need to download the scripts and data in this repository. No GitHub account should needed for this tutorial as this is a public repository. You can do this by opening the Terminal and running:
```{bash}
cd ~/Documents
git checkout https://github.com/Danderson123/primers-containers
cd primers-containers
```
3) Singularity and Docker should already be installed on your virtual machine. Please test these using:
```{bash}
singularity --version
docker --version
```

## Task

### Background

[Minimap2](https://github.com/lh3/minimap2) is a fast sequence alignment tool to align DNA or mRNA sequences to a reference genome, built in C. Minimap2 has precompiled binaries available but what if we want to stick to a specific version and save our users from having to find the correct binary? We can install the binary in a singularity container!

### Defining the singularity container

The software we want to include in our container is often specified by writing a singularity recipe as a module definition file ending with `.def`. [task/minimap2.def](https://github.com/Danderson123/primers-containers/blob/master/task/minimap2-recipe.def) is a recipe I have written to install minimap2 in a singularity container. Singularity containers are essentially completely bare instances of the specified operating system, so we have to specify every programme that we will need to install the software we want to put in our container. We will be downloading `curl` and `tar` so we can download and expand the binary.

### Building the singularity container

We can then build the singularity container using:
```{bash}
singularity build --fakeroot --force task/minimap2.img task/minimap2.def
```
Fortunately, we no longer need root privileges to build the container so can do this on the cluster. If you already had a container built locally though, you can transfer it to the cluster using:
```{bash}
scp -r ~/Documents/primers-containers/task/minimap2.img <EBI USERNAME>@codon-login.ebi.ac.uk:~
```

### Executing the container

Now we have built the container, let's test it out by mapping some SARS-CoV-2 Nanopore reads to a reference sequence and output a SAM alignment file. We can do this like so:
```{bash}
singularity run task/minimap2.img task/MN908947.3.fasta task/ERR5729799.fastq.gz > task/ERR5729799_mapped.sam
```

## Bonus: Converting a Dockerfile to a singularity container

Now that we have seen how useful singularity containers can be in helping us make software containers without root privileges, let's see how we can convert a Docker container, which does need root privileges to build and run, to a singularity container.

Unfortunately we cannot build a singularity container from a Dockerfile directly, but we can build one from a running docker container if we have root privileges like this:

```{bash}
docker build -t local/minimap2:latest bonus
singularity build --force bonus/minimap2.img docker-daemon://local/minimap2:latest
```

