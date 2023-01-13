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

## Task 1

### Background

[Minimap2](https://github.com/lh3/minimap2) is a fast sequence alignment tool to align DNA or mRNA sequences to a reference genome, built in C. Minimap2 has precompiled binaries available but what if we want to stick to a specific version and save our users from having to find the correct binary? We can install the binary in a singularity container!

### Defining the singularity container

The software we want to include in our container is often specified by writing a singularity recipe as a module definition file ending with `.def`. [task1/minimap2.def](https://github.com/Danderson123/primers-containers/blob/master/task1/minimap2-recipe.def) is a recipe I have written to install minimap2 in a singularity container. Singularity containers are essentially completely bare instances of the specified operating system, so we have to specify every programme that we will need to install the software we want to put in our container. We will be downloading `curl` and `tar` so we can download and expand the binary.

### Building the singularity container

We can then build the singularity container using:
```{bash}
singularity build --fakeroot --force task1/minimap2.img task1/minimap2.def
```
Fortunately, we no longer need root privileges to build the container so can do this on the cluster. If you already had a container built locally though, you can transfer it to the cluster using:
```{bash}
scp -r task1/minimap2.img <EBI USERNAME>@codon-login.ebi.ac.uk:~
```

### Executing the container

Now we have built the container, let's test it out by mapping some SARS-CoV-2 Nanopore reads to a reference sequence and output a SAM alignment file. We can do this like so:
```{bash}
singularity run task1/minimap2.img -a data/MN908947.3.fasta data/ERR5729799.fastq.gz > task1/ERR5729799_mapped.sam
```

## Task 2

### Background

Containers become even more convenient when you realise we can install multiple pieces of software in a single container. To demonstrate this I have written a second singularity recipe that installs minimap2 like before, but now also installs Samtools and Artemis via Conda, 2 tools that we can use to visualise how well the Nanopore reads map to the reference. Samtools is a package we need to convert our SAM file to a BAM file that we can examine with Artemis.

### Building the singularity container

Artemis is available through conda so we will install this in our container first, then use conda to install Artemis.

Let's build the singularity container like this:
```{bash}
singularity build --fakeroot --force task2/minimap2_and_artemis.img task2/minimap2_and_artemis.def
```

### Executing the container

To map the Nanopore reads we can run:
```{bash}
singularity run task2/minimap2_and_artemis.img minimap2 -a data/MN908947.3.fasta data/ERR5729799.fastq.gz > task2/ERR5729799_mapped.sam
```
Then to run Samtools to sort the SAM and make the BAM file:
```{bash}
singularity run task2/minimap2_and_artemis.img samtools sort task2/ERR5729799_mapped.sam > task2/ERR5729799_mapped.bam
```
We now need to index the BAM with Samtools:
```{bash}
singularity run task2/minimap2_and_artemis.img samtools index task2/ERR5729799_mapped.bam
```
Then to run artemis:
```{bash}
singularity run task2/minimap2_and_artemis.img art
```

## Bonus: Converting a Dockerfile to a singularity container

Now that we have seen how useful singularity containers can be in helping us make software containers without root privileges, let's see how we can convert a Docker container, which does need root privileges to build and run, to a singularity container.

Unfortunately we cannot build a singularity container from a Dockerfile directly, but we can build one from a running docker container if we have root privileges like this:

```{bash}
docker build -t local/minimap2:latest bonus
sudo singularity build --force bonus/minimap2.img docker-daemon://local/minimap2:latest
```
