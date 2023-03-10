# Containers Tutorial

Task for the Primers for Predocs containers session on 18th January 2023.

## Singularity

Singularity is a tool that allows us to create and run software containers and greatly improves the portability and reproducibility of our analyses.

## Setup

1) Log into the virtual machine as described by the course instructions.
2) The scripts and files for this session should already be available via Penelopeprime on your virtual machine. You will need to copy these to your machine using the following command:
```{bash}
cd ~/Documents
cp -r /media/penelopeprime/Primer4PredocsJan23/primers-containers .
cd primers-containers
```
If this is not the case or you want to try this tutorial on your own machine, you can get the data and scripts from the GitHub repository. No GitHub account should needed for this as it is a public repository. You can do this by opening the Terminal and running:
```{bash}
cd ~/Documents
git clone https://github.com/Danderson123/primers-containers
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

The software we want to include in our container is often specified by writing a singularity recipe as a module definition file ending with `.def`. [task1/minimap2.def](https://github.com/Danderson123/primers-containers/blob/master/task1/minimap2.def) is a recipe I have written to install minimap2 in a singularity image. Singularity containers are essentially completely bare instances of the specified operating system, so we have to specify every programme that we will need to install the software we want to put in our container. We will be downloading `curl` and `tar` so we can download and expand the binary.

### Building the singularity image

We can then build the singularity container using:
```{bash}
singularity build --fakeroot --force task1/minimap2.simg task1/minimap2.def
```
Fortunately, we no longer need root privileges to build the container so can do this on the cluster. If you already had a container built locally though, you can transfer it to the cluster using:
```{bash}
scp -r task1/minimap2.simg <EBI username>@codon-login.ebi.ac.uk:~
```

### Executing the container

Now we have built the image, let's test it out by mapping some SARS-CoV-2 Nanopore reads to a reference sequence and output a SAM alignment file. We can do this like so:
```{bash}
singularity run task1/minimap2.simg -a data/MN908947.3.fasta data/ERR5729799.fastq.gz > task1/ERR5729799_mapped.sam
```

## Task 2

### Background

Containers become even more convenient when you realise we can run multiple pieces of software from a single container. To demonstrate this I have written a second singularity recipe that installs minimap2 like before, but now also installs Samtools and Artemis via Conda, 2 tools that we can use to visualise how well the Nanopore reads map to the reference. Samtools is a package we need to convert our SAM file to a BAM file that we can examine with Artemis.

### Building the singularity image

Samtools is available through conda so we will install this in our image first, then use conda to install Artemis.

Let's build the singularity image like this:
```{bash}
singularity build --fakeroot task2/msa.simg task2/msa.def
```

### Executing the container

To map the Nanopore reads we can run:
```{bash}
singularity run task2/msa.simg minimap2 -a data/MN908947.3.fasta data/ERR5729799.fastq.gz > task2/ERR5729799_mapped.sam
```
Then to run Samtools to sort the SAM and make the BAM file:
```{bash}
singularity run task2/msa.simg samtools sort task2/ERR5729799_mapped.sam > task2/ERR5729799_mapped.bam
```
We now need to index the BAM with Samtools:
```{bash}
singularity run task2/msa.simg samtools index task2/ERR5729799_mapped.bam
```
Then to run artemis:
```{bash}
singularity run task2/msa.simg art
```

## Task 3

### Background

We do not actually need to repeat the building process everytime we want to build an image for a new tool available with a Singularity recipe. Singularity and Docker images are available, shareable and runnable from container registries. Container registries include:
* [DockerHub](https://hub.docker.com/)
* [SingularityHub](https://singularityhub.com/)
* [quay.io](https://quay.io/)
* [GitHub Container Registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
* [Biocontainers](https://biocontainers.pro/registry)<sup>[1](https://academic.oup.com/bioinformatics/article/33/16/2580/3096437?login=true)</sup>: Builds and hosts a large number of containers for Bioinformatics tools. All tools available via the conda Bioconda channel will have a Biocontainer.

## Running a remote container

Let's run a container for Samtools v1.14 available from the Quay container registry:
```{bash}
img="docker://quay.io/biocontainers/samtools:1.14--hb421002_0"
singularity -s exec "$img" samtools view -h task2/ERR5729799_mapped.bam
```

## Sandbox development of singularity images

Developing Singularity containers can be frustrating at times as we have to rebuild a container every time we run into an error. To make our lives a bit easier, we can use Singularity's `--sandbox` option to speed up development. The Singularity sandbox mimics a container environment so that we can test commands we want to add to a container, preventing us from having to rebuild a container and saving us a lot of time. We can develop containers using sandbox using the commands below:
```{bash}
mkdir sandbox-dev
cd sandbox-dev
singularity build --fakeroot --sandbox <container name> ../data/template.def
```
To start an interactive shell within the container we run:
```{bash}
singularity shell --writable <container name>
```
The `--writable` option allows us to install software and save our changes to the container.

We should now see our terminal update to:
```{bash}
Singularity playground:~>
```
If we want to make a recipe (rather than just a container itself), we will need to manually copy and paste successful sandbox commands to a separate document.

## Bonus: Converting a Dockerfile to a singularity container

Now that we have seen how useful singularity containers can be in helping us make software containers without root privileges, let's see how we can convert a Docker container, which does need root privileges to build and run, to a singularity container.

Unfortunately we cannot build a singularity container from a Dockerfile directly, but we can build one from a running docker container if we have root privileges like this:

```{bash}
docker build -t local/minimap2:latest bonus
sudo singularity build --force bonus/minimap2.simg docker-daemon://local/minimap2:latest
```

## Additional resources

For a much more in-depth tutorial, see Michael Hall's 2019 Singularity tutorial: https://github.com/mbhall88/eipp-2019-singularity#run-and-serving-applications
