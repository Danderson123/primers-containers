# Containers Tutorial

Task for the Primers for Predocs containers session on 18th January 2023.

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

[Minimap2](https://github.com/lh3/minimap2) is a fast sequence alignment tool to align DNA or mRNA sequences to a reference genome, built in C. Minimap2 has precompiled binaries available but what if we want to stick to a specific version and save our users from having to find the correct binary? We can install the binary in a container!