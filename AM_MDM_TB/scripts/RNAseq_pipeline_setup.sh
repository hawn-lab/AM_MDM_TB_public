#!/bin/bash

##### Installation and setup to run RNAseq_pipeline.sh ######
##### Amazon Linux AMI #####
# Increase base EBS to 32 GB
# Add 2 TB addtl data EBS

#Note: using Anaconda build comes with pre-installed conda but this clutters the 
#environment and makes subsequent package installs VERY slow. Better to do basic
#Amazon AMI and install conda

########################################
## Update OS
########################################
sudo yum update
sudo yum upgrade

########################################
## Setup AWS CLI
########################################
aws configure
### AKIA2SL7UQV42G5ZXOQL
### LZQTQ8yDyBZJEkDgU352sF5cGVdMLP39TwI+b+Sr
### us-west-2
### text

########################################
## Setup addtl volumes
########################################
# List all volumes
lsblk

### Format volume for data
sudo mkfs -t ext4 /dev/nvme1n1
sudo mkdir -p ~/project1
sudo mount /dev/nvme1n1 ~/project1/
### Change permissions
sudo chmod 777 -R ~/project1/

# Check
lsblk

########################################
## Python 3
########################################
sudo yum install python3

########################################
## Conda
########################################
# Install conda
sudo mkdir ~/project1/applications
sudo chmod 777 -R ~/project1/applications
cd ~/project1/applications
# Change to correct URL if not using Linux 64-bit
sudo curl -O https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
sudo bash Anaconda3-2019.10-Linux-x86_64.sh
### Save to /home/ec2-user/project1/applications/anaconda3

eval "$(/home/ec2-user/project1/applications/anaconda3/bin/conda shell.bash hook)"
conda init

# Close and reopen shell for changes to take effect

########################################
## Conda packages
########################################
# Install packages in conda
### Configure channel priority
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority false
conda config --set allow_conda_downgrades true

### Install packages
### Copy and run each line one at a time or will get write error
conda install -c conda-forge fastqc samtools
conda install adapterremoval bedtools
conda install -c bioconda/label/cf201901 picard star subread
conda update --all

#Check installs
conda list

########################################
## git
########################################
#sudo yum install git git-all

########################################
## git packages
########################################
# Install ENSEMBL

#cd ~/project1/applications
#git clone https://github.com/Ensembl/ensembl-vep.git
#cd 

########################################
## Link S3 storage
########################################
### Install fuse
### One line at a time
sudo amazon-linux-extras install epel
sudo yum install s3fs-fuse

### Setup key
echo AKIA2SL7UQV42G5ZXOQL:LZQTQ8yDyBZJEkDgU352sF5cGVdMLP39TwI+b+Sr > ~/.passwd-s3fs
chmod 600  ~/.passwd-s3fs

################# END ##################