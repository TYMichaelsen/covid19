#!/bin/bash
set -e
cd ~
# install basic tools
sudo apt-get update
sudo apt-get install -y \
  tmux nano git curl

#installs dependencies required for compiling singularity from source
#from: https://sylabs.io/guides/3.5/admin-guide/installation.html
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    uuid-dev \
    uidmap \
    libssl-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup-bin

#install nvidia-drivers
sudo add-apt-repository -y ppa:graphics-drivers/ppa
sudo apt-get update
sudo apt-get install -y nvidia-driver-440

# Install CUDA 10.2 host driver for base image to use.
#wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
#sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600
#sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
#sudo add-apt-repository "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/ /"
#sudo apt-get update
#sudo apt-get -y install cuda

# Install R
export R_VERSION=3.6.3
echo "deb http://cran.r-project.org/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/r.list
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
apt-get update
apt-get install -y --no-install-recommends \
  r-base=${R_VERSION}* \
  r-base-core=${R_VERSION}* \
  r-base-dev=${R_VERSION}* \
  r-recommended=${R_VERSION}* \
  r-base-html=${R_VERSION}* \
  r-doc-html=${R_VERSION}* \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libcairo2-dev \
  libxt-dev

#install RStudio
wget -q https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.2.5033-amd64.deb
sudo dpkg -i rstudio-1.2.5033-amd64.deb
rm rstudio-1.2.5033-amd64.deb