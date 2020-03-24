#!/bin/bash
set -e
#installs dependencies required for compiling singularity from source
#from: https://sylabs.io/guides/3.5/admin-guide/installation.html
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    uuid-dev \
    uidmap \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup-bin 

