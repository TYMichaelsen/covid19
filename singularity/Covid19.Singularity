BootStrap: shub
From: biocyberman/singularity-r:3.6.3

%labels
  Maintainer Jeremy Nicklas
  RStudio_Version 1.2.5033

%help
  This will run r and applications for covid19 analysis without rstudio. This image has GPU/CUDA support

%apprun covidan 
  exec covidan "${@}"

%runscript
  exec covidan "${@}"

%environment
  export PATH=/opt/covid19/bin:${PATH}

%setup
   install -Dv \
    install_cuda.sh \
    ${SINGULARITY_ROOTFS}/usr/local/bin/install_cuda


%post
  # Software versions
  export R_VERSION=3.6.3
  export LANG=C.UTF-8 LC_ALL=C.UTF-8
  export DEBIAN_FRONTEND=noninteractive

 #general requirements
 cd /opt
 apt-get update && \
   apt-get install -y --no-install-recommends --no-install-suggests \
   nano less git wget make g++ ca-certificates curl libidn11 \
   parallel=20161222* software-properties-common build-essential \
   gawk
 
 #install GUPPY (3.4.3) (with GPU support) !!remember to add the correct path to $PATH
 wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.4.3_linux64.tar.gz && \
   tar -zxf ont-guppy_3.4.3_linux64.tar.gz -C /opt && \
   rm ont-guppy_3.4.3_linux64.tar.gz
 echo 'export PATH="/opt/ont-guppy/bin:${PATH}"' >> $SINGULARITY_ENVIRONMENT
 
 #install miniconda3 in silent mode (4.7.12.1)
 wget -q https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O miniconda.sh && \
   bash miniconda.sh -b -p /opt/miniconda && \
   rm miniconda.sh && \
   ln -s /opt/miniconda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
   echo ". /opt/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
 #echo "conda activate base" >> ~/.bashrc
 echo 'export PATH="/opt/miniconda/bin:${PATH}"' >> $SINGULARITY_ENVIRONMENT
 export PATH="/opt/miniconda/bin:${PATH}"
 
 #install minimap2 (v2.17) (SKIPPED)
 #/opt/miniconda/bin/conda install -c bioconda minimap2=2.17 -y

 #install cutadapt
 conda install -c bioconda cutadapt=1.18 -y
  
 #install artic-ncov2019-medaka conda environment
 git clone --recursive https://github.com/artic-network/artic-ncov2019.git && \
   conda env create -f artic-ncov2019/environment-medaka.yml
 #install artic-ncov2019-medaka conda environment (SKIPPED)
 #conda env create -f artic-ncov2019/environment.yml
 
 #install nextstrain+augur+auspice conda environment
 mkdir -p /opt/nextstrain && \
   cd /opt/nextstrain && \
   curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml && \
   conda env create -f nextstrain.yml && \
   git clone https://github.com/nextstrain/auspice.git && \
   cd /opt/nextstrain/auspice && \
   git checkout -q v2.11.3 && \
   conda run -n nextstrain npm update && \
   conda run -n nextstrain npm install && \
   conda run -n nextstrain npm run build && \
   conda run -n nextstrain npm link

 #install primer prospector
 conda create --name primerprospector -c bioconda/label/cf201901 primerprospector
 
 #install bedtools (binary)
 mkdir -p /opt/bin
 wget -q https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary && \
   mv bedtools.static.binary /opt/bin/bedtools && \
   chmod +x /opt/bin/bedtools
 echo 'export PATH="/opt/bin:${PATH}"' >> $SINGULARITY_ENVIRONMENT
 
 #install MinKnow
 wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add - && \
   echo "deb http://mirror.oxfordnanoportal.com/apt bionic-stable non-free" | \
   tee /etc/apt/sources.list.d/nanoporetech.sources.list && \
   apt-get update && \
   apt-get install -y --no-install-recommends --no-install-suggests \
   minion-nc
 
 # Clean up
 apt-get clean
 apt-get autoclean
 apt-get autoremove
 rm -rf /var/lib/apt/lists/*
 
 unset DEBIAN_FRONTEND
 
 
