FROM nvidia/cuda:10.2-base-ubuntu18.04
WORKDIR /opt

#environment
ARG DEBIAN_FRONTEND=noninteractive
ENV R_VERSION 3.6.0
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

#general requirements
RUN apt-get update && \
  apt-get install -y --no-install-recommends --no-install-suggests \
  nano less git wget make g++ ca-certificates curl libidn11 \
  parallel=20161222* software-properties-common build-essential

#install GUPPY (3.4.3) (with GPU support) !!remember to add the correct path to $PATH
RUN wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.4.3_linux64.tar.gz && \
  tar -zxf ont-guppy_3.4.3_linux64.tar.gz && \
  rm /opt/ont-guppy_3.4.3_linux64.tar.gz
ENV PATH="/opt/ont-guppy/bin:${PATH}"

#install miniconda3 in silent mode (4.7.12.1)
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O miniconda.sh && \
  bash miniconda.sh -b -p /opt/miniconda && \
  rm /opt/miniconda.sh && \
  ln -s /opt/miniconda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
  echo ". /opt/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
#RUN echo "conda activate base" >> ~/.bashrc
ENV PATH="/opt/miniconda/bin:${PATH}"

#install minimap2 (v2.17) (SKIPPED)
#RUN /opt/miniconda/bin/conda install -c bioconda minimap2=2.17 -y

#install artic-ncov2019-medaka conda environment
RUN git clone --recursive https://github.com/artic-network/artic-ncov2019.git && \
  conda env create -f artic-ncov2019/environment-medaka.yml
#install artic-ncov2019-medaka conda environment (SKIPPED)
#conda env create -f artic-ncov2019/environment.yml

#install nextstrain+augur+auspice conda environment
RUN mkdir -p /opt/nextstrain && \
  cd /opt/nextstrain && \
  curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml && \
  conda env create -f nextstrain.yml && \
  git clone https://github.com/nextstrain/auspice.git && \
  cd /opt/nextstrain/auspice && \
  git checkout -q v2.10.0
#install auspice with a separate command to catch error
#inspect manually to verify it's only the postinstall script that fails
RUN conda run -n nextstrain npm install --global /opt/nextstrain/auspice \
  || echo "*** Auspice install returned an error, please inspect ***"

#install R and R packages
#specific version of R to be installed
ENV R_BASE_VERSION 3.6.3
#add R CRAN mirror to apt sources to install R using APT
RUN apt-key adv --keyserver keyserver.ubuntu.com \
  --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
  echo "deb https://mirrors.dotsrc.org/cran/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/cran.list && \
  apt-get update && \
  apt-get -y install --no-install-recommends --no-install-suggests \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    r-base=${R_BASE_VERSION}* \
    r-base-core=${R_BASE_VERSION}* \
    r-base-dev=${R_BASE_VERSION}* \
    r-base-html=${R_BASE_VERSION}* \
    r-recommended=${R_BASE_VERSION}*
RUN R -e "install.packages(c('data.table', 'tidyverse', 'rmarkdown'), Ncpus = 100)"

#install RStudio Server (SKIPPED)
#RUN apt-get -y install --no-install-recommends --no-install-suggests \
#  gdebi-core
#RUN wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.2.5033-amd64.deb
#RUN gdebi rstudio-server-1.2.5033-amd64.deb

#install primer prospector
RUN conda create --name primerprospector -c bioconda/label/cf201901 primerprospector

#install bedtools (binary)
RUN mkdir -p /opt/bin
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary && \
  mv bedtools.static.binary /opt/bin/bedtools && \
  chmod +x /opt/bin/bedtools
ENV PATH="/opt/bin:${PATH}"

#install MinKnow
RUN wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add - && \
  echo "deb http://mirror.oxfordnanoportal.com/apt bionic-stable non-free" | \
  tee /etc/apt/sources.list.d/nanoporetech.sources.list && \
  apt-get update && \
  apt-get install -y --no-install-recommends --no-install-suggests \
  minion-nc

WORKDIR /home/user
