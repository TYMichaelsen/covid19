FROM nvidia/cuda:10.2-base-ubuntu18.04
WORKDIR /opt

#environment
ARG DEBIAN_FRONTEND=noninteractive
ENV R_VERSION 3.6.0
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

#general requirements
RUN apt-get update
RUN apt-get install -y \
  --no-install-recommends --no-install-suggests \
  nano less git wget make g++ ca-certificates

#install GUPPY (3.4.3) (with GPU support) !!remember to add the correct path to $PATH
RUN wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.4.3_linux64.tar.gz
RUN tar -zxf ont-guppy_3.4.3_linux64.tar.gz
RUN rm /opt/ont-guppy_3.4.3_linux64.tar.gz
ENV PATH="/opt/ont-guppy/bin:${PATH}"

#install miniconda3 in silent mode (4.7.12.1)
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p /opt/miniconda
RUN rm /opt/miniconda.sh
RUN ln -s /opt/miniconda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
RUN echo ". /opt/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate base" >> ~/.bashrc
ENV PATH="/opt/miniconda/bin:${PATH}"

#install minimap2 (v2.17) (SKIPPED)
#RUN /opt/miniconda/bin/conda install -c bioconda minimap2=2.17 -y

#install artic-ncov2019 conda environment
RUN git clone --recursive https://github.com/artic-network/artic-ncov2019.git
RUN conda env create -f artic-ncov2019/environment.yml

WORKDIR /home/user