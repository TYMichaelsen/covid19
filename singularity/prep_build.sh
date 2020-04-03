#!/usr/bin/env bash
set -ex
mkdir -p build

if [ ! -f ./build/ont-guppy_3.4.3_linux64.tar.gz ] ; then
	wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.4.3_linux64.tar.gz -O build/ont-guppy_3.4.3_linux64.tar.gz
fi
if [ ! -d ./build/artic-ncov2019 ]; then
	git clone --recursive https://github.com/artic-network/artic-ncov2019.git build/artic-ncov2019
	cd build/artic-ncov2019
	git checkout 2127aef
	cd -
fi
if [ ! -f ./build/miniconda.sh ] ; then
	wget -q https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O build/miniconda.sh
fi
if [ ! -d ./build/auspice ]; then
	git clone https://github.com/nextstrain/auspice.git build/auspice
	cd build/auspice
	git checkout -q v2.11.3
	cd -
fi
