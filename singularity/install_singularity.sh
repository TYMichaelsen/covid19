#!/bin/bash
set -e
if [ ! $1 ]
then
	echo "Usage:"
	echo "bash $(basename "$0") path/to/softwaredir "
	echo "(note \"/singularity\" will be appended at the end)"
	exit 1
fi
export usersoftwaredir=$1
mkdir -p $usersoftwaredir
cd $usersoftwaredir #or somewhere else
echo "downloading go..."
wget https://dl.google.com/go/go1.14.1.linux-amd64.tar.gz
echo "unpacking..."
tar -C $usersoftwaredir -zxvf go1.14.1.linux-amd64.tar.gz
rm go1.14.1.linux-amd64.tar.gz

#go is only needed to build singularity from source
export GOPATH=${HOME}/go
export PATH=${usersoftwaredir}/go/bin:${PATH}:${GOPATH}/bin

echo "downloading singularity..."
git clone https://github.com/sylabs/singularity.git
cd singularity
git checkout v3.5.3
echo "installing singularity into ${usersoftwaredir}/singularity..."
./mconfig --without-suid --prefix=${usersoftwaredir}/singularity
make -j -C ./builddir
make -j -C ./builddir install #remove sudo if without setuid
#rm -rf $GOPATH
echo "export PATH=${usersoftwaredir}/singularity/bin:$PATH" >> ${HOME}/.bashrc
echo ". ${usersoftwaredir}/singularity/etc/bash_completion.d/singularity" >> ${HOME}/.bashrc
source ~/.bashrc
unset GOPATH

