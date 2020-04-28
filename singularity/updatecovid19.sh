#!/bin/bash
set -e
dir=/opt/covid19
gitdir=${dir}/covid19
mkdir -p $gitdir
rm -rf $gitdir

git clone -b workstations https://github.com/kasperskytte/covid19.git $gitdir
if [ ! -f ${dir}/covid19_1.5.sif ]
then
  singularity pull -F covid19_1.5.sif library://kasperskytte/default/covid19:1.5
fi
cd ${dir}
singularity build -F covid19_latest.sif ${gitdir}/singularity/covid19.singularity
