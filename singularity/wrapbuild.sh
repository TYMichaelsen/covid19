#!/usr/bin/env bash

BASEIMG="covid19_1.5.sif"
FINALIMG="covid19_latest.sif"

THISDIR=$PWD
REPONAME=${THISDIR##*/}
REPODIR=$(dirname ${THISDIR})

buildImage(){
        make
}

installImage(){
    INSTDIR="/srv/rbd/thecontainer" 
    if [[ ! -f ${THISDIR}/${FINALIMG} ]]; then
        buildImage
    fi
    local TimeStamp=$(date +%y%m%d_%H%M)
    cp ${THISDIR}/${FINALIMG} ${INSTDIR}/covid19_${TimeStamp}.sif
    cd ${INSTDIR}
    ln -snf covid19_${TimeStamp}.sif covid19_latest.sif
    echo Image installed and linked:
    ls -lah ${INSTDIR}/covid19_latest.sif
}

case $1 in
    'build')
        buildImage
        ;;
    'install')
        installImage
        ;;
    'clean')
        rm $BASEIMG
        rm $FINALIMG
        ;;
    *)
        echo "Unknown command $1"
	      echo "Usage: $0 [install|build|clean]"
	      exit 1
        ;;
esac

