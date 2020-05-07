#!/usr/bin/env bash
if [ "$1" != "install" ] && [ "$1" != "build" ]
then
	echo "Usage: $0 [install|build]"
	exit 1
fi

BUILDDIR=/tmp/$UID
mkdir -p $BUILDDIR

BASEIMG="covid19_1.5.sif"
FINALIMG="covid19_latest.sif"

THISDIR=$PWD
REPONAME=${THISDIR##*/}
REPODIR=$(dirname ${THISDIR})

buildImage(){
    # Handle build under /srv/rbd
    if [[ "$THISDIR" == "/srv/rbd"* ]]; then
        echo Cleaning up $BUILDDIR ...
        rm -rf ${BUILDDIR}
        echo Copying files to $BUILDDIR ...
        cp -a ${REPODIR} ${BUILDDIR}/ 
        cd ${BUILDDIR}/${REPONAME}
        echo Working under $PWD ...
        make
        if [ $? -eq 0 ]; then
            mv $BASEIMG $FINALIMG  ${THISDIR}/
            rm -rf ${BUILDDIR}
        else
            echo "Build failed. Check your defition file or build command"
        fi

    else # build under $HOME
        make
    fi
}

installImage(){
    INSTDIR="/srv/rbd/thecontainer" 
    if [[ ! -f ${THISDIR}/${FINALIMG} ]]; then
        buildImage
    else
    	local TimeStamp=$(date +%Y%m%d)
        cp ${THISDIR}/${FINALIMG} ${INSTDIR}/covid19_${TimeStamp}.sif
        cd ${INSTDIR}
        ln -snf covid19_${TimeStamp}.sif covid19_latest.sif
        echo Image installed and linked:
        ls -lah ${INSTDIR}/covid19_latest.sif
    fi
}

if [ "$1" == 'build' ]; then buildImage; fi
if [ "$1" == 'install' ]; then installImage; fi
