#!/usr/bin/env bash
# Usage: $0 <image_file_name.sif> <Singularity_defition_file>

if [ $# -lt 2 ]; then
 echo "Usage: $0 <image_file_name.sif> <Singularity_defition_file>"
 exit 1
fi

BUILDDIR=/tmp/$UID
mkdir -p $BUILDDIR

THISDIR=$PWD

DEFFILE="$2"
IMGFILE="$1"

cp "$DEFFILE" $BUILDDIR/

cd $BUILDDIR

echo -e "command: singularity build --fakeroot $IMGFILE ${DEFFILE} "
singularity build --fakeroot "$IMGFILE" ${DEFFILE} 

if [ $? -eq 0 ]; then
  mv $IMGFILE $THISDIR/
  rm ${BUILDDIR}/${DEFFILE}
else
 echo "Build failed. Check your defition file or build command"
fi
