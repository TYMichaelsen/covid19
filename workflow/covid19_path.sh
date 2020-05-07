#!/bin/sh
# Prepare PATH for running scripts of COVID19 workflow
# This should be sourced in ~/bash_profile ~/bashrc, or /etc/profile.d
# Written by Vang Le-Quy
# Date: 20200506

SCRIPTDIR="/srv/rbd/covid19/git/covid19/workflow"

if [ -d $SCRIPTDIR ] && [[ "$PATH" != *"$SCRIPTDIR"* ]]; then
    export PATH=$PATH:$SCRIPTDIR
fi
