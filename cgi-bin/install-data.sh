#!/bin/sh
# Installation script of PipelineMiRNA data
# Takes as argument the path of Mirdup zip file
# and the path in which to install

if [ -z "$1" ]
  then
    echo "No argument supplied"
    echo "Usage: `basename $0` <ZIP_PATH> <DEST_PATH>"
    exit 1
fi
ZIP_PATH=$1

if [ -z "$2" ]
  then
    echo "No argument supplied"
    echo "Usage: `basename $0` <ZIP_PATH> <DEST_PATH>"
    exit 1
fi
DEST_PATH=$2

cd $DEST_PATH
MIRDUP=$DEST_PATH/mirdup
if [ ! -f $MIRDUP/plant.model ]
then
    echo "Deploying Mirdup data"
    mkdir $MIRDUP
    unzip -d $MIRDUP $ZIP_PATH/model_plants.zip
else
    echo "Mirdup data already there"
fi

