#!/bin/sh
# Installation script of PipelineMiRNA data
# Takes as argument the path in which to install

if [ -z "$1" ]
  then
    echo "No argument supplied"
	echo "Usage: `basename $0` <PATH>"
	exit 1
fi
ROOT_PATH=$(readlink -f $1)

cd $ROOT_PATH
MIRDUP=$ROOT_PATH/mirdup
if [ ! -f $MIRDUP/plant.model ]
then
    echo "Deploying Mirdup data"
    mkdir $MIRDUP
    cd $MIRDUP
    wget http://www.cs.mcgill.ca/~blanchem/mirdup/model_plants.zip
    unzip model_plants.zip
    rm model_plants.zip
    cd $ROOT_PATH
else
    echo "Mirdup data already there"
fi

