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

echo "Deploying Mirdup data"
MIRDUP=$ROOT_PATH/mirdup
mkdir $MIRDUP
cd $MIRDUP
wget http://www.cs.mcgill.ca/~blanchem/mirdup/model_plants.zip
unzip model_plants.zip
rm model_plants.zip
cd $ROOT_PATH

