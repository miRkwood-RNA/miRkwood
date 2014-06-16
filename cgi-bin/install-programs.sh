#!/bin/sh
# Installation script of PipelineMiRNA programs
# Takes as argument the path in which to install

if [ -z "$1" ]
  then
    echo "No argument supplied"
	echo "Usage: `basename $0` <PATH>"
	exit 1
fi
ROOT_PATH=$(readlink -f $1)
cd $ROOT_PATH

VARNA=$ROOT_PATH/VARNA.jar
if [ ! -f $VARNA ]
then
    echo "Installing VARNA"
    wget http://varna.lri.fr/bin/VARNAv3-9.jar -O $VARNA
    cd $ROOT_PATH
else
    echo "VARNA already installed"
fi

MIRDUP=$ROOT_PATH/miRdup-1.2
if [ ! -e $MIRDUP ]
then
    echo "Installing MiRdup 1.2"
    mkdir $MIRDUP
    cd $MIRDUP
    wget http://www.cs.mcgill.ca/~blanchem/mirdup/miRdup_1.2.zip
    unzip miRdup_1.2.zip
    cd $ROOT_PATH
else
    echo "MiRdup 1.2 already installed"
fi

rm -f *.gz *.zip
