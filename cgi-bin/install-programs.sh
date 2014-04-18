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

#echo "Installing ViennaRNA 2.1.2"
#wget http://www.tbi.univie.ac.at/~ronny/RNA/ViennaRNA-2.1.2.tar.gz
#tar xf ViennaRNA-2.1.2.tar.gz
#cd ViennaRNA-2.1.2
#./configure
#make
#cd $ROOT_PATH

echo "Installing VARNA 3.9"
wget http://varna.lri.fr/bin/VARNAv3-9.jar
cd $ROOT_PATH

echo "Installing MiRdup 1.2"
mkdir miRdup_1.2
cd miRdup_1.2
wget http://www.cs.mcgill.ca/~blanchem/mirdup/miRdup_1.2.zip
unzip miRdup_1.2.zip
cd $ROOT_PATH

rm *.gz *.zip
