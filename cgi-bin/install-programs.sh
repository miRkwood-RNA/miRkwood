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

echo "Installing ViennaRNA 2.1.2"
wget http://www.tbi.univie.ac.at/~ronny/RNA/ViennaRNA-2.1.2.tar.gz
tar xf ViennaRNA-2.1.2.tar.gz
cd ViennaRNA-2.1.2
./configure
make
cd $ROOT_PATH

echo "Installing Randfold 2.0"

echo "Installing Squid"
wget http://selab.janelia.org/software/squid/squid.tar.gz
tar xf squid.tar.gz
cd squid-1.9g/
./configure
make install
cd $ROOT_PATH

wget http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003/downloads/randfold-2.0.tar.gz
tar xf randfold-2.0.tar.gz
cd randfold-2.0
make
cd $ROOT_PATH

echo "Installing selfcontain"
wget http://kim.bio.upenn.edu/software/rna-sc/selfcontain_unix.zip
unzip selfcontain_unix.zip
cd $ROOT_PATH

echo "Installing VARNA 3.9"
wget http://varna.lri.fr/bin/VARNAv3-9.jar
cd $ROOT_PATH

echo "Installing Exonerate"
#wget https://www.ebi.ac.uk/%7Eguy/exonerate/exonerate-2.2.0.tar.gz
#tar xf exonerate-2.2.0.tar.gz
#cd exonerate-2.2.0
#./configure
#make
wget https://www.ebi.ac.uk/%7Eguy/exonerate/exonerate-2.2.0-i386.tar.gz
tar xf exonerate-2.2.0-i386.tar.gz
cd $ROOT_PATH

echo "Installing MiRdup 1.2"
mkdir miRdup_1.2
cd miRdup_1.2
wget http://www.cs.mcgill.ca/~blanchem/mirdup/miRdup_1.2.zip
unzip miRdup_1.2.zip
cd $ROOT_PATH

echo "Installing NCBI Blast 2.2.28"
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.28+-src.tar.gz
tar xf ncbi-blast-2.2.28+-src.tar.gz
cd ncbi-blast-2.2.28+-src/c++/
./configure
make
cd $ROOT_PATH

rm *.gz *.zip
