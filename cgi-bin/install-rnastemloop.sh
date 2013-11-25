#!/bin/sh
# Installation script of RNAstemloop
# Takes as argument the path in which to install

if [ -z "$1" ]
  then
    echo "No argument supplied"
	echo "Usage: `basename $0` <PATH>"
	exit 1
fi
ROOT_PATH=$(readlink -f $1)
cd $ROOT_PATH

echo "Installing RNAstemloop"
RNASTEMLOOP=$ROOT_PATH/extractPutativeStemloop
svn checkout svn+ssh://scm.gforge.inria.fr/svnroot/sequoia/pipelineMiRNA/software/extractPutativeStemloop
cd $RNASTEMLOOP
make
mv bin/RNAstemloop $ROOT_PATH/
rm -rf $RNASTEMLOOP
