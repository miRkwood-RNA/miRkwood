#!/bin/sh
# Script that runs the test coverage.
# Takes as argument the path in which to run
RESOURCE_PATH=$(dirname $(readlink -f $0))
if [ -z "$1" ]
  then
    echo "No argument supplied"
	echo "Usage: `basename $0` <PATH>"
	exit 1
fi
EXECUTION_PATH=$(readlink -f $1)
cd $EXECUTION_PATH
dzil cover -select_re=lib/miRkwood/* -outputdir $EXECUTION_PATH/cover_db
#cover -delete && HARNESS_PERL_SWITCHES=-MDevel::Cover=+ignore,\.t,\.pl prove --lib t && cover
