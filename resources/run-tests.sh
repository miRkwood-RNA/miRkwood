#!/bin/sh
# Script that runs the tests.
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
prove --r --lib --timer --formatter=TAP::Formatter::JUnit t/ | tee test-results.xml
