#!/bin/sh
BASEDIR=$(dirname $0)
sh $BASEDIR/use-locallib.sh
prove --r --lib --timer --formatter=TAP::Formatter::JUnit t/ | tee test-results.xml
sh t/functional/run-functional-tests.sh | tee functional-test-results.txt
