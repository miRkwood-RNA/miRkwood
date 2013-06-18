#!/bin/bash
BASEDIR=$(dirname $0)
perlcritic -profile resources/perlcritic.rc lib/ scripts/ cgi-bin/ | tee perlcritic.txt
