#!/bin/sh
BASEDIR=$(dirname $0)
perl --version
cpanm --version
cpanm --local-lib=~/perl5 local::lib && sh $BASEDIR/use-locallib.sh
cpanm --quiet --notest --skip-satisfied Dist::Zilla
dzil authordeps | cpanm --quiet --notest --skip-satisfied
dzil listdeps | cpanm --quiet

