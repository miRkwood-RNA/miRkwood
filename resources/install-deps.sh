#!/bin/sh
BASEDIR=$(dirname $0)
perl --version
cpanm --version
cpanm --local-lib=~/perl5 local::lib && sh $BASEDIR/use-locallib.sh
cpanm --quiet --notest --skip-satisfied Dist::Zilla Test::Compile
dzil authordeps | cpanm --quiet --notest --skip-satisfied --force
cpanm Test::Compile Dist::Zilla::Plugin::Test::Compile
dzil listdeps | cpanm

