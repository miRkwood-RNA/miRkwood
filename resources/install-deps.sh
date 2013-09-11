#!/bin/sh
perl --version
cpanm --version
cpanm --local-lib=~/perl5 local::lib && sh use-locallib.sh
cpanm --quiet --notest --skip-satisfied Dist::Zilla
dzil authordeps | cpanm --quiet --notest --skip-satisfied
dzil listdeps | cpanm --quiet

