#!/bin/sh
perl --version
cpanm --version
cpanm --local-lib=~/perl local::lib && eval `perl -I ~/perl/lib/perl5/ -Mlocal::lib=~/perl`
cpanm --quiet --notest --skip-satisfied Dist::Zilla
dzil authordeps | cpanm --quiet --notest --skip-satisfied
dzil listdeps | cpanm --quiet
