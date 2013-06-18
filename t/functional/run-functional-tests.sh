#!/bin/sh

# Importing TAP library
if [ -z "$LIBTAP_SH_HOME" ]
then
    if [ ! -f libtap.sh ]
    then
       wget --quiet "http://git.eyrie.org/?p=devel/c-tap-harness.git;a=blob_plain;f=tests/tap/libtap.sh;hb=HEAD" -O libtap.sh
    fi
    LIBTAP="./libtap.sh"
else
    LIBTAP="$LIBTAP_SH_HOME/libtap.sh"
fi
. "$LIBTAP"

# Initialisation
SCRIPT=$(readlink -f $0)
BASEDIR=$(dirname $SCRIPT)
ROOTDIR=$BASEDIR/../..

# Testing
plan 1

### Testing script of the PipelineMiRNA
rm -rf $BASEDIR/output/sequenceSomething/ && mkdir -p $BASEDIR/output/sequenceSomething/ && cp $BASEDIR/data/sequenceSomething.fas $BASEDIR/output/sequenceSomething/sequenceUpload.fas
perl -I$ROOTDIR/lib $ROOTDIR/scripts/execute_scripts.pl unChecked  mfeiChecked randfoldChecked UNSCChecked UNalignChecked $BASEDIR/output/sequenceSomething/
ok 'First test' [ `diff --exclude=pvalue.txt -qr $BASEDIR/output/sequenceSomething/ $BASEDIR/expected/sequenceSomething/ | wc -l` -lt 1 ]
