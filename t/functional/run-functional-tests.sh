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
plan 3

### Testing script of the PipelineMiRNA
rm -rf $BASEDIR/output/fullpipeline1/ && mkdir -p $BASEDIR/output/fullpipeline1/ && cp $BASEDIR/data/sequenceSomething.fas $BASEDIR/output/fullpipeline1/sequenceUpload.fas
perl -I$ROOTDIR/lib $ROOTDIR/scripts/execute_scripts.pl unChecked  mfeiChecked randfoldChecked UNSCChecked UNalignChecked $BASEDIR/output/fullpipeline1/ ATpepTAIR10
ok 'Full pipeline' [ `diff --exclude=.svn --exclude=pvalue.txt -qr $BASEDIR/output/fullpipeline1/ $BASEDIR/expected/fullpipeline1/ | wc -l` -eq 0 ]

rm -rf $BASEDIR/output/fullpipeline2/ && mkdir -p $BASEDIR/output/fullpipeline2/ && cp $BASEDIR/data/filtercds_in.fas $BASEDIR/output/fullpipeline2/sequenceUpload.fas
perl -I$ROOTDIR/lib $ROOTDIR/scripts/execute_scripts.pl checked  mfeiChecked randfoldChecked UNSCChecked UNalignChecked $BASEDIR/output/fullpipeline2/ ATpepTAIR10
ok 'Full pipeline with FilteringCDS' [ `diff --exclude=.svn --exclude=pvalue.txt -qr $BASEDIR/output/fullpipeline2/ $BASEDIR/expected/fullpipeline2/ | wc -l` -eq 0 ]

rm -rf $BASEDIR/output/filtercds/ && mkdir -p $BASEDIR/output/filtercds/ && cp $BASEDIR/data/filtercds_in.fas $BASEDIR/output/filtercds/sequenceUpload.fas
perl -I$ROOTDIR/lib $ROOTDIR/scripts/filterCDS.pl $ROOTDIR/data/ $BASEDIR/output/filtercds/ ATpepTAIR10
ok 'FilterCDS' [ `diff --exclude=.svn --exclude=outBlast.txt -qr $BASEDIR/expected/filtercds/ $BASEDIR/output/filtercds/ | wc -l` -eq 0 ]
