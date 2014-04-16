#!/bin/sh

# Initialisation
SCRIPT=$(readlink -f $0)
BASEDIR=$(dirname $SCRIPT)
ROOTDIR=$BASEDIR/../..

# Importing TAP library
if [ -z "$LIBTAP_SH_HOME" ]
then
    if [ ! -f $BASEDIR/libtap.sh ]
    then
       wget --quiet "http://git.eyrie.org/?p=devel/c-tap-harness.git;a=blob_plain;f=tests/tap/libtap.sh;hb=HEAD" -O libtap.sh
    fi
    LIBTAP="$BASEDIR/libtap.sh"
else
    LIBTAP="$LIBTAP_SH_HOME/libtap.sh"
fi
. "$LIBTAP"



# Testing
plan 3

### Testing script of the PipelineMiRNA
EXCLUDES="--exclude=.svn --exclude=pvalue.txt --exclude=outBlast.txt --exclude=*miRdupOutput.txt --exclude=*.log --exclude=*.cfg"

rm -rf $BASEDIR/output/filtercds/ && mkdir -p $BASEDIR/output/filtercds/ && cp $BASEDIR/data/filtercds_in.fas $BASEDIR/output/filtercds/sequenceUpload.fas
perl -I$ROOTDIR/lib $ROOTDIR/scripts/filterCDS.pl $ROOTDIR/data/ $BASEDIR/output/filtercds/ ATpepTAIR10
ok 'FilterCDS' [ `diff $EXCLUDES -qr $BASEDIR/expected/filtercds/ $BASEDIR/output/filtercds/ | wc -l` -eq 0 ]

rm -rf $BASEDIR/output/fullpipeline1/ && mkdir -p $BASEDIR/output/fullpipeline1/ && cp $BASEDIR/data/sequenceSomething.fas $BASEDIR/output/fullpipeline1/sequenceUpload.fas
perl -I$ROOTDIR/lib $ROOTDIR/scripts/execute_scripts.pl unChecked strandChecked mfeiChecked randfoldChecked UNalignChecked $BASEDIR/output/fullpipeline1/ ATpepTAIR10
DIFF=$(diff $EXCLUDES $IGNORES -qr $BASEDIR/output/fullpipeline1/ $BASEDIR/expected/fullpipeline1/ | wc -l)
ok 'Full pipeline' [ $DIFF -eq 0 ]

rm -rf $BASEDIR/output/fullpipeline2/ && mkdir -p $BASEDIR/output/fullpipeline2/ && cp $BASEDIR/data/filtercds_in.fas $BASEDIR/output/fullpipeline2/sequenceUpload.fas
perl -I$ROOTDIR/lib $ROOTDIR/scripts/execute_scripts.pl checked  strandChecked mfeiChecked randfoldChecked UNalignChecked $BASEDIR/output/fullpipeline2/ ATpepTAIR10
DIFF=$(diff $EXCLUDES -I 'fullpipeline' -qr $BASEDIR/output/fullpipeline2/ $BASEDIR/expected/fullpipeline2/ | wc -l)
ok 'Full pipeline with FilteringCDS' [ $DIFF -eq 0 ]

#rm -rf $BASEDIR/output/fullpipeline3/ && mkdir -p $BASEDIR/output/fullpipeline3/ && cp $BASEDIR/data/sequenceSomething.fas $BASEDIR/output/fullpipeline3/sequenceUpload.fas
#perl -I$ROOTDIR/lib $ROOTDIR/scripts/execute_scripts.pl UNchecked  mfeiChecked randfoldChecked alignChecked $BASEDIR/output/fullpipeline3/ ATpepTAIR10
#ok 'Full pipeline with alignment' [ `diff $EXCLUDES -qr $BASEDIR/output/fullpipeline3/ $BASEDIR/expected/fullpipeline3/ | wc -l` -eq 0 ]
