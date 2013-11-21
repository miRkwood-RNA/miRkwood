# Initialisation
SCRIPT=$(readlink -f $0)
BASEDIR=$(dirname $SCRIPT)
ROOTDIR=$BASEDIR/../..

TESTFOLDER=$BASEDIR/output/quicktest

rm -rf $TESTFOLDER && mkdir -p $TESTFOLDER && cp $BASEDIR/data/sequenceSomething.fas $TESTFOLDER/sequenceUpload.fas
perl -I$ROOTDIR/lib $ROOTDIR/scripts/execute_scripts.pl 0 0 1 1 $TESTFOLDER ATpepTAIR10
