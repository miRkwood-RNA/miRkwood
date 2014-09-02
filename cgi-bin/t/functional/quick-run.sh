# Initialisation
SCRIPT=$(readlink -f $0)
BASEDIR=$(dirname $SCRIPT)
ROOTDIR=$BASEDIR/../..

TESTFOLDER=$BASEDIR/output/quicktest

rm -rf $TESTFOLDER && mkdir -p $TESTFOLDER
perl -I$ROOTDIR/lib $ROOTDIR/bin/mirkwood.pl --output $TESTFOLDER $BASEDIR/data/sequenceSomething.fas --no-process
