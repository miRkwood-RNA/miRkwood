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
       wget --quiet "http://git.eyrie.org/?p=devel/c-tap-harness.git;a=blob_plain;f=tests/tap/libtap.sh;hb=HEAD" -O $BASEDIR/libtap.sh
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

rm -rf $BASEDIR/output/fullpipeline1/
perl -I$ROOTDIR/lib $ROOTDIR/bin/mirkwood.pl --output $BASEDIR/output/fullpipeline1/ $BASEDIR/data/sequenceSomething.fas --align --no-process
DIFF=$(diff $EXCLUDES -I 'fullpipeline' -qr $BASEDIR/output/fullpipeline1/ $BASEDIR/expected/fullpipeline1/ | wc -l)
ok 'Full pipeline' [ $DIFF -eq 0 ]

rm -rf $BASEDIR/output/fullpipeline2/
perl -I$ROOTDIR/lib $ROOTDIR/bin/mirkwood.pl --output $BASEDIR/output/fullpipeline2/ $BASEDIR/data/filtercds_in.fas --align --no-process --species-mask Arabidopsis_thaliana
DIFF=$(diff $EXCLUDES -I 'fullpipeline' -qr $BASEDIR/output/fullpipeline2/ $BASEDIR/expected/fullpipeline2/ | wc -l)
ok 'Full pipeline with FilteringCDS' [ $DIFF -eq 0 ]

rm -rf $BASEDIR/output/fullpipeline-bam/
perl -I$ROOTDIR/lib $ROOTDIR/bin/mirkwood-bam.pl --output $BASEDIR/output/fullpipeline-bam/ $BASEDIR/../data/Clusters.reads-Athaliana_167-ChrC.bam --genome $BASEDIR/../data/Clusters.Athaliana_167-ChrC.fa  --no-process
DIFF=$(diff $EXCLUDES -I 'fullpipeline' -qr $BASEDIR/output/fullpipeline-bam/ $BASEDIR/expected/fullpipeline-bam/ | wc -l)
ok 'Full BAM pipeline' [ $DIFF -eq 0 ]
