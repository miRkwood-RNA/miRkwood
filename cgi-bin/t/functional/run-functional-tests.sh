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
plan 4

if [ ! -d "$BASEDIR/output/" ]; then
    mkdir "$BASEDIR/output/"
fi

### Testing script of miRkwood
EXCLUDES="--exclude=.svn --exclude=pvalue.txt --exclude=outBlast.txt --exclude=*miRdupOutput.txt --exclude=*.log --exclude=*.cfg --exclude=*.png --exclude=*.yml --exclude=*.html"

rm -rf $BASEDIR/output/fullpipeline1/
perl -I$ROOTDIR/lib $ROOTDIR/bin/mirkwood.pl --output $BASEDIR/output/fullpipeline1/ $BASEDIR/data/sequenceSomething.fas --align
perl compare_results.pl $BASEDIR/output/fullpipeline1/ $BASEDIR/expected/fullpipeline1/ > $BASEDIR/output/diff_fullpipeline1
DIFF=$(perl compare_results.pl $BASEDIR/output/fullpipeline1/ $BASEDIR/expected/fullpipeline1/ | wc -l)
ok 'Full pipeline' [ $DIFF -eq 0 ]


rm -rf $BASEDIR/output/fullpipeline2/
perl -I$ROOTDIR/lib $ROOTDIR/bin/mirkwood.pl --output $BASEDIR/output/fullpipeline2/ $BASEDIR/data/filtercds_in.fas --align --species-mask Arabidopsis_thaliana
perl compare_results.pl $BASEDIR/output/fullpipeline2/ $BASEDIR/expected/fullpipeline2/ > $BASEDIR/output/diff_fullpipeline2
DIFF=$(perl compare_results.pl $BASEDIR/output/fullpipeline2/ $BASEDIR/expected/fullpipeline2/ | wc -l)
ok 'Full pipeline with coding region masking (using BLAST)' [ $DIFF -eq 0 ]


rm -rf $BASEDIR/output/fullpipeline-bam/
perl -I$ROOTDIR/lib $ROOTDIR/bin/mirkwood-bam.pl --output $BASEDIR/output/fullpipeline-bam/ $BASEDIR/../data/Clusters.reads-Athaliana_167-ChrC.bam --genome $BASEDIR/../data/Clusters.Athaliana_167-ChrC.fa
perl compare_results.pl $BASEDIR/output/fullpipeline-bam/ $BASEDIR/expected/fullpipeline-bam/ > $BASEDIR/output/diff_fullpipeline-bam
DIFF=$(perl compare_results.pl $BASEDIR/output/fullpipeline-bam/ $BASEDIR/expected/fullpipeline-bam/ | wc -l)
ok 'Full BAM pipeline' [ $DIFF -eq 0 ]


rm -rf $BASEDIR/output/fullpipeline-bed/
perl -I$ROOTDIR/lib $ROOTDIR/bin/mirkwood-bed.pl --align --mfei --filter_otherRNA --filter_CDS --filter_multimapped --output $BASEDIR/output/fullpipeline-bed/ --genome $BASEDIR/../../data/genomes/Arabidopsis_thaliana.fasta $BASEDIR/../../data/SRR051927_test.bed
perl compare_results.pl $BASEDIR/output/fullpipeline-bed/ $BASEDIR/expected/fullpipeline-bed/ > $BASEDIR/output/diff_fullpipeline-bed
DIFF=$(perl compare_results.pl $BASEDIR/output/fullpipeline-bed/ $BASEDIR/expected/fullpipeline-bed/ | wc -l)
ok 'Full BED pipeline' [ $DIFF -eq 0 ]
