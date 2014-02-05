#*/bin/sh
ROOT_FOLDER=results
WORK_FOLDER=~/Tuile/pipelineMiRNA-web/R
for THRESHOLD in 5 10 20 100 200 500
do
	for CHR in ChrM ChrC
	do
#		touch data/$CHR.bed
#		bedtools coverage -d -abam tutorial.bam -b $CHR.bed > $CHR.bedgraph
		FOLDER=$ROOT_FOLDER/Threshold$THRESHOLD
		mkdir -p $FOLDER
		Rscript $WORK_FOLDER/testlocalscoremiRNA-revamped.R $FOLDER $WORK_FOLDER/data/$CHR.bedgraph $CHR $THRESHOLD
	done
done
