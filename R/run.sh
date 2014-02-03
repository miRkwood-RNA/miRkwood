#*/bin/sh
ROOT_FOLDER=results
cd data
for CHR in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC
do
#	touch data/$CHR.bed
#	bedtools coverage -d -abam tutorial.bam -b $CHR.bed > $CHR.bedgraph
	FOLDER=$ROOT_FOLDER/$CHR
	mkdir -f $FOLDER
	Rscript ~/Tuile/pipelineMiRNA-web/R/testlocalscoremiRNA-revamped.R $FOLDER data/$CHR.bedgraph
done
