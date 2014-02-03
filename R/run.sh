#*/bin/sh
ROOT_FOLDER=results
for CHR in ChrM ChrC Chr1 Chr2 Chr3 Chr4 Chr5
do
#	touch data/$CHR.bed
#	bedtools coverage -d -abam tutorial.bam -b $CHR.bed > $CHR.bedgraph
	FOLDER=$ROOT_FOLDER/$CHR
	mkdir -p $FOLDER
	Rscript ~/Tuile/pipelineMiRNA-web/R/testlocalscoremiRNA-revamped.R $FOLDER /home/jeanfred/Tuile/pipelineMiRNA-web/R/data/$CHR.bedgraph $CHR
done
