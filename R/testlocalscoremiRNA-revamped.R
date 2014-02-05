#!/usr/bin/Rscript
#Usage:
#Rscript testlocalscoremiRNA-test.R <working_dir> <bedfile>

# Loading library
library(pvalues)

# Reading arguments
args<-commandArgs(TRUE)
working_directory = args[1]
bedfile = args[2]
name = args[3]
seuil = as.integer(args[4])

# Setting global variables
posdep = 0
nbpermutations = 80

#
calcfinalres=function(x, start_position, minadjtags=2, permutations_number=nbpermutations){
  score=localScore(x)
  if(nrow(score)==0){resScLoc=NULL}else{
    lengthsig=length(x)-c(0,cumsum(score[,"length"])[-nrow(score)])
    subseqMatrix=cbind("nbmeas"=score[,"length"], "lengthsig"=lengthsig, "localscore"=score[,"score"])
    keeppos=which(subseqMatrix[,"nbmeas"]>=minadjtags)
    resScLoc=data.frame(score[keeppos,],pvalueResampling(subseqMatrix[keeppos,], permutations_number, x, nbThreads=1, seed=123))
  }
  signifreg=resScLoc[which(resScLoc$pvalue<=0.05),]
  finalres=data.frame("startreg"=signifreg$start+start_position,"endreg"=signifreg$end+start_position,signifreg)
  finalres
}

processChr = function(counts, start_position){
    xint=counts-seuil
    finalres = calcfinalres(xint, start_position=start_position)
    finalres[,c("startreg","endreg")]
}


main = function(file_path, threshold, start_position, chr_name){
    signal = read.table(file_path)
    setwd(working_directory)
    start_file = paste(chr_name, '-begin.txt', sep="")
    cat("Start", file=start_file, sep="\n")

    counts_chr = signal[,5]
    res_chr = processChr(counts_chr, start_position)
    end_file = paste(chr_name, '-done.txt', sep="")
    cat("Done", file=end_file, sep="\n")

    resTogive = data.frame("chr"=chr_name, res_chr)
    output_name = paste(chr_name, '-peaks.txt', sep="")
    write.table(resTogive, output_name, row.names=FALSE)
    resTogive
}

result = main(bedfile, seuil, posdep, name)

