setwd("/home/jeanfred/Workbench/)

library(pvalues)

calcfinalres=function(x,posdep,minadjtags=2,B=500){
  score=localScore(x)
  if(nrow(score)==0){resScLoc=NULL}else{
    lengthsig=length(x)-c(0,cumsum(score[,"length"])[-nrow(score)])
    subseqMatrix=cbind("nbmeas"=score[,"length"], "lengthsig"=lengthsig, "localscore"=score[,"score"])
    keeppos=which(subseqMatrix[,"nbmeas"]>=minadjtags)
    resScLoc=data.frame(score[keeppos,],pvalueResampling(subseqMatrix[keeppos,],B, x, nbThreads=1, seed=123))
  }
  signifreg=resScLoc[which(resScLoc$pvalue<=0.05),]
  finalres=data.frame("startreg"=signifreg$start+posdep,"endreg"=signifreg$end+posdep,signifreg)
  finalres
}

signal=read.table("countschr2.bedgraphBaseByBase")
#On a deux régions dans notre BED
countsRegInt=signal[which(signal[,2]==19150000),5]
countsBruit=signal[which(signal[,2]==200000),5]

#tests seuils 3 a 8
llseuil=lapply(3:8,FUN=function(seuil){
  xint=countsRegInt-seuil
  xb=countsBruit-seuil
  finalresint=calcfinalres(xint,posdep=19150000)
  finalresb=calcfinalres(xb,posdep=200000)
  res=list(finalresint,finalresb)
  res  
})
# llseuil est une liste

allfinalresint=lapply(llseuil,FUN=function(x) x[[1]])
allfinalresb=lapply(llseuil,FUN=function(x) x[[2]])
#On récupère le 1 et le 2 élément de la liste


#filtre pour background

filterreg=lapply(allfinalresb,FUN=function(x) x[which(x$length<=150),])
filterreg

##########graphes

seuil=3
#les graphes sont centres sur le milieu du pic et font tous font la meme taille, possible grace au filtre
setwd(paste(c("seuil",seuil,"/filter/"),collapse=""))
finalresb=filterreg[[seuil]]
decalmidreg=200
posdep=200000
midreg=round((finalresb[,"start"]+finalresb[,"end"])/2)
for (i in 1:nrow(finalresb)){
  jpeg(paste(c("backgroundPeak",finalresb$startreg[i],".jpg"),collapse=""))
  xlimmin=max(1,(midreg[i]-decalmidreg))
  xlimmax=min(length(countsBruit),(midreg[i]+decalmidreg))
  plot((xlimmin+posdep):(xlimmax+posdep),countsBruit[xlimmin:xlimmax],ylim=c(0,1000),
       xlab="position",ylab="counts",main="Peak in background",type='h')
  abline(v=finalresb[i,"startreg"]-0.5,col="blue")
  abline(v=finalresb[i,"endreg"]+0.5,col="blue")
  dev.off()
}

posdep=19150000
finalresint=allfinalresint[[seuil]]
midreg=round((finalresint[,"start"]+finalresint[,"end"])/2)
for (i in 1:nrow(finalresint)){
  jpeg(paste(c("peakInterest",finalresint$startreg[i],".jpg"),collapse=""))
  xlimmin=max(1,(midreg[i]-decalmidreg))
  xlimmax=min(length(countsRegInt),(midreg[i]+decalmidreg))
  plot((xlimmin+posdep):(xlimmax+posdep),countsRegInt[xlimmin:xlimmax],ylim=c(0,1000),
       xlab="position",ylab="counts",main="Peak in region of interest",type='h')
  abline(v=finalresint[i,"startreg"]-0.5,col="blue")
  abline(v=finalresint[i,"endreg"]+0.5,col="blue")
  dev.off()
}

setwd("../..")

resTogive=rbind(data.frame("chr"="chr2",allfinalresint[[3]][,c("startreg","endreg")],"region"="interest"),
                data.frame("chr"="chr2",filterreg[[3]][,c("startreg","endreg")],"region"="background"))

write.table(resTogive[order(resTogive[,"startreg"]),],"peaksfound.txt",row.names=FALSE)



