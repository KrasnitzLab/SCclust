# Please edit the filenames so you can read and save
# prostate <- "nyu003.benign.1"

# to compare
# diff breakPointPins/Gleason6.1/averageNoEndsNewest/GL6.1simP.txt out/Gleason6.1/GL6.1simP.txt

prostate <- "nyu003.benign.1"
savedir <- "./out_many/nyu003/"
rdir <- "./out/nyu003/"

truefile <- paste(savedir,prostate,"trueP.txt",sep="")
simfile <- paste(savedir,prostate,"simP.txt",sep="")
source("./simFisher.R")
pinmat <- as.matrix(read.table(paste(rdir,prostate,"smear1bpPinMat.txt",sep=""),header=T))
pins<-as.matrix(read.table(paste(rdir,prostate,"smear1bpPins.txt",sep=""),header=T))
eviltwins <- c("CJA1593","CJA2772","CJA3438","CJA0951")
# eviltwins <- c("CJA1593","CJA2772","CJA0951")
pinmat <- pinmat[,setdiff(dimnames(pinmat)[[2]],eviltwins),drop=F]
m<-vector(mode="list",length=length(unique(pins[,"sign"])))
for(i in 1:length(unique(pins[,"sign"])))m[[i]]<-pinmat[pins[,"sign"]==unique(pins[,"sign"])[i],,drop=F]
require(parallel)
usecores <- 20 # detectCores()-4
vtrue <- simFisher(m,nsim=1,nsweep=0,seedme=123,distrib="Rparallel",njobs=usecores,combo="Fisher", savedir=savedir)
write(vtrue,file=truefile)
msim <- simFisher(m,nsim=500,nsweep=200,seedme=123,distrib="Rparallel",njobs=usecores,combo="Fisher", savedir=savedir)
write(msim,file=simfile)
# quit(save="no")

