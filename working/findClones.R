#Please edit the filenames so you can read and save
hcdir<-"/Volumes/user/krasnitz/prostateSingleCell/breakPointPins/NYU002HGPIN/averageNoEnds/"
hcname<-"NYU002smear1bpLog10FisherHC"
datadir<-"/Volumes/user/krasnitz/prostateSingleCell/NYU002HGPIN/new/"
goodslicefile<-"NYU003slices20kGood.txt"
shortfile<-"NYU003short20k.txt"
prostate<-"NYU002"
savedir<-"/Volumes/user/krasnitz/prostateSingleCell/breakPointPins/NYU002HGPIN/averageNoEnds/"
pcut<-1e-3
rdir<-"/Volumes/user/krasnitz/prostateSingleCell/breakPointPins/FisherP/"
fdrthresh<-(-3)
truefile<-paste(rdir,prostate,"trueP.txt",sep="")
simfile<-paste(rdir,prostate,"simP.txt",sep="")
source("/Volumes/user/krasnitz/prostateSingleCell/Rtools/TreePy.R")
cellnames<-
	dimnames(read.table(paste(savedir,prostate,"smear1bpPinMat.txt",sep=""),header=T))[[2]]
eviltwins<-c("CJA1593","CJA2772","CJA3438","CJA0951")
jguide<-read.table("/Volumes/user/krasnitz/prostateSingleCell/annot/joan02.guide111114.txt",header=T,as.is=T,sep="\t",comment.char="",fill=T)
log10data<-F
nsim<-500
hcmethod<-"average"
vtrue<-scan(truefile)
msim<-scan(simfile)
smsim<-sort(msim)
rm(msim)
gc()
usmsim<-unique(smsim)
csim<-tapply(match(smsim,usmsim),match(smsim,usmsim),length)
rm(smsim)
gc()
if(log10data)usmsim<-10^usmsim
if(log10data)vtrue<-10^vtrue
pos<-(1:length(vtrue))[order(vtrue)]
svtrue<-sort(vtrue)
usvtrue<-unique(svtrue)
ctrue<-tapply(match(svtrue,usvtrue),match(svtrue,usvtrue),length)
zc<-csim/cumsum(csim)
zc[usmsim<pcut]<-0
mylm<-lm(log(cumsum(csim[(cumsum(csim)/sum(csim))<0.1])/sum(csim))~
	poly(log(usmsim[(cumsum(csim)/sum(csim))<0.1]),degree=2,raw=T))
#	log(usmsim[(cumsum(csim)/sum(csim))<0.1])+I(log(usmsim[(cumsum(csim)/sum(csim))<0.1])^2))
#mylm<-lm(log(cumsum(csim[(cumsum(csim)/sum(csim))<0.1])/sum(csim))~
#	log(usmsim[(cumsum(csim)/sum(csim))<0.1]))
#mylm<-lm(log(cumsum(csim[usmsim<1])/sum(csim))~log(usmsim[usmsim<1]))
#mylm<-lm(log(cumsum(csim[usmsim<usmsim[which.max(zc)]])/sum(csim))~
#	log(usmsim[usmsim<usmsim[which.max(zc)]]))
#logfdr<-mylm$coefficients[2]*log(usvtrue)+mylm$coefficients[1]-log(cumsum(ctrue)/sum(ctrue))
logfdr<-mylm$coefficients[3]*log(usvtrue)^2+mylm$coefficients[2]*log(usvtrue)+
	mylm$coefficients[1]-log(cumsum(ctrue)/sum(ctrue))
#logfdr<-mylm$coefficients[2]*log(usvtrue[usvtrue<1])+
#	max(log(cumsum(ctrue[usvtrue<1])/sum(ctrue)))-max(log(cumsum(csim[usmsim<1])/sum(csim)))+
#	mylm$coefficients[1]-log(cumsum(ctrue[usvtrue<1])/sum(ctrue))
logfdrlong<-logfdr[match(vtrue,usvtrue)]
mdist<-matrix(ncol=(1+sqrt(1+8*length(vtrue)))/2,nrow=(1+sqrt(1+8*length(vtrue)))/2,data=0)
mdist[upper.tri(mdist)]<-log10(vtrue)
mdist<-pmin(mdist,t(mdist))
dimnames(mdist)<-list(cellnames,cellnames)
mfdr<-matrix(ncol=(1+sqrt(1+8*length(vtrue)))/2,nrow=(1+sqrt(1+8*length(vtrue)))/2,data=0)
mfdr[upper.tri(mfdr)]<-logfdrlong/log(10)
mfdr<-pmin(mfdr,t(mfdr))
dimnames(mfdr)<-list(cellnames,cellnames)
newcellnames<-setdiff(cellnames,eviltwins)
mdist<-mdist[newcellnames,newcellnames]
mfdr<-mfdr[newcellnames,newcellnames]
hc<-hclust(as.dist(mdist),method=hcmethod)
leaflist<-vector(mode="list",length=nrow(hc$merge))
mergefdr<-rep(NA,nrow(hc$merge))
nodesize<-rep(NA,nrow(hc$merge))
for(i in 1:nrow(hc$merge)){
	if(hc$merge[i,1]<0)leaflist[[i]]<-(-hc$merge[i,1])
	else leaflist[[i]]<-leaflist[[hc$merge[i,1]]]
	if(hc$merge[i,2]<0)leaflist[[i]]<-c(leaflist[[i]],(-hc$merge[i,2]))
	else leaflist[[i]]<-c(leaflist[[i]],leaflist[[hc$merge[i,2]]])
	nodesize[i]<-length(leaflist[[i]])
	mergefdr[i]<-
		max(mfdr[leaflist[[i]],leaflist[[i]]][upper.tri(mfdr[leaflist[[i]],leaflist[[i]]])])
}
hc$mergefdr<-mergefdr
hc$nodesize<-nodesize
hc$leaflist<-leaflist
rm(mergefdr,nodesize,leaflist)
fdrthresh<-(-3)
clonemerge<-hc$merge[hc$mergefdr<fdrthresh,]
clonenodes<-(1:nrow(hc$merge))[hc$mergefdr<fdrthresh]
clonenodes<-clonenodes[!(clonenodes%in%c(clonemerge[,1],clonemerge[,2]))]
hc$fdrthresh<-fdrthresh
hc$clonenodes<-clonenodes
rm(clonemerge,clonenodes,fdrthresh)
#curve(exp(mylm$coefficients[1]+mylm$coefficients[2]*log(x)),from=min(usvtrue),to=max(usvtrue),
#log="xy")
curve(exp(mylm$coefficients[1]+mylm$coefficients[2]*log(x)+mylm$coefficients[3]*log(x)^2),
from=min(usvtrue),to=max(usvtrue),log="xy")
points(usmsim,cumsum(csim)/sum(csim))
points(usvtrue,cumsum(ctrue)/sum(ctrue),col="red")
#lines(usmsim,exp(mylm$coefficients[1]+mylm$coefficients[2]*log(usmsim)))
for(co in 2:3){
hcut<-hc
#co<-2
hcut$height<-hcut$height[hcut$nodesize>=co]
hcut$mergefdr<-hcut$mergefdr[hcut$nodesize>=co]
hcut$merge<-hcut$merge[hcut$nodesize>=co,]
hcut$leaflist<-hcut$leaflist[hcut$nodesize>=co]
hcut$nodesize<-hcut$nodesize[hcut$nodesize>=co]
plot(10^(hcut$mergefdr[order(hcut$mergefdr)]),cummax(hcut$nodesize[order(hcut$mergefdr)]),
	log="xy",ylab="max. cohesive node size",xlab="FDR",ylim=c(1,max(hcut$nodesize)))
plot(10^(hcut$mergefdr[order(hcut$mergefdr)]),sapply(1:nrow(hcut$merge),
	function(j,x)length(unique(unlist(x[1:j]))),x=hcut$leaflist[order(hcut$mergefdr)]),log="xy",
	ylab="# of cells in cohesive nodes",xlab="FDR",ylim=c(1,max(hcut$nodesize)))
coreid<-jguide[match(newcellnames,jguide[,"seq.unit.id"]),"sector"]
ucoreid<-unique(coreid)
for(uc in ucoreid){
	coreleaves<-lapply(hcut$leaflist,intersect,y=(1:length(newcellnames))[coreid==uc])
	points(10^(hcut$mergefdr[order(hcut$mergefdr)]),sapply(1:nrow(hcut$merge),
  	function(j,x)length(unique(unlist(x[1:j]))),x=coreleaves[order(hcut$mergefdr)]))
	text(x=10^(hcut$mergefdr[order(hcut$mergefdr)]),sapply(1:nrow(hcut$merge),
		function(j,x)length(unique(unlist(x[1:j]))),x=coreleaves[order(hcut$mergefdr)]),cex=0.7,
		labels=substring(uc,first=1,last=1),offset=0)
}
fdrseries<-c(1e-15,1e-12,1e-9,1e-6,1e-3)
coresNclones<-matrix(ncol=length(fdrseries)+1,nrow=length(ucoreid),dimnames=list(ucoreid,
	c(as.character(fdrseries),"Ncells")))
for(uc in ucoreid){
	coreleaves<-lapply(hcut$leaflist,intersect,y=(1:length(newcellnames))[coreid==uc])
	z<-cbind(c(10^(hcut$mergefdr[order(hcut$mergefdr)]),fdrseries),c(sapply(1:nrow(hcut$merge),
		function(j,x)length(unique(unlist(x[1:j]))),x=coreleaves[order(hcut$mergefdr)]),
		rep(-1,length(fdrseries))))
	z<-z[order(z[,1]),]
	coresNclones[uc,]<-c(pmax(cummax(z[,2]),0)[z[,2]==-1],length(unique(unlist(coreleaves))))
}
write.table(coresNclones,
	paste(savedir,prostate,"smear1bpCoresNclonesCo",as.character(co),".txt", sep=""),
	col.names=T,row.names=T,sep="\t",quote=F)
}
pytable<-TreePy(data=as.dist(mdist),method="average")
pytable<-cbind(pytable,hc$mergefdr)
dimnames(pytable)[[2]][ncol(pytable)]<-"log10fdr"
hcname<-paste(prostate,"smear1bpLog10FisherHC",sep="")
assign(hcname,hc)
save(list=hcname,file=paste(savedir,hcname,".rda",sep=""))
write.table(mdist,paste(savedir,prostate,"smear1bpLog10FisherP.txt",sep=""),col.names=T,
	row.names=T,sep="\t",quote=F)
write.table(mfdr,paste(savedir,prostate,"smear1bpLog10FisherFDR.txt",sep=""),
	col.names=T,row.names=T,sep="\t",quote=F)
write.table(pytable,paste(savedir,prostate,"smear1bpFisherTreePy.txt", sep=""),col.names=T,
	row.names=F,sep="\t",quote=F)
savehistory(paste(savedir,"Rhistory",Sys.Date(),".txt",sep=""))
quit(save="no")
load(paste(hcdir,hcname,".rda",sep=""))
hcorig<-get(hcname)
clonemerge<-hcorig$merge[hcorig$mergefdr<fdrthresh,]
clonenodes<-(1:nrow(hcorig$merge))[hcorig$mergefdr<fdrthresh]
clonenodes<-clonenodes[!(clonenodes%in%c(clonemerge[,1],clonemerge[,2]))]
hcorig$fdrthresh<-fdrthresh
hcorig$clonenodes<-clonenodes
rm(clonenodes)
tshortAll<-read.table(paste(datadir,shortfile,sep=""),header=T,as.is=T)
truefile<-paste(rdir,prostate,"trueP.txt",sep="")
simfile<-paste(rdir,prostate,"simP.txt",sep="")
source("/Volumes/user/krasnitz/prostateSingleCell/Rtools/TreePy.R")
cellnames<-
	dimnames(read.table(paste(savedir,prostate,"smear1bpPinMat.txt",sep=""),header=T))[[2]]
eviltwins<-c("CJA1593","CJA2772","CJA3438","CJA0951")
jguide<-read.table("/Volumes/user/krasnitz/prostateSingleCell/annot/joan02.guide111114.txt",header=T,as.is=T,sep="\t",comment.char="",fill=T)
log10data<-F
nsim<-500
hcmethod<-"average"
vtrue<-scan(truefile)
msim<-scan(simfile)
smsim<-sort(msim)
rm(msim)
gc()
usmsim<-unique(smsim)
csim<-tapply(match(smsim,usmsim),match(smsim,usmsim),length)
rm(smsim)
gc()
if(log10data)usmsim<-10^usmsim
if(log10data)vtrue<-10^vtrue
pos<-(1:length(vtrue))[order(vtrue)]
svtrue<-sort(vtrue)
usvtrue<-unique(svtrue)
ctrue<-tapply(match(svtrue,usvtrue),match(svtrue,usvtrue),length)
zc<-csim/cumsum(csim)
zc[usmsim<pcut]<-0
mylm<-lm(log(cumsum(csim[(cumsum(csim)/sum(csim))<0.1])/sum(csim))~
	poly(log(usmsim[(cumsum(csim)/sum(csim))<0.1]),degree=2,raw=T))
#	log(usmsim[(cumsum(csim)/sum(csim))<0.1])+I(log(usmsim[(cumsum(csim)/sum(csim))<0.1])^2))
#mylm<-lm(log(cumsum(csim[(cumsum(csim)/sum(csim))<0.1])/sum(csim))~
#	log(usmsim[(cumsum(csim)/sum(csim))<0.1]))
#mylm<-lm(log(cumsum(csim[usmsim<1])/sum(csim))~log(usmsim[usmsim<1]))
#mylm<-lm(log(cumsum(csim[usmsim<usmsim[which.max(zc)]])/sum(csim))~
#	log(usmsim[usmsim<usmsim[which.max(zc)]]))
#logfdr<-mylm$coefficients[2]*log(usvtrue)+mylm$coefficients[1]-log(cumsum(ctrue)/sum(ctrue))
logfdr<-mylm$coefficients[3]*log(usvtrue)^2+mylm$coefficients[2]*log(usvtrue)+
	mylm$coefficients[1]-log(cumsum(ctrue)/sum(ctrue))
#logfdr<-mylm$coefficients[2]*log(usvtrue[usvtrue<1])+
#	max(log(cumsum(ctrue[usvtrue<1])/sum(ctrue)))-max(log(cumsum(csim[usmsim<1])/sum(csim)))+
#	mylm$coefficients[1]-log(cumsum(ctrue[usvtrue<1])/sum(ctrue))
logfdrlong<-logfdr[match(vtrue,usvtrue)]
mdist<-matrix(ncol=(1+sqrt(1+8*length(vtrue)))/2,nrow=(1+sqrt(1+8*length(vtrue)))/2,data=0)
mdist[upper.tri(mdist)]<-log10(vtrue)
mdist<-pmin(mdist,t(mdist))
dimnames(mdist)<-list(cellnames,cellnames)
mfdr<-matrix(ncol=(1+sqrt(1+8*length(vtrue)))/2,nrow=(1+sqrt(1+8*length(vtrue)))/2,data=0)
mfdr[upper.tri(mfdr)]<-logfdrlong/log(10)
mfdr<-pmin(mfdr,t(mfdr))
dimnames(mfdr)<-list(cellnames,cellnames)
newcellnames<-setdiff(cellnames,eviltwins)
mdist<-mdist[newcellnames,newcellnames]
mfdr<-mfdr[newcellnames,newcellnames]
hc<-hclust(as.dist(mdist),method=hcmethod)
leaflist<-vector(mode="list",length=nrow(hc$merge))
mergefdr<-rep(NA,nrow(hc$merge))
nodesize<-rep(NA,nrow(hc$merge))
for(i in 1:nrow(hc$merge)){
	if(hc$merge[i,1]<0)leaflist[[i]]<-(-hc$merge[i,1])
	else leaflist[[i]]<-leaflist[[hc$merge[i,1]]]
	if(hc$merge[i,2]<0)leaflist[[i]]<-c(leaflist[[i]],(-hc$merge[i,2]))
	else leaflist[[i]]<-c(leaflist[[i]],leaflist[[hc$merge[i,2]]])
	nodesize[i]<-length(leaflist[[i]])
	mergefdr[i]<-
		max(mfdr[leaflist[[i]],leaflist[[i]]][upper.tri(mfdr[leaflist[[i]],leaflist[[i]]])])
}
hc$mergefdr<-mergefdr
hc$nodesize<-nodesize
hc$leaflist<-leaflist
rm(mergefdr,nodesize,leaflist)
fdrthresh<-(-6)
clonemerge<-hc$merge[hc$mergefdr<fdrthresh,]
clonenodes<-(1:nrow(hc$merge))[hc$mergefdr<fdrthresh]
clonenodes<-clonenodes[!(clonenodes%in%c(clonemerge[,1],clonemerge[,2]))]
hc$fdrthresh<-fdrthresh
hc$clonenodes<-clonenodes
rm(clonemerge,clonenodes,fdrthresh)
#curve(exp(mylm$coefficients[1]+mylm$coefficients[2]*log(x)),from=min(usvtrue),to=max(usvtrue),
#log="xy")
curve(exp(mylm$coefficients[1]+mylm$coefficients[2]*log(x)+mylm$coefficients[3]*log(x)^2),
from=min(usvtrue),to=max(usvtrue),log="xy")
points(usmsim,cumsum(csim)/sum(csim))
points(usvtrue,cumsum(ctrue)/sum(ctrue),col="red")
#lines(usmsim,exp(mylm$coefficients[1]+mylm$coefficients[2]*log(usmsim)))
for(co in 2:3){
hcut<-hc
#co<-2
hcut$height<-hcut$height[hcut$nodesize>=co]
hcut$mergefdr<-hcut$mergefdr[hcut$nodesize>=co]
hcut$merge<-hcut$merge[hcut$nodesize>=co,]
hcut$leaflist<-hcut$leaflist[hcut$nodesize>=co]
hcut$nodesize<-hcut$nodesize[hcut$nodesize>=co]
plot(10^(hcut$mergefdr[order(hcut$mergefdr)]),cummax(hcut$nodesize[order(hcut$mergefdr)]),
	log="xy",ylab="max. cohesive node size",xlab="FDR",ylim=c(1,max(hcut$nodesize)))
plot(10^(hcut$mergefdr[order(hcut$mergefdr)]),sapply(1:nrow(hcut$merge),
	function(j,x)length(unique(unlist(x[1:j]))),x=hcut$leaflist[order(hcut$mergefdr)]),log="xy",
	ylab="# of cells in cohesive nodes",xlab="FDR",ylim=c(1,max(hcut$nodesize)))
coreid<-jguide[match(newcellnames,jguide[,"seq.unit.id"]),"sector"]
ucoreid<-unique(coreid)
for(uc in ucoreid){
	coreleaves<-lapply(hcut$leaflist,intersect,y=(1:length(newcellnames))[coreid==uc])
	points(10^(hcut$mergefdr[order(hcut$mergefdr)]),sapply(1:nrow(hcut$merge),
  	function(j,x)length(unique(unlist(x[1:j]))),x=coreleaves[order(hcut$mergefdr)]))
	text(x=10^(hcut$mergefdr[order(hcut$mergefdr)]),sapply(1:nrow(hcut$merge),
		function(j,x)length(unique(unlist(x[1:j]))),x=coreleaves[order(hcut$mergefdr)]),cex=0.7,
		labels=substring(uc,first=1,last=1),offset=0)
}
fdrseries<-c(1e-15,1e-12,1e-9,1e-6,1e-3)
coresNclones<-matrix(ncol=length(fdrseries)+1,nrow=length(ucoreid),dimnames=list(ucoreid,
	c(as.character(fdrseries),"Ncells")))
for(uc in ucoreid){
	coreleaves<-lapply(hcut$leaflist,intersect,y=(1:length(newcellnames))[coreid==uc])
	z<-cbind(c(10^(hcut$mergefdr[order(hcut$mergefdr)]),fdrseries),c(sapply(1:nrow(hcut$merge),
		function(j,x)length(unique(unlist(x[1:j]))),x=coreleaves[order(hcut$mergefdr)]),
		rep(-1,length(fdrseries))))
	z<-z[order(z[,1]),]
	coresNclones[uc,]<-c(pmax(cummax(z[,2]),0)[z[,2]==-1],length(unique(unlist(coreleaves))))
}
write.table(coresNclones,
	paste(savedir,prostate,"smear1bpCoresNclonesCo",as.character(co),".txt", sep=""),
	col.names=T,row.names=T,sep="\t",quote=F)
}
pytable<-TreePy(data=as.dist(mdist),method="average")
pytable<-cbind(pytable,hc$mergefdr)
dimnames(pytable)[[2]][ncol(pytable)]<-"log10fdr"
hcname<-paste(prostate,"smear1bpLog10FisherHC",sep="")
assign(hcname,hc)
save(list=hcname,file=paste(savedir,hcname,".rda",sep=""))
write.table(mdist,paste(savedir,prostate,"smear1bpLog10FisherP.txt",sep=""),col.names=T,
	row.names=T,sep="\t",quote=F)
write.table(mfdr,paste(savedir,prostate,"smear1bpLog10FisherFDR.txt",sep=""),
	col.names=T,row.names=T,sep="\t",quote=F)
write.table(pytable,paste(savedir,prostate,"smear1bpFisherTreePy.txt", sep=""),col.names=T,
	row.names=F,sep="\t",quote=F)
savehistory(paste(savedir,"Rhistory",Sys.Date(),".txt",sep=""))
quit(save="no")
