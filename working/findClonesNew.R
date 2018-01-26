#Please edit the filenames so you can read and save
hcdir<-"/Volumes/user/krasnitz/prostateSingleCell/breakPointPins/NYU011Gl7.5/averageNoEndsNewer/"
#hcname<-"GL7.2smear1bpLog10FisherHC"
#datadir<-"/Volumes/user/krasnitz/prostateSingleCell/Gleason9.2/newer/"
#goodslicefile<-"GL6.1slices20kGood.txt"
#shortfile<-"GL6.1short20k.txt"
prostate<-"nyu011.GL7.5"
#savedir<-"/Volumes/user/krasnitz/prostateSingleCell/breakPointPins/NYU009B9/averageNoEndsNew/"
#rdir<-"/Volumes/user/krasnitz/prostateSingleCell/breakPointPins/FisherP/"
savedir<-hcdir
rdir<-savedir
fdrthresh<-(-2)	#FDR criterion for clone nodes
sharemin<-0.90	#A feature is considered shared if present in sharemin fraction of leaves in a node
nshare<-3	#Minimal number of shared features in a clone node
bymax<-T	#Use maximal of mean FDR for the node to find clones?
lmmax<-0.001 #A parameter in a linear fit to empirical null distribution of Fisher p-values
graphic<-F	#To plot or not
truefile<-paste(rdir,prostate,"trueP.txt",sep="")
simfile<-paste(rdir,prostate,"simP.txt",sep="")
source("/Volumes/user/krasnitz/prostateSingleCell/Rtools/TreePy.R")
cellnames<-
	dimnames(read.table(paste(hcdir,prostate,"smear1bpPinMat.txt",sep=""),header=T))[[2]]
#There is another twin filter upstream in the pipeline, but this does not hurt as long as we keep
#track of the actual cell IDs in the tree.
eviltwins<-c("CJA1593","CJA2772","CJA3438","CJA0951")	#1 of each known twin pair
jguide<-read.table("/Volumes/user/krasnitz/prostateSingleCell/annot/joan02.guide111114.txt",
	header=T,as.is=T,sep="\t",comment.char="",fill=T)
pinmat<-read.table(paste(hcdir,prostate,"smear1bpPinMat.txt",sep=""),header=T,as.is=T) #Incidence
log10data<-F	#I.e.,expect to read p-values, not log p-values
nsim<-500	#The size of a sample from the null
hcmethod<-"average"
vtrue<-scan(truefile)	#True pairwise Fisher p-values
msim<-scan(simfile)	#nsim sets of Fisher p-values sampled from the null
smsim<-sort(msim)	#Sort the true and the null data and get counts for each unique value
rm(msim)
gc()
usmsim<-unique(smsim)
csim<-tapply(match(smsim,usmsim),match(smsim,usmsim),length)
rm(smsim)
gc()
if(log10data)usmsim<-10^usmsim
if(log10data)vtrue<-10^vtrue
svtrue<-sort(vtrue)
usvtrue<-unique(svtrue)
ctrue<-tapply(match(svtrue,usvtrue),match(svtrue,usvtrue),length)
#Observed p-values are often far lower than any p-value sampled from the null; to determine FDR in
#such cases get a power-law fit to the low-p tail of the null CDF and use it to extrapolate to very 
#low p-values. Use the actual null CDF to estimate FDR for higher p-values.
mylm<-lm(log(cumsum(csim[(cumsum(csim)/sum(csim))<lmmax])/sum(csim))~
	log(usmsim[(cumsum(csim)/sum(csim))<lmmax]))
if(exists("mylm"))logfdrmod<-
	mylm$coefficients[2]*log(usvtrue)+mylm$coefficients[1]-log(cumsum(ctrue)/sum(ctrue))
if(graphic){
	curve(exp(mylm$coefficients[1]+mylm$coefficients[2]*log(x)),from=min(usvtrue),to=max(usvtrue),
		log="xy")
	points(usmsim,cumsum(csim)/sum(csim))
	points(usvtrue,cumsum(ctrue)/sum(ctrue),col="red")
}
#FDR is computed by comparing true to simulated CDF. However, empirical true CDF is only defined 
#for usvtrue, while empirical simulated CDF only for usmsim. Where possible, find empirical CDF for
#usvtrue by linear interpolation from the flanking values of usmsim.
z<-cbind(c(log(usvtrue),log(usmsim)),c(log(cumsum(ctrue)/sum(ctrue)),log(cumsum(csim)/sum(csim))),
	c(rep(0,length(usvtrue)),rep(1,length(usmsim))))
z<-z[order(z[,1]),]
simlow<-cumsum(z[,3])[z[,3]==0]
x1pos<-match(simlow,cumsum(z[,3]))[simlow>0]
x2pos<-match(simlow+1,cumsum(z[,3]))[simlow>0]
x1<-z[x1pos,1]
x2<-z[x2pos,1]
y1<-z[x1pos,2]
y2<-z[x2pos,2]
logfdrinterp<-rep(0,length(usvtrue))
logfdrinterp[simlow>0]<-(y2-y1)*log(usvtrue)[simlow>0]/(x2-x1)+(y1*x2-y2*x1)/(x2-x1)-
	log(cumsum(ctrue)/sum(ctrue))[simlow>0]
logfdr<-logfdrinterp
lmu<-max(usmsim[(cumsum(csim)/sum(csim))<lmmax])
if(is.finite(lmu)&min(usvtrue)<min(usmsim)){
	logfdr[usvtrue<lmu&usvtrue>min(usmsim)]<-
		(logfdrmod[usvtrue<lmu&usvtrue>min(usmsim)]*
		(log(lmu)-log(usvtrue[usvtrue<lmu&usvtrue>min(usmsim)]))-
		logfdrinterp[usvtrue<lmu&usvtrue>min(usmsim)]*
		(log(min(usmsim))-log(usvtrue[usvtrue<lmu&usvtrue>min(usmsim)])))/(log(lmu)-log(min(usmsim)))
	logfdr[usvtrue<min(usmsim)]<-logfdrmod[usvtrue<min(usmsim)]
}
logfdr<-cummax(logfdr)
logfdr[logfdr>0]<-0
if(graphic)plot(usvtrue,exp(logfdr),log="xy")
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
#Grow a tree and add multiple items to the standard hclust object
hc<-hclust(as.dist(mdist),method=hcmethod)
#Leaf indices for each node, in the order of the original labels
leaflist<-vector(mode="list",length=nrow(hc$merge))
#Maximal pairwise FDR anywhere in the node
mergefdr<-rep(NA,nrow(hc$merge))
#Mean FDR for the node
meanfdr<-rep(NA,nrow(hc$merge))
#Number of leaves in the node
nodesize<-rep(NA,nrow(hc$merge))
#Leaf lables for the node
labellist<-vector(mode="list",length=nrow(hc$merge))
#For each node and each feature determine the fraction of leaves in the node with the feature
sharing<-matrix(ncol=nrow(hc$merge),nrow=nrow(pinmat))
complexity<-rep(NA,nrow(hc$merge)) #Mean number of features per leaf in a node
for(i in 1:nrow(hc$merge)){
	if(hc$merge[i,1]<0)leaflist[[i]]<-(-hc$merge[i,1])
	else leaflist[[i]]<-leaflist[[hc$merge[i,1]]]
	if(hc$merge[i,2]<0)leaflist[[i]]<-c(leaflist[[i]],(-hc$merge[i,2]))
	else leaflist[[i]]<-c(leaflist[[i]],leaflist[[hc$merge[i,2]]])
	labellist[[i]]<-hc$labels[leaflist[[i]]]
	complexity[i]<-mean(colSums(pinmat[,labellist[[i]]]))
	sharing[,i]<-rowMeans(pinmat[,labellist[[i]]])
	nodesize[i]<-length(leaflist[[i]])
	mergefdr[i]<-
		max(mfdr[leaflist[[i]],leaflist[[i]]][upper.tri(mfdr[leaflist[[i]],leaflist[[i]]])])
	meanfdr[i]<-
		mean(mfdr[leaflist[[i]],leaflist[[i]]][upper.tri(mfdr[leaflist[[i]],leaflist[[i]]])])
}
hc$mergefdr<-mergefdr
hc$meanfdr<-meanfdr
hc$nodesize<-nodesize
hc$leaflist<-leaflist
hc$labellist<-labellist
hc$sharing<-sharing
hc$complexity<-complexity
rm(mergefdr,meanfdr,nodesize,leaflist,labellist,sharing,complexity)
shareacross<-colSums(hc$sharing>sharemin)	#Number of features (approximately) shared across the node
#A node is considered compliant if FDR is below and its sharing across above threshod for the node
#and all its descendants
if(bymax)compliant<-(hc$mergefdr<fdrthresh&(shareacross-shareacross[nrow(hc$merge)])>nshare)
if(!bymax)compliant<-(hc$meanfdr<fdrthresh&(shareacross-shareacross[nrow(hc$merge)])>nshare)
leftchild<-(hc$merge[,1]<0)
leftchild[hc$merge[,1]>0]<-compliant[hc$merge[hc$merge[,1]>0,1]]
rightchild<-(hc$merge[,2]<0)
rightchild[hc$merge[,2]>0]<-compliant[hc$merge[hc$merge[,2]>0,2]]
newcompliant<-compliant&leftchild&rightchild
while(!all(newcompliant==compliant)){
	compliant<-newcompliant
	leftchild<-(hc$merge[,1]<0)
	leftchild[hc$merge[,1]>0]<-compliant[hc$merge[hc$merge[,1]>0,1]]
	rightchild<-(hc$merge[,2]<0)
	rightchild[hc$merge[,2]>0]<-compliant[hc$merge[hc$merge[,2]>0,2]]
	newcompliant<-compliant&leftchild&rightchild
}
#Clone nodes are maximum compliant nodes
clonenodes<-setdiff((1:nrow(hc$merge))[compliant],c(hc$merge[compliant,1],hc$merge[compliant,2]))
hc$fdrthresh<-fdrthresh
hc$clonenodes<-clonenodes
hc$bymax<-bymax
hc$shareacross<-shareacross
hc$sharemin<-sharemin
hc$nshare<-nshare
rm(clonenodes,shareacross)
cloneleaves<-unique(unlist(hc$leaflist[hc$clonenodes]))
coreid<-jguide[match(newcellnames,jguide[,"seq.unit.id"]),"sector"][cloneleaves]
hc$coresNclones<-table(coreid,coreid)
if(graphic){
	plot(hc,labels=F)
	abline(h=hc$height[hc$clonenodes],lty=2)
}
hcf<-hclust(as.dist(mfdr),method=hcmethod)
leaflist<-vector(mode="list",length=nrow(hcf$merge))
nodesize<-rep(NA,nrow(hcf$merge))
labellist<-vector(mode="list",length=nrow(hcf$merge))
sharing<-matrix(ncol=nrow(hcf$merge),nrow=nrow(pinmat))
complexity<-rep(NA,nrow(hcf$merge))
for(i in 1:nrow(hcf$merge)){
	if(hcf$merge[i,1]<0)leaflist[[i]]<-(-hcf$merge[i,1])
	else leaflist[[i]]<-leaflist[[hcf$merge[i,1]]]
	if(hcf$merge[i,2]<0)leaflist[[i]]<-c(leaflist[[i]],(-hcf$merge[i,2]))
	else leaflist[[i]]<-c(leaflist[[i]],leaflist[[hcf$merge[i,2]]])
	labellist[[i]]<-hcf$labels[leaflist[[i]]]
	complexity[i]<-mean(colSums(pinmat[,labellist[[i]]]))
	sharing[,i]<-rowMeans(pinmat[,labellist[[i]]])
	nodesize[i]<-length(leaflist[[i]])
}
hcf$nodesize<-nodesize
hcf$leaflist<-leaflist
hcf$labellist<-labellist
hcf$sharing<-sharing
hcf$complexity<-complexity
rm(nodesize,leaflist,labellist,sharing,complexity)
shareacross<-colSums(hcf$sharing>sharemin)
compliant<-(hcf$height<fdrthresh&shareacross>nshare)
leftchild<-(hcf$merge[,1]<0)
leftchild[hcf$merge[,1]>0]<-compliant[hcf$merge[hcf$merge[,1]>0,1]]
rightchild<-(hc$merge[,2]<0)
rightchild[hcf$merge[,2]>0]<-compliant[hcf$merge[hcf$merge[,2]>0,2]]
newcompliant<-compliant&leftchild&rightchild
while(!all(newcompliant==compliant)){
	compliant<-newcompliant
	leftchild<-(hcf$merge[,1]<0)
	leftchild[hcf$merge[,1]>0]<-compliant[hcf$merge[hcf$merge[,1]>0,1]]
	rightchild<-(hcf$merge[,2]<0)
	rightchild[hcf$merge[,2]>0]<-compliant[hcf$merge[hcf$merge[,2]>0,2]]
	newcompliant<-compliant&leftchild&rightchild
}
clonenodes<-setdiff((1:nrow(hcf$merge))[compliant],c(hcf$merge[compliant,1],
	hcf$merge[compliant,2]))
hcf$fdrthresh<-fdrthresh
hcf$sharemin<-sharemin
hcf$sharemin<-nshare
hcf$clonenodes<-clonenodes
hang<-rep(0,length(hcf$height))
hang[hcf$merge[hcf$merge[,1]>0,1]]<-
	hcf$height[hcf$merge[,1]>0]-hcf$height[hcf$merge[hcf$merge[,1]>0,1]]
hang[hcf$merge[hcf$merge[,2]>0,2]]<-
	hcf$height[hcf$merge[,2]>0]-hcf$height[hcf$merge[hcf$merge[,2]>0,2]]
hcf$hang<-hang
hcf$shareacross<-shareacross
rm(clonenodes,hang)
cloneleaves<-unique(unlist(hcf$leaflist[hcf$clonenodes]))
coreid<-jguide[match(newcellnames,jguide[,"seq.unit.id"]),"sector"][cloneleaves]
hcf$coresNclones<-table(coreid,coreid)
if(graphic){
	plot(hcf,labels=F)
	abline(h=hcf$height[hcf$clonenodes],lty=2)
}
pytableP<-TreePy(data=as.dist(mdist),method="average")
pytableP<-cbind(pytableP,hc$mergefdr)
dimnames(pytableP)[[2]][ncol(pytableP)]<-"log10fdr"
pytableF<-TreePy(data=as.dist(mfdr),method="average")
hcname<-paste(prostate,"smear1bpLog10FisherHCP",sep="")
assign(hcname,hc)
save(list=hcname,file=paste(savedir,hcname,".rda",sep=""))
hcfname<-paste(prostate,"smear1bpLog10FisherHCF",sep="")
assign(hcfname,hcf)
save(list=hcfname,file=paste(savedir,hcfname,".rda",sep=""))
write.table(mdist,paste(savedir,prostate,"smear1bpLog10FisherP.txt",sep=""),col.names=T,
	row.names=T,sep="\t",quote=F)
write.table(mfdr,paste(savedir,prostate,"smear1bpLog10FisherFDR.txt",sep=""),
	col.names=T,row.names=T,sep="\t",quote=F)
write.table(pytableP,paste(savedir,prostate,"smear1bpFisherTreePyP.txt", sep=""),col.names=T,
	row.names=F,sep="\t",quote=F)
write.table(pytableF,paste(savedir,prostate,"smear1bpFisherTreePyF.txt", sep=""),col.names=T,
	row.names=F,sep="\t",quote=F)
#savehistory(paste(savedir,"RhistoryClones.",Sys.Date(),".",
#  strsplit(as.character(Sys.time()),split=" ")[[1]][2],".txt",sep=""))
quit(save="no")
