#The overall goal of this script is to determine whether clones identifiend by
#makeClonesSoft contain subclones. This is done, for each clone, by setting up
#a sub-matrix of "pinmat", whose columns are cells in the clone. We remove
#from the sub-matrix rows that are all 0s or all 1s., then analyze it +-
#exactly as we did the "pinmat" of the entire "prostate". For each clone, the 
#results are saved in a separate sub-directory with a long, informative name.
#
#Note that this is a script with a lot of hacks, not a properly structured
#module or function. I was running out of time and wished to analyze multiple 
#patient cases in one run.
#
#load a tree, its pin table and dump pin sub-tables and trees for larger clones
require(parallel)	#included in R core
usecores<-detectCores()-4	#wigclust nodes have 32 cores; leave 4 free
#hcClimb finds "soft" clones, starting with "hard" clones
source("/Volumes/user/krasnitz/prostateSingleCell/Rtools/hcClimbNew.R")
#it should not be too difficult to get rid of TBEST
require(TBEST)
#code for computing p-values on randomized "pin" matrices
source("/Volumes/user/krasnitz/prostateSingleCell/Rtools/simFisher.R")
#code for conversion from R to Python tree objects
source("/Volumes/user/krasnitz/prostateSingleCell/Rtools/TreePy.R")
#
#casemat is set up to loop over multiple patient cases; I don't think this
#loop should be part of the pipeline
casemat<-matrix(ncol=2,data=c("Gleason9.1","Gleason9.2","NYU001Gl7","NYU004Gl7.4","NYU005Gl6.2",
"NYU007Gl7.2","NYU010Gl7.3","NYU011Gl7.5","GL9.1","GL9.2","nyu001.GL7.1","nyu004.GL7.4","nyu005.GL6.2","nyu007.GL7.2","nyu010.GL7.3","nyu011.GL7.5"))
#
#loop over patient cases begins here
#
for(casenums in 1:nrow(casemat)){
#find and read all we need for the patient case
hcdir<-paste("/Volumes/user/krasnitz/prostateSingleCell/breakPointPins/",casemat[casenums,1],
	"/averageNoEndsNewest/",sep="")
pindir<-hcdir
hctype<-"P"	#i.e., the tree will be grown using Fisher log p-values as
						#dissimilarities
prostate<-casemat[casenums,2]
nmin<-6	#do not look for sub-clones in clones smaller than this
nsim<-500	#the number of randomizations to sample from the null distribution 
					#for Fisher's p-value
#
hcmethod<-"average"	#use average linkage to grow the tree
log10data<-F	#the input for clone analysis will be p-values, not their logs
fdrthresh<-(-2)	#here through "graphic" the vars mean the same thing as in
								#findClonesSoft
sharemin<-0.85
baseshare<-3	#cells in a clone share more features than cells in the whole
							#prostate
bymax<-T
lmmax<-0.001
climbfromsize<-2
climbtoshare<-3
clonedef<-"soft"
graphic<-F
#load the R tree object for the case
load(paste(hcdir,prostate,"smear1bpLog10FisherHC",hctype,".rda",sep=""))
#give it a short name
hcorig<-get(paste(prostate,"smear1bpLog10FisherHC",hctype,sep=""))
rm(list=paste(prostate,"smear1bpLog10FisherHC",hctype,sep=""))
#read "pinmat" and "pins" for the case
pinmat<-
	as.matrix(read.table(paste(pindir,prostate,"smear1bpPinMat.txt",sep=""),header=T,as.is=T))
pins<-read.table(paste(pindir,prostate,"smear1bpPins.txt",sep=""),header=T,as.is=T)
gc()	#I tend to gc a lot; does not hurt, I suppose
#identify clones with at least nmin cells
bigclones<-
	unique(hcorig$softclones[clonedef,])[hcorig$nodesize[unique(hcorig$softclones[clonedef,])]>=nmin]
#loop over big clones
if(length(bigclones)>0)for(clone in bigclones){
#find the clone sub-matrix of the "pinmat' and the corresponding pins
	clonemat<-pinmat[,dimnames(pinmat)[[2]]%in%hcorig$labellist[[clone]],drop=F]
	clonepins<-pins[(rowSums(clonemat)>0)&(rowSums(clonemat)<ncol(clonemat)),,drop=F]
	clonemat<-clonemat[(rowSums(clonemat)>0)&(rowSums(clonemat)<ncol(clonemat)),,drop=F]
	nshare<-baseshare+sum(rowSums(clonemat)>(sharemin*ncol(clonemat)))
	cellnames<-dimnames(clonemat)[[2]]
	#descriptive name for the clone directory; make it and write into it
	clonedir<-paste(hcdir,"clone",hctype,as.character(clone),".",Sys.Date(),".",
		strsplit(as.character(Sys.time()),split=" ")[[1]][2],"/",sep="")
	system(paste("mkdir ",clonedir,sep=""))
	write.table(clonemat,paste(clonedir,prostate,hctype,as.character(clone),
		"smear1bpPinMat.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
	write.table(clonepins,paste(clonedir,prostate,hctype,as.character(clone),
		"smear1bpPins.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
	#prepare input and run simFisher for the clone, exactly as for the prostate
	m<-vector(mode="list",length=length(unique(clonepins[,"sign"])))
	for(i in 1:length(unique(clonepins[,"sign"])))m[[i]]<-
		clonemat[clonepins[,"sign"]==unique(clonepins[,"sign"])[i],,drop=F]
	vtrue<-simFisher(m,nsim=1,nsweep=0,seedme=123,distrib="Rparallel",njobs=usecores,
		combo="Fisher")
	write(vtrue,file=paste(clonedir,prostate,hctype,as.character(clone),"trueP.txt",sep=""))
	msim<-simFisher(m,nsim=nsim,nsweep=200,seedme=123,distrib="Rparallel",njobs=usecores,
		combo="Fisher")
	write(msim,file=paste(clonedir,prostate,hctype,as.character(clone),"simP.txt",sep=""))
	rm(clonepins,m)
	gc()
	#Now we have both simulated and observed p-values; proceed as in findClones
	#Sort Fisher p's sampled from the null
	smsim<-sort(msim)
	rm(msim)
	gc()
	usmsim<-unique(smsim)
	#for each of the unique p-values in smsim, find its count in smsim
	csim<-tapply(match(smsim,usmsim),match(smsim,usmsim),length)
	rm(smsim)
	gc()
	if(log10data)usmsim<-10^usmsim
	if(log10data)vtrue<-10^vtrue
	pos<-(1:length(vtrue))[order(vtrue)]	#unused?
	svtrue<-sort(vtrue)
	usvtrue<-unique(svtrue)
	#ctrue are counts of observed unique Fisher p's
	ctrue<-tapply(match(svtrue,usvtrue),match(svtrue,usvtrue),length)
	#Here we model the log of the cumulative null distribution of Fisher p 
	#as linear in log(Fisher p) in the range where the cumulative distribution
	# is < lmmmax. This will fail if the count of lowest three unique Fisher p's
	# is > lmmax of the total sample (we need at least 3 points for a linear
	# fit).
	mylm<-lm(log(cumsum(csim[(cumsum(csim)/sum(csim))<lmmax])/sum(csim))~
  	log(usmsim[(cumsum(csim)/sum(csim))<lmmax]))
	#If the model holds in the entire range of observed Fisher p's, the FDR
	#for any observed p is the ratio of the model prediction at that p to the
	#observed cumlattive distribution at that p
	if(exists("mylm"))logfdrmod<-
  	mylm$coefficients[2]*log(usvtrue)+mylm$coefficients[1]-log(cumsum(ctrue)/sum(ctrue))
	if(graphic){
  	curve(exp(mylm$coefficients[1]+mylm$coefficients[2]*log(x)),from=min(c(usvtrue,usmsim)),
			to=max(c(usvtrue,usmsim)),log="xy")
  	points(usmsim,cumsum(csim)/sum(csim))
  	points(usvtrue,cumsum(ctrue)/sum(ctrue),col="red")
	}
	#If an observed p is in the range of p's sampled from the null, we will 
	#compute FDR by interpolation; if it is smaller than any p sampled from the
	#null, we will use the linear model to extrapolate.
	#For each observed p, find the bracketing values of the sampled p
	z<-cbind(c(log(usvtrue),log(usmsim)),c(log(cumsum(ctrue)/sum(ctrue)),
		log(cumsum(csim)/sum(csim))),c(rep(0,length(usvtrue)),rep(1,length(usmsim))))
	z<-z[order(z[,1]),]
	simlow<-cumsum(z[,3])[z[,3]==0]
	#Logical: these values of the observed p can be interpolated
	valid<-(simlow>0)&(simlow<sum(z[,3]))
	#Left bracketing sampled p index
	x1pos<-match(simlow,cumsum(z[,3]))[valid]
	#Right bracketing sampled p index
	x2pos<-match(simlow+1,cumsum(z[,3]))[valid]
	x1<-z[x1pos,1]	#Left bracketing sampled p value
	x2<-z[x2pos,1]	#Right bracketing sampled p value
	y1<-z[x1pos,2]	#Cumulative null at left bracket
	y2<-z[x2pos,2]	#Cumulative null at right bracket
	#Compute linear interpolation for log(FDR), i.e.
	#log(FDR) = [ interpolated log(cumulative null) - log(cumulative observed) ],
	#where valid (where it is possible to interpolate).
	logfdrinterp<-rep(0,length(usvtrue))
	logfdrinterp[valid]<-(y2-y1)*log(usvtrue)[valid]/(x2-x1)+(y1*x2-y2*x1)/(x2-x1)-
  log(cumsum(ctrue)/sum(ctrue))[valid]
	logfdr<-logfdrinterp
	#The largest p sampled from the null for which cumulative null < lmmax
	#Gets -Infinity is there is no such p in the sample
	lmu<-max(usmsim[(cumsum(csim)/sum(csim))<lmmax])
	#We wish a smooth transtition between the extrapolated and the interpolated 
	#FDR; we therefore interpolate between the linear model prediction for FDR
	#and its interpolated value in the region 
	#min(p from the null) < (observed p) < lmu
	if(is.finite(lmu)&min(usvtrue)<min(usmsim)){
  	logfdr[usvtrue<lmu&usvtrue>min(usmsim)]<-
			(logfdrmod[usvtrue<lmu&usvtrue>min(usmsim)]*
			(log(lmu)-log(usvtrue[usvtrue<lmu&usvtrue>min(usmsim)]))-
			logfdrinterp[usvtrue<lmu&usvtrue>min(usmsim)]*
			(log(min(usmsim))-log(usvtrue[usvtrue<lmu&usvtrue>min(usmsim)])))/
			(log(lmu)-log(min(usmsim)))
		logfdr[usvtrue<min(usmsim)]<-logfdrmod[usvtrue<min(usmsim)]
	}
	#Make sure that FDR never decrease...
	logfdr<-cummax(logfdr)
	#... and never exceed 1
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
	write.table(mdist,paste(clonedir,prostate,hctype,as.character(clone),
		"smear1bpLog10FisherP.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
	write.table(mfdr,paste(clonedir,prostate,hctype,as.character(clone),
		"smear1bpLog10FisherFDR.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
	hc<-hclust(as.dist(mdist),method=hcmethod)
	leaflist<-vector(mode="list",length=nrow(hc$merge))
	mergefdr<-rep(NA,nrow(hc$merge))
	meanfdr<-rep(NA,nrow(hc$merge))
	nodesize<-rep(NA,nrow(hc$merge))
	labellist<-vector(mode="list",length=nrow(hc$merge))
	sharing<-matrix(ncol=nrow(hc$merge),nrow=nrow(pinmat))
	complexity<-rep(NA,nrow(hc$merge))
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
	shareacross<-colSums(hc$sharing>sharemin)
	if(bymax)compliant<-(hc$mergefdr<fdrthresh&shareacross>nshare)
	if(!bymax)compliant<-(hc$meanfdr<fdrthresh&shareacross>nshare)
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
	clonenodes<-setdiff((1:nrow(hc$merge))[compliant],
		c(hc$merge[compliant,1],hc$merge[compliant,2]))
	hc$fdrthresh<-fdrthresh
	hc$clonenodes<-clonenodes
	hc$bymax<-bymax
	hc$shareacross<-shareacross
	hc$sharemin<-sharemin
	hc$nshare<-nshare
	rm(clonenodes,shareacross)
	if(!is.null(hc$clonenodes))hc$softclones<-hcClimb(hc,minsize=climbfromsize,
		minshare=climbtoshare+hc$shareacross[nrow(hc$merge)])
	if(graphic){
		plot(hc,labels=F)
		abline(h=hc$height[hc$clonenodes],lty=2)
	}
	pytableP<-TreePy(data=as.dist(mdist),method="average")
	pytableP<-cbind(pytableP,hc$mergefdr)
	write.table(pytableP,paste(clonedir,prostate,hctype,as.character(clone),
		"smear1bpFisherTreePyP.txt", sep=""),col.names=T,row.names=F,sep="\t",quote=F)
	hcname<-paste(prostate,hctype,as.character(clone),"smear1bpLog10FisherHCP",sep="")
	assign(hcname,hc)
	save(list=hcname,file=paste(clonedir,hcname,".rda",sep=""))
}	
}
quit("no")
