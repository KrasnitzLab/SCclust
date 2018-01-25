clonetype<-"soft"
subcloneTooBig<-0.8
hcfile<-system("ls|grep smear1bpLog10FisherHCP.rda",intern=T)
prostate<-substring(hcfile,first=1,
	last=nchar(hcfile)-nchar("smear1bpLog10FisherHCP.rda"))
load(hcfile)
hc<-get(substring(hcfile,first=1,last=nchar(hcfile)-4))
rm(list=substring(hcfile,first=1,last=nchar(hcfile)-4))
clonetable<-data.frame(hc$labels,rep(0,length(hc$labels)),
	rep(0,length(hc$labels)),stringsAsFactors=F)
dimnames(clonetable)[[2]]<-c("ID","clone","subclone")
for(nodes in unique(hc$softclones[clonetype,]))
	clonetable[hc$leaflist[[nodes]],"clone"]<-nodes
clonedirs<-system("ls|grep cloneP",intern=T)
if(length(clonedirs)>0)for(dirs in clonedirs){
	setwd(dirs)
	hcfile<-system("ls|grep smear1bpLog10FisherHCP.rda",intern=T)
	load(hcfile)
	hc<-get(substring(hcfile,first=1,last=nchar(hcfile)-4))
	rm(list=substring(hcfile,first=1,last=nchar(hcfile)-4))
	clunique<-unique(hc$softclones[clonetype,])
	if(length(clunique)>1)for(nodes in clunique)
		clonetable[match(hc$labellist[[nodes]],clonetable[,1]),"subclone"]<-nodes
	if(length(clunique)==1)
		if(hc$nodesize[clunique]<(subcloneTooBig*max(hc$nodesize)))
			clonetable[match(hc$labellist[[clunique]],clonetable[,1]),"subclone"]<-
				clunique
	setwd("..")
}
write.table(clonetable,paste(prostate,"smear1bpFisherPcloneTracks.txt",sep=""),
	col.names=T,row.names=F,sep="\t",quote=F)
