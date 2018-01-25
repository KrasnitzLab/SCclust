hcClimb<-function(hc,minsize,minshare,sectors=NULL){
	if(is.null(sectors))sectors<-rep(1,length(hc$labels))
	else if(length(sectors)!=length(hc$labels))
		stop("Numbers of leaves and sector assignments do not match")
	hard2soft<-matrix(nrow=2,ncol=0,dimnames=list(c("hard","soft"),NULL))
	hardsectors<-unique(sectors[unique(
		unlist(hc$leaflist[hc$clonenodes[hc$nodesize[hc$clonenodes]>minsize]]))])
	for(cnode in hc$clonenodes)if(hc$nodesize[cnode]>minsize){
		nodenow<-cnode
		ancestor<-cnode
		#while(nodenow<nrow(hc$merge)&hc$shareacross[row(hc$merge)[hc$merge==nodenow]]>=minshare)
		while(nodenow<nrow(hc$merge)){
			nodenow<-row(hc$merge)[hc$merge==nodenow]
			if(hc$shareacross[nodenow]>=minshare&
				all(sectors[hc$leaflist[[nodenow]]]%in%hardsectors))ancestor<-nodenow
		}
		hard2soft<-cbind(hard2soft,c(cnode,ancestor))
	}
	hard2soft
}
