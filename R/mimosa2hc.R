#' This function returns an hclust object coresponding to the partition tree
#' computed by mimain and given by the argument mymi of class "mimosa". The node
#' heights are computed recursively, by taking the maximum of the heights of the
#' two merging nodes and adding to those the mymi$height of the parent node.
mimosa2hc<-function(mymi){
	sortby<-mymi$pathcode
	rawone<-as.raw(T)
	raweighty<-rawShift(rawone,7)
	alldone<-(sortby&raweighty)==raweighty
	while(!all(alldone)){
		sortby[!alldone]<-rawShift(sortby[!alldone],1)
		alldone<-(sortby&raweighty)==raweighty
	}
	sortby<-as.integer(sortby)
	icode<-as.integer(mymi$pathcode)
	iorder<-order(sortby,decreasing=T)
	isort<-icode[iorder]
	usort<-unique(isort)
	mergemat<-matrix(ncol=2,nrow=length(isort)-1)
	hivec<-rep(0,length(isort)-1)
	starthere<-1
	#nvals<-rev(tapply(match(isort,usort),match(isort,usort),length))
	nvals<-tapply(match(isort,usort),match(isort,usort),length)
	noderow<-(-match(usort,icode))
	for(ind in 1:length(usort)){
		if(nvals[ind]>1){
			myord<-iorder[isort==usort[ind]]
			mergemat[starthere,]<-(-myord[1:2])
			hivec[starthere]<-0
			if(nvals[ind]>2){
				mergemat[(starthere+1):(starthere+nvals[ind]-2),1]<-(-myord[-(1:2)])
				mergemat[(starthere+1):(starthere+nvals[ind]-2),2]<-
					starthere:(starthere+nvals[ind]-3)
				hivec[(starthere+1):(starthere+nvals[ind]-2)]<-0
			}
			noderow[ind]<-starthere+nvals[ind]-2
			starthere<-starthere+nvals[ind]-1
		}
	}
	while(length(usort)>1){
		parents<-usort%/%2
		maxpar<-max(parents)
		mergeus<-which(parents==maxpar)
		mergemat[starthere,]<-sort(noderow[mergeus])
		if(all(noderow[mergeus]<0))mergemat[starthere,]<-rev(mergemat[starthere,])
		hiadd<-0
		if(mergemat[starthere,1]>0)hiadd<-hivec[mergemat[starthere,1]]
		if(mergemat[starthere,2]>0)hiadd<-max(hivec[mergemat[starthere,2]],hiadd)
		hivec[starthere]<-hiadd+mymi$height[match(maxpar,as.integer(mymi$upath))]
		noderow<-c(noderow[-mergeus],starthere)
		usort<-c(usort[-mergeus],maxpar)
		starthere<-starthere+1
	}
	mihc<-list(merge=mergemat,height=hivec,order=iorder,labels=mymi$leafnames)
	class(mihc)<-"hclust"
	mihc
}
