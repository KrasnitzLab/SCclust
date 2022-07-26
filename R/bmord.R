bmord<-function(m){	#Lexicographic column order of a binary matrix a dumb way
	if(all(bm[1,]))bm
	bm<-bm[,order(bm[1,]),drop=F]
	bm[-1,bm[1,]==0]<-bmord(bm[-1,bm[1,]==0,drop=F])
	bm[-1,bm[1,]==1]<-bmord(bm[-1,bm[1,]==1,drop=F])
	bm
}
