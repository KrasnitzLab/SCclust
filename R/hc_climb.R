

hc_climb<-function(hc, minsize, minshare){
  hard2soft<-matrix(nrow=2,ncol=0,dimnames=list(c("hard","soft"),NULL))
  for(cnode in hc$clonenodes) {
    if(hc$nodesize[cnode]>minsize){
          nodenow<-cnode
          ancestor<-cnode
          while(nodenow<nrow(hc$merge)){
              nodenow<-row(hc$merge)[hc$merge==nodenow]
              if(hc$count_pins_share[nodenow]>=minshare)ancestor<-nodenow
          }
          
          hard2soft<-cbind(hard2soft,c(cnode,ancestor))
      }
  }
  return(hard2soft)
}