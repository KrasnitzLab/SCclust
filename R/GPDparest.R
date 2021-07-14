#Given a vector "exc" of exceedances, estimate the parameters of the 
#generalized Pareto distribution using the method of moments.
#The alpha and k parameter estimates are according to 
#J.R.M. Hosking and J.R. Wallis, Technometrics 29(3) (1987) 339-349.
#Value: a list containing a moments estimate of the GPD parameters alpha
#(scale) and k (shape), and large-sample approximation for their 
#covariance matrix. Note that the covariance matrix does not exist for 
#k<= -0.25, although its elements can still be formally computed.
GPDparest<-function(exc){
	empmean<-mean(exc)
	empvar<-var(exc)
	estalpha<-0.5*empmean*(1+empmean^2/empvar)
	estk=0.5*(-1+empmean^2/empvar)
	if(estk<=-0.25)warning("GPD k<=-1/4, unreliable estimate!")
	infront<-(1+estk)^2/(1+2*estk)/(1+3*estk)/(1+4*estk)/length(exc)
	offdiag<-estalpha*(1+2*estk)*(1+4*estk+12*estk^2)
	covmat<-matrix(nrow=2,data=infront*c(2*estalpha^2*(1+6*estk+12*estk^2),
		offdiag,offdiag,(1+2*estk)^2*(1+estk+6*estk^2)))
	list(estalpha=estalpha,estk=estk,covmat=covmat)
}
