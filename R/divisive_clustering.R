
squeeze_vector <- function(v){
    if(!is.raw(v))
        v <- as.raw(v)
    return(
        packBits(
            c(v, rep(as.raw(F), (8 - length(v) %% 8) %% 8)), type="raw"))
}

showBits <- function(r) {
    stats::symnum(as.logical(rawToBits(r)))
}


build_incidence_table <- function(m) {
    assertthat::assert_that(class(m) == "matrix")
    
    flog.warn("elements of incidence will be coerced to raw")
    incidence <- list(
            apply(m, 2, squeeze_vector),
            apply(m, 1, squeeze_vector))
    
    if(is.null(dim(incidence[[1]]))) {
        dim(incidence[[1]]) <- c(1, length(incidence[[1]]))
    }
    if(is.null(dim(incidence[[2]]))) {
        dim(incidence[[2]]) <- c(1, length(incidence[[2]]))
    }

    class(incidence)<-"incidencetable"

    return(incidence)    
}

load_incidence_table <- function(filename) {
    m <- load_matrix(filename)
    incidence <- build_incidence_table(m)

    return(incidence)
}


replicate_incidence <- function(incidence, from) {
    assertthat::assert_that(from == 1 | from == 2)
    if (from == 1) 
        to <- 2
    else {
       to <- 1
    }

    rawone <- as.raw(T)

	inblocks  <- ((1 : ncol(incidence[[to]])) - 1) %/% 8 + 1
	whichbits <- ((1 : ncol(incidence[[to]])) - 1) %% 8 + 1

    # func <- function(v) {
    #     return( squeeze_vector(rawShift(incidence[[from]][v[1],], 1 - v[2]) & rawone) )
    # }

	incidence[[to]] <- apply(
        cbind(inblocks, whichbits), 1, 
        function(v) squeeze_vector(rawShift(incidence[[from]][v[1],], 1 - v[2]) & rawone)
    )

    return(incidence)
}

consolidate_incidence <- function(incidence, from) {
    assertthat::assert_that(from == 1 | from == 2)
    if (from == 1) 
        to <- 2
    else {
       to <- 1
    }

	allup<-c(rep(!as.raw(F), ncol(incidence[[to]])%/%8),
		squeeze_vector(rep(as.raw(T),ncol(incidence[[to]])%%8)))
	alldown<-rep(as.raw(F),nrow(incidence[[from]]))

	incidence[[from]]<-incidence[[from]][ ,
        (colSums(incidence[[from]]==allup)!=nrow(incidence[[from]]))&
		(colSums(incidence[[from]]==alldown)!=nrow(incidence[[from]])),drop=F]

	return(replicate_incidence(incidence, from))
}


contingencies<-function(incidence, partition){

	pinc <- colSums(
        matrix(
            ncol=ncol(incidence[[2]]),
            data=sapply(partition & incidence[[2]], setBits)))
	rsi <- colSums(
        matrix(
            ncol=ncol(incidence[[2]]),
		    data=sapply(incidence[[2]], setBits)))

	sp  <- sum(sapply(partition, setBits))
	rsp <- rep(sp, ncol(incidence[[2]]))
	lp  <- ncol(incidence[[1]])
	rlp <- rep(lp, ncol(incidence[[2]]))

    contables  <- cbind(pinc, rsi - pinc, rsp - pinc, rlp - rsi - rsp + pinc)
    pmarginals <- c(sp, lp - sp)

    return(list(contables=contables, pmarginals=pmarginals))
}

misum <- function(contingencies) {
    contribs<-(allcounts$contables!=0)
    result <- sum(allcounts$contables[contribs]*log(allcounts$contables[contribs]))-
        nrow(allcounts$contables)*sum(allcounts$pmarginals[allcounts$pmarginals!=0]*
        log(allcounts$pmarginals[allcounts$pmarginals!=0]))
    return(result)
}
