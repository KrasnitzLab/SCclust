#' Pack the elements of a (0,1)-valued integer, numeric, raw or (F,T)-valued 
#'logical vector 
#' yesno as bits into a vector of type raw. If the length of yesno is not a 
#' multiple of 8 bits of type into, it is padded by the raw 00 to the 
#' nearest such multiple.
squeeze_vector <- function(v){
    if(!is.raw(v))
        v <- as.raw(v)
    return(
        packBits(
            c(v, rep(as.raw(F), (8 - length(v) %% 8) %% 8)), type="raw"))
}

showBits <- function(r) {
    paste(stats::symnum(as.logical(rawToBits(r))), collapse="")
}


build_incidence_table <- function(m) {
	if(data.class(m) == "data.frame") {
		m <- data.matrix(m);
	}

    assertthat::assert_that(data.class(m) == "matrix")
    
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

incidence_table2matrix <- function(incidence) {
	cols <- ncol(incidence[[1]])
	rows <- ncol(incidence[[2]])

	inblocks  <- ((1 : ncol(incidence[[2]])) - 1) %/% 8 + 1
	whichbits <- ((1 : ncol(incidence[[2]])) - 1) %% 8 + 1

    rawone <- as.raw(T)

	r <- matrix(rep(0, rows * cols), ncol=cols)
	indexes <- cbind(inblocks, whichbits)

	for(row in 1:rows) {
		index <- indexes[row, ]
		r[row, ] <- as.integer(rawShift(incidence[[1]][index[1],], 1 - index[2]) & rawone)
	}
	# assertthat::assert_that(ncol(r) == cols)
	# assertthat::assert_that(nrow(r) == rows)

	return(r)
}


#' incidence[[1]] is the original (F,T)-valued incidence matrix, packed column by
#' column into raws. incidence[[2]] is the same matrix, packed row by row (or the 
#' transposed matrix, packed column by column) into raws. Re-create
#' incidence[[3-from]] from incidence[[from]]. This is useful if, e.g., all-0
#' and/or all-1 columns are dropped from incidence[[from]]. It is simpler to
#' restore consistency by re-creating the [[3-from]] from the [[from]] 
#' than by modifying the rows in the [[3-from]].
replicate_incidence <- function(incidence, from) {
    assertthat::assert_that(from == 1 | from == 2)

	assertthat::assert_that(!is.null(dim(incidence[[1]])))
	assertthat::assert_that(!is.null(dim(incidence[[1]])))

    if (from == 1) 
        to <- 2
    else {
       to <- 1
    }

		if(ncol(incidence[[from]])==0){
			incidence[[to]]<-matrix(ncol=ncol(incidence[[to]]),nrow=0)
			return(incidence)
		}

    rawone <- as.raw(T)

	inblocks  <- ((1 : ncol(incidence[[to]])) - 1) %/% 8 + 1
	whichbits <- ((1 : ncol(incidence[[to]])) - 1) %% 8 + 1

	incidence[[to]] <- apply(
        cbind(inblocks, whichbits), 1, 
        function(v) squeeze_vector(rawShift(incidence[[from]][v[1],], 1 - v[2]) & rawone)
    )

    if(is.null(dim(incidence[[1]]))) {
        dim(incidence[[1]]) <- c(1, length(incidence[[1]]))
    }
    if(is.null(dim(incidence[[2]]))) {
        dim(incidence[[2]]) <- c(1, length(incidence[[2]]))
    }
  	assertthat::assert_that(all(t(apply(incidence[[2]],2,
			rawToBits))[1:ncol(incidence[[2]]),1:ncol(incidence[[1]])]==
			apply(incidence[[1]],2,rawToBits)[1:ncol(incidence[[2]]),
			1:ncol(incidence[[1]])]))

    return(incidence)
}


#' This function removes uninformative all-0 and all-1 columns or rows from the
#' incidence table. We do the removal in incidence[[1]] or incidence[[2]], 
#' depending on whether columns or rows are removed. Then we call 
#' replicateIncidence to restore consistency.
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
		(colSums(incidence[[from]]==alldown)!=nrow(incidence[[from]])), drop=F]

	return(replicate_incidence(incidence, from))
}


#' incidence[[1]] is the original (F,T)-valued incidence matrix, packed column by
#' column into raws. incidence[[2]] is the same matrix, packed row by row (or the 
#' transposed matrix, packed column by column) into raws. partition is an
#' (F,T)-valued vector of the same length as the column number of  the original 
#' incidence matrix, packed into raws. Return a list with 2 items:
#' contables - matrix of 4 columns, where each row
#'	    contains a contingency table between partition and a row of incidence.
#' pmarginals - a vector of length 2 of marginals for partition
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

#' conti is a list of contingencies and marginals, as defined for the value of
#' function contingencies. It contains contingency tables and their marginals.
#' Compute mutual information in each contingency table and sum over all tables.
misum <- function(conti) {
    contribs<-(conti$contables!=0)
    result <- sum(conti$contables[contribs]*log(conti$contables[contribs]))-
        nrow(conti$contables)*sum(conti$pmarginals[conti$pmarginals!=0]*
        log(conti$pmarginals[conti$pmarginals!=0]))
    return(result)
}


#' The lists saspars and swappars contain, respectively, parameters of the 
#' simulated annealing schedule and of the incidence matrix randomization. The
#' meaning is as follows
#' saspars
#' 	 restarts: how many times SA must be restarted from a random partition
#' 	 cooler: reduce the temperature by this factor after every cycle
#' 	 sweepspercycle: a sweep means that a partition reassgnment is proposed once
#' 	 for each column in the incidence matrix. Perform this number of sweeps per 
#' 	 cycle
#' 	 maxcycles: the maximal number of cycles to execute
#' 	 stopatfreezeout: stop if the objective function is the same at the end of the
#' 	 current cycle as it was at the end of the previous one
#' 	 epsilon: ignore relative changes of the objective function smaller than this
#' swappars
#' 	 configs: generate this many randomized incidence matrices and, for each,
#' 	 compute the optimal objective function, to create a sampling from the null
#' 	 distribution
#' 	 burnin: perform this many swaps to get the 1st randomized incidence matrix
#' 	 permeas: perform this many swaps between subsequent objective function 
#' 	 computations
#' 	 for each swap, either a pair of columns or a pair of rows is chosen. The
#' 	 choice between columns and rows is made at random, and this is the probability
#' 	 of choosing columns
default_saspars <- list(
    restarts=10,
    cooler=1.12,
    acceptance=0.234,
    sweepspercycle=10,
	maxcycles=20,
    stopatfreezeout=T,
    epsilon=0.0001)

default_swappars <- list(
    configs=500,
    burnin=1000,
    permeas=500,
    choosemargin=0.5)

# #' pathcode is a vector of raws. Each element corresponds to a leaf in a node of a
# #' tree and encodes the path from the root to the leaf. depth is an integer 
# #' length of the path.
# nodesplit<-function(parition, pathcode) {
# 	rawToBits(partition)[1:length(pathcode)]&rawShift(pathcode,1)
# 	depth <- depth+1
# )



#' recompute mutual information with the column "updateme" of the incidence 
#' matrix moved from its current to the other part.
#' We use the incidence table packed as bits into raws row by row,
#' i.e., incidence[[2]], for this purpose.
alist.miupdate<-alist(updateme=,
{
	inblocks<-(updateme-1)%/%8+1
	whichbits<-(updateme-1)%%8+1
	myp<-as.logical(rawone&rawShift(partition[inblocks],1-whichbits))
	din<-1-2*myp
	newpmarginals<<-allcounts$pmarginals+c(din,-din)
	mybits<-as.logical(rawone&rawShift(incidence[[2]][inblocks,],1-whichbits))

	newct<<-allcounts$contables

	dp<-(!myp)-myp
	dcont<-mybits*dp
	newct[,1]<<-newct[,1]+dcont
	newct[,2]<<-newct[,2]-dcont

	dcont<-(!mybits)*dp
	newct[,3]<<-newct[,3]+dcont
	newct[,4]<<-newct[,4]-dcont
	
	contribs<-(newct!=0)
	sum(newct[contribs]*log(newct[contribs]))-
		nrow(newct)*sum(newpmarginals[newpmarginals!=0]*
		log(newpmarginals[newpmarginals!=0]))
}
)


#' for a given incidence matrix (in the calling environment) find a column
#' partition maximizing the objective function, defined as mutual 
#' information MI between the partition and each row, summed over all rows
mimax <- function(incidence, saspars=default_saspars) {
    rawone<-as.raw(T)
	
	partition <- squeeze_vector(
        sample(as.raw(c(T,F)), size=ncol(incidence[[1]]), replace=T))
	allcounts <- contingencies(incidence, partition)
	bestmi <- misum(allcounts)
	bestpartition <- partition

	miupdate<-as.function(alist.miupdate)

	for(newstart in 1 : saspars$restarts){

		partition <- squeeze_vector(
            sample(as.raw(c(T,F)), size=ncol(incidence[[1]]), replace=T))

	    allcounts <- contingencies(incidence, partition)
		minow <- misum(allcounts)
		newpmarginals <- allcounts$pmarginals
		newct <- allcounts$contables

		deltami <- sapply(1:ncol(incidence[[1]]), miupdate) - minow
		beta <- log(saspars$acceptance)/mean(-abs(deltami))/saspars$cooler

		for(cycles in 1 : saspars$maxcycles){
			beta <- beta * saspars$cooler
			miin <- minow
			for(sweeps in 1 : saspars$sweepspercycle){
				for(updateme in 1 : ncol(incidence[[1]])){
					newmi <- miupdate(updateme)
    				if((exp(beta*(newmi-minow))) > runif(1)){
						allcounts$contables <- newct
						allcounts$pmarginals <- newpmarginals
						inblocks <- (updateme-1)%/%8+1
						whichbits <- (updateme-1)%%8+1
						partition[inblocks]<-
							xor(partition[inblocks], rawShift(rawone,whichbits-1))
						minow <- newmi
						if(minow > bestmi){
							bestmi <- minow
							bestpartition <- partition
						}
					}
				}
			}
			if(saspars$stopatfreezeout & (abs(minow - miin) < saspars$epsilon))
                break
		}
	}
	partition <- bestpartition
	allcounts <- contingencies(incidence, partition)
	result <- list(
        partition=partition,
        contables=allcounts$contables,
		pmarginals=allcounts$pmarginals,
        beta=beta,
        cycles=cycles,
        mi=bestmi)
    return(result)
}


#' incidence[[1]] is the original (F,T)-valued incidence matrix, packed column by
#' column into raws. incidence[[2]] is the same matrix, packed row by row (or the 
#' transposed matrix, packed column by column) into raws. The following, if 
#' mymargin = 1, encodes choosing a random pair of columns in the original matrix 
#' and performing a maximal-possible raw- and column-sum preserving swap of Ts and
#' Fs between these two columns. The columns and rows are interchanged if 
#' mymargin = 2. The argument niter is the total number of swaps to be performed.
#' The swaps are applied to the incidence matrix in the calling environment.
mpshuffle<- function(incidence, niter, choosemargin=default_swappars$choosemargin) {
	for(iter in 1:niter){
		rawone <- as.raw(1)
		mymargin <- 1 + (runif(1) > choosemargin)
        assertthat::assert_that(mymargin == 1 | mymargin == 2)

		if(ncol(incidence[[3-mymargin]]) < 2) {
			next
		}

		weswap <- sample(ncol(incidence[[3-mymargin]]), size=2)

		inblocks <- (weswap - 1) %/% 8 + 1
		whichbits <- (weswap - 1) %% 8 + 1
    
		oneup <- rawone&rawShift(incidence[[mymargin]][inblocks[1],],1-whichbits[1])
		twoup <- rawone&rawShift(incidence[[mymargin]][inblocks[2],],1-whichbits[2])

		twonotone <- which((twoup&!oneup)==rawone)
		onenottwo <- which((oneup&!twoup)==rawone)

        if(length(onenottwo) == 0 | length(twonotone) == 0) {
            next
		}

		nswap <- min(length(onenottwo), length(twonotone))
		swapus <- sort(
            c(onenottwo[sample(length(onenottwo),size=nswap)],
			twonotone[sample(length(twonotone),size=nswap)]))

		incidence[[mymargin]][inblocks[1],swapus]<-
			xor(incidence[[mymargin]][inblocks[1],swapus], 
			rawShift(rawone,whichbits[1]-1))

		incidence[[mymargin]][inblocks[2],swapus]<-
			xor(incidence[[mymargin]][inblocks[2], swapus],
			rawShift(rawone,whichbits[2] - 1))
		swapusrows <- (swapus - 1) %/% 8 + 1
		swapusbits <- (swapus - 1) %% 8 + 1
		swapmask <- as.vector(tapply(swapusbits, swapusrows,
			function(v)Reduce("|", Map(rawShift, rawone,v-1))))
		urows<-unique(swapusrows)
		incidence[[3 - mymargin]][urows, weswap]<-
			xor(incidence[[3 - mymargin]][urows, weswap], swapmask)
	}
  	assertthat::assert_that(all(t(apply(incidence[[2]],2,
		rawToBits))[1:ncol(incidence[[2]]),1:ncol(incidence[[1]])]==
		apply(incidence[[1]],2,rawToBits)[1:ncol(incidence[[2]]),
		1:ncol(incidence[[1]])]))
    return(incidence)
}


#' repeatedly randomize the incidence table (incidence) and, for each randomized
#' table, maximize mutual information (MI) over all possible partitions using 
#' simulated annealing (SA). SA parameters are given by saspars and randomization
#' parameters by swappars. Return a vector of best MI values.
randomimax<-function(incidence, saspars=default_saspars, swappars=default_swappars){
	bestmirand<-rep(NA, swappars$configs)

	for(config in 1:swappars$configs){
        niter <- swappars$burnin * (config==1) + swappars$permeas * (config>1)
		# flog.debug("randomimax config=%s; niter=%s", config, niter)
		incidence <- mpshuffle(
            incidence,
            niter,
            choosemargin=swappars$choosemargin)
		bestmirand[config] <- mimax(incidence, saspars=saspars)$mi
	}
	return(bestmirand)
}


initial_pathcode <- function(incidence, maxgens=7) {
	if(maxgens > 7){
		flog.warn("initial_pathcode maxgens (%s) > 7; reseting to 7", maxgens)
		maxgens <- 7
	}

	rawone <- as.raw(T)	#least significant bit set
	rawzero <- as.raw(F)	#no bits set
	rawff <- !rawzero	#all bits set
	raweighty <- rawShift(rawone,7)	#most significant bit set
	firstpath <- rawzero

	for(pushme in 0:(7-maxgens))
		firstpath <- rawShift(firstpath,1) | rawone

	pathcode <- rep(firstpath, ncol(incidence[[1]]))

	return(pathcode)
}


#' perform divisive hierarchical clustering (recursive binary partitioning)
#' of data represented by a binary incidence matrix incidence. The leaves of
#' the hierarchical tree correspond to the columns of incidence. Other arguments
#' are maxgens, the maximal number of edges from the root to a leaf; maxempv,
#' the maximal empirical p-value below which the node split is accepted; saspars
#' and swappars lists as explained separately
#' incidence can be a data frame, a matrix or an object of class incidencetable,
#' i.e., a list of two items, each representing the incidence matrix packed into
#' raws. A data frame will be coerced to matrix, whose entries will then be 
#' coerced to raw. An object of class incidencetable will be generated from this
#' matrix.
#' The value is an object of class "mimosa", a list with items
#' incidence, the input incidence matrix as incidencetable; pathcode,a raw 
#' vector encoding, for each leaf, its path to the root; upath, a raw vector
#' encoding, for each node, its path to the root; height, the "height" of each
#' node above its 2 descendants, as given by -log(empirical p-value for the node
#' split); maxgens, the maxgens argument above; leafnames, the column names
#' of the input incidence matrix; call, the function call giving this value.
mimain<-function(incidence, 
	maxgens=7, 
	maxempv=0.05,
	saspars=default_saspars, 
	swappars=default_swappars){

	upath <- NULL
	height <- NULL

	#' Recursively partition items represented as columns of a binary incidence 
	#' table (IT). Thus, each row of the table represents a feature. Find a partition
	#' whose mutual information (MI) with the features is maximal. If the observed
	#' maximal MI is an outlier of a null distribution obtained by sampling from 
	#' randomized ITs, the optimal partition is accepted, and partitioning of each of 
	#' the two parts is attempted by calling this function recursively. 
	#' The items in the IT are interpreted as leaves of a node in a tree. The function
	#' call generates a branch of the tree with the current node as root.
	#' arguments: 
	#' 	incidence, a list with two items, IT packed into raws by the column
	#' 	and by the row, respectively.
	#' 	pathcode, a vector of as many raws as there are columns in the IT, encoding
	#' 	the path from the tree root to the current node.
	#' value: a raw vector, encoding, for each leaf, the path from the root to the 
	#' leaf. In addition, height and upath are updated in the calling environment.
	minode <- function(
		incidence, pathcode, 
		maxgens=7, maxempv=0.05,
		saspars=default_saspars, swappars=default_swappars) {

		incidence <- replicate_incidence(incidence, from=1)
		incidence <- consolidate_incidence(incidence, from=2)
		flog.debug("minode: incidence dim(%s)", ncol(incidence[[1]]))

		if((ncol(incidence[[1]])*nrow(incidence[[1]])*ncol(incidence[[2]])*nrow(incidence[[2]]))!=0){

			bestsplit <- mimax(incidence, saspars=saspars)
			nullmi <- randomimax(incidence, saspars=saspars, swappars=swappars)

			empv <- (sum(nullmi > bestsplit$mi) + 1) / (length(nullmi) + 2)

			flog.debug("minode: empv=%s; maxempv=%s", empv, maxempv)
			#cat("empv\t",empv,"\tmaxempv\t",maxempv,"\n")

			rawone <- as.raw(T)
			rawzero <- as.raw(F)

			if(empv < maxempv){
				#cat("Looking for upath, pathcode[1]=",pathcode[1],"\n")
				if(is.null(upath))
					upath <<-pathcode[1]
				else
					upath <<- c(upath,pathcode[1])
				#cat("upath set\t",upath,"\n")
				if(is.null(height))
					height<<-(-log(empv))
				else
					height <<- c(height,-log(empv))
				#cat("height set\t",height,"\n")
				longpartition <- rawToBits(bestsplit$partition)[1: ncol(incidence[[1]])]
				#	cat("longpartition\t",longpartition,"\npathcode before update\t",
				#		pathcode[1],"\t",rawToBits(pathcode[1]),"\n")
				pathcode <- longpartition | rawShift(pathcode,1)
				#	cat("pathcode after update \nfirst \t",
				#		pathcode[1],"\t",rawToBits(pathcode[1]),"\nlast\t",
				#		pathcode[length(pathcode)],"\t",
				#		rawToBits(pathcode[length(pathcode)]),"\n")
				raweighty <- rawShift(rawone, 7)

				if((pathcode[1] & raweighty) == raweighty)
					return(pathcode)

				subincidence <- list(
					incidence[[1]][,longpartition==rawone, drop=F], incidence[[2]])
				pathcode[longpartition==rawone] <-
					minode(
						subincidence, pathcode[longpartition==rawone],
						maxgens=maxgens, maxempv=maxempv,
						saspars=saspars, swappars=swappars)

				subincidence<-list(incidence[[1]][,longpartition==rawzero, drop=F],incidence[[2]])
				pathcode[longpartition==rawzero] <-
					minode(
						subincidence, pathcode[longpartition==rawzero],
						maxgens=maxgens, maxempv=maxempv,
						saspars=saspars, swappars=swappars)
			}
		}
		return(pathcode)
	}

	assertthat::assert_that(class(incidence) == "incidencetable")

	if(maxgens>7){
		warning("Input maxgens > 7, reset to 7")
		maxgens<-7
	}

	pathcode <- initial_pathcode(incidence, maxgens=maxgens)
	#cat("Initial pathcode [1]\t",rawToBits(pathcode[1]),"\n")
	upath<-NULL
	height<-NULL

	pathcode <- minode(
		incidence, pathcode,
		maxgens=maxgens, maxempv=maxempv,
		saspars=saspars, swappars=swappars)
	
	result <- list(
		incidence=incidence,
		pathcode=pathcode,
		upath=upath,
		height=height,
		maxgens=maxgens,
		call=match.call())
	class(result)<-"mimosa"

	return(result)
}
