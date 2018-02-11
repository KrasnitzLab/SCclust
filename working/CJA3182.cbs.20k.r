devtools::load_all()

data_dir <- Sys.getenv("SGAINS_DATA")
assertthat::assert_that(file.exists(data_dir))

library("DNAcopy")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

remove.segment <- function( rsShort, rsSegnum, ratioData, sd.undo ) {

	appendLeft <- TRUE
	checkSdundo <- FALSE

	if (rsSegnum == 1) {
		appendLeft <- FALSE
	} else {
	if (rsSegnum == nrow(rsShort)) {
		appendLeft <- TRUE
	} else {
		rightIndex <- rsSegnum + 1
		leftIndex <- rsSegnum - 1

		if (rsShort[rightIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
			appendLeft <- TRUE
		} else {
		if (rsShort[leftIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
			appendLeft <- FALSE
		} else {
		if (abs(rsShort[leftIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"]) < abs(rsShort[rightIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"])) {
			appendLeft <- TRUE
			checkSdundo <- TRUE
		} else {
			appendLeft <- FALSE
			checkSdundo <- TRUE
		}}}
	}}

	appendIndex <- 99999999
	if (appendLeft) {
		appendIndex <- rsSegnum - 1
	} else {
		appendIndex <- rsSegnum + 1
	}

	tempShort <- rsShort
	newLocStart <- -1
	newLocEnd <- -1
	if (appendLeft) {
		tempShort[appendIndex, "loc.end"] <- tempShort[rsSegnum, "loc.end"]
		tempShort[appendIndex, "seg.end"] <- tempShort[rsSegnum, "seg.end"]
	} else {
		tempShort[appendIndex, "loc.start"] <- tempShort[rsSegnum, "loc.start"]
		tempShort[appendIndex, "seg.start"] <- tempShort[rsSegnum, "seg.start"]
	}

	tempShort[appendIndex, "num.mark"] <- tempShort[appendIndex, "num.mark"] + tempShort[rsSegnum, "num.mark"]
	tempShort[appendIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[appendIndex, "seg.start"]:tempShort[appendIndex, "seg.end"]], base=2))

	cat("append", tempShort[appendIndex, "chrom"], tempShort[appendIndex, "loc.start"], tempShort[appendIndex, "loc.end"], tempShort[appendIndex, "num.mark"], tempShort[appendIndex, "seg.mean"], tempShort[appendIndex, "seg.start"], tempShort[appendIndex, "seg.end"], "\n")

	tempShort <- tempShort[-rsSegnum, ]
	tempShort$segnum <- seq(1:nrow(tempShort))

	if (checkSdundo) {
		thisSd <- -1
		if (appendLeft) {
			leftIndex <- appendIndex
			rightIndex <- appendIndex + 1
		} else {
			leftIndex <- appendIndex - 2
			rightIndex <- appendIndex - 1
		}
		#thisSd <- sd(ratioData[tempShort$seg.start[leftIndex]:tempShort$seg.start[rightIndex], "lowratio"])
		thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)

		if (abs(tempShort$seg.mean[leftIndex] - tempShort$seg.mean[rightIndex]) < (sd.undo * thisSd) ) {

			cat("left", tempShort[leftIndex, "chrom"], tempShort[leftIndex, "loc.start"], tempShort[leftIndex, "loc.end"], tempShort[leftIndex, "num.mark"], tempShort[leftIndex, "seg.mean"], tempShort[leftIndex, "seg.start"], tempShort[leftIndex, "seg.end"], "\n")
			cat("right", tempShort[rightIndex, "chrom"], tempShort[rightIndex, "loc.start"], tempShort[rightIndex, "loc.end"], tempShort[rightIndex, "num.mark"], tempShort[rightIndex, "seg.mean"], tempShort[rightIndex, "seg.start"], tempShort[rightIndex, "seg.end"], "\n")

			##  remove breakpoint
			tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
			tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
			tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
			tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base=2))
			tempShort <- tempShort[-rightIndex, ]
			tempShort$segnum <- seq(1:nrow(tempShort))
		}
	}

	return(tempShort)
}


sdundo.all <- function (sdShort, ratioData, sd.undo) {

	tempShort <- sdShort
	thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)

	while ( TRUE ) {

		chrom <- tempShort$chrom
		chrom.shift <- c(tempShort$chrom[-1], tempShort$chrom[1])

		breakpoints <- which(chrom == chrom.shift)
		cat("sdundo.all intrachrom breakpoints", length(breakpoints), "\n")

		if (length(breakpoints) < 1) {
			break
		}

		breakpoints.shift <- breakpoints + 1

		undo.breakpoints <- breakpoints[which(abs(tempShort$seg.mean[breakpoints] - tempShort$seg.mean[breakpoints.shift]) < thisSd * sd.undo)]

		cat("sdundo.all undo breakpoints", length(undo.breakpoints), "\n")

		if (length(undo.breakpoints) < 1) {
			break
		}

		undo.breakpoints.shift <- undo.breakpoints + 1

		undo.df <- tempShort[undo.breakpoints, ]
		undo.df$seg.mean.diff <- abs(tempShort$seg.mean[undo.breakpoints] - tempShort$seg.mean[undo.breakpoints.shift])

		min.index <- which.min(undo.df$seg.mean.diff)

		leftIndex <- undo.df$segnum[min.index]
		rightIndex <- leftIndex + 1

		cat("sdundo.all left", tempShort[leftIndex, "chrom"], tempShort[leftIndex, "loc.start"], tempShort[leftIndex, "loc.end"], tempShort[leftIndex, "num.mark"], tempShort[leftIndex, "seg.mean"], tempShort[leftIndex, "seg.start"], tempShort[leftIndex, "seg.end"], "\n")
		cat("sdundo.all right", tempShort[rightIndex, "chrom"], tempShort[rightIndex, "loc.start"], tempShort[rightIndex, "loc.end"], tempShort[rightIndex, "num.mark"], tempShort[rightIndex, "seg.mean"], tempShort[rightIndex, "seg.start"], tempShort[rightIndex, "seg.end"], "\n")

		tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
		tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
		tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
		tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base=2))
		tempShort <- tempShort[-rightIndex, ]
		tempShort$segnum <- seq(1:nrow(tempShort))

	}

	return(tempShort)

}


cbs.segment01 <- function(indir, outdir, varbin.gc, varbin.data, sample.name, alt.sample.name, alpha, nperm, undo.SD, min.width) {
	gc <- read.table(varbin.gc, header=T)

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisRatio <- read.table(paste(indir, varbin.data, sep="/"), header=F)
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
	a <- thisRatio$bincount + 1
	thisRatio$ratio <- a / mean(a)
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)
	#thisRatio$log2lowratio <- log(thisRatio$lowratio, base=2)

	set.seed(25)
	CNA.object <- CNA(log(thisRatio$lowratio, base=2), thisRatio$chrom, thisRatio$chrompos, data.type="logratio", sampleid=sample.name)
	smoothed.CNA.object <- smooth.CNA(CNA.object)
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2)
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	cbs.long <- m[, 1]

	#####  NEW STUFF  also check min.width=2 above

	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.varbin.short.cbs.txt", sep=""), quote=F, row.names=F)

	workShort <- thisShort
	workShort$segnum <- 0
	workShort$seg.start <- 0
	workShort$seg.end <- 0
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		workShort$seg.start[i] <- thisStart
		workShort$seg.end[i] <- thisEnd
		workShort$segnum[i] <- i
		prevEnd = thisEnd
	}

	discardSegments <- TRUE
	while (discardSegments) {
		orderShort <- workShort[order(workShort$num.mark, abs(workShort$seg.mean)), ]
		if (orderShort[1, "num.mark"] < min.width) {
			workShort <- remove.segment(workShort, orderShort[1, "segnum"], thisRatio, undo.SD)
		} else {
			discardSegments <- FALSE
		}
	}

	workShort <- sdundo.all(workShort, thisRatio, undo.SD)
	thisShort <- workShort

	#####  END NEW STUFF

	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	thisRatio$seg.mean.LOWESS <- m[, 1]

	thisGrid <- seq(1.5, 5.5, by=0.05)
	thisOuter <- thisRatio$seg.mean.LOWESS %o% thisGrid
	thisOuterRound <- round(thisOuter)
	thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
	thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
	thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
	thisError <- min(thisOuterColsums)
	thisShredded <- length(which(thisRatio$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

	thisRatio$ratio.quantal <- thisRatio$lowratio * thisMultiplier
	thisRatio$seg.quantal <- thisRatio$seg.mean.LOWESS * thisMultiplier

	thisRatio$cbs.seg <- cbs.long
	thisRatio$cbs.seg.quantal <- cbs.long * thisMultiplier

	thisQuantalStats <- data.frame("ploidy"=thisMultiplier, "error"=thisError, "shredded"=thisShredded)

	write.table(thisQuantalStats, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.varbin.quantal.stats.txt", sep=""), quote=F, row.names=F)
	write.table(thisRatio, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.varbin.data.txt", sep=""), quote=F, row.names=F)
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.varbin.short.txt", sep=""), quote=F, row.names=F)

	return(thisRatio)
}

indir = file.path(data_dir, "nyu003/CJA3182/processed")
assertthat::assert_that(file.exists(indir))
outdir = "./out/nyu003/"
assertthat::assert_that(file.exists(indir))
varbin_gc_file = file.path(data_dir, "varbin_orig/varbin.gc.content.20k.bowtie.k50.hg19.txt.gz")
assertthat::assert_that(file.exists(varbin_gc_file))

ratio_data <-cbs.segment01(
  indir=indir,
  outdir=outdir,
  varbin.gc=varbin_gc_file,
  varbin.data="CJA3182.varbin.20k.txt",
  sample.name="CJA3182",
  alt.sample.name="FC64BEMAAXX lane 8 FA014 A1 bc1 NYU_003_5_PBXW0032 2C",
  alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5)

gc <- load_table(varbin_gc_file)
varbin_file <- file.path(indir, "CJA3182.varbin.20k.txt")
assertthat::assert_that(file.exists(varbin_file))

varbin_data <- read.table(varbin_file, header=F, as.is=T)
names(varbin_data) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
varbin_data$gc.content <- gc$gc.content

out_data <- cbs.segment_1(varbin_data, gc, alpha=0.05, nperm=1000, undo.SD=1.0, min.width = 5)
head(out_data)
head(ratio_data)

assertthat::assert_that(all(ratio_data$ratio.quantal - out_data$ratio.quantal < 1e-9))
assertthat::assert_that(all(ratio_data$seg.quantal - out_data$seg.quantal <1e-9))

cytofile<-"./annot/cytoBandHG19.txt"
centromere<-c("p11","q11")
chromrange<-1:22

cyto<-read.table(cytofile,header=F,as.is=T)
cyto[cyto[,1]=="chrX",1]<-"chr23"
cyto[cyto[,1]=="chrY",1]<-"chr24"
cyto[,1]<-as.numeric(substring(cyto[,1],first=4,last=nchar(cyto[,1])))
cyto<-cyto[order(cyto[,1]),]
centroleft<-cyto[grep(centromere[1],cyto[,4]),]
centroright<-cyto[grep(centromere[2],cyto[,4]),]
centroleft<-centroleft[match(unique(centroleft[,1]),centroleft[,1]),]
centroright<-centroright[nrow(centroright):1,]
centroright<-centroright[match(unique(centroright[,1]),centroright[,1]),]
centroright<-centroright[nrow(centroright):1,]
dropareas<-cbind(centroleft[,c(1,2)],centroright[,3])
dimnames(dropareas)[[2]]<-c("chrom","from","to")

varbin_data$chrom <- numeric_chrom(varbin_data$chrom)
head(varbin_data)
head(dropareas)
