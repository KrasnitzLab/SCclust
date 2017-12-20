library("DNAcopy", lib.loc="/mnt/wigclust5/data/safe/kendall/DNAcopy_1.26.0")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
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
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

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
	
	thisQuantalStats <- data.frame("ploidy"=thisMultiplier, "error"=thisError, "shredded"=thisShredded)
	
	chr <- thisRatio$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

	png(paste(outdir, "/", sample.name, ".20k.wg.png", sep=""), height=800, width=1200)
	plot(x=thisRatio$abspos, y=thisRatio$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatio$abspos, y=thisRatio$lowratio, col="#CCCCCC")
	points(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA")
	lines(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	hlines <- c(1, 2, 3, 4, 5, 6)

	png(paste(outdir, "/", sample.name, ".20k.wg.quantal.png", sep=""), height=800, width=1200)
	plot(x=thisRatio$abspos, y=thisRatio$ratio.quantal, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatio$abspos, y=thisRatio$ratio.quantal, col="#CCCCCC")
	points(x=thisRatio$abspos, y=thisRatio$seg.quantal, col="#0000AA")
	lines(x=thisRatio$abspos, y=thisRatio$seg.quantal, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()

	write.table(thisQuantalStats, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.varbin.quantal.stats.txt", sep=""), quote=F, row.names=F) 
	write.table(thisRatio, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.varbin.short.txt", sep=""), quote=F, row.names=F) 


	bad <- c(845, 846, 847, 848, 849, 850, 2208, 2209, 2210, 2211, 2212, 5021, 5022, 5023, 5024, 5025, 11642, 11643, 11644, 11645, 11646, 16232, 16233, 16234, 16235, 16236, 17781, 17782, 17783, 17801, 18420, 18421, 18422, 18423, 19939, 19940, 19941, 19942, 19943)
	thisRatioNobig <- thisRatio[-bad, ]

	set.seed(25) 
	CNA.object <- CNA(log(thisRatioNobig$lowratio, base=2), thisRatioNobig$chrom, thisRatioNobig$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatioNobig), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatioNobig$seg.mean.LOWESS <- m[, 1]

	thisGrid <- seq(1.5, 5.5, by=0.05)
	thisOuter <- thisRatioNobig$seg.mean.LOWESS %o% thisGrid
	thisOuterRound <- round(thisOuter)
	thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
	thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
	thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
	thisError <- min(thisOuterColsums)
	thisShredded <- length(which(thisRatioNobig$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

	thisRatioNobig$ratio.quantal <- thisRatioNobig$lowratio * thisMultiplier
	thisRatioNobig$seg.quantal <- thisRatioNobig$seg.mean.LOWESS * thisMultiplier
	
	thisQuantalStatsNobig <- data.frame("ploidy"=thisMultiplier, "error"=thisError, "shredded"=thisShredded)
	
	chr <- thisRatioNobig$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	# Use the vlines abspos positions from above.  Because these start at
	# the third bin on the acrocentric chromosomes the vlines end up to
	# the right of the centromere rather than the left which is wrong.
	#vlines <- c(1, thisRatioNobig$abspos[which(chr != chr.shift) + 1], thisRatioNobig$abspos[nrow(thisRatioNobig)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

	png(paste(outdir, "/", sample.name, ".20k.wg.nobad.png", sep=""), height=800, width=1200)
	plot(x=thisRatioNobig$abspos, y=thisRatioNobig$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC", cex=1)
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$lowratio, col="#CCCCCC", lwd=1.0)
	points(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.mean.LOWESS, col="#0000AA", cex=1)
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.mean.LOWESS, col="#0000AA", lwd=1.0)
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	hlines <- c(1, 2, 3, 4, 5, 6)

	png(paste(outdir, "/", sample.name, ".20k.wg.nobad.quantal.png", sep=""), height=800, width=1200)
	plot(x=thisRatioNobig$abspos, y=thisRatioNobig$ratio.quantal, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$ratio.quantal, col="#CCCCCC")
	points(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.quantal, col="#0000AA")
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.quantal, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()

	write.table(thisQuantalStatsNobig, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.nobad.varbin.quantal.stats.txt", sep=""), quote=F, row.names=F) 
	write.table(thisRatioNobig, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.nobad.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.20k.k50.nobad.varbin.short.txt", sep=""), quote=F, row.names=F) 

}

