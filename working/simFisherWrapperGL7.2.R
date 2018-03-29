# Please edit the filenames so you can read and save
# prostate <- "nyu003.benign.1"

# to compare
# diff breakPointPins/Gleason6.1/averageNoEndsNewest/GL6.1simP.txt out/Gleason6.1/GL6.1simP.txt
devtools::load_all()

data_dir <- Sys.getenv("SGAINS_DATA")
assertthat::assert_that(file.exists(data_dir))

prostate <- "GL7.2"
savedir <- file.path(
    data_dir,
    "gleason7.2/out_many/")

rdir <- file.path(
    data_dir,
    "gleason7.2/results")
assertthat::assert_that(file.exists(savedir))
assertthat::assert_that(file.exists(rdir))

pinmat_filename <- file.path(
    rdir,
    "nyu007.GL7.2smear1bpPinMat.txt")
pins_filename <- file.path(
    rdir,
    "nyu007.GL7.2smear1bpPins.txt")
assertthat::assert_that(file.exists(pinmat_filename))
assertthat::assert_that(file.exists(pins_filename))

truefile <- paste(savedir,prostate,"trueP.txt",sep="")
simfile <- paste(savedir,prostate,"simP.txt",sep="")

print(truefile)
print(simfile)

source("./simFisher.R")
pinmat <- as.matrix(read.table(pinmat_filename,header=T))
pins<-as.matrix(read.table(pins_filename,header=T))
# eviltwins <- c("CJA1593","CJA2772","CJA3438","CJA0951")
eviltwins <- c("CJA1593","CJA2772","CJA0951")

pinmat <- pinmat[,setdiff(dimnames(pinmat)[[2]],eviltwins),drop=F]
m<-vector(mode="list",length=length(unique(pins[,"sign"])))
for(i in 1:length(unique(pins[,"sign"])))m[[i]]<-pinmat[pins[,"sign"]==unique(pins[,"sign"])[i],,drop=F]
require(parallel)
usecores <- 20 # detectCores()-4
vtrue <- simFisher(
    m,nsim=1,nsweep=0,seedme=123,
    distrib="Rparallel",njobs=usecores,combo="Fisher", savedir=savedir)
write(vtrue,file=truefile)
msim <- simFisher(
    m,nsim=500,nsweep=200,seedme=123,
    distrib="Rparallel",njobs=usecores,combo="Fisher", savedir=savedir)
write(msim,file=simfile)
# quit(save="no")

