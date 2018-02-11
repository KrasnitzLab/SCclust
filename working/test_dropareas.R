devtools::load_all()


data_dir <- Sys.getenv("SGAINS_DATA")
assertthat::assert_that(file.exists(data_dir))

nyu003_scgv_dir <- file.path(data_dir, "nyu003/scgv")
segment_df <- load_table(file.path(nyu003_scgv_dir, "uber.hg19.nyu003.benign.1.20k.seg.quantal.R.seg.txt"))
assertthat::are_equal(nrow(segment_df), 19943)

a<-round(as.matrix(segment_df[,-(1:3)]))
head(a)
dim(a)
dim(segment_df)
head(segment_df[,(1:3)])
b<-a-rbind(matrix(nrow=1,ncol=ncol(a),data=0),a[-nrow(a),])

head(rbind(matrix(nrow=1,ncol=ncol(a),data=0),a[-nrow(a),]))
