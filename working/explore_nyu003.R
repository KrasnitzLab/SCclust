devtools::load_all()

output_dir <- "/home/lubo/Work/SCclust/working/breakPointPins/NYU003B9/scgv"
casename <- "nyu003"
filenames <- case_filenames(output_dir, casename)

uberfile <- file.path(output_dir, "uber.hg19.nyu003.benign.1.20k.lowratio.quantal.R.ratio.txt")
uber_df <- load_table(uberfile)
cells_df <- uber_cells(uber_df, skip=3)
dim(cells_df)

pinmatfile <- file.path(output_dir, "nyu003.benign.1smear1bpPinMat.featuremat.txt")
pinmat_df <- load_table(pinmatfile)
cells2_df <- uber_cells(pinmat_df, skip=0)
dim(cells2_df)

write.table(cells2_df, file=filenames$cells, row.names = F, col.names = T, quote = F)
