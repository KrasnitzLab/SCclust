devtools::load_all()

output_dir <- "/home/lubo/Work/SCclust/working/breakPointPins/NYU011Gl7.5/scgv"
casename <- "nyu011"
filenames <- case_filenames(output_dir, casename)

uberfile <- file.path(output_dir, "uber.hg19.nyu011.GL7.5.20k.seg.quantal.R.seg.txt")
uber_df <- load_table(uberfile)
cells_df <- uber_cells(uber_df, skip=3)
dim(cells_df)

pinmatfile <- file.path(output_dir, "nyu011.GL7.5smear1bpPinMat.featuremat.txt")
pinmat_df <- load_table(pinmatfile)
cells2_df <- uber_cells(pinmat_df, skip=0)
dim(cells2_df)

filenames$cells
write.table(cells2_df, file=filenames$cells, row.names = F, col.names = T, quote = F)
