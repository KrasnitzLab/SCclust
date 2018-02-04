context("file utils work as expected")

test_that("varbin input file are constructed", {

  input_dir <- "/home/lubo/Work/sgains/data/bowtie-0.12.7/nyu003_R50_B20k_b0.12.7/varbin"
  suffix_pattern <- "\\.varbin.20k.txt$"

  files_list = varbin_input_files(input_dir, suffix_pattern)

  print(head(files_list))
  print(dim(files_list))

})
