

varbin_input_files <- function(input_file_dir, suffix_pattern="") {
  assertthat::assert_that(file.exists(input_file_dir))

  # pattern.file <- paste("\\.varbin.", Nk, ".txt$", sep = "")
  names <- list.files(path = input_file_dir, pattern = suffix_pattern)
  paths <- lapply(names, function(fn) file.path(input_file_dir, fn))
  cells <- sub("\\..*","", lapply(names, basename))

  names <- c(unlist(names))
  paths <- c(unlist(paths))
  cells <- c(unlist(cells))

  files <- data.frame(names, cells, paths)
  return(files)
}


load_table <- function(filename) {
  assertthat::assert_that(file.exists(filename))

  df <- read.table(filename, header=T, as.is=T)
  return(df)
}

uber_cells <- function(df, skip=3) {
  cell <- colnames(df)[c((skip+1):ncol(df))]
  return(data.frame(cell))
}

cells_filename <- function(output_dir, casename) {
  assertthat::assert_that(file.exists(output_dir))

  filename <- paste(casename, ".cells.txt")
  outname <- file.path(output_dir, filename)
  return(outname)
}

ratio_filename <- function(output_dir, casename) {

}

case_filenames <- function(output_dir, casename) {
  suffixes <- list(
    'ratio', 'clone', 'tree', 'seg',
    'featuremat', 'features',
    'guide', 'genome', 'cells')

  filenames <- lapply(
    suffixes,
    function(s) file.path(output_dir, paste(casename, ".", s, ".txt", sep="")))
  names(filenames) <- suffixes
  return(filenames)
}
