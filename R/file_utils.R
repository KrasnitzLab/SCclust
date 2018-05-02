

varbin_input_files <- function(input_file_dir, suffix_pattern="") {
  assertthat::assert_that(file.exists(input_file_dir))

  # pattern.file <- paste("\\.varbin.", Nk, ".txt$", sep = "")
  names <- list.files(path = input_file_dir, pattern = suffix_pattern)
  filenames <- lapply(names, function(fn) file.path(input_file_dir, fn))
  cells <- sub("\\..*","", lapply(names, basename))

  names <- c(unlist(names))
  filenames <- c(unlist(filenames))
  cells <- c(unlist(cells))

  return(data.frame(names, cells, paths=filenames, stringsAsFactors=FALSE))
}


load_table <- function(filename) {
  flog.debug("loading table: %s", filename)
  assertthat::assert_that(file.exists(filename))

  df <- read.table(filename, header=T, as.is=T)
  return(df)
}

load_cytobands <- function(filename) {
  flog.debug("loading cytobands: %s", filename)
  assertthat::assert_that(file.exists(filename))
  
  return(read.table(cyto_file, header=F,as.is=T))
}

save_table <- function(filename, df) {
  flog.debug("writing table: %s", filename)
  write.table(df, filename, col.names=T, row.names=F, sep="\t", quote=F)
}

save_mat <- function(filename, mat) {
  flog.debug("writing matrix: %s", filename)
  write(mat, file=filename)
}

load_mat <- function(filename) {
  flog.debug("loading matrix: %s", filename)
  return(as.matrix(scan(filename)))
}

uber_cells <- function(df, skip=3) {
  cells <- colnames(df)[c((skip+1):ncol(df))]
  return(data.frame(cells, stringsAsFactors=F))
}

cells_filename <- function(output_dir, casename) {
  assertthat::assert_that(file.exists(output_dir))

  filename <- paste(casename, ".cells.txt")
  outname <- file.path(output_dir, filename)
  return(outname)
}

case_filenames <- function(output_dir, casename) {
  suffixes <- list(
    'ratio', 'clone', 'tree', 'seg',
    'featuremat', 'features',
    'guide', 'genome', 'cells',
    'sim_pv', 'true_pv')

  filenames <- lapply(
    suffixes,
    function(s) file.path(output_dir, paste(casename, ".", s, ".txt", sep="")))
  names(filenames) <- suffixes
  return(filenames)
}
