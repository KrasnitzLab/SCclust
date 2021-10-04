
#' Collects all bin count files from given directory
#' @param input_file_dir directory to scan for bin count files
#' @param suffix_pattern suffix to select files from input directory
#' @return data frame with filenames, cell names and file basenames.
#' @export 
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


#' @export
load_table <- function(filename) {
  assertthat::assert_that(file.exists(filename))

  df <- read.table(filename, header=T, as.is=T)
  return(df)
}

#' @export
load_cytobands <- function(filename) {
  assertthat::assert_that(file.exists(filename))
  
  return(read.table(cyto_file, header=F,as.is=T))
}

#' @export
save_table <- function(filename, df) {
  write.table(df, filename, col.names=T, row.names=F, sep="\t", quote=F)
}

#' @export
save_mat <- function(filename, mat) {
  write(mat, file=filename)
}

#' @export
load_mat <- function(filename) {
  return(as.matrix(scan(filename)))
}

#' @export
load_matrix <- function(filename) {
    m <- as.matrix(read.table(filename, sep="\t", header=F))
    return(m)
}

#' @export
uber_cells <- function(df, skip=3) {
  cells <- colnames(df)[c((skip+1):ncol(df))]
  return(data.frame(cells, stringsAsFactors=F))
}

#' @export
cells_filename <- function(output_dir, casename) {
  assertthat::assert_that(file.exists(output_dir))

  filename <- paste(casename, ".cells.txt")
  outname <- file.path(output_dir, filename)
  return(outname)
}

#' Constructs names for various output files based on `output_dir` and
#' `casename`
#' 
#' @export 
case_filenames <- function(output_dir, casename) {
  suffixes <- list(
    'ratio', 'clone', 'tree', 'seg',
    'featuremat', 'features', 'ploidies',
    'guide', 'genome', 'cells',
    'sim_pv', 'true_pv')

  filenames <- lapply(
    suffixes,
    function(s) file.path(output_dir, paste(casename, ".", s, ".txt", sep="")))
  names(filenames) <- suffixes
  return(filenames)
}
