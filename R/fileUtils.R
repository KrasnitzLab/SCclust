

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
