# Load helper functions ---------------------------------------------------

get_data <- function(path) {
  if(!is.null(CSV_ROOT_DIR)) {
    paste0(CSV_ROOT_DIR, path)
  } else {
    stop("Must specify location of CSV directory!")
  }
}

specd <- function(x, k) trimws(format(round(x, k), nsmall=k))