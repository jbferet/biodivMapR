#' Reads ENVI hdr file
#'
#' @param hdr_path character. Path of the hdr file
#'
#' @return list of the content of the hdr file
#' @export

read_envi_header <- function(hdr_path) {
  # header <- paste(header, collapse = "\n")
  if (!grepl("\\.hdr$", hdr_path) & !grepl("\\.HDR$", hdr_path))
    stop("File extension should be .hdr or .HDR")

  hdr <- readLines(hdr_path)
  ## check ENVI at beginning of file
  if (!grepl("ENVI", hdr[1])) {
    stop("Not an ENVI header (ENVI keyword missing)")
  } else {
    hdr <- hdr [-1]
  }
  ## remove curly braces and put multi-line key-value-pairs into one line
  hdr <- gsub("\\{([^}]*)\\}", "\\1", hdr)
  l <- grep("\\{", hdr)
  r <- grep("\\}", hdr)
  if (length(l) != length(r))
    stop("Error matching curly braces in header (differing numbers).")
  if (any(r <= l))
    stop("Mismatch of curly braces in header.")
  hdr[l] <- sub("\\{", "", hdr[l])
  hdr[r] <- sub("\\}", "", hdr[r])
  for (i in rev(seq_along(l)))
    hdr <- c(hdr [seq_len(l [i] - 1)],
             paste(hdr [l [i]:r [i]], collapse = "\n"),
             hdr [-seq_len(r [i])])

  ## split key = value constructs into list with keys as names
  hdr <- sapply(hdr, split_line, "=", USE.NAMES = FALSE)
  names(hdr) <- tolower(names(hdr))
  ## process numeric values
  tmp <- names(hdr) %in% c(
    "samples", "lines", "bands", "header offset", "data type",
    "byte order", "default bands", "data ignore value",
    "wavelength", "fwhm", "data gain values", "band names"
  )
  hdr[tmp] <- lapply(hdr[tmp], function(x) {
    tryCatch({as.numeric(unlist(strsplit(x, ","))) },
             warning = function(w){ x })
  })
  return(hdr)
}


