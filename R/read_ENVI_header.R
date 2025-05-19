#' Reads ENVI hdr file
#'
#' @param HDRpath Path of the hdr file
#'
#' @return list of the content of the hdr file
#' @export

read_ENVI_header <- function(HDRpath) {
  # header <- paste(header, collapse = "\n")
  if (!grepl(".hdr$", HDRpath) & !grepl(".HDR$", HDRpath))
    stop("File extension should be .hdr or .HDR")

  HDR <- readLines(HDRpath)
  ## check ENVI at beginning of file
  if (!grepl("ENVI", HDR[1])) {
    stop("Not an ENVI header (ENVI keyword missing)")
  } else {
    HDR <- HDR [-1]
  }
  ## remove curly braces and put multi-line key-value-pairs into one line
  HDR <- gsub("\\{([^}]*)\\}", "\\1", HDR)
  l <- grep("\\{", HDR)
  r <- grep("\\}", HDR)
  if (length(l) != length(r))
    stop("Error matching curly braces in header (differing numbers).")
  if (any(r <= l))
    stop("Mismatch of curly braces in header.")
  HDR[l] <- sub("\\{", "", HDR[l])
  HDR[r] <- sub("\\}", "", HDR[r])
  for (i in rev(seq_along(l)))
    HDR <- c(HDR [seq_len(l [i] - 1)],
             paste(HDR [l [i]:r [i]], collapse = "\n"),
             HDR [-seq_len(r [i])])

  ## split key = value constructs into list with keys as names
  HDR <- sapply(HDR, split_line, "=", USE.NAMES = FALSE)
  names(HDR) <- tolower(names(HDR))
  ## process numeric values
  tmp <- names(HDR) %in% c(
    "samples", "lines", "bands", "header offset", "data type",
    "byte order", "default bands", "data ignore value",
    "wavelength", "fwhm", "data gain values"
  )
  HDR [tmp] <- lapply(HDR [tmp], function(x) {
    as.numeric(unlist(strsplit(x, ",")))
  })
  return(HDR)
}


