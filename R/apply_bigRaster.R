#' This function performs big raster processing
#' It requires user to write a function expecting vectors / matrices / dataframes
#'  corresponding to raster data read sequentially, a list of input files, a list
#'  of output files, and it is good to go.
#'
#' @param funct function. function to be applied to process chunks of big rasters
#' @param input_rasters list. path for input files to be processed
#' @param input_args list. additional arguments sent to funct
#' @param output_rasters list. path for output file to be written
#' @param output_lyrs numeric. number of bands for output raster
#' @param filetype character. name of GDAL driver
#' @param bandNames list. list of band names corresponding to output files
#' @param unit_nrows numeric. define to read a number of lines multiple of unit_nrows
#' @param maxRows numeric. maximumnumber of rows to process at once
#'
#' @return output_rasters list. path for output files
#' @importFrom terra rast blocks readStart writeStart readValues writeValues writeStop
#' @import tools
#' @export

apply_bigRaster <- function(funct, input_rasters, input_args = NULL,
                            output_rasters, output_lyrs = NULL,
                            filetype = 'EHdr', bandNames = NULL,
                            unit_nrows = NULL, maxRows = NULL){

  # how many layers will be read & written decides on amount of data to be read
  nbLayers <- 0
  # initiate reading each input file
  if (is.null(names(input_rasters)))
    names(input_rasters) <- seq_len(length(input_rasters))
  r_in <- list()
  for (fid in names(input_rasters)){
    if (!is.null(input_rasters[[fid]])){
      if (inherits(input_rasters[[fid]], what = 'SpatRaster')){
        r_in[[fid]] <- input_rasters[[fid]]
      } else if (inherits(input_rasters[[fid]], what = 'list') &
                 all(unlist(lapply(input_rasters[[fid]], file.exists)))){
        r_in[[fid]] <- terra::rast(lapply(input_rasters[[fid]], terra::rast))
        input_rasters[[fid]] <- unlist(input_rasters[[fid]])
      } else if (inherits(input_rasters[[fid]], what = 'character') &
                 file.exists(input_rasters[[fid]])){
        r_in[[fid]] <- terra::rast(input_rasters[[fid]])
      }
      nbLayers <- nbLayers + dim(r_in[[fid]])[3]
    } else if (is.null(input_rasters[[fid]])){
      input_rasters[[fid]] <- NULL
    }
  }

  # specify number of bands for output rasters
  if (is.null(output_lyrs))
    output_lyrs <- dim(r_in[[fid]])[3]
  if (length(output_lyrs)==1)
    output_lyrs <- as.list(rep(output_lyrs,length(output_rasters)))
  names(output_lyrs) <- names(output_rasters)
  nbLayers <- nbLayers + sum(unlist(output_lyrs))

  # Adjust size of individual chunks
  brast <- terra::rast(input_rasters[[1]])
  blk <- terra::blocks(brast, n = nbLayers)
  # adjust number of lines to be read each time
  if (!is.null(maxRows)){
    if (blk$nrows[1]>maxRows){
      nrows_indiv <- maxRows
      blkud <- list()
      blkud$nrows <- rep(nrows_indiv,floor(sum(blk$nrows)/nrows_indiv))
      if (sum(blkud$nrows)<sum(blk$nrows)){
        blkud$nrows <- c(blkud$nrows, sum(blk$nrows)-sum(blkud$nrows))
      }
      blkud$row <- NULL
      blkud$n <- length(blkud$nrows)
      blkud$row <- c(1,cumsum(blkud$nrows)+1)[1:blkud$n]
      blk$row <- blkud$row
      blk$nrows <- blkud$nrows
      blk$n <- blkud$n
    }
  }

  if (!is.null(unit_nrows)){
    if (length(blk$nrows) >1){
      nrows_indiv <- blk$nrows[1] - (blk$nrows[1] %% unit_nrows)
      if (nrows_indiv==0) nrows_indiv <- unit_nrows
      blkud <- list()
      blkud$nrows <- rep(nrows_indiv,floor(sum(blk$nrows)/nrows_indiv))
      if (sum(blkud$nrows)<sum(blk$nrows)){
        blkud$nrows <- c(blkud$nrows, sum(blk$nrows)-sum(blkud$nrows))
      }
      blkud$row <- NULL
      blkud$n <- length(blkud$nrows)
      blkud$row <- c(1,cumsum(blkud$nrows)+1)[1:blkud$n]
      blk$row <- blkud$row
      blk$nrows <- blkud$nrows
      blk$n <- blkud$n
    }
  }
  for (fid in names(input_rasters))
    terra::readStart(r_in[[fid]])

  # initiate writing each output file
  if (is.null(names(output_rasters)))
    names(output_rasters) <- seq_len(length(output_rasters))
  r_out <- r_outbs <- list()
  for (fid in names(output_rasters)){
    # r_out[[fid]] <- terra::rast(input_rasters[[1]], nlyrs = min(output_lyrs[[fid]], dim(brast)[3]))
    # if (output_lyrs[[fid]] > dim(brast)[3]){
    #   dim(r_out[[fid]])[3] <- output_lyrs[[fid]]
    # }
    if (!is.null(input_args$datatype)){
      datatype <- input_args$datatype
    } else {
      datatype <- 'FLT4S'
    }
    r_out[[fid]] <- terra::rast(input_rasters[[1]])
    dim(r_out[[fid]])[3] <- output_lyrs[[fid]]
    r_outbs[[fid]] <- terra::writeStart(r_out[[fid]],
                                        filename = output_rasters[[fid]],
                                        overwrite = TRUE,
                                        filetype = filetype,
                                        datatype = datatype)
  }

  # init progressbar
  pgbarlength <- blk$n
  pb <- progress::progress_bar$new(
    format = "Processing raster [:bar] :percent in :elapsedfull , estimated time remaining :eta",
    total = pgbarlength, clear = FALSE, width= 100)
  # loop over blocks
  for (i in seq_along(blk$row)) {
    # read input files
    input_data <- list()
    for (fid in names(input_rasters)){
      input_data[[fid]] <- terra::readValues(r_in[[fid]], row = blk$row[i],
                                             nrows = blk$nrows[i],
                                             dataframe = TRUE)}
    # Apply Function
    out <- funct(input_data = input_data, input_args = input_args)
    # Write output files
    for (fid in names(output_rasters)){
      terra::writeValues(x = r_out[[fid]], v = out[[fid]],
                         blk$row[i], nrows = blk$nrows[i])}
    pb$tick()
  }
  # close output files
  # for (fid in names(output_rasters)){
  #   r_out[[fid]]  < - terra::writeStop(r_out[[fid]])
  # }
  for (fid in names(output_rasters)){
    r_out[[fid]]  < - terra::writeStop(r_out[[fid]])
    if (filetype %in% c('EHdr', 'ENVI'))
      harmonize_envi_hdr(output_raster = output_rasters[[fid]],
                         filetype = filetype,
                         bandNames = bandNames[[fid]])
  }
  rm(list=setdiff(ls(), "output_rasters"))
  gc()
  return(output_rasters)
}
