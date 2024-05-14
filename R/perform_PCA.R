#' Performs PCA for all images and create PCA file with either all or a selection of PCs
#'
#' @param input_raster_path character. path for image to be processed
#' @param output_dir character. Path for output directory
#' @param input_rast_wl numeric. spectral bands corresponding to input_raster_path
#' @param input_mask_path character. path for mask corresponding to the image
#' @param Continuum_Removal boolean. Set to TRUE if continuum removal should be applied
#' @param TypePCA character. Type of PCA: choose either "PCA" or "SPCA"
#' @param NbPCs_To_Keep numeric. number of components to ke saved in the PCA file. default = 30 if set to FALSE (or nb PC if <30)
#' @param Excluded_WL numeric. Water Vapor Absorption domains (in nanometers, min and max WL).
#' Can also be used to exclude spectific domains.
#' dims = N x 2 (N = number of domains to be eliminated)
#' @param nbPix_PCA numeric. number of pixels to use to compute PCA
#' @param nbIter numeric. nb of iterations averaged to compute diversity indices
#' @param maxRows numeric. max number of rows in each block
#' @param filetype character. gdal driver
#'
#' @return list of paths corresponding to resulting PCA files
#' @importFrom bigRaster apply_bigRaster
#' @importFrom terra rast
#' @export

perform_PCA  <- function(input_raster_path, output_dir, input_rast_wl = NULL,
                         input_mask_path = NULL, Continuum_Removal = TRUE,
                         TypePCA = 'SPCA', NbPCs_To_Keep = 30,
                         Excluded_WL = NULL, nbPix_PCA = 1e6,
                         nbIter = 20, maxRows = 100, filetype = 'COG') {
  # check if format of raster data is as expected
  input_rast <- terra::rast(input_raster_path)
  if (!is.null(input_rast_wl))
  {names(input_rast) <- input_rast_wl
  check_data(input_data = input_rast, arguments = 'input_rast')
  }
  input_mask <- NULL
  if (!is.null(input_mask_path)) {
    input_mask <- terra::rast(input_mask_path)
    check_data(input_data = input_mask, arguments = 'input_mask')
  }
  # Identify water vapor absorption bands in image and possibly other spectral domains to discard
  SpectralFilter <- exclude_spectral_domains(input_rast = input_rast,
                                             Excluded_WL = Excluded_WL)
  # Extract valid data subset and check validity
  print("Extract pixels & perform PCA")
  # define number of pixels to be extracted from the image for each iteration
  Pix_Per_Partition <- define_pixels_per_iter(input_rast = input_rast,
                                              input_mask = input_mask,
                                              nbPix = nbPix_PCA,
                                              nbIter = nbIter)
  nbSamples <- nbIter * Pix_Per_Partition
  # extract a random selection of pixels from image
  if (!TypePCA=='MNF'){
    extent_area <- get_raster_extent(input_rast[[1]])
    Subset <- sample_from_raster(extent_area = extent_area,
                                 nbSamples = nbSamples,
                                 input_rast = input_rast,
                                 input_mask = input_mask)
    # Subset <- sample_exact_raster(extent_area = extent_area,
    #                               nbSamples = nbSamples,
    #                               input_rast = input_rast,
    #                               input_mask = input_mask)
    Subset$ID <- NULL
  }
  # if needed, apply continuum removal
  if (Continuum_Removal == TRUE) {
    Subset <- apply_continuum_removal(Spectral_Data = Subset,
                                      Spectral = SpectralFilter)
  } else {
    if (length(SpectralFilter$WaterVapor)>0) Subset <- Subset[, -SpectralFilter$WaterVapor]
  }
  # if number of pixels available inferior number initial sample size
  if (nrow(Subset) < nbSamples) {
    nbSamples <- nrow(Subset)
    nbIter <- ceiling(nbSamples / Pix_Per_Partition)
    Pix_Per_Partition <- floor(nbSamples / nbIter)
    nbSamples <- nbIter * Pix_Per_Partition
  }
  # clean reflectance data from inf and constant values
  CleanData <- rm_invariant_bands(Subset, SpectralFilter)
  Subset <- CleanData$DataMatrix
  SpectralFilter <- CleanData$Spectral
  # Compute PCA #1 on Subset
  print(paste('perform',TypePCA,'on image subset'))
  if (TypePCA == "PCA" | TypePCA == "SPCA") {
    PCA_model <- pca(Subset, TypePCA)
    # } else if(TypePCA=="MNF"){
    #   PCA_model <- mnf(Subset, coordPix)
  }
  # Number of PCs computed and written in the PCA file: 30 if hyperspectral
  Nb_PCs <- dim(PCA_model$x)[2]
  if (Nb_PCs > NbPCs_To_Keep) Nb_PCs <- NbPCs_To_Keep
  PCA_model$Nb_PCs <- Nb_PCs
  # CREATE PCA FILE CONTAINING ONLY SELECTED PCs
  PCA_Files <- list('PCA'= file.path(output_dir,
                                     paste0('OutputPCA_', Nb_PCs, '_PCs')))
  if (filetype %in% c('COG', 'tif', 'geoTiff', 'GeoTIFF', 'GTiff')){
    PCA_Files$PCA <- paste0(PCA_Files$PCA,'.tif')
  }
  funct <- wrapperBig_PCA
  input_PCA <- list('main' = input_raster_path,
                    'mask' = input_mask_path)
  input_args <- list('CR' = Continuum_Removal, 'Spectral' = CleanData$Spectral,
                     'PCA_model' = PCA_model, 'Nb_PCs' = Nb_PCs)
  bigRaster::apply_bigRaster(funct = funct, filetype = filetype,
                             input_rasters = input_PCA,
                             input_args = input_args,
                             output_lyrs = Nb_PCs,
                             bandNames = list('PCA' = paste0('PC#',seq_len(Nb_PCs))),
                             output_rasters = PCA_Files,
                             maxRows = maxRows)
  # save workspace for this stage
  WS_Save <- file.path(output_dir, "PCA_info.RData")
  my_list <- list("PCA_Files" = PCA_Files, "Pix_Per_Partition" = Pix_Per_Partition,
                  "nbIter" = nbIter, "MaskPath" = input_mask_path,
                  "PCA_model" = PCA_model, "SpectralFilter" = SpectralFilter,
                  "TypePCA" = TypePCA)
  save(PCA_Files, Pix_Per_Partition, nbIter, input_mask_path,
       PCA_model, SpectralFilter, TypePCA, file = WS_Save)
  return(my_list)
}
