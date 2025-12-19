#' Performs PCA for all images and create PCA file with either all or a selection of PCs
#'
#' @param input_raster_path character. path for image to be processed
#' @param output_dir character. Path for output directory
#' @param input_rast_wl numeric. spectral bands corresponding to input_raster_path
#' @param input_mask_path character. path for mask corresponding to the image
#' @param Continuum_Removal boolean. Set to TRUE if continuum removal should be applied
#' @param TypePCA character. Type of PCA: choose either "PCA" or "SPCA"
#' @param nb_pcs_to_keep numeric. number of components to ke saved in the PCA file. default = 30 if set to FALSE (or nb PC if <30)
#' @param Excluded_WL numeric. Water Vapor Absorption domains (in nanometers, min and max WL).
#' Can also be used to exclude spectific domains.
#' dims = N x 2 (N = number of domains to be eliminated)
#' @param nb_pix_pca numeric. number of pixels to use to compute PCA
#' @param nb_iter numeric. nb of iterations averaged to compute diversity indices
#' @param maxRows numeric. max number of rows in each block
#' @param filetype character. gdal driver
#'
#' @return list of paths corresponding to resulting PCA files
#' @importFrom bigRaster apply_bigRaster
#' @importFrom terra rast
#' @export

perform_PCA  <- function(input_raster_path, output_dir, input_rast_wl = NULL,
                         input_mask_path = NULL, Continuum_Removal = TRUE,
                         TypePCA = 'SPCA', nb_pcs_to_keep = 30,
                         Excluded_WL = NULL, nb_pix_pca = 1e6,
                         nb_iter = 10, maxRows = 100, filetype = 'GTiff') {
  # check if format of raster data is as expected
  input_rast <- terra::rast(input_raster_path)
  if (!is.null(input_rast_wl)){
    names(input_rast) <- input_rast_wl
    check_data(input_data = input_rast, arguments = 'input_rast')
  }
  input_mask <- NULL
  if (!is.null(input_mask_path)) {
    input_mask <- terra::rast(input_mask_path)
    check_data(input_data = input_mask, arguments = 'input_mask')
  }
  if (!is.null(input_rast_wl)){
    # Identify water vapor absorption bands in image and possibly other spectral domains to discard
    spectral_filter <- exclude_spectral_domains(input_rast = input_rast,
                                                Excluded_WL = Excluded_WL)
    # Extract valid data subset and check validity
    print("Extract pixels & perform PCA")
    # define number of pixels to be extracted from the image for each iteration
    pix_per_partition <- define_pixels_per_iter(input_rast = input_rast,
                                                input_mask = input_mask,
                                                nb_pix = nb_pix_pca,
                                                nb_iter = nb_iter)
    nb_samples <- nb_iter * pix_per_partition
  } else {
    nb_samples <- nb_pix_pca
    spectral_filter <- nb_iter <- pix_per_partition <- NULL
    Continuum_Removal <- FALSE
  }
  # extract a random selection of pixels from image
  if (!TypePCA=='MNF'){
    extent_area <- get_raster_extent(input_rast[[1]])
    img_subset <- sample_from_raster(extent_area = extent_area,
                                     nb_samples = nb_samples,
                                     input_rast = input_rast,
                                     input_mask = input_mask)
    img_subset$ID <- NULL
  }
  # if needed, apply continuum removal
  if (Continuum_Removal == TRUE) {
    img_subset <- apply_continuum_removal(spectral_data = img_subset,
                                          spectral = spectral_filter)
  } else {
    if (length(spectral_filter$WaterVapor)>0)
      img_subset <- img_subset[, -spectral_filter$WaterVapor]
  }
  # if number of pixels available inferior number initial sample size
  if (nrow(img_subset) < nb_samples) {
    nb_samples <- nrow(img_subset)
    if (!is.null(pix_per_partition)){
      nb_iter <- ceiling(nb_samples / pix_per_partition)
      pix_per_partition <- floor(nb_samples / nb_iter)
      nb_samples <- nb_iter * pix_per_partition
    }
  }
  if (!is.null(spectral_filter)){
    # clean reflectance data from inf and constant values
    CleanData <- rm_invariant_bands(img_subset, spectral_filter)
    img_subset <- CleanData$DataMatrix
    spectral_filter <- CleanData$Spectral
  } else {
    CleanData <- list()
    CleanData$Spectral <- list('Bands2Keep'= seq_len(dim(img_subset)[2]))
  }
  # Compute PCA #1 on img_subset
  print(paste('perform',TypePCA,'on image subset'))
  if (TypePCA == "PCA" | TypePCA == "SPCA") {
    PCA_model <- pca(img_subset, TypePCA)
    # } else if(TypePCA=="MNF"){
    #   PCA_model <- mnf(img_subset, coordPix)
  }
  # Number of PCs computed and written in the PCA file: 30 if hyperspectral
  Nb_PCs <- dim(PCA_model$x)[2]
  if (Nb_PCs > nb_pcs_to_keep)
    Nb_PCs <- nb_pcs_to_keep
  PCA_model$Nb_PCs <- Nb_PCs
  # CREATE PCA FILE CONTAINING ONLY SELECTED PCs
  PCA_Files <- list('PCA'= file.path(output_dir,
                                     paste0('OutputPCA_', Nb_PCs, '_PCs')))
  if (filetype %in% c('COG', 'tif', 'geoTiff', 'GeoTIFF', 'GTiff'))
    PCA_Files$PCA <- paste0(PCA_Files$PCA,'.tiff')
  funct <- wrapperBig_PCA

  if (length(input_raster_path)>1){
    mainlist <- as.list(input_raster_path)
    names(mainlist) <- file_path_sans_ext(basename(input_raster_path))
  } else {
    mainlist <- input_raster_path
  }
  input_PCA <- list('main' = mainlist,
                    'mask' = input_mask_path)
  input_args <- list('CR' = Continuum_Removal, 'Spectral' = CleanData$Spectral,
                     'PCA_model' = PCA_model, 'Nb_PCs' = Nb_PCs)
  bigRaster::apply_bigRaster(funct = funct, filetype = filetype,
                             input_rasters = input_PCA,
                             input_args = input_args,
                             output_lyrs = Nb_PCs,
                             bandNames = list('PCA' = paste0('PC#',
                                                             seq_len(Nb_PCs))),
                             output_rasters = PCA_Files,
                             maxRows = maxRows)
  # save workspace for this stage
  WS_Save <- file.path(output_dir, "PCA_info.RData")
  my_list <- list("PCA_Files" = PCA_Files,
                  "pix_per_partition" = pix_per_partition,
                  "nb_iter" = nb_iter, "MaskPath" = input_mask_path,
                  "PCA_model" = PCA_model,
                  "spectral_filter" = spectral_filter,
                  "TypePCA" = TypePCA)
  save(PCA_Files, pix_per_partition, nb_iter, input_mask_path,
       PCA_model, spectral_filter, TypePCA, file = WS_Save)
  return(my_list)
}
