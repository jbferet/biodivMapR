# ==============================================================================
# biodivMapR
# Lib_MapAlphaDiversity.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This Library produces maps of alpha diversity indicators (Shannon, Simpson,
# Fischer...) based on spectral species file
# ==============================================================================

#' maps alpha diversity indicators  based on prior selection of PCs
#'
#' @param Input_Image_File character. Path and name of the image to be processed.
#' @param Input_Mask_File character. Path and name of the mask corresponding to the image to be processed.
#' @param Output_Dir character. Output directory.
#' @param window_size numeric. Size of spatial units (in pixels) to compute diversity.
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...).
#' @param nbclusters numeric. Number of clusters defined in k-Means.
#' @param MinSun numeric. Minimum proportion of sunlit pixels required to consider plot.
#' @param pcelim numeric. Minimum contribution (in \%) required for a spectral species.
#' @param Index_Alpha character. Either 'Shannon', 'Simpson' or 'Fisher'.
#' @param FullRes boolean. Full resolution.
#' @param LowRes boolean. Low resolution.
#' @param MapSTD boolean. map of standard deviation of the alpha diversity map (over repetitions)
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param ClassifMap character. If FALSE, perform standard biodivMapR based on SpectralSpecies.
#'                              else corresponds to path for a classification map.
#'
#' @return None
#' @importFrom stars read_stars write_stars
#' @export
map_alpha_div <- function(Input_Image_File = FALSE,
                          Input_Mask_File = FALSE,
                          Output_Dir = '',
                          window_size = 10,
                          TypePCA = "SPCA",
                          nbclusters = 50,
                          MinSun = 0.25,
                          pcelim = 0.02,
                          Index_Alpha = "Shannon",
                          FullRes = FALSE, LowRes = TRUE, MapSTD = TRUE,
                          nbCPU = 1, MaxRAM = 0.25, ClassifMap = FALSE) {

  # 1- get path for spectral species path, and possibly update Input_Image_File
  # and nbclusters if using classification map as input data
  SSDpathlist <- get_SSpath(Output_Dir, Input_Image_File, TypePCA, ClassifMap, nbclusters)
  Spectral_Species_Path <- SSDpathlist$Spectral_Species_Path
  SSD_Dir <- SSDpathlist$SSD_Dir
  Input_Image_File <- SSDpathlist$Input_Image_File
  nbclusters <- SSDpathlist$nbclusters

  # 2- COMPUTE ALPHA DIVERSITY
  ALPHA <- compute_alpha_metrics(Spectral_Species_Path = Spectral_Species_Path,
                                 SSD_Dir = SSD_Dir,
                                 window_size = window_size,
                                 Input_Mask_File = Input_Mask_File,
                                 nbclusters = nbclusters,
                                 MinSun = MinSun,
                                 pcelim = pcelim,
                                 nbCPU = nbCPU,
                                 MaxRAM = MaxRAM,
                                 Index_Alpha = Index_Alpha)

  # 3- SAVE ALPHA DIVERSITY MAPS
  print("Write alpha diversity maps")
  # which spectral indices will be computed
  Shannon <- Simpson <- Fisher <- FALSE
  if (length((grep("Shannon", Index_Alpha))) > 0) Shannon <- TRUE
  if (length((grep("Simpson", Index_Alpha))) > 0) Simpson <- TRUE
  if (length((grep("Fisher", Index_Alpha))) > 0) Fisher <- TRUE

  Output_Dir_Alpha <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "ALPHA")
  HDR <- ALPHA$HDR
  if (Shannon == TRUE) {
    Index <- "Shannon"
    HDR$`band names` <- Index
    Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
    write_raster(Image = ALPHA$Shannon, HDR = HDR, ImagePath = Alpha_Path,
                 window_size = window_size, FullRes = FullRes, LowRes = LowRes,
                 SmoothImage = TRUE)
    if (MapSTD == TRUE) {
      Index <- "Shannon_SD"
      HDR$`band names` <- Index
      Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
      write_raster(ALPHA$Shannon_SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    }
  }

  if (Fisher == TRUE) {
    Index <- "Fisher"
    HDR$`band names` <- Index
    Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
    write_raster(ALPHA$Fisher, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    if (MapSTD == TRUE) {
      Index <- "Fisher_SD"
      HDR$`band names` <- Index
      Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
      write_raster(ALPHA$Fisher_SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    }
  }

  if (Simpson == TRUE) {
    Index <- "Simpson"
    HDR$`band names` <- Index
    Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
    write_raster(ALPHA$Simpson, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    if (MapSTD == TRUE) {
      Index <- "Simpson_SD"
      HDR$`band names` <- Index
      Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
      write_raster(ALPHA$Simpson_SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    }
  }
  return(invisible())
}


#' compute alpha diversity from spectral species computed for a plot
#' expecting a matrix of spectral species (n pixels x p repetitions)
#'
#' @param SpectralSpecies_Plot numeric. matrix of spectral species
#' @param pcelim minimum contribution of spectral species to estimation of diversity
#' each spectral species with a proprtion < pcelim is eliminated before computation of diversity
#
#' @return list of alpha diversity metrics
#' @export

compute_ALPHA_FromPlot <- function(SpectralSpecies_Plot,pcelim = 0.02){

  nb_partitions <- dim(SpectralSpecies_Plot)[2]
  Richness.tmp <- Shannon.tmp <- Fisher.tmp <- Simpson.tmp <- vector(length = nb_partitions)
  for (i in seq_len(nb_partitions)){
    # compute distribution of spectral species
    Distritab <- table(SpectralSpecies_Plot[,i])
    # compute distribution of spectral species
    Pixel.Inventory <- as.data.frame(Distritab)
    SumPix <- sum(Pixel.Inventory$Freq)
    ThreshElim <- pcelim*SumPix
    ElimZeros <- which(Pixel.Inventory$Freq<ThreshElim)
    if (length(ElimZeros)>=1){
      Pixel.Inventory <- Pixel.Inventory[-ElimZeros,]
    }
    if (length(which(Pixel.Inventory$Var1==0))==1){
      Pixel.Inventory <- Pixel.Inventory[-which(Pixel.Inventory$Var1==0),]
    }
    Alpha <- get_alpha_metrics(Pixel.Inventory$Freq)
    # Alpha diversity
    Richness.tmp[i] <- as.numeric(Alpha$Richness)
    Fisher.tmp[i] <- Alpha$fisher
    Shannon.tmp[i] <- Alpha$Shannon
    Simpson.tmp[i] <- Alpha$Simpson
  }
  Richness <- mean(Richness.tmp)
  Fisher <- mean(Fisher.tmp)
  Shannon <- mean(Shannon.tmp)
  Simpson <- mean(Simpson.tmp)
  alpha <- list('Richness' = Richness, 'Fisher' = Fisher,
                'Shannon' = Shannon, 'Simpson' = Simpson,
                'Richness_ALL' = Richness.tmp, 'Fisher_ALL' = Fisher.tmp,
                'Shannon_ALL' = Shannon.tmp, 'Simpson_ALL' = Simpson.tmp)
  return(alpha)
}

#' Computes alpha diversity metrics based on spectral species
#
#' @param Spectral_Species_Path character. path for spectral species file
#' @param SSD_Dir character. path for spectral species distribution file to be written
#' @param window_size numeric. size of spatial units (in pixels) to compute diversity
#' @param Input_Mask_File character. Path and name of the mask corresponding to the image to be processed.
#' @param nbclusters numeric. number of clusters defined in k-Means
#' @param MinSun numeric. minimum proportion of sunlit pixels required to consider plot
#' @param pcelim numeric. percentage of occurence of a cluster below which cluster is eliminated
#' @param nbCPU numeric. number of CPUs available
#' @param MaxRAM numeric. maximum RAM available
#' @param Index_Alpha list. list of alpha diversity indices to be computed from spectral species
#' @param progressbar boolean. should progress bar be displayed (set to TRUE only if no conflict of parallel process)
#
#' @return list of mean and SD of alpha diversity metrics
#' @import cli
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor handlers with_progress
#' @importFrom stats sd
#' @export

compute_alpha_metrics <- function(Spectral_Species_Path,
                                  SSD_Dir,
                                  window_size,
                                  Input_Mask_File = FALSE,
                                  nbclusters,
                                  MinSun = 0.25,
                                  pcelim = 0.02,
                                  nbCPU = 1,
                                  MaxRAM = 0.25,
                                  Index_Alpha = 'Shannon', progressbar = FALSE) {

  # Prepare files to read and write: spectral species (R), spectral species
  # distribution (W) and sunlit proportion per window (W)
  # optional mask file
  SS_HDR <- get_HDR_name(Spectral_Species_Path)
  HDR_SS <- read_ENVI_header(SS_HDR)
  nbPieces <- split_image(HDR_SS, MaxRAM)
  # prepare for reading spectral species file
  SeqRead_SS <- where_to_read_kernel(HDR = HDR_SS,
                                     nbPieces = nbPieces,
                                     SE_Size = window_size)

  if (Input_Mask_File==FALSE){
    SeqRead_Mask <- FALSE
    HDR_Mask <- FALSE
  } else {
    Mask_HDR <- get_HDR_name(Input_Mask_File)
    HDR_Mask <- read_ENVI_header(Mask_HDR)
    # prepare for reading mask file
    SeqRead_Mask <- where_to_read_kernel(HDR_Mask, nbPieces, window_size)
  }

  # prepare for writing spectral species distribution file
  SSD_Path <- file.path(SSD_Dir, "SpectralSpecies_Distribution")
  HDR_SSD <- prepare_HDR_SSD(HDR_SS, SSD_Path, nbclusters, window_size = window_size)
  SeqWrite_SSD <- where_to_write_kernel(HDR_SS = HDR_SS,
                                        HDR_SSD = HDR_SSD,
                                        nbPieces = nbPieces,
                                        SE_Size = window_size)
  # prepare for writing sunlit proportion file
  Sunlit_Path <- file.path(SSD_Dir, "SpectralSpecies_Distribution_Sunlit")
  HDR_Sunlit <- prepare_HDR_Sunlit(HDR_SSD, Sunlit_Path)
  SeqWrite_Sunlit <- where_to_write_kernel(HDR_SS, HDR_Sunlit, nbPieces, window_size)

  # for each piece of image
  ReadWrite <- RW_bytes_all(SeqRead_SS, SeqWrite_SSD, SeqWrite_Sunlit, SeqRead_Mask)

  if (nbPieces>1 & progressbar == TRUE){
    handlers(global = TRUE)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = nbPieces)
      ALPHA <- lapply(ReadWrite,
                      FUN = convert_PCA_to_SSD,
                      Spectral_Species_Path = Spectral_Species_Path,
                      HDR_SS = HDR_SS, HDR_SSD = HDR_SSD, HDR_Sunlit = HDR_Sunlit,
                      window_size = window_size, nbclusters = nbclusters,
                      MinSun = MinSun, pcelim = pcelim, Index_Alpha = Index_Alpha,
                      SSD_Path = SSD_Path,
                      Input_Mask_File = Input_Mask_File, HDR_Mask = HDR_Mask,
                      Sunlit_Path = Sunlit_Path,
                      p= p, nbCPU = nbCPU, nbPieces = nbPieces)
    })
  } else {
    ALPHA <- lapply(ReadWrite,
                    FUN = convert_PCA_to_SSD,
                    Spectral_Species_Path = Spectral_Species_Path,
                    HDR_SS = HDR_SS, HDR_SSD = HDR_SSD, HDR_Sunlit = HDR_Sunlit,
                    window_size = window_size, nbclusters = nbclusters,
                    Input_Mask_File = Input_Mask_File, HDR_Mask = HDR_Mask,
                    MinSun = MinSun, pcelim = pcelim, Index_Alpha = Index_Alpha,
                    SSD_Path = SSD_Path,
                    Sunlit_Path = Sunlit_Path,
                    p= NULL, nbCPU = nbCPU, nbPieces = nbPieces)
  }

  # create ful map from chunks
  Shannon_Mean_Chunk <- Fisher_Mean_Chunk <- Simpson_Mean_Chunk <- list()
  Shannon_SD_Chunk <- Fisher_SD_Chunk <- Simpson_SD_Chunk <- list()
  for (i in seq_len(length(ALPHA))) {
    Shannon_Mean_Chunk[[i]] <- ALPHA[[i]]$Shannon
    Fisher_Mean_Chunk[[i]] <- ALPHA[[i]]$Fisher
    Simpson_Mean_Chunk[[i]] <- ALPHA[[i]]$Simpson
    Shannon_SD_Chunk[[i]] <- ALPHA[[i]]$Shannon_SD
    Fisher_SD_Chunk[[i]] <- ALPHA[[i]]$Fisher_SD
    Simpson_SD_Chunk[[i]] <- ALPHA[[i]]$Simpson_SD
  }
  Shannon.Mean <- do.call(rbind, Shannon_Mean_Chunk)
  Fisher.Mean <- do.call(rbind, Fisher_Mean_Chunk)
  Simpson.Mean <- do.call(rbind, Simpson_Mean_Chunk)
  Shannon_SD <- do.call(rbind, Shannon_SD_Chunk)
  Fisher_SD <- do.call(rbind, Fisher_SD_Chunk)
  Simpson_SD <- do.call(rbind, Simpson_SD_Chunk)
  # prepare HDR for alpha diversity
  HDR <- HDR_SSD
  HDR$bands <- 1
  HDR$`data type` <- 4
  HDR$lines <- dim(Shannon.Mean)[1]
  HDR$samples <- dim(Shannon.Mean)[2]
  my_list <- list(
    "Shannon" = Shannon.Mean, "Fisher" = Fisher.Mean, "Simpson" = Simpson.Mean,
    "Shannon_SD" = Shannon_SD, "Fisher_SD" = Fisher_SD, "Simpson_SD" = Simpson_SD, "HDR" = HDR
  )
  return(my_list)
}

#' Convert PCA into SSD based on previous clustering
#'
#' @param ReadWrite list. includes byte coordinates where to read and write from files
#' @param Spectral_Species_Path character. path for Spectral_Species raster file
#' @param HDR_SS list. HDR file for spectral species
#' @param HDR_SSD list. HDR file for spectral species distribution
#' @param HDR_Sunlit list. HDR file for spectral species distribution
#' @param window_size numeric. window size
#' @param nbclusters numeric. number of clusters
#' @param Input_Mask_File character. input mask file (useful when using classification map)
#' @param HDR_Mask list. HDR file for mask
#' @param MinSun numeric. Minimum proportion of sunlit pixels required to consider plot.
#' @param pcelim numeric. Minimum contribution (in \%) required for a spectral species.
#' @param Index_Alpha list. list of spectral metrics to be computed
#' @param SSD_Path character. path for spectral species distribution file
#' @param Sunlit_Path character. path for sunlit proportion file
#' @param p list. progressor object for progress bar
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param nbPieces numeric. number of pieces to split file read into
#' @param progressbar boolean. should progress bar be displayed (set to TRUE only if no conflict of parallel process)
#'
#' @return list of mean and SD of alpha diversity metrics
#' @export

convert_PCA_to_SSD <- function(ReadWrite, Spectral_Species_Path,
                               HDR_SS, HDR_SSD, HDR_Sunlit,
                               window_size, nbclusters,
                               Input_Mask_File = FALSE, HDR_Mask = FALSE,
                               MinSun = 0.25, pcelim = 0.02,
                               Index_Alpha = c('Shannon'),
                               SSD_Path, Sunlit_Path,
                               p = NULL, nbCPU = 1, nbPieces= 1,
                               progressbar = FALSE) {

  ImgFormat <- "3D"
  SSD_Format <- ENVI_type2bytes(HDR_SSD)
  SS_Format <- ENVI_type2bytes(HDR_SS)
  Sunlit_Format <- ENVI_type2bytes(HDR_Sunlit)
  if (typeof(HDR_Mask)=='list'){
    Mask_Format <- ENVI_type2bytes(HDR_Mask)
  }

  # read image chunk
  SS_Chunk <- read_BIL_image_subset(Spectral_Species_Path, HDR_SS,
                                    ReadWrite$RW_SS$Byte_Start, ReadWrite$RW_SS$lenBin,
                                    ReadWrite$RW_SS$nbLines, SS_Format, ImgFormat)

  if (!Input_Mask_File==FALSE){
    Mask_Chunk <- read_BIL_image_subset(Input_Mask_File, HDR_Mask,
                                        ReadWrite$RW_Mask$Byte_Start, ReadWrite$RW_Mask$lenBin,
                                        ReadWrite$RW_Mask$nbLines, Mask_Format, ImgFormat)
  } else {
    Mask_Chunk <- FALSE
  }

  # compute SSD from each window of the image chunk
  SSD_Alpha <- compute_SSD(Image_Chunk = SS_Chunk, window_size = window_size,
                           nbclusters = nbclusters, MinSun = MinSun, pcelim = pcelim,
                           Index_Alpha = Index_Alpha, nbCPU = nbCPU, nbPieces = nbPieces,
                           Mask_Chunk = Mask_Chunk, progressbar = progressbar)

  # write spectral Species ditribution file
  fidSSD <- file(description = SSD_Path, open = "r+b", blocking = TRUE,
                 encoding = getOption("encoding"), raw = FALSE)
  if (!ReadWrite$RW_SSD$Byte_Start == 1) {
    seek(fidSSD, where = ReadWrite$RW_SSD$Byte_Start - 1, origin = "start", rw = "write")
  }
  SSD_Chunk <- aperm(array(SSD_Alpha$SSD, c(ReadWrite$RW_SSD$nbLines, HDR_SSD$samples, HDR_SSD$bands)), c(2, 3, 1))
  writeBin(c(SSD_Chunk), fidSSD, size = SSD_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
  close(fidSSD)

  # write PCsunlit pixels corresponding to SSD file
  fidSunlit <- file(description = Sunlit_Path, open = "r+b", blocking = TRUE,
                    encoding = getOption("encoding"), raw = FALSE)
  if (!ReadWrite$RW_Sunlit$Byte_Start == 1) {
    seek(fidSunlit, where = ReadWrite$RW_Sunlit$Byte_Start - 1, origin = "start", rw = "write")
  }
  Sunlit.Chunk <- t(SSD_Alpha$PCsun)
  writeBin(c(Sunlit.Chunk), fidSunlit, size = Sunlit_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
  close(fidSunlit)

  #if progress bar called
  if (!is.null(p)){p()}

  # a bit of cleaning in the environment
  rm(SSD_Chunk)
  rm(Sunlit.Chunk)
  gc()

  # gt mean value and standard deviation for alpha diversity metrics
  Shannon_Mean_Chunk <- apply(SSD_Alpha$Shannon, 1:2, mean)
  Fisher_Mean_Chunk <- apply(SSD_Alpha$Fisher, 1:2, mean)
  Simpson_Mean_Chunk <- apply(SSD_Alpha$Simpson, 1:2, mean)
  Shannon_SD_Chunk <- apply(SSD_Alpha$Shannon, 1:2, sd)
  Fisher_SD_Chunk <- apply(SSD_Alpha$Fisher, 1:2, sd)
  Simpson_SD_Chunk <- apply(SSD_Alpha$Simpson, 1:2, sd)
  rm(SSD_Alpha)
  gc()
  my_list <- list("Shannon" = Shannon_Mean_Chunk, "Fisher" = Fisher_Mean_Chunk,
                  "Simpson" = Simpson_Mean_Chunk, "Shannon_SD" = Shannon_SD_Chunk,
                  "Fisher_SD" = Fisher_SD_Chunk, "Simpson_SD" = Simpson_SD_Chunk)
  return(my_list)
}

#' compute spectral species distribution from original spectral species map
#'
#' @param Image_Chunk numeric. 3D image chunk of spectral species
#' @param window_size numeric. size of spatial units (in pixels) to compute diversity
#' @param nbclusters number of clusters defined in k-Means
#' @param MinSun minimum proportion of sunlit pixels required to consider plot
#' @param pcelim minimum proportion for a spectral species to be included
#' @param Index_Alpha list. list of alpha diversity indices to be computed from spectral species
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param nbPieces numeric. number of pieces in which original image is split
#' @param Mask_Chunk numeric. 3D image chunk of mask (optional)
#' @param progressbar boolean. should progress bar be displayed (set to TRUE only if no conflict of parallel process)
#
#' @return list of alpha diversity metrics for each iteration
#' @importFrom vegan fisher.alpha
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom snow splitList
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export

compute_SSD <- function(Image_Chunk, window_size, nbclusters,
                        MinSun = 0.25, pcelim = 0.02,
                        Index_Alpha = "Shannon", nbCPU = 1, nbPieces = 1,
                        Mask_Chunk = FALSE, progressbar = FALSE) {

  nbi <- ceiling(dim(Image_Chunk)[1] / window_size)
  nbj <- ceiling(dim(Image_Chunk)[2] / window_size)
  # nbi <- round(dim(Image_Chunk)[1] / window_size)
  # nbj <- round(dim(Image_Chunk)[2] / window_size)
  nb_partitions <- dim(Image_Chunk)[3]
  SSDMap <- array(NA, c(nbi, nbj, nb_partitions * nbclusters))
  shannonIter <- FisherAlpha <- SimpsonAlpha <- array(NA, dim = c(nbi, nbj, nb_partitions))
  PCsun <- matrix(NA, nrow = nbi, ncol = nbj)

  alphaIdx <- list()
  alphaIdx$Shannon <- alphaIdx$Simpson <- alphaIdx$Fisher <- FALSE
  if (length((grep("Shannon", Index_Alpha))) > 0) alphaIdx$Shannon <- TRUE
  if (length((grep("Simpson", Index_Alpha))) > 0) alphaIdx$Simpson <- TRUE
  if (length((grep("Fisher", Index_Alpha))) > 0) alphaIdx$Fisher <- TRUE

  # split the image chunk into a list of windows
  listAlpha <- listIJ <- list()
  ij <- 0
  for (ii in seq_len(nbi)) {
    # for each window in the column
    for (jj in seq_len(nbj)) {
      li <- ((ii - 1) * window_size) + 1
      ui <- min(c(ii * window_size,dim(Image_Chunk)[1]))
      lj <- ((jj - 1) * window_size) + 1
      uj <- min(c(jj * window_size,dim(Image_Chunk)[2]))
      # put all iterations in a 2D matrix shape
      ijit <- t(matrix(Image_Chunk[li:ui, lj:uj, ], ncol = nb_partitions))
      # keep non zero values
      if (typeof(Mask_Chunk)=='integer'){
        ijit_mask <- t(matrix(Mask_Chunk[li:ui, lj:uj, ], ncol = 1))
        ijit <- matrix(ijit[, which(!ijit_mask[1, ] == 0)], nrow = nb_partitions)
      } else {
        ijit <- matrix(ijit[, which(!ijit[1, ] == 0)], nrow = nb_partitions)

      }
      nbPix_Sunlit <- dim(ijit)[2]
      PCsun[ii, jj] <- nbPix_Sunlit / window_size**2
      ij <- ij+1
      listAlpha[[ij]] <- listIJ[[ij]] <- list()
      listIJ[[ij]]$i <- ii
      listIJ[[ij]]$j <- jj
      listAlpha[[ij]]$nbPix_Sunlit <- nbPix_Sunlit
      listAlpha[[ij]]$data <- ijit
      listAlpha[[ij]]$PCsun <- PCsun[ii, jj]
    }
  }

  # tictoc::tic()
  # alphaSSD <- pbapply::pblapply(listAlpha, compute_ALPHA_SSD_per_window, MinSun, nb_partitions, alphaIdx)
  # tictoc::toc()

  if (nbCPU > 1){
    listAlpha <- snow::splitList(listAlpha,ncl = nbCPU)
    # plan(multisession, workers = nbCPU)
    cl <- parallel::makeCluster(nbCPU)
    plan("cluster", workers = cl)  ## same as plan(multisession, workers = nbCPU)
    # add a progress bar if the image was read in one piece only
    if (nbPieces ==1){
      if (progressbar==TRUE){
        handlers(global = TRUE)
        handlers("cli")
        with_progress({
          p <- progressr::progressor(steps = nbCPU)
          alphaSSD <- future.apply::future_lapply(X = listAlpha,
                                                  FUN = compute_ALPHA_SSD_per_window_list,
                                                  nb_partitions = nb_partitions,
                                                  nbclusters = nbclusters,
                                                  alphaIdx = alphaIdx,
                                                  MinSun = MinSun, pcelim = pcelim, p = p)
        })
      } else {
        alphaSSD <- future.apply::future_lapply(X = listAlpha,
                                                FUN = compute_ALPHA_SSD_per_window_list,
                                                nb_partitions = nb_partitions,
                                                nbclusters = nbclusters,
                                                alphaIdx = alphaIdx,
                                                MinSun = MinSun, pcelim = pcelim)
      }
    } else {
      alphaSSD <- future.apply::future_lapply(X = listAlpha,
                                              FUN = compute_ALPHA_SSD_per_window_list,
                                              nb_partitions = nb_partitions,
                                              nbclusters = nbclusters,
                                              alphaIdx = alphaIdx,
                                              MinSun = MinSun, pcelim = pcelim, p = NULL)
    }
    parallel::stopCluster(cl)
    plan(sequential)
    shannonIter_list <- SimpsonAlpha_list <- FisherAlpha_list <- SSDMap_list <- list()
    for (i in seq_len(nbCPU)){
      shannonIter_list <- append(shannonIter_list,alphaSSD[[i]]$shannonIter)
      SimpsonAlpha_list <- append(SimpsonAlpha_list,alphaSSD[[i]]$SimpsonAlpha  )
      FisherAlpha_list <- append(FisherAlpha_list,alphaSSD[[i]]$FisherAlpha  )
      SSDMap_list <- append(SSDMap_list,alphaSSD[[i]]$SSDMap)
    }

  } else {
    alphaSSD <- lapply(X = listAlpha, FUN = compute_ALPHA_SSD_per_window,
                       nb_partitions = nb_partitions,
                       nbclusters = nbclusters,
                       alphaIdx = alphaIdx,
                       MinSun = MinSun,
                       pcelim = pcelim)
    shannonIter_list <- lapply(alphaSSD,'[[', 'shannonIter')
    SimpsonAlpha_list <- lapply(alphaSSD,'[[', 'SimpsonAlpha')
    FisherAlpha_list <- lapply(alphaSSD,'[[', 'FisherAlpha')
    SSDMap_list <- lapply(alphaSSD,'[[', 'SSDMap')
  }
  for (i in seq_len(ij)){
    ii <- listIJ[[i]]$i
    jj <- listIJ[[i]]$j
    SSDMap[ii, jj,] <- SSDMap_list[[i]]
    shannonIter[ii, jj, ] <- shannonIter_list[[i]]
    SimpsonAlpha[ii, jj, ] <- SimpsonAlpha_list[[i]]
    FisherAlpha[ii, jj, ] <- FisherAlpha_list[[i]]
  }
  my_list <- list("Shannon" = shannonIter, "Fisher" = FisherAlpha,
                  "Simpson" = SimpsonAlpha, "SSD" = SSDMap, "PCsun" = PCsun)
  return(my_list)
}


#' defines path for the spectral species raster file. Depends on preprocessing:
#' - if standard dimensionality reduction (PCA, SPCA, MNF) was used, the
#' corresponding directory is pointed
#' - if no dimensionality reduction, or an alternative dimensionality reduction
#' (spectral indices, tasseled cap... ) was used, user defined subdirectory is used
#' - Diversity metrics can also be computed from a classification map
#'
#' @param Output_Dir character. Output directory.
#' @param Input_Image_File character. Path and name of the image to be processed.
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...).
#' @param ClassifMap character. If FALSE, perform standard biodivMapR based on SpectralSpecies.
#'                              else corresponds to path for a classification map.
#' @param nbclusters numeric. number of clusters (possibly updated)
#'
#' @return list. includes Spectral_Species_Path
#' @importFrom stars read_stars write_stars
#' @importFrom tools file_path_sans_ext
#' @export

get_SSpath <- function(Output_Dir, Input_Image_File, TypePCA, ClassifMap, nbclusters){
  # if 'standard use' of biodivMapR
  if (ClassifMap == FALSE){
    Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "SpectralSpecies")
    Spectral_Species_Path <- file.path(Output_Dir_SS, "SpectralSpecies")
  } else {
    # if biodivMapR used with classification map
    message("Classification Map will be used instead of SpectralSpecies")
    message("Classes are expected to be integer values")
    message("conversion to ENVI File if not the format of the original classification map")
    message("Class '0' will be considered as No Data and excluded")
    if (! file.exists(ClassifMap)){
      message("classification map is not found:")
      print(ClassifMap)
      stop()
    } else {
      message("updating nbclusters based on number of classes")
      ClassifRaster <- stars::read_stars(ClassifMap,proxy = FALSE)
      Classif_Values <- ClassifRaster[[1]]
      CheckClasses <- tryCatch(
        {
          CheckClasses <- max(Classif_Values,na.rm = TRUE)
        },
        error = function(cond) {
          CheckClasses <- length(unique(Classif_Values))
          return(CheckClasses)
        },
        finally = {
        }
      )
      nbclusters <- CheckClasses
      # nbclusters <- max(Classif_Values,na.rm = TRUE)
      # nbclusters <- length(unique(Classif_Values))
      message(paste("Number of classes : "),nbclusters)
    }
    # save classification map in proper format in output directory
    # if not expected file format for Spectral Species map
    info <- get_gdal_info(ClassifMap)
    driver <- info$driverShortName
    df <- unique(info$bands$type)
    if (driver=='ENVI' & df =='Byte'){
      if (Input_Image_File==FALSE){
        Input_Image_File <- tools::file_path_sans_ext(basename(ClassifMap))
      }
      Spectral_Species_Path <- ClassifMap
    } else {
      if (Input_Image_File==FALSE){
        Input_Image_File <- tools::file_path_sans_ext(basename(ClassifMap))
      }
      Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "UserClassification")
      Spectral_Species_Path <- file.path(Output_Dir_SS, "UserClassification")
      if (! file.exists(Spectral_Species_Path)){
        stars::write_stars(ClassifRaster, Spectral_Species_Path, driver =  "ENVI",type='Byte')
      } else {
        message("This already existing classification map will be used")
        print(Spectral_Species_Path)
        message("Please delete it and re-run if you updated classification since last run")
      }
    }
  }
  Output_Dir_SSD <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "SpectralSpecies")
  res <- list('Spectral_Species_Path' = Spectral_Species_Path,
              'SSD_Dir' = Output_Dir_SSD,
              'Input_Image_File' = Input_Image_File,
              'nbclusters' = nbclusters)
  return(res)
}

#' computes shannon index from a distribution
#' (faster than version implemented in vegan package)
#'
#' @param Distrib Distribution
#'
#' @return Shannon index correspnding to the distribution
#' @export

get_Shannon <- function(Distrib) {
  Distrib <- Distrib / sum(Distrib, na.rm = TRUE)
  Distrib <- Distrib[which(!Distrib == 0)]
  shannon <- -1 * sum(Distrib * log(Distrib), na.rm = TRUE)
  return(shannon)
}

#' computes Simpson index from a distribution
#' (faster than version implemented in vegan package)
#'
#' @param Distrib Distribution
#'
#' @return Simpson index correspnding to the distribution
#' @export

get_Simpson <- function(Distrib) {
  Distrib <- Distrib / sum(Distrib, na.rm = TRUE)
  Simpson <- 1 - sum(Distrib * Distrib, na.rm = TRUE)
  return(Simpson)
}

#' computes alpha diversity metrics and spectral species distribution for a list
#' of windows (individual spatial unit) defined in biodivmapR
#'
#' @param listAlpha list.
#' @param nb_partitions numeric. number of iterations performed with biodivMapR
#' @param nbclusters numeric. number of clusters defined in k-Means
#' @param alphaIdx list.
#' @param MinSun minimum proportion of sunlit pixels required to consider plot
#' @param pcelim minimum proportion for a spectral species to be included
#' @param p list. progressor object for progress bar


#' @return list_assd list of alpha metrics and spectral species distribution
#' corresponding to the list of windows
#' @export

compute_ALPHA_SSD_per_window_list <- function(listAlpha, nb_partitions, nbclusters,
                                              alphaIdx, MinSun = 0.25, pcelim = 0.02,
                                              p = NULL) {
  alphaSSD <- lapply(X = listAlpha,
                     FUN = compute_ALPHA_SSD_per_window,
                     nb_partitions = nb_partitions, nbclusters = nbclusters,
                     alphaIdx = alphaIdx, MinSun = MinSun, pcelim = pcelim)

  shannonIter <- lapply(alphaSSD,'[[','shannonIter')
  SimpsonAlpha <- lapply(alphaSSD,'[[','SimpsonAlpha')
  FisherAlpha <- lapply(alphaSSD,'[[','FisherAlpha')
  SSDMap <- lapply(alphaSSD,'[[',4)
  list_assd <- list('shannonIter' = shannonIter, 'SimpsonAlpha' = SimpsonAlpha,
                    'FisherAlpha' = FisherAlpha, 'SSDMap' = SSDMap)
  if (!is.null(p)){p()}
  return(list_assd)
}

#' computes alpha diversity metrics and spectral species distribution for a
#' window (individual spatial unit) defined in biodivmapR
#'
#' @param listAlpha list.
#' @param nb_partitions numeric. number of iterations performed with biodivMapR
#' @param nbclusters numeric. number of clusters defined in k-Means
#' @param alphaIdx list.
#' @param MinSun minimum proportion of sunlit pixels required to consider plot
#' @param pcelim minimum proportion for a spectral species to be included

#' @importFrom snow splitRows
#' @return res list of alpha metrics and spectral species distribution
#' corresponding to the window
#' @export

compute_ALPHA_SSD_per_window <- function(listAlpha, nb_partitions, nbclusters,
                                         alphaIdx, MinSun = 0.25, pcelim = 0.02) {

  if (listAlpha$PCsun > MinSun) {
    # compute diversity metrics and spectral species distribution for each iteration
    SSD <- apply(X = listAlpha$data,MARGIN = 1,FUN = table)
    # unexplainablebug occuring sometimes and resulting in impossibility to apply
    # table function to each row
    if (!typeof(SSD)=='list'){
      SSD <- lapply(X = snow::splitRows(listAlpha$data, ncl = nrow(listAlpha$data)), FUN = table)
    }

    # future optimization using data.table to be implemented beforehands
    # data2 = matrix(sample(seq(0,50), 2000000,replace = T),nrow=20000)
    # tictoc::tic()
    # SSD = (data2 %>% t() %>% as.data.table %>% melt(measure.vars=paste0('V',1:20000)))[,.(count=.N), by=c('variable', 'value')]
    # tictoc::toc()


    alphawin <- lapply(SSD, FUN = compute_ALPHA_per_window,
                       nbPix_Sunlit = listAlpha$nbPix_Sunlit,
                       alphaIdx = alphaIdx,
                       nbclusters = nbclusters,
                       pcelim = pcelim)
    shannonIter <- unlist(lapply(alphawin,'[[',1))
    SimpsonAlpha <- unlist(lapply(alphawin,'[[',2))
    FisherAlpha <- unlist(lapply(alphawin,'[[',3))
    SSDMap <- unlist(lapply(alphawin,'[[',4))

  } else {
    shannonIter <- NA*vector(length = nb_partitions)
    FisherAlpha <- NA*vector(length = nb_partitions)
    SimpsonAlpha <- NA*vector(length = nb_partitions)
    SSDMap <- NA*vector(length = nb_partitions*nbclusters)
  }
  res <- list('shannonIter' = shannonIter, 'SimpsonAlpha' = SimpsonAlpha,
              'FisherAlpha' = FisherAlpha, 'SSDMap' = SSDMap)
  return(res)
}

#' computes alpha diversity metrics and update spectral species distribution for
#' a window defined in biodivmapR
#'
#' @param SSD numeric.
#' @param nbPix_Sunlit numeric.
#' @param alphaIdx list. list of diversity metrics to be computed
#' @param nbclusters numeric. number of clusters defined in k-Means
#' @param pcelim minimum proportion for a spectral species to be included

#' @return shannon, Simpson, Fisher, SSDMap
#' @export

compute_ALPHA_per_window <- function(SSD, nbPix_Sunlit, alphaIdx, nbclusters,
                                     pcelim = 0.02){

  ClusterID <- as.numeric(names(SSD))
  if (pcelim > 0) {
    KeepSS <- which(SSD >= pcelim * nbPix_Sunlit)
    ClusterID <- ClusterID[KeepSS]
    SSD <- SSD[KeepSS]
  }
  SSDMap0 <- NA*vector(length = nbclusters)
  if (length(ClusterID)>0){
    SSDMap0[ClusterID] <- SSD
  }
  SimpsonAlpha <- shannonIter <- FisherAlpha <- NA
  if (alphaIdx$Shannon == TRUE) { shannonIter <- get_Shannon(as.numeric(SSD)) }
  if (alphaIdx$Simpson == TRUE) { SimpsonAlpha <- get_Simpson(as.numeric(SSD)) }
  if (alphaIdx$Fisher == TRUE) {
    if (length(SSD) > 2) {
      FisherAlpha <- fisher.alpha(as.numeric(SSD))
    }
  } else {
    FisherAlpha <- 0
  }
  res <- list('shannon' = shannonIter,
              'simpson' = SimpsonAlpha,
              'fisher' = FisherAlpha,
              'SSD' = SSDMap0)
  return(res)
}

#' prepares Spectral species distribution (SSD) file and header
#'
#' @param HDR_SS list. header of the spectral species file
#' @param SSD_Path character. path for spectral species file
#' @param nbclusters numeric. number of clusters
#' @param window_size numeric. window size
#'
#' @return HDR_SSD
#' @export

prepare_HDR_SSD <- function(HDR_SS, SSD_Path, nbclusters, window_size){
  # prepare SS distribution map
  HDR_SSD <- HDR_SS
  # define number of bands
  HDR_SSD$bands <- HDR_SS$bands * nbclusters
  # define image size
  HDR_SSD$samples <- ceiling(HDR_SS$samples / window_size)
  HDR_SSD$lines <- ceiling(HDR_SS$lines / window_size)
  # HDR_SSD$samples <- round(HDR_SS$samples / window_size)
  # HDR_SSD$lines <- round(HDR_SS$lines / window_size)
  # HDR_SSD$samples <- floor(HDR_SS$samples / window_size)
  # HDR_SSD$lines <- floor(HDR_SS$lines / window_size)
  HDR_SSD$interleave <- 'bil'
  HDR_SSD$`file type` <- NULL
  HDR_SSD$classes <- NULL
  HDR_SSD$`class names` <- NULL
  HDR_SSD$`class lookup` <- NULL
  # change resolution
  HDR_SSD <- change_resolution_HDR(HDR_SSD, window_size)
  HDR_SSD$`band names` <- NULL
  # create SSD file
  fidSSD <- file(
    description = SSD_Path, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidSSD)
  headerFpath <- paste(SSD_Path, ".hdr", sep = "")
  write_ENVI_header(HDR_SSD, headerFpath)
  return(HDR_SSD)
}

#' prepare for writing sunlit proportion file
#'
#' @param HDR_SSD list. header of the spectral species distribution file
#' @param Sunlit_Path character. path for spectral species file
#'
#' @return HDR_Sunlit
#' @export

prepare_HDR_Sunlit <- function(HDR_SSD, Sunlit_Path){
  # prepare proportion of sunlit pixels from each spatial unit
  HDR_Sunlit <- HDR_SSD
  # define number of bands
  HDR_Sunlit$bands <- 1
  # define number of bands
  HDR_Sunlit$`data type` <- 4
  # create SSD Sunlit mask
  fidSunlit <- file(
    description = Sunlit_Path, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidSunlit)
  headerFpath <- paste(Sunlit_Path, ".hdr", sep = "")
  write_ENVI_header(HDR_Sunlit, headerFpath)
  return(HDR_Sunlit)
}

#' identifies bytes where to read and write for each piece of image and all files
#'
#' @param SeqRead_SS list. coordinates corresponding to spectral species file
#' @param SeqWrite_SSD list. coordinates corresponding to spectral species distribution file
#' @param SeqWrite_Sunlit list. coordinates corresponding to sunlit proportion
#' @param SeqRead_Mask list. coordinates corresponding to optional mask file
#'
#' @return ReadWrite list of bytes coordinates for each read and write
#' @export

RW_bytes_all <- function(SeqRead_SS, SeqWrite_SSD, SeqWrite_Sunlit,
                         SeqRead_Mask = FALSE){
  nbPieces <- length(SeqRead_SS$ReadByte_Start)
  ReadWrite <- list()
  for (i in seq_len(nbPieces)) {
    ReadWrite[[i]] <- list()
    # ReadWrite[[i]]$RW_SS <- ReadWrite[[i]]$RW_SSD <- ReadWrite[[i]]$RW_Sunlit <- list()
    ReadWrite[[i]]$RW_SS <- RW_bytes(SeqRead_SS, i)
    ReadWrite[[i]]$RW_SSD <- RW_bytes(SeqWrite_SSD, i)
    ReadWrite[[i]]$RW_Sunlit <- RW_bytes(SeqWrite_Sunlit, i)
    if (typeof(SeqRead_Mask)=='list'){
      ReadWrite[[i]]$RW_Mask <- RW_bytes(SeqRead_Mask, i)
    }
  }
  return(ReadWrite)
}

#' identifies bytes where to read or write for each piece of an image in a given file
#'
#' @param SeqRW list. coordinates corresponding to spectral species file
#' @param piece numeric. identifier of the piece of image
#'
#' @return RW0 list of bytes coordinates for each read and write
#' @export

RW_bytes <- function(SeqRW, piece){
  RW0 <- list()
  RW0$Byte_Start <- SeqRW$ReadByte_Start[piece]
  RW0$nbLines <- SeqRW$Lines_Per_Chunk[piece]
  RW0$lenBin <- SeqRW$ReadByte_End[piece] - SeqRW$ReadByte_Start[piece] + 1
  return(RW0)
}
