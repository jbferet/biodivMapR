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
map_alpha_div <- function(Input_Image_File=FALSE, Output_Dir='', window_size=10,
                          TypePCA = "SPCA", nbclusters = 50,
                          MinSun = 0.25, pcelim = 0.02,
                          Index_Alpha = "Shannon", FullRes = TRUE,
                          LowRes = FALSE, MapSTD = FALSE,
                          nbCPU = FALSE, MaxRAM = FALSE,
                          ClassifMap = FALSE) {

  # if 'standard use' of biodivMapR
  if (ClassifMap == FALSE){
    Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "SpectralSpecies")
    Spectral_Species_Path <- file.path(Output_Dir_SS, "SpectralSpecies")
  } else {
    message("Classification Map will be used instead of SpectralSpecies")
    message("Classes are expected to be integer values")
    message("conversion to ENVI File if not the format of the original classification map")
    message("Class '0' will be considered as No Data and excluded")
    if (! file.exists(ClassifMap)){
      message("classification map is not found:")
      print(ClassifMap)
    } else {
      message("updating nbclusters based on number of classes")
      ClassifRaster <- stars::read_stars(ClassifMap,proxy = FALSE)
      Classif_Values <- ClassifRaster[[1]]
      nbclusters <- max(Classif_Values,na.rm = TRUE)
      message(paste("Number of classes : "),nbclusters)

      # save classification map in proper format in output directory
      # if not expected file format for Spectral Species map
      driver <- attr(rgdal::GDALinfo(ClassifMap,returnStats = FALSE), 'driver')
      df <- unique(attr(rgdal::GDALinfo(ClassifMap,returnStats = FALSE),"df")$GDType)
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
  }
  # 1- COMPUTE ALPHA DIVERSITY
  ALPHA <- compute_alpha_metrics(Spectral_Species_Path = Spectral_Species_Path, window_size = window_size,
                                 nbclusters = nbclusters, MinSun = MinSun, pcelim = pcelim,
                                 nbCPU = nbCPU, MaxRAM = MaxRAM, Index_Alpha = Index_Alpha)
  # 2- SAVE ALPHA DIVERSITY MAPS
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
      write_raster(ALPHA$Shannon.SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
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
      write_raster(ALPHA$Fisher.SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
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
      write_raster(ALPHA$Simpson.SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
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
  for (i in 1:nb_partitions){
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
  alpha <- list('Richness' = Richness, 'Fisher' = Fisher, 'Shannon' = Shannon, 'Simpson' = Simpson,
                'Richness_ALL' = Richness.tmp, 'Fisher_ALL' = Fisher.tmp, 'Shannon_ALL' = Shannon.tmp, 'Simpson_ALL' = Simpson.tmp)
  return(alpha)
}


# Map alpha diversity metrics based on spectral species
#
# @param Spectral_Species_Path character. path for spectral species file to be written
# @param window_size numeric. size of spatial units (in pixels) to compute diversity
# @param nbclusters numeric. number of clusters defined in k-Means
# @param pcelim numeric. percentage of occurence of a cluster below which cluster is eliminated
# @param nbCPU numeric. number of CPUs available
# @param MaxRAM numeric. maximum RAM available
# @param Index_Alpha list. list of alpha diversity indices to be computed from spectral species
# @param MinSun numeric. minimum proportion of sunlit pixels required to consider plot
#
# @return list of mean and SD of alpha diversity metrics
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats sd
compute_alpha_metrics <- function(Spectral_Species_Path, window_size, nbclusters,
                                  MinSun, pcelim, nbCPU = FALSE, MaxRAM = FALSE, Index_Alpha = "Shannon") {
  ##      read SpectralSpecies file and write distribution per spatial unit   ##
  SS_HDR <- get_HDR_name(Spectral_Species_Path)
  HDR_SS <- read_ENVI_header(SS_HDR)
  if (MaxRAM == FALSE) {
    MaxRAM <- 0.25
  }
  nbPieces_Min <- split_image(HDR_SS, MaxRAM)
  if (nbCPU == FALSE) {
    nbCPU <- 1
  }
  if (nbPieces_Min < nbCPU) {
    nbPieces_Min <- nbCPU
  }
  SeqRead.SS <- where_to_read_kernel(HDR_SS, nbPieces_Min, window_size)

  ##          prepare SS distribution map and corresponding sunlit map        ##
  # prepare SS distribution map
  SSD_Path <- paste(Spectral_Species_Path, "_Distribution", sep = "")
  HDR_SSD <- HDR_SS
  # define number of bands
  HDR_SSD$bands <- HDR_SS$bands * nbclusters
  # define image size
  HDR_SSD$samples <- round(HDR_SS$samples / window_size)
  HDR_SSD$lines <- round(HDR_SS$lines / window_size)
  HDR_SSD$interleave <- 'bil'
  HDR_SSD$`file type` <- NULL
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
  SeqWrite.SSD <- where_to_write_kernel(HDR_SS, HDR_SSD, nbPieces_Min, window_size)

  # prepare proportion of sunlit pixels from each spatial unit
  Sunlit_Path <- paste(SSD_Path, "_Sunlit", sep = "")
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
  SeqWrite.Sunlit <- where_to_write_kernel(HDR_SS, HDR_Sunlit, nbPieces_Min, window_size)

  # for each piece of image
  ReadWrite <- list()
  for (i in 1:nbPieces_Min) {
    ReadWrite[[i]] <- list()
    ReadWrite[[i]]$RW_SS <- ReadWrite[[i]]$RW_SSD <- ReadWrite[[i]]$RW_Sunlit <- list()
    ReadWrite[[i]]$RW_SS$Byte_Start <- SeqRead.SS$ReadByte_Start[i]
    ReadWrite[[i]]$RW_SS$nbLines <- SeqRead.SS$Lines_Per_Chunk[i]
    ReadWrite[[i]]$RW_SS$lenBin <- SeqRead.SS$ReadByte_End[i] - SeqRead.SS$ReadByte_Start[i] + 1

    ReadWrite[[i]]$RW_SSD$Byte_Start <- SeqWrite.SSD$ReadByte_Start[i]
    ReadWrite[[i]]$RW_SSD$nbLines <- SeqWrite.SSD$Lines_Per_Chunk[i]
    ReadWrite[[i]]$RW_SSD$lenBin <- SeqWrite.SSD$ReadByte_End[i] - SeqWrite.SSD$ReadByte_Start[i] + 1

    ReadWrite[[i]]$RW_Sunlit$Byte_Start <- SeqWrite.Sunlit$ReadByte_Start[i]
    ReadWrite[[i]]$RW_Sunlit$nbLines <- SeqWrite.Sunlit$Lines_Per_Chunk[i]
    ReadWrite[[i]]$RW_Sunlit$lenBin <- SeqWrite.Sunlit$ReadByte_End[i] - SeqWrite.Sunlit$ReadByte_Start[i] + 1
  }
  ImgFormat <- "3D"
  SSD_Format <- ENVI_type2bytes(HDR_SSD)
  SS_Format <- ENVI_type2bytes(HDR_SS)
  Sunlit_Format <- ENVI_type2bytes(HDR_Sunlit)

  # multiprocess of spectral species distribution and alpha diversity metrics
  plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
  Schedule_Per_Thread <- ceiling(nbPieces_Min / nbCPU)
  ALPHA <- future_lapply(ReadWrite,
    FUN = convert_PCA_to_SSD, Spectral_Species_Path = Spectral_Species_Path,
    HDR_SS = HDR_SS, HDR_SSD = HDR_SSD, SS_Format = SS_Format, SSD_Format = SSD_Format,
    ImgFormat = ImgFormat, window_size = window_size, nbclusters = nbclusters, MinSun = MinSun,
    pcelim = pcelim, Index_Alpha = Index_Alpha, SSD_Path = SSD_Path, Sunlit_Path = Sunlit_Path,
    Sunlit_Format = Sunlit_Format, future.scheduling = Schedule_Per_Thread
  )
  plan(sequential)
  # create ful map from chunks
  Shannon_Mean_Chunk <- Fisher_Mean_Chunk <- Simpson_Mean_Chunk <- list()
  Shannon_SD_Chunk <- Fisher_SD_Chunk <- Simpson_SD_Chunk <- list()
  for (i in 1:length(ALPHA)) {
    Shannon_Mean_Chunk[[i]] <- ALPHA[[i]]$Shannon
    Fisher_Mean_Chunk[[i]] <- ALPHA[[i]]$Fisher
    Simpson_Mean_Chunk[[i]] <- ALPHA[[i]]$Simpson
    Shannon_SD_Chunk[[i]] <- ALPHA[[i]]$Shannon.SD
    Fisher_SD_Chunk[[i]] <- ALPHA[[i]]$Fisher.SD
    Simpson_SD_Chunk[[i]] <- ALPHA[[i]]$Simpson.SD
  }
  Shannon.Mean <- do.call(rbind, Shannon_Mean_Chunk)
  Fisher.Mean <- do.call(rbind, Fisher_Mean_Chunk)
  Simpson.Mean <- do.call(rbind, Simpson_Mean_Chunk)
  Shannon.SD <- do.call(rbind, Shannon_SD_Chunk)
  Fisher.SD <- do.call(rbind, Fisher_SD_Chunk)
  Simpson.SD <- do.call(rbind, Simpson_SD_Chunk)
  # prepare HDR for alpha diversity
  HDR <- HDR_SSD
  HDR$bands <- 1
  HDR$`data type` <- 4
  HDR$lines <- dim(Shannon.Mean)[1]
  HDR$samples <- dim(Shannon.Mean)[2]
  my_list <- list(
    "Shannon" = Shannon.Mean, "Fisher" = Fisher.Mean, "Simpson" = Simpson.Mean,
    "Shannon.SD" = Shannon.SD, "Fisher.SD" = Fisher.SD, "Simpson.SD" = Simpson.SD, "HDR" = HDR
  )
  return(my_list)
}

# Convert PCA into SSD based on previous clustering
#
# @param ReadWrite
# @param Spectral_Species_Path
# @param HDR_SS
# @param HDR_SSD
# @param SS_Format
# @param SSD_Format
# @param ImgFormat
# @param window_size
# @param nbclusters
# @param MinSun
# @param pcelim
# @param Index_Alpha
# @param SSD_Path
# @param Sunlit_Path
# @param Sunlit_Format
#
# @param
# @param
# @return
convert_PCA_to_SSD <- function(ReadWrite, Spectral_Species_Path, HDR_SS, HDR_SSD,
                               SS_Format, SSD_Format, ImgFormat, window_size, nbclusters,
                               MinSun, pcelim, Index_Alpha, SSD_Path, Sunlit_Path, Sunlit_Format) {
  SS_Chunk <- read_BIL_image_subset(
    Spectral_Species_Path, HDR_SS,
    ReadWrite$RW_SS$Byte_Start, ReadWrite$RW_SS$lenBin,
    ReadWrite$RW_SS$nbLines, SS_Format, ImgFormat
  )
  SSD_Alpha <- compute_SSD(Image_Chunk = SS_Chunk, window_size = window_size,
                           nbclusters = nbclusters, MinSun = MinSun, pcelim = pcelim,
                           Index_Alpha = Index_Alpha)
  # write spectral Species ditribution file
  fidSSD <- file(
    description = SSD_Path, open = "r+b", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  if (!ReadWrite$RW_SSD$Byte_Start == 1) {
    seek(fidSSD, where = ReadWrite$RW_SSD$Byte_Start - 1, origin = "start", rw = "write")
  }
  SSD_Chunk <- aperm(array(SSD_Alpha$SSD, c(ReadWrite$RW_SSD$nbLines, HDR_SSD$samples, HDR_SSD$bands)), c(2, 3, 1))
  writeBin(c(SSD_Chunk), fidSSD, size = SSD_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
  close(fidSSD)

  # write PCsunlit pixels corresponding to SSD file
  fidSunlit <- file(
    description = Sunlit_Path, open = "r+b", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  if (!ReadWrite$RW_Sunlit$Byte_Start == 1) {
    seek(fidSunlit, where = ReadWrite$RW_Sunlit$Byte_Start - 1, origin = "start", rw = "write")
  }
  Sunlit.Chunk <- t(SSD_Alpha$PCsun)
  writeBin(c(Sunlit.Chunk), fidSunlit, size = Sunlit_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
  close(fidSunlit)

  rm(SSD_Chunk)
  rm(Sunlit.Chunk)
  gc()
  Shannon_Mean_Chunk <- apply(SSD_Alpha$Shannon, 1:2, mean)
  Fisher_Mean_Chunk <- apply(SSD_Alpha$Fisher, 1:2, mean)
  Simpson_Mean_Chunk <- apply(SSD_Alpha$Simpson, 1:2, mean)
  Shannon_SD_Chunk <- apply(SSD_Alpha$Shannon, 1:2, sd)
  Fisher_SD_Chunk <- apply(SSD_Alpha$Fisher, 1:2, sd)
  Simpson_SD_Chunk <- apply(SSD_Alpha$Simpson, 1:2, sd)
  rm(SSD_Alpha)
  gc()
  my_list <- list(
    "Shannon" = Shannon_Mean_Chunk, "Fisher" = Fisher_Mean_Chunk, "Simpson" = Simpson_Mean_Chunk,
    "Shannon.SD" = Shannon_SD_Chunk, "Fisher.SD" = Fisher_SD_Chunk, "Simpson.SD" = Simpson_SD_Chunk
  )
  return(my_list)
}

# compute spectral species distribution from original spectral species map
#
# @param Image_Chunk 3D image chunk of spectral species
# @param window_size size of spatial units (in pixels) to compute diversity
# @param nbclusters number of clusters defined in k-Means
# @param MinSun minimum proportion of sunlit pixels required to consider plot
# @param Index_Alpha list. list of alpha diversity indices to be computed from spectral species
# @param pcelim minimum proportion for a spectral species to be included
#
# @return list of alpha diversity metrics for each iteration
#' @importFrom vegan fisher.alpha
compute_SSD <- function(Image_Chunk, window_size, nbclusters, MinSun, pcelim, Index_Alpha = "Shannon") {
  nbi <- round(dim(Image_Chunk)[1] / window_size)
  nbj <- round(dim(Image_Chunk)[2] / window_size)
  nb_partitions <- dim(Image_Chunk)[3]
  SSDMap <- array(NA, c(nbi, nbj, nb_partitions * nbclusters))
  shannonIter <- FisherAlpha <- SimpsonAlpha <- array(NA, dim = c(nbi, nbj, nb_partitions))
  PCsun <- matrix(NA, nrow = nbi, ncol = nbj)

  # which spectral indices will be computed
  Shannon <- Simpson <- Fisher <- FALSE
  if (length((grep("Shannon", Index_Alpha))) > 0) Shannon <- TRUE
  if (length((grep("Simpson", Index_Alpha))) > 0) Simpson <- TRUE
  if (length((grep("Fisher", Index_Alpha))) > 0) Fisher <- TRUE

  # for each kernel in the line
  for (ii in 1:nbi) {
    # for each kernel in the column
    for (jj in 1:nbj) {
      li <- ((ii - 1) * window_size) + 1
      ui <- min(c(ii * window_size,dim(Image_Chunk)[1]))
      lj <- ((jj - 1) * window_size) + 1
      uj <- min(c(jj * window_size,dim(Image_Chunk)[2]))
      # put all iterations in a 2D matrix shape
      ijit <- t(matrix(Image_Chunk[li:ui, lj:uj, ], ncol = nb_partitions))
      # keep non zero values
      ijit <- matrix(ijit[, which(!ijit[1, ] == 0)], nrow = nb_partitions)
      nbPix_Sunlit <- dim(ijit)[2]
      PCsun[ii, jj] <- nbPix_Sunlit / window_size**2
      if (PCsun[ii, jj] > MinSun) {
        # for each iteration
        for (it in 1:nb_partitions) {
          lbk <- (it - 1) * nbclusters
          SSD <- as.vector(table(ijit[it, ]))
          ClusterID <- sort(unique(ijit[it, ]))
          if (pcelim > 0) {
            KeepSS <- which(SSD >= pcelim * nbPix_Sunlit)
            ClusterID <- ClusterID[KeepSS]
            SSD <- SSD[KeepSS]
          }
          SSDMap[ii, jj, (lbk + ClusterID)] <- SSD
          if (Shannon == TRUE) {
            shannonIter[ii, jj, it] <- get_Shannon(SSD)
          }
          if (Simpson == TRUE) {
            SimpsonAlpha[ii, jj, it] <- get_Simpson(SSD)
          }
          if (Fisher == TRUE) {
            if (length(SSD) > 2) {
              FisherAlpha[ii, jj, it] <- fisher.alpha(SSD)
            } else {
              FisherAlpha[ii, jj, it] <- 0
            }
          }
        }
      } else {
        shannonIter[ii, jj, ] <- NA
        FisherAlpha[ii, jj, ] <- NA
        SimpsonAlpha[ii, jj, ] <- NA
      }
    }
  }
  my_list <- list("Shannon" = shannonIter, "Fisher" = FisherAlpha, "Simpson" = SimpsonAlpha, "SSD" = SSDMap, "PCsun" = PCsun)
  return(my_list)
}

# computes shannon index from a distribution
#
# @param Distrib Distribution
#
# @return Shannon index correspnding to the distribution
get_Shannon <- function(Distrib) {
  Distrib <- Distrib / sum(Distrib, na.rm = TRUE)
  Distrib <- Distrib[which(!Distrib == 0)]
  shannon <- -1 * sum(Distrib * log(Distrib), na.rm = TRUE)
  return(shannon)
}

# computes Simpson index from a distribution
#
# @param Distrib Distribution
#
# @return Simpson index correspnding to the distribution
get_Simpson <- function(Distrib) {
  Distrib <- Distrib / sum(Distrib, na.rm = TRUE)
  Simpson <- 1 - sum(Distrib * Distrib, na.rm = TRUE)
  return(Simpson)
}
