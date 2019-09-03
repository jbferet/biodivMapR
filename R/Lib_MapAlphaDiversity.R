# ==============================================================================
# biodivMapR
# Lib_MapAlphaDiversity.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
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
#' @param FullRes boolean. Full resolution.
#' @param LowRes boolean. Low resolution.
#' @param MapSTD boolean. map of standard deviation of the alpha diversity map (over repetitions)
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param Index_Alpha character. Either 'Shannon', 'Simpson' or 'Fisher'.
#'
#' @export
map_alpha_div <- function(Input_Image_File, Output_Dir, window_size,
                                TypePCA = "SPCA", nbclusters = 50,
                                MinSun = 0.25, pcelim = 0.02,
                                Index_Alpha = "Shannon", FullRes = TRUE,
                                LowRes = FALSE, MapSTD = FALSE,
                                nbCPU = FALSE, MaxRAM = FALSE) {
  Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "SpectralSpecies")
  Output_Dir_PCA <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "PCA")
  Spectral_Species_Path <- paste(Output_Dir_SS, "SpectralSpecies", sep = "")
  # 1- COMPUTE ALPHA DIVERSITY
  ALPHA <- compute_alpha_metrics(Spectral_Species_Path, window_size, nbclusters, MinSun, pcelim, nbCPU = nbCPU, MaxRAM = MaxRAM, Index_Alpha = Index_Alpha)
  # 2- SAVE ALPHA DIVERSITY MAPS
  print("Write alpha diversity maps")
  # which spectral indices will be computed
  Shannon <- Simpson <- Fisher <- FALSE
  if (length((grep("Shannon", Index_Alpha))) > 0) Shannon <- TRUE
  if (length((grep("Simpson", Index_Alpha))) > 0) Simpson <- TRUE
  if (length((grep("Fisher", Index_Alpha))) > 0) Fisher <- TRUE

  Output_Dir_Alpha <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "ALPHA")
  if (Shannon == TRUE) {
    Index <- "Shannon"
    Alpha_Path <- paste(Output_Dir_Alpha, Index, "_", window_size, sep = "")
    write_raster_alpha(ALPHA$Shannon, ALPHA$HDR, Alpha_Path, window_size, Index, FullRes = FullRes, LowRes = LowRes)
    if (MapSTD == TRUE) {
      Index <- "Shannon_SD"
      Alpha_Path <- paste(Output_Dir_Alpha, Index, "_", window_size, sep = "")
      write_raster_alpha(ALPHA$Shannon.SD, ALPHA$HDR, Alpha_Path, window_size, Index, FullRes = FullRes, LowRes = LowRes)
    }
  }

  if (Fisher == TRUE) {
    Index <- "Fisher"
    Alpha_Path <- paste(Output_Dir_Alpha, Index, "_", window_size, sep = "")
    write_raster_alpha(ALPHA$Fisher, ALPHA$HDR, Alpha_Path, window_size, Index, FullRes = FullRes, LowRes = LowRes)
    if (MapSTD == TRUE) {
      Index <- "Fisher_SD"
      Alpha_Path <- paste(Output_Dir_Alpha, Index, "_", window_size, sep = "")
      write_raster_alpha(ALPHA$Fisher.SD, ALPHA$HDR, Alpha_Path, window_size, Index, FullRes = FullRes, LowRes = LowRes)
    }
  }

  if (Simpson == TRUE) {
    Index <- "Simpson"
    Alpha_Path <- paste(Output_Dir_Alpha, Index, "_", window_size, sep = "")
    write_raster_alpha(ALPHA$Simpson, ALPHA$HDR, Alpha_Path, window_size, Index, FullRes = FullRes, LowRes = LowRes)
    if (MapSTD == TRUE) {
      Index <- "Simpson_SD"
      Alpha_Path <- paste(Output_Dir_Alpha, Index, "_", window_size, sep = "")
      write_raster_alpha(ALPHA$Simpson.SD, ALPHA$HDR, Alpha_Path, window_size, Index, FullRes = FullRes, LowRes = LowRes)
    }
  }
  return()
}

# Map alpha diversity metrics based on spectral species
#
# @param Spectral_Species_Path path for spectral species file to be written
# @param window_size size of spatial units (in pixels) to compute diversity
# @param nbclusters number of clusters defined in k-Means
# @param pcelim
# @param nbCPU
# @param MaxRAM
# @param Index_Alpha
# @param MinSun minimum proportion of sunlit pixels required to consider plot
#
# @return list of mean and SD of alpha diversity metrics
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats sd
compute_alpha_metrics <- function(Spectral_Species_Path, window_size, nbclusters, MinSun, pcelim, nbCPU = FALSE, MaxRAM = FALSE, Index_Alpha = "Shannon") {
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
  HDR_SSD$samples <- floor(HDR_SS$samples / window_size)
  HDR_SSD$lines <- floor(HDR_SS$lines / window_size)
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
    ReadWrite[[i]]$RW.SS <- ReadWrite[[i]]$RW.SSD <- ReadWrite[[i]]$RW.Sunlit <- list()
    ReadWrite[[i]]$RW.SS$Byte_Start <- SeqRead.SS$ReadByte_Start[i]
    ReadWrite[[i]]$RW.SS$nbLines <- SeqRead.SS$Lines_Per_Chunk[i]
    ReadWrite[[i]]$RW.SS$lenBin <- SeqRead.SS$ReadByte_End[i] - SeqRead.SS$ReadByte_Start[i] + 1

    ReadWrite[[i]]$RW.SSD$Byte_Start <- SeqWrite.SSD$ReadByte_Start[i]
    ReadWrite[[i]]$RW.SSD$nbLines <- SeqWrite.SSD$Lines_Per_Chunk[i]
    ReadWrite[[i]]$RW.SSD$lenBin <- SeqWrite.SSD$ReadByte_End[i] - SeqWrite.SSD$ReadByte_Start[i] + 1

    ReadWrite[[i]]$RW.Sunlit$Byte_Start <- SeqWrite.Sunlit$ReadByte_Start[i]
    ReadWrite[[i]]$RW.Sunlit$nbLines <- SeqWrite.Sunlit$Lines_Per_Chunk[i]
    ReadWrite[[i]]$RW.Sunlit$lenBin <- SeqWrite.Sunlit$ReadByte_End[i] - SeqWrite.Sunlit$ReadByte_Start[i] + 1
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
  my_list <- list(
    "Shannon" = Shannon.Mean, "Fisher" = Fisher.Mean, "Simpson" = Simpson.Mean,
    "Shannon.SD" = Shannon.SD, "Fisher.SD" = Fisher.SD, "Simpson.SD" = Simpson.SD, "HDR" = HDR_SSD
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
  SS_Chunk <- read_image_subset(
    Spectral_Species_Path, HDR_SS,
    ReadWrite$RW.SS$Byte_Start, ReadWrite$RW.SS$lenBin,
    ReadWrite$RW.SS$nbLines, SS_Format, ImgFormat
  )
  SSD_Alpha <- compute_SSD(SS_Chunk, window_size, nbclusters, MinSun, pcelim, Index_Alpha = Index_Alpha)
  # write spectral Species ditribution file
  fidSSD <- file(
    description = SSD_Path, open = "r+b", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  if (!ReadWrite$RW.SSD$Byte_Start == 1) {
    seek(fidSSD, where = ReadWrite$RW.SSD$Byte_Start - 1, origin = "start", rw = "write")
  }
  SSD_Chunk <- aperm(array(SSD_Alpha$SSD, c(ReadWrite$RW.SSD$nbLines, HDR_SSD$samples, HDR_SSD$bands)), c(2, 3, 1))
  writeBin(c(SSD_Chunk), fidSSD, size = SSD_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
  close(fidSSD)

  # write PCsunlit pixels corresponding to SSD file
  fidSunlit <- file(
    description = Sunlit_Path, open = "r+b", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  if (!ReadWrite$RW.Sunlit$Byte_Start == 1) {
    seek(fidSunlit, where = ReadWrite$RW.Sunlit$Byte_Start - 1, origin = "start", rw = "write")
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
# @param Index_Alpha
# @param pcelim minimum proportion for a spectral species to be included
#
# @return list of alpha diversity metrics for each iteration
#' @importFrom vegan fisher.alpha
compute_SSD <- function(Image_Chunk, window_size, nbclusters, MinSun, pcelim, Index_Alpha = "Shannon") {
  nbi <- floor(dim(Image_Chunk)[1] / window_size)
  nbj <- floor(dim(Image_Chunk)[2] / window_size)
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
      ui <- ii * window_size
      lj <- ((jj - 1) * window_size) + 1
      uj <- jj * window_size
      # put all iterations in a 2D matrix shape
      ijit <- t(matrix(Image_Chunk[li:ui, lj:uj, ], ncol = nb_partitions))
      # keep non zero values
      ijit <- matrix(ijit[, which(!ijit[1, ] == 0)], nrow = nb_partitions)
      nb.Pix.Sunlit <- dim(ijit)[2]
      PCsun[ii, jj] <- nb.Pix.Sunlit / window_size**2
      if (PCsun[ii, jj] > MinSun) {
        # for each iteration
        for (it in 1:nb_partitions) {
          lbk <- (it - 1) * nbclusters
          SSD <- as.vector(table(ijit[it, ]))
          ClusterID <- sort(unique(ijit[it, ]))
          if (pcelim > 0) {
            KeepSS <- which(SSD >= pcelim * nb.Pix.Sunlit)
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

# Writes image of alpha diversity indicator (1 band) and smoothed alpha diversity
#
# @param Image 2D matrix of image to be written
# @param HDR_SSD hdr template derived from SSD to modify
# @param ImagePath path of image file to be written
# @param window_size spatial units use dto compute diversiy (in pixels)
# @param Index name of the index (eg. Shannon)
# @param FullRes should full resolution image be written (original pixel size)
# @param LowRes should low resolution image be written (one value per spatial unit)
#
# @return
write_raster_alpha <- function(Image, HDR_SSD, ImagePath, window_size, Index, FullRes = TRUE, LowRes = FALSE) {

  # Write image with resolution corresponding to window_size
  HDR_Alpha <- HDR_SSD
  HDR_Alpha$bands <- 1
  HDR_Alpha$`data type` <- 4
  HDR_Alpha$`band names` <- Index
  Image_Format <- ENVI_type2bytes(HDR_Alpha)
  if (LowRes == TRUE) {
    headerFpath <- paste(ImagePath, ".hdr", sep = "")
    write_ENVI_header(HDR_Alpha, headerFpath)
    ImgWrite <- array(Image, c(HDR_Alpha$lines, HDR_Alpha$samples, 1))
    ImgWrite <- aperm(ImgWrite, c(2, 3, 1))
    fidOUT <- file(
      description = ImagePath, open = "wb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    writeBin(c(ImgWrite), fidOUT, size = Image_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
  }
  if (FullRes == TRUE) {
    # Write image with Full native resolution
    HDR_Full <- HDR_Alpha
    HDR_Full$samples <- HDR_Alpha$samples * window_size
    HDR_Full$lines <- HDR_Alpha$lines * window_size
    HDR_Full <- revert_resolution_HDR(HDR_Full, window_size)
    ImagePath_FullRes <- paste(ImagePath, "_Fullres", sep = "")
    headerFpath <- paste(ImagePath_FullRes, ".hdr", sep = "")
    write_ENVI_header(HDR_Full, headerFpath)
    Image_Format <- ENVI_type2bytes(HDR_Full)
    Image_FullRes <- matrix(NA, ncol = HDR_Full$samples, nrow = HDR_Full$lines)
    for (i in 1:HDR_SSD$lines) {
      for (j in 1:HDR_SSD$samples) {
        Image_FullRes[((i - 1) * window_size + 1):(i * window_size), ((j - 1) * window_size + 1):(j * window_size)] <- Image[i, j]
      }
    }
    ImgWrite <- array(Image_FullRes, c(HDR_Full$lines, HDR_Full$samples, 1))
    ImgWrite <- aperm(ImgWrite, c(2, 3, 1))
    fidOUT <- file(
      description = ImagePath_FullRes, open = "wb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    writeBin(c(ImgWrite), fidOUT, size = Image_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
    # zip resulting file
    ZipFile(ImagePath_FullRes)
  }
  # write smoothed map
  SizeFilt <- 1
  nbi <- dim(Image)[1]
  nbj <- dim(Image)[2]
  if (min(c(nbi, nbj)) > (2 * SizeFilt + 1)) {
    Image_Smooth <- mean_filter(Image, nbi, nbj, SizeFilt)
    Image_Smooth[which(is.na(Image))] <- NA
    ImagePath.Smooth <- paste(ImagePath, "_MeanFilter", sep = "")
    headerFpath <- paste(ImagePath.Smooth, ".hdr", sep = "")
    Image_Format <- ENVI_type2bytes(HDR_Alpha)
    if (LowRes == TRUE) {
      write_ENVI_header(HDR_Alpha, headerFpath)
      ImgWrite <- array(Image_Smooth, c(HDR_Alpha$lines, HDR_Alpha$samples, 1))
      ImgWrite <- aperm(ImgWrite, c(2, 3, 1))
      fidOUT <- file(
        description = ImagePath.Smooth, open = "wb", blocking = TRUE,
        encoding = getOption("encoding"), raw = FALSE
      )
      writeBin(c(ImgWrite), fidOUT, size = Image_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
      close(fidOUT)
    }
    if (FullRes == TRUE) {
      # Write image with Full native resolution
      ImagePath_FullRes <- paste(ImagePath.Smooth, "_Fullres", sep = "")
      headerFpath <- paste(ImagePath_FullRes, ".hdr", sep = "")
      write_ENVI_header(HDR_Full, headerFpath)
      Image_Format <- ENVI_type2bytes(HDR_Full)
      Image_FullRes <- matrix(NA, ncol = HDR_Full$samples, nrow = HDR_Full$lines)
      for (i in 1:HDR_SSD$lines) {
        for (j in 1:HDR_SSD$samples) {
          Image_FullRes[((i - 1) * window_size + 1):(i * window_size), ((j - 1) * window_size + 1):(j * window_size)] <- Image_Smooth[i, j]
        }
      }
      ImgWrite <- array(Image_FullRes, c(HDR_Full$lines, HDR_Full$samples, 1))
      ImgWrite <- aperm(ImgWrite, c(2, 3, 1))
      fidOUT <- file(
        description = ImagePath_FullRes, open = "wb", blocking = TRUE,
        encoding = getOption("encoding"), raw = FALSE
      )
      writeBin(c(ImgWrite), fidOUT, size = Image_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
      close(fidOUT)
      # zip resulting file
      ZipFile(ImagePath_FullRes)
    }
  }
  return("")
}
