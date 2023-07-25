# ==============================================================================
# biodivMapR
# Lib_MapBetaDiversity.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This Library computes bray curtis dissimilarity among spatial units based on
# spectral species distribution file and generates a RGB representation with NMDS
# ==============================================================================

#' maps beta diversity indicator based on spectral species distribution
#'
#' @param Input_Image_File character. Path and name of the image to be processed.
#' @param Output_Dir character. Output directory.
#' @param window_size numeric. Dimensions of the spatial unit.
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...).
#' @param nbclusters numeric. Number of clusters defined in k-Means.
#' @param Nb_Units_Ordin numeric. Maximum number of spatial units to be processed in NMDS.
#' --> 1000 will be fast but may not capture important patterns if large area
#' --> 4000 will be slow but may show better ability to capture landscape patterns
#' @param MinSun numeric. Minimum proportion of sunlit pixels required to consider plot.
#' @param pcelim numeric. Minimum contribution in percent required for a spectral species.
#' @param scaling character. PCO or NMDS
#' @param dimMDS numeric. number of dimensions for the scaling.
#' @param FullRes boolean.
#' @param LowRes boolean.
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param ClassifMap character. If FALSE, perform standard biodivMapR based on SpectralSpecies.
#'                              else corresponds to path for a classification map.
#'
#' @return PCoA_model
#' @export

map_beta_div <- function(Input_Image_File = FALSE,
                         Output_Dir = '',
                         window_size = 10,
                         TypePCA = 'SPCA',
                         nbclusters = 50,
                         Nb_Units_Ordin = 2000,
                         MinSun = 0.25, pcelim = 0.02, scaling = 'PCO', dimMDS = 3,
                         FullRes = FALSE, LowRes = TRUE,
                         nbCPU = 1, MaxRAM = 0.25, ClassifMap = FALSE) {

  # 1- get path for spectral species path, and possibly update Input_Image_File
  # and nbclusters if using classification map as input data
  SSDpathlist <- get_SSpath(Output_Dir, Input_Image_File, TypePCA, ClassifMap, nbclusters)
  Spectral_Species_Path <- SSDpathlist$Spectral_Species_Path
  Input_Image_File <- SSDpathlist$Input_Image_File
  nbclusters <- SSDpathlist$nbclusters
  SSD_Dir <- SSDpathlist$SSD_Dir
  HDR_SS <- read_ENVI_header(get_HDR_name(Spectral_Species_Path))

  # 2- compute beta diversity
  Beta <- compute_beta_metrics(ClusterMap_Path = Spectral_Species_Path,
                               SSD_Dir = SSD_Dir,
                               MinSun = MinSun,
                               Nb_Units_Ordin = Nb_Units_Ordin,
                               nb_partitions = HDR_SS$bands,
                               nbclusters = nbclusters, pcelim = pcelim,
                               scaling = scaling, dimMDS = dimMDS,
                               nbCPU = nbCPU, MaxRAM = MaxRAM)

  # Create images corresponding to Beta-diversity
  print("Write beta diversity maps")
  Index <- paste("BetaDiversity_BCdiss_", scaling, sep = "")
  Output_Dir_BETA <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, 'BETA')
  Beta.Path <- file.path(Output_Dir_BETA, paste(Index, "_", window_size, sep = ""))
  write_raster(Image = Beta$BetaDiversity, HDR = Beta$HDR, ImagePath = Beta.Path,
               window_size = window_size, FullRes = FullRes, LowRes = TRUE,
               SmoothImage = FALSE)

  PCoA_model <- Beta$PCoA_model
  return(PCoA_model)
}

#' #' computes NMDS
#' #
#' #' @param MatBCdist BC dissimilarity matrix
#' #' @param dimMDS numeric. number of dimensions of the NMDS
#' #
#' #' @return BetaNMDS_sel
#' #' @importFrom future plan multisession sequential
#' #' @importFrom future.apply future_lapply
#' #' @importFrom ecodist nmds
#' #' @importFrom utils find
#' #' @export
#'
#' compute_NMDS <- function(MatBCdist,dimMDS=3) {
#'   nbiterNMDS <- 4
#'   if (Sys.info()["sysname"] == "Windows") {
#'     nbCoresNMDS <- 2
#'   } else if (Sys.info()["sysname"] == "Linux") {
#'     nbCoresNMDS <- 4
#'   }
#'   # multiprocess of spectral species distribution and alpha diversity metrics
#'   # plan(multisession, workers = nbCoresNMDS) ## Parallelize using four cores
#'   plan(multisession, workers = nbCoresNMDS) ## Parallelize using four cores
#'   BetaNMDS <- future_lapply(MatBCdist,
#'                             FUN = nmds,
#'                             mindim = dimMDS, maxdim = dimMDS, nits = 1,
#'                             future.packages = c("ecodist"))
#'   plan(sequential)
#'   # find iteration with minimum stress
#'   Stress <- vector(length = nbiterNMDS)
#'   for (i in 1:nbiterNMDS) {
#'     Stress[i] <- BetaNMDS[[i]]$stress
#'   }
#'   print("Stress obtained for NMDS iterations:")
#'   print(Stress)
#'   print("Rule of thumb")
#'   print("stress < 0.05 provides an excellent represention in reduced dimensions")
#'   print("stress < 0.1 is great")
#'   print("stress < 0.2 is good")
#'   print("stress > 0.3 provides a poor representation")
#'   MinStress <- find(Stress == min(Stress))
#'   BetaNMDS_sel <- BetaNMDS[[MinStress]]$conf
#'   BetaNMDS_sel <- data.frame(BetaNMDS_sel[[1]])
#'   return(BetaNMDS_sel)
#' }

#' Identifies ordination coordinates based on nearest neighbors
#'
#' @param SSD_subset numeric. matrix corresponding to Spectral species distribution for a set of windows
#' @param Beta_Ordination_sel ordination of dissimilarity matrix for a selection of spatial units
#' @param Sample_Sel numeric. Samples selected during ordination
#' @param nb_partitions number of k-means then averaged
#' @param nbclusters number of clusters
#' @param pcelim numeric. Minimum contribution in percent required for a spectral species
#' @param nbCPU numeric. number of CPUs available
#' @param SamplesPerThread numeric. number of samples to be processed per thread for dissimilarity matrix
#'
#' @return Ordination_est coordinates of each spatial unit in ordination space
#' @import cli
#' @importFrom snow splitRows
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor handlers with_progress
#' @export

ordination_to_NN <- function(SSD_subset,
                             Beta_Ordination_sel,
                             Sample_Sel,
                             nb_partitions, nbclusters,
                             pcelim = 0.02, nbCPU = 1,
                             SamplesPerThread = 2000) {

  # define number of samples to be sampled each time during parallel processing
  # (set to 2000 to keep reasonable size of dissimilarity matrix)
  nb_subsets <- round(dim(SSD_subset)[1]/SamplesPerThread)
  if (nb_subsets == 0) nb_subsets <- 1
  SSD_subset <- snow::splitRows(SSD_subset, ncl = nb_subsets)
  # compute ordination coordinates from each SSD_subset
  if (nbCPU>1){
    plan(multisession, workers = nbCPU) ## Parallelize using four cores
    handlers(global = TRUE)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = nb_subsets)
      OutPut <- future_lapply(SSD_subset,
                              FUN = ordination_to_NN_list,
                              Sample_Sel = Sample_Sel,
                              Beta_Ordination_sel = Beta_Ordination_sel,
                              nb_partitions = nb_partitions,
                              nbclusters = nbclusters, pcelim = pcelim, p = p,
                              future.packages = c("vegan", "dissUtils", "tools",
                                                  "snow", "matlab"))
      })
    plan(sequential)
  } else {
    handlers(global = TRUE)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = nb_subsets)
      OutPut <- lapply(SSD_subset,
                       FUN = ordination_to_NN_list,
                       Sample_Sel = Sample_Sel,
                       Beta_Ordination_sel = Beta_Ordination_sel,
                       nb_partitions = nb_partitions,
                       nbclusters = nbclusters, pcelim = pcelim, p = p)
    })
  }
  Ordination_est <- do.call("rbind", OutPut)
  gc()
  return(Ordination_est)
}
#' computes nearest neighbors of a SSD subset with samples used for the adjustment
#' of the ordination, then gets the coordinates of each window in the ordination space
#
#' @param SSD_subset numeric. subset of spectral species distribution file
#' @param Sample_Sel numeric. Samples selected during ordination
#' @param Beta_Ordination_sel numeric. ordination of dissimilarity matrix for a selection of spatial units
#' @param nb_partitions numeric. Number of partitions (repetitions) to be computed then averaged.
#' @param nbclusters numeric. Number of clusters defined in k-Means.
#' @param pcelim numeric. Minimum contribution in percents required for a spectral species
#' @param p function.
#
#' @return OutPut list of mean and SD of alpha diversity metrics
#' @importFrom snow splitCols
#' @export

ordination_to_NN_list <- function(SSD_subset,
                                  Sample_Sel,
                                  Beta_Ordination_sel,
                                  nb_partitions,
                                  nbclusters,
                                  pcelim, p = NULL) {

  # compute the mean BC dissimilarity sequentially for each iteration
  SSD_subsub <- snow::splitCols(x = SSD_subset, ncl = nb_partitions)
  Sample_Sel2 <- snow::splitCols(x = Sample_Sel, ncl = nb_partitions)
  MatBCtmp <- list()
  for (i in 1:nb_partitions){
    MatBCtmp[[i]] <- list()
    MatBCtmp[[i]]$mat1 <- SSD_subsub[[i]]
    MatBCtmp[[i]]$mat2 <- Sample_Sel2[[i]]
  }
  MatBCtmp0 <- lapply(X = MatBCtmp, FUN = compute_BCdiss, pcelim)
  MatBCtmp <- Reduce('+', MatBCtmp0)/nb_partitions
  rm(MatBCtmp0)
  gc()
  # get the knn closest neighbors from each kernel
  knn <- 3
  OutPut <- compute_NN_from_ordination(MatBCtmp, knn, Beta_Ordination_sel)
  if (!is.null(p)){p()}
  return(OutPut)
}


#' compute beta diversity from spectral species computed for a plot
#' expecting a matrix of spectral species (n pixels x p repetitions)
#'
#' @param SpectralSpecies_Plots list. list of matrices of spectral species corresponding to plots
#' @param nbclusters numeric. number of clusters
#' @param Hellinger boolean. set TRUE to compute Hellinger distance matrices based on Euclidean distances
#' @param pcelim numeric. Minimum contribution (in %) required for a spectral species
#' each spectral species with a proprtion < pcelim is eliminated before computation of diversity
#
#' @return Mean bray curtis dissimilarity matrix for all plots, and individual BC matrices corresponding to each repetitions
#' @importFrom vegan vegdist
#' @export

compute_BETA_FromPlots <- function(SpectralSpecies_Plots,nbclusters,Hellinger = FALSE, pcelim = 0.02){

  nbPolygons<- length(SpectralSpecies_Plots)
  Pixel_Inventory_All <- Pixel_Hellinger_All <- list()
  nb_partitions <- dim(SpectralSpecies_Plots[[1]])[2]
  # for each plot
  for (plot in 1:nbPolygons){
    # for each repetition
    Pixel_Inventory <- Pixel_Hellinger <- list()
    for (i in 1:nb_partitions){
      # compute distribution of spectral species
      Distritab <- table(SpectralSpecies_Plots[[plot]][,i])
      # compute distribution of spectral species
      Pixel_Inventory[[i]] <- as.data.frame(Distritab)
      SumPix <- sum(Pixel_Inventory[[i]]$Freq)
      ThreshElim <- pcelim*SumPix
      ElimZeros <- which(Pixel_Inventory[[i]]$Freq<ThreshElim)
      if (length(ElimZeros)>=1){
        Pixel_Inventory[[i]] <- Pixel_Inventory[[i]][-ElimZeros,]
      }
      if (length(which(Pixel_Inventory[[i]]$Var1==0))==1){
        Pixel_Inventory[[i]] <- Pixel_Inventory[[i]][-which(Pixel_Inventory[[i]]$Var1==0),]
      }
    }
    Pixel_Inventory_All[[plot]] <- Pixel_Inventory
    if (Hellinger == TRUE){
      for (i in 1:nb_partitions){
        # compute Hellinger distance
        Pixel_Hellinger[[i]] <- Pixel_Inventory[[i]]
        Pixel_Hellinger[[i]]$Freq <- sqrt(Pixel_Hellinger[[i]]$Freq/sum(Pixel_Hellinger[[i]]$Freq))
      }
      Pixel_Hellinger_All[[plot]] <- Pixel_Hellinger
    }
  }

  # for each pair of plot, compute beta diversity indices
  BC <- list()
  for(i in 1:nb_partitions){
    MergeDiversity <- matrix(0,nrow = nbclusters,ncol = nbPolygons)
    for(j in 1:nbPolygons){
      SelSpectralSpecies <- as.numeric(as.vector(Pixel_Inventory_All[[j]][[i]]$Var1))
      SelFrequency <- Pixel_Inventory_All[[j]][[i]]$Freq
      MergeDiversity[SelSpectralSpecies,j] = SelFrequency
    }
    BC[[i]] <- vegan::vegdist(t(MergeDiversity),method="bray")
  }
  BC_mean <- 0*BC[[1]]
  for(i in 1:nb_partitions){
    BC_mean <- BC_mean+BC[[i]]
  }
  BC_mean <- BC_mean/nb_partitions

  # Hellinger
  if (Hellinger==TRUE){
    # for each pair of plot, compute Euclidean distance on Hellinger
    Hellmat <- list()
    for(i in 1:nb_partitions){
      MergeDiversity <- matrix(0,nrow = nbclusters,ncol = nbPolygons)
      for(j in 1:nbPolygons){
        SelSpectralSpecies <- as.numeric(as.vector(Pixel_Hellinger_All[[j]][[i]]$Var1))
        SelFrequency <- Pixel_Hellinger_All[[j]][[i]]$Freq
        MergeDiversity[SelSpectralSpecies,j] = SelFrequency
      }
      Hellmat[[i]] <- vegan::vegdist(t(MergeDiversity),method="euclidean")
    }
    Hellinger_mean <- 0*Hellmat[[1]]
    for(i in 1:nb_partitions){
      Hellinger_mean <- Hellinger_mean+Hellmat[[i]]
    }
    Hellinger_mean <- Hellinger_mean/nb_partitions
  }
  beta <- list('BrayCurtis' = BC_mean, 'BrayCurtis_ALL' = BC,
               'Hellinger' = Hellinger_mean, 'Hellinger_ALL' = Hellmat)
  return(beta)
}

#' computes beta diversity metrics
#'
#' @param ClusterMap_Path character. File containing spectral species or classes from prior classification
#' @param SSD_Dir character. path for spectral species distribution file to be written
#' @param MinSun numeric. minimum proportion of sunlit pixels required to consider plot
#' @param Nb_Units_Ordin numeric. maximum number of spatial units to be processed in Ordination
#' @param nb_partitions numeric. number of iterations
#' @param nbclusters numeric. number of clusters defined in k-Means
#' @param pcelim numeric. Minimum contribution in percents required for a spectral species
#' @param scaling character. scaling method, PCO or NMDS
#' @param dimMDS numeric. number of dimensions for the scaling.
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param SamplesPerThread numeric. number of samples to be processed per thread for dissimilarity matrix
#'
#' @return my_list a list of information including results of 3 dimensional ordination with PCoA, PCoA model, map of sunlit windows
#' @importFrom labdsv pco
#' @importFrom stats as.dist
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @export

compute_beta_metrics <- function(ClusterMap_Path,
                                 SSD_Dir,
                                 MinSun,
                                 Nb_Units_Ordin,
                                 nb_partitions,
                                 nbclusters, pcelim,
                                 scaling = 'PCO', dimMDS=3,
                                 nbCPU = 1, MaxRAM = 0.25,
                                 SamplesPerThread = 2000) {

  # 1- perform analysis for a subset of windows
  # Define path for images to be used
  SSD_Path <- file.path(SSD_Dir, 'SpectralSpecies_Distribution')
  ImPathSunlit <- file.path(SSD_Dir, 'SpectralSpecies_Distribution_Sunlit')

  # Get illuminated pixels based on  SpectralSpecies_Distribution_Sunlit
  Sunlit_Pixels <- get_sunlit_pixels(ImPathSunlit, MinSun)
  Select_Sunlit <- Sunlit_Pixels$Select_Sunlit
  nb_Sunlit <- length(Select_Sunlit)

  # Define spatial units processed through ordination and those processed through
  # Nearest neighbor based on the first ordination batch
  print("subset Spectral Species distribution")
  RandPermKernels <- sample(seq(1, nb_Sunlit, by = 1))
  if (Nb_Units_Ordin <= nb_Sunlit) {
    Kernels_NN <- RandPermKernels[(Nb_Units_Ordin + 1):nb_Sunlit]
  } else {
    Nb_Units_Ordin <- nb_Sunlit
    Kernels_NN <- c()
  }

  # read samples from spectral species distribution file
  Kernels_Ordination <- RandPermKernels[1:Nb_Units_Ordin]
  Sample_Sel <- extract_samples_from_image(SSD_Path, Sunlit_Pixels$coordTotSort[Kernels_Ordination,])
  gc()

  # create a Bray curtis dissimilarity matrix for each iteration
  print("compute BC dissimilarity for selected kernels")
  Sample_Sel_list <- snow::splitCols(x = Sample_Sel,ncl = nb_partitions)
  plan(multisession, workers = nbCPU)
  handlers(global = TRUE)
  handlers("cli")
  with_progress({
    p <- progressr::progressor(steps = nb_partitions)
    MatBCtmp0 <- future.apply::future_lapply(Sample_Sel_list,
                                             FUN = getBCdiss,
                                             pcelim = pcelim, p = p)
  })
  plan(sequential)
  MatBC <- Reduce('+', MatBCtmp0)/nb_partitions
  rm(MatBCtmp0)
  gc()

  # Perform Ordination based on BC dissimilarity matrix
  print("perform Ordination on the BC dissimilarity matrix averaged from all iterations")
  # parallel computing of Ordination can be run on 2 cores on Windows.
  # core management seems better on linux --> 4 cores possible
  MatBCdist <- as.dist(MatBC, diag = FALSE, upper = FALSE)
  BetaPCO <- NULL
  # if (scaling == 'NMDS') {
  #   Beta_Ordination_sel <- compute_NMDS(MatBCdist)
  #   PCname <- 'NMDS'
  # } else if (scaling == 'PCO') {
    BetaPCO <- labdsv::pco(MatBCdist, k = dimMDS)
    Beta_Ordination_sel <- BetaPCO$points
    PCname <- 'PCoA'
  # }

  # 2- perform analysis on all the image, using nearest neighbor between
  # subset previously analysed and spatial units excluded from Ordination
  print("BC dissimilarity between samples selected for Ordination and remaining")
  # first define in how many pieces the image will be split
  SSD_HDR <- get_HDR_name(SSD_Path)
  HDR_SSD <- read_ENVI_header(SSD_HDR)
  nbPieces <- split_image(HDR_SSD, MaxRAM)

  coordTotSort <- snow::splitRows(x = Sunlit_Pixels$coordTotSort, ncl = nbPieces)
  Ordination_est <- list()
  for (i in 1:nbPieces){
    message(paste('Computing beta diversity for image subset #',i,' / ',nbPieces))
    # read SSD raster
    SSD_subset <- extract_samples_from_image(SSD_Path, coordTotSort[[i]])
    # perform ordination for the subset
    Ordination_est[[i]] <- ordination_to_NN(SSD_subset = SSD_subset,
                                            Beta_Ordination_sel = Beta_Ordination_sel,
                                            Sample_Sel = Sample_Sel,
                                            nb_partitions = nb_partitions, nbclusters = nbclusters,
                                            pcelim = pcelim, nbCPU = nbCPU,
                                            SamplesPerThread = SamplesPerThread)
  }
  Ordination_est  <- do.call('rbind', Ordination_est)

  # Reconstuct spatialized beta diversity map from previous computation
  Sunlit_HDR <- get_HDR_name(ImPathSunlit)
  HDR_Sunlit <- read_ENVI_header(Sunlit_HDR)
  BetaDiversity <- as.matrix(Ordination_est, ncol = dimMDS)
  BetaDiversityRGB <- array(NA, c(as.double(HDR_Sunlit$lines), as.double(HDR_Sunlit$samples), dimMDS))
  BetaTmp <- matrix(NA, nrow = as.double(HDR_Sunlit$lines), ncol = as.double(HDR_Sunlit$samples))
  for (i in 1:dimMDS) {
    BetaTmp[Select_Sunlit] <- BetaDiversity[, i]
    BetaDiversityRGB[, , i] <- BetaTmp
  }
  listvar <- ls()
  # update HDR file
  HDR_Beta <- HDR_Sunlit
  HDR_Beta$bands <- dimMDS
  HDR_Beta$`data type` <- 4
  PCs <- list()
  for (i in 1:dimMDS) {
    PCs <- c(PCs, paste(PCname,'#', i,sep = ''))
  }
  PCs <- paste(PCs, collapse = ', ')
  HDR_Beta$`band names` <- PCs

  rm(list = listvar[-which(listvar == 'BetaDiversityRGB' | listvar == 'Select_Sunlit' |
                             listvar == 'HDR_Beta' | listvar == 'BetaPCO')])
  gc()
  my_list <- list('BetaDiversity' = BetaDiversityRGB, 'Select_Sunlit' = Select_Sunlit,
                  'PCoA_model' = BetaPCO, 'HDR' = HDR_Beta)
  return(my_list)
}

#' compute the bray curtis dissimilarity matrix corresponding to a list of kernels
#' (rows) defined by their spectral species (columns)
#' SSDList is a list containing spectral species distribution for two sets of kernels ([[1]] and [[2]])
#' pcelim is the threshold for minimum contributin of a spctral species to be kept
#
#' @param SSDList list. list of 2 groups to compute BC dissimilarity from
#' @param pcelim numeric. minimum proportion required for a species to be included (currently deactivated)
#
#' @return MatBC matrix of bray curtis dissimilarity matrix
#' @importFrom dissUtils diss
#' @export

compute_BCdiss <- function(SSDList, pcelim=0.02) {
  # compute the proportion of each spectral species
  # Here, the proportion is computed with regards to the total number of sunlit pixels
  # One may want to determine if the results are similar when the proportion is computed
  # with regards to the total number of pixels (se*se)
  # however it would increase dissimilarity between kernels with different number of sunlit pixels

  # SSD <- list()
  # for (i in 1:length(SSDList)) {
  #   # get the total number of sunlit pixels in spatial unit
  #   SumSpecies <- rowSums(SSDList[[i]])
  #   elim <- which(SumSpecies == 0)
  #   if (length(elim) > 0) {
  #     SumSpecies[elim] <- 1
  #     SSDList[[i]][elim, ] <- 0
  #   }
  #   SSD[[i]] <- apply(SSDList[[i]], 2, function(x, c1) x / c1, 'c1' = SumSpecies)
  #   SSD[[i]][which(SSD[[i]] < pcelim)] <- 0
  # }

  SSD <- lapply(SSDList,FUN = normalize_SSD, pcelim = pcelim)
  # matrix of bray curtis dissimilarity (size = nb kernels x nb kernels)
  # Here use the package "dissUtils" to compute dissimilarity matrices sequentially
  MatBC <- dissUtils::diss(SSD[[1]], SSD[[2]], method = 'braycurtis')
  # EDIT 06-Feb-2019: added this to fix problem when empty kernels occur, leading to NA BC value
  BCNA <- which(is.na(MatBC) == TRUE)
  if (length(BCNA) > 0) {
    MatBC[BCNA] <- 1
  }
  return(MatBC)
}


#' Compute the weighted coordinates of a spatial unit based on nearest neighbors used during PCoA
#
#' @param NN list. coordinates and ID of nearest neighbors from BetaDiversity0 matrix
#' @param knn number of neighbors
#' @param BetaDiversity0 PCoA coordinates of reference samples
#
#' @return estimated NMDS position based on nearest neighbors from NMDS
#' @export

WeightedCoordsNN <- function(NN, knn, BetaDiversity0) {

  # get distance and ID of NN samples
  DistNN <- NN$x[1:knn]
  IdNN <- NN$ix[1:knn]
  # final location weighted by location of NN
  # if exact same location as nearest neighbor
  if (DistNN[1]==0){
    MDSpos <- BetaDiversity0[IdNN[1], ]
  } else {
    # total dissimilarity from k nearest neighbors to weight
    Dist_Tot <- 1/sum(1/DistNN)
    if (ncol(BetaDiversity0)>1){
      MDSpos <- colSums((Dist_Tot/DistNN) * BetaDiversity0[IdNN, ])
    } else {
      MDSpos <- sum((Dist_Tot/DistNN) * BetaDiversity0[IdNN, ])
    }
  }
  Ordin_est <- MDSpos
  return(Ordin_est)
}

#' compute the nearest neighbors among kernels used in NMDS
#
#' @param MatBC3 matrix of BC dissimilarity between the kernels excluded from Ordination (rows)
#' @param knn numeric. number of neighbors
#' @param BetaDiversity0 numeric. matrix of ordinated coordinates computed from dissimilarity matrix
#
#' @return Ordin_est estimated NMDS position based on nearest neighbors from NMDS
#' @export

compute_NN_from_ordination <- function(MatBC3, knn, BetaDiversity0) {

  dimMDS <- ncol(BetaDiversity0)
  # get nearest neighbors (coordinates and ID)
  NN <- apply(MatBC3, 1, FUN = sort,index.return = TRUE)
  Ordin_est <- lapply(X = NN,
                      FUN = WeightedCoordsNN,
                      knn = knn, BetaDiversity0 = BetaDiversity0)
  Ordin_est <- do.call(rbind,Ordin_est)
  return(Ordin_est)
}

#' Gets sunlit pixels from SpectralSpecies_Distribution_Sunlit
#'
#' @param ImPathSunlit path for SpectralSpecies_Distribution_Sunlit
#' @param MinSun minimum proportion of sunlit pixels required
#'
#' @return list of sunlit pixels
#' @export

get_sunlit_pixels <- function(ImPathSunlit, MinSun) {

  # Filter out spatial units showing poor illumination
  Sunlit_HDR <- get_HDR_name(ImPathSunlit)
  HDR_Sunlit <- read_ENVI_header(Sunlit_HDR)
  nbpix <- as.double(HDR_Sunlit$lines) * as.double(HDR_Sunlit$samples)
  fid <- file(
    description = ImPathSunlit, open = 'rb', blocking = TRUE,
    encoding = getOption('encoding'), raw = FALSE
  )
  Sunlit <- readBin(fid, numeric(), n = nbpix, size = 4)
  close(fid)
  Sunlit <- aperm(array(Sunlit, dim = c(HDR_Sunlit$samples, HDR_Sunlit$lines)))
  # define sunlit spatial units
  Select_Sunlit <- which(Sunlit > MinSun)
  # define where to extract each datapoint in the file
  coordi <- ((Select_Sunlit - 1) %% HDR_Sunlit$lines) + 1
  coordj <- floor((Select_Sunlit - 1) / HDR_Sunlit$lines) + 1
  # sort based on line and col (important for optimal scan of large files)
  coordTot <- cbind(coordi, coordj)
  # sort samples from top to bottom in order to optimize read/write of the image
  # indxTot saves the order of the data for reconstruction later
  indxTot <- order(coordTot[, 1], coordTot[, 2], na.last = FALSE)
  coordTotSort <- coordTot[indxTot, ]
  Select_Sunlit <- Select_Sunlit[indxTot]

  my_list <- list('Select_Sunlit' = Select_Sunlit, 'coordTotSort' = coordTotSort)
  return(my_list)
}

#' Computes BC dissimilarity for a given matrix
#'
#' @param Mat1 numeric. matrix of spectral species distribution
#' @param pcelim numeric. Minimum contribution in percent required for a spectral species.
#' @param p list. progressor object for progress bar
#'
#' @return BCtmp numeric. BC dissimilarity matrix corresponding to Mat1 and Mat2
#' @export
getBCdiss <- function(Mat1, pcelim = 0.02, p = NULL){
  SSDList <- list()
  SSDList[[1]] <- SSDList[[2]] <- Mat1
  BCtmp <- compute_BCdiss(SSDList, pcelim)
  #if progress bar called
  if (!is.null(p)){p()}
  return(BCtmp)
}

#' Performs normalization of spectral species distribution based on cumulative
#' sum of spectral species
#'
#' @param SSDList numeric. matrix corresponding to element of a list of SSD
#' @param pcelim numeric. Minimum contribution in percent required for a spectral species.
#'
#' @return SSD numeric. matrix corresponding to normalized spectral species distribution
#' @export
normalize_SSD <- function(SSDList, pcelim = 0.02){

  SumSpecies <- rowSums(SSDList)
  elim <- which(SumSpecies == 0)
  if (length(elim) > 0) {
    SumSpecies[elim] <- 1
    SSDList[elim, ] <- 0
  }
  SSD <- apply(SSDList, 2, function(x, c1) x / c1, 'c1' = SumSpecies)
  SSD[which(SSD < pcelim)] <- 0
  return(SSD)
}

#' #' Computes BC dissimilarity between distributions defined by a subset of columns
#' #' of two matrices
#' #'
#' #' @param lub list. lower and upper values of the columns used in the computation
#' #' @param Mat1 numeric. matrix of spectral species distribution #1
#' #' @param Mat2 numeric. matrix of spectral species distribution #2
#' #' @param pcelim numeric. Minimum contribution in percent required for a spectral species.
#' #' @param p list. progressor object for progress bar
#' #'
#' #' @return BCtmp numeric. BC dissimilarity matrix corresponding to Mat1 and Mat2
#' #' @export
#' getBCdiss <- function(lub, Mat1, Mat2, pcelim = 0.02, p = NULL){
#'   SSDList <- list()
#'   SSDList[[1]] <- Mat1[, lub$lb:lub$ub]
#'   SSDList[[2]] <- Mat2[, lub$lb:lub$ub]
#'   BCtmp <- compute_BCdiss(SSDList, pcelim)
#'   #if progress bar called
#'   if (!is.null(p)){p()}
#'   return(BCtmp)
#' }

