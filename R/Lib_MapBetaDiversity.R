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
#' @param nb_partitions numeric. Number of partitions (repetitions) to be computed then averaged.
#' @param nbclusters numeric. Number of clusters defined in k-Means.
#' @param Nb_Units_Ordin numeric. Maximum number of spatial units to be processed in NMDS.
#' --> 1000 will be fast but may not capture important patterns if large area
#' --> 4000 will be slow but may show better ability to capture landscape patterns
#' @param MinSun numeric. Minimum proportion of sunlit pixels required to consider plot.
#' @param pcelim numeric. Minimum contribution in percent required for a spectral species.
#' @param scaling character. PCO or NMDS
#' @param FullRes boolean.
#' @param LowRes boolean.
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param ClassifMap character. If FALSE, perform standard biodivMapR based on SpectralSpecies.
#'                              else corresponds to path for a classification map.
#' @param dimMDS numeric. number of dimensions for the scaling.
#'
#' @return PCoA_model
#' @export

map_beta_div <- function(Input_Image_File=FALSE, Output_Dir='', window_size=10,
                         TypePCA = 'SPCA', nb_partitions = 20,nbclusters = 50,
                         Nb_Units_Ordin = 2000, MinSun = 0.25,
                         pcelim = 0.00, scaling = 'PCO', FullRes = TRUE,
                         LowRes = FALSE, nbCPU = 1, MaxRAM = 0.25,
                         ClassifMap = FALSE, dimMDS=3) {

  if (ClassifMap == FALSE){
    Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, 'SpectralSpecies')
    Spectral_Species_Path <- file.path(Output_Dir_SS, 'SpectralSpecies')
  } else {
    message("Classification Map will be used instead of SpectralSpecies")
    message("Classes are expected to be integer values")
    message("conversion to ENVI File if not the format of the original classification map")
    message("Class '0' will be considered as No Data and excluded")
    if (! file.exists(ClassifMap)){
      message("classification map is not found:")
      print(ClassifMap)
    } else {
      Spectral_Species_Path <- ClassifMap
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
  Output_Dir_BETA <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, 'BETA')

  Beta <- compute_beta_metrics(ClusterMap_Path = Spectral_Species_Path, MinSun = MinSun,
                               Nb_Units_Ordin = Nb_Units_Ordin, nb_partitions = nb_partitions,
                               nbclusters = nbclusters, pcelim = pcelim, scaling = scaling,
                               nbCPU = nbCPU, MaxRAM = MaxRAM, dimMDS=dimMDS)
  # Create images corresponding to Beta-diversity
  print("Write beta diversity maps")
  Index <- paste("BetaDiversity_BCdiss_", scaling, sep = "")
  Beta.Path <- file.path(Output_Dir_BETA, paste(Index, "_", window_size, sep = ""))
  write_raster(Image = Beta$BetaDiversity, HDR = Beta$HDR, ImagePath = Beta.Path,
               window_size = window_size, FullRes = FullRes, LowRes = TRUE,
               SmoothImage = FALSE)

  PCoA_model <- Beta$PCoA_model
  return(PCoA_model)
}

#' computes NMDS
#
#' @param MatBCdist BC dissimilarity matrix
#' @param dimMDS numeric. number of dimensions of the NMDS
#
#' @return BetaNMDS_sel
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @importFrom ecodist nmds
#' @importFrom utils find
#' @export

compute_NMDS <- function(MatBCdist,dimMDS=3) {
  nbiterNMDS <- 4
  if (Sys.info()["sysname"] == "Windows") {
    nbCoresNMDS <- 2
  } else if (Sys.info()["sysname"] == "Linux") {
    nbCoresNMDS <- 4
  }
  # multiprocess of spectral species distribution and alpha diversity metrics
  plan(multiprocess, workers = nbCoresNMDS) ## Parallelize using four cores
  BetaNMDS <- future_lapply(MatBCdist, FUN = nmds, mindim = dimMDS, maxdim = dimMDS, nits = 1, future.packages = c("ecodist"))
  plan(sequential)
  # find iteration with minimum stress
  Stress <- vector(length = nbiterNMDS)
  for (i in 1:nbiterNMDS) {
    Stress[i] <- BetaNMDS[[i]]$stress
  }
  print("Stress obtained for NMDS iterations:")
  print(Stress)
  print("Rule of thumb")
  print("stress < 0.05 provides an excellent represention in reduced dimensions")
  print("stress < 0.1 is great")
  print("stress < 0.2 is good")
  print("stress > 0.3 provides a poor representation")
  MinStress <- find(Stress == min(Stress))
  BetaNMDS_sel <- BetaNMDS[[MinStress]]$conf
  BetaNMDS_sel <- data.frame(BetaNMDS_sel[[1]])
  return(BetaNMDS_sel)
}

#' Identifies ordination coordinates based on nearest neighbors
#'
#' @param Beta_Ordination_sel ordination of dissimilarity matrix for a selection of spatial units
#' @param SSD_Path character. path for spectral species distribution file
#' @param Sample_Sel numeric. Samples selected during ordination
#' @param coordTotSort numeric. coordinates of sunlit spatial units
#' @param nb_partitions number of k-means then averaged
#' @param nbclusters number of clusters
#' @param pcelim numeric. Minimum contribution in percent required for a spectral species
#' @param nbCPU numeric. number of CPUs available
#'
#' @return Ordination_est coordinates of each spatial unit in ordination space
#' @importFrom snow splitRows
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @export

ordination_to_NN <- function(Beta_Ordination_sel, SSD_Path, Sample_Sel, coordTotSort,
                             nb_partitions, nbclusters, pcelim, nbCPU = 1) {
  nb_Sunlit <- dim(coordTotSort)[1]
  # define number of samples to be sampled each time during paralle processing
  nb_samples_per_sub <- round(1e7 / dim(Sample_Sel)[1])
  # number of paralle processes to run
  nb.sub <- round(nb_Sunlit / nb_samples_per_sub)
  if (nb.sub == 0) nb.sub <- 1
  id.sub <- splitRows(as.matrix(seq(1, nb_Sunlit, by = 1), ncol = 1), ncl = nb.sub)
  # compute ordination coordinates from each subpart
  Nb_Units_Ordin <- dim(Sample_Sel)[1]
  if (nbCPU>1){
    plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
    Schedule_Per_Thread <- ceiling(nb.sub / nbCPU)
    OutPut <- future_lapply(id.sub,
                            FUN = ordination_parallel, coordTotSort = coordTotSort, SSD_Path = SSD_Path,
                            Sample_Sel = Sample_Sel, Beta_Ordination_sel = Beta_Ordination_sel, Nb_Units_Ordin = Nb_Units_Ordin,
                            nb_partitions = nb_partitions, nbclusters = nbclusters, pcelim = pcelim,
                            future.scheduling = Schedule_Per_Thread,
                            future.packages = c("vegan", "dissUtils", "R.utils", "tools", "snow", "matlab")
    )
    plan(sequential)
  } else {
    OutPut <- lapply(id.sub, FUN = ordination_parallel, coordTotSort = coordTotSort,
                     SSD_Path = SSD_Path, Sample_Sel = Sample_Sel,
                     Beta_Ordination_sel = Beta_Ordination_sel, Nb_Units_Ordin = Nb_Units_Ordin,
                     nb_partitions = nb_partitions, nbclusters = nbclusters, pcelim = pcelim)
  }
  Ordination_est <- do.call("rbind", OutPut)
  gc()
  return(Ordination_est)
}

#' applies results of ordination to full image based on nearest neighbors
#
#' @param id.sub numeric. ID of a subset of sunlit spatial units
#' @param coordTotSort numeric. coordinates of sunlit spatial units
#' @param SSD_Path character. path for spectral species distribution file
#' @param Sample_Sel numeric. Samples selected during ordination
#' @param Beta_Ordination_sel numeric. ordination of dissimilarity matrix for a selection of spatial units
#' @param Nb_Units_Ordin numeric. Maximum number of spatial units to be processed in NMDS.
#' @param nb_partitions numeric. Number of partitions (repetitions) to be computed then averaged.
#' @param nbclusters numeric. Number of clusters defined in k-Means.
#' @param pcelim numeric. Minimum contribution in percents required for a spectral species
#
#' @return OutPut list of mean and SD of alpha diversity metrics
#' @export

ordination_parallel <- function(id.sub, coordTotSort, SSD_Path, Sample_Sel, Beta_Ordination_sel,
                                Nb_Units_Ordin, nb_partitions, nbclusters, pcelim) {

  # get Spectral species distribution
  coordPix <- coordTotSort[id.sub, ]
  SSD_NN <- extract_samples_from_image(SSD_Path, coordPix)
  # compute the mean BC dissimilarity sequentially for each iteration
  MatBCtmp <- matrix(0, nrow = nrow(id.sub), ncol = Nb_Units_Ordin)
  SSDList <- list()
  for (i in 1:nb_partitions) {
    lb <- 1 + (i - 1) * nbclusters
    ub <- i * nbclusters
    SSDList[[1]] <- SSD_NN[, lb:ub]
    SSDList[[2]] <- Sample_Sel[, lb:ub]
    MatBCtmp0 <- compute_BCdiss(SSDList, pcelim)
    MatBCtmp <- MatBCtmp + MatBCtmp0
  }
  MatBCtmp <- MatBCtmp / nb_partitions
  # get the knn closest neighbors from each kernel
  knn <- 3
  OutPut <- compute_NN_from_ordination(MatBCtmp, knn, Beta_Ordination_sel)
  return(OutPut)
}


#' compute beta diversity from spectral species computed for a plot
#' expecting a matrix of spectral species (n pixels x p repetitions)
#'
#' @param SpectralSpecies_Plots list. list of matrices of spectral species corresponding to plots
#' @param nbclusters numeric. number of clusters
#' @param pcelim numeric. Minimum contribution (in %) required for a spectral species
#' each spectral species with a proprtion < pcelim is eliminated before computation of diversity
#
#' @return Mean bray curtis dissimilarity matrix for all plots, and individual BC matrices corresponding to each repetitions
#' @importFrom vegan vegdist
#' @export

compute_BETA_FromPlots <- function(SpectralSpecies_Plots,nbclusters,pcelim = 0.02){

  nbPolygons<- length(SpectralSpecies_Plots)
  Pixel.Inventory.All <- list()
  nb_partitions <- dim(SpectralSpecies_Plots[[1]])[2]
  # for each plot
  for (plot in 1:nbPolygons){
    # for each repetition
    Pixel.Inventory <- list()
    for (i in 1:nb_partitions){
      # compute distribution of spectral species
      Distritab <- table(SpectralSpecies_Plots[[plot]][,i])
      # compute distribution of spectral species
      Pixel.Inventory[[i]] <- as.data.frame(Distritab)
      SumPix <- sum(Pixel.Inventory[[i]]$Freq)
      ThreshElim <- pcelim*SumPix
      ElimZeros <- which(Pixel.Inventory[[i]]$Freq<ThreshElim)
      if (length(ElimZeros)>=1){
        Pixel.Inventory[[i]] <- Pixel.Inventory[[i]][-ElimZeros,]
      }
      if (length(which(Pixel.Inventory[[i]]$Var1==0))==1){
        Pixel.Inventory[[i]] <- Pixel.Inventory[[i]][-which(Pixel.Inventory[[i]]$Var1==0),]
      }
    }
    Pixel.Inventory.All[[plot]] <- Pixel.Inventory
  }

  # for each pair of plot, compute beta diversity indices
  BC <- list()
  for(i in 1:nb_partitions){
    MergeDiversity <- matrix(0,nrow = nbclusters,ncol = nbPolygons)
    for(j in 1:nbPolygons){
      SelSpectralSpecies <- as.numeric(as.vector(Pixel.Inventory.All[[j]][[i]]$Var1))
      SelFrequency <- Pixel.Inventory.All[[j]][[i]]$Freq
      MergeDiversity[SelSpectralSpecies,j] = SelFrequency
    }
    BC[[i]] <- vegan::vegdist(t(MergeDiversity),method="bray")
  }
  BC_mean <- 0*BC[[1]]
  for(i in 1:nb_partitions){
    BC_mean <- BC_mean+BC[[i]]
  }
  BC_mean <- BC_mean/nb_partitions
  beta <- list('BrayCurtis' = BC_mean, 'BrayCurtis_ALL' = BC)
  return(beta)
}



#' computes beta diversity metrics
#'
#' @param ClusterMap_Path character. File containing spectral species or classes from prior classification
#' @param MinSun numeric. minimum proportion of sunlit pixels required to consider plot
#' @param Nb_Units_Ordin numeric. maximum number of spatial units to be processed in Ordination
#' @param nb_partitions numeric. number of iterations
#' @param nbclusters numeric. number of clusters defined in k-Means
#' @param pcelim numeric. Minimum contribution in percents required for a spectral species
#' @param scaling character. scaling method, PCO or NMDS
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param dimMDS numeric. number of dimensions for the scaling.
#'
#' @return
#' @importFrom labdsv pco
#' @importFrom stats as.dist
#' @export

compute_beta_metrics <- function(ClusterMap_Path, MinSun, Nb_Units_Ordin, nb_partitions,
                                 nbclusters, pcelim, scaling = 'PCO', nbCPU = FALSE,
                                 MaxRAM = FALSE, dimMDS=3) {
  # Define path for images to be used
  SSD_Path <- paste(ClusterMap_Path, '_Distribution', sep = '')
  ImPathSunlit <- paste(ClusterMap_Path, '_Distribution_Sunlit', sep = '')
  # Get illuminated pixels based on  SpectralSpecies_Distribution_Sunlit
  Sunlit_Pixels <- get_sunlit_pixels(ImPathSunlit, MinSun)
  Select_Sunlit <- Sunlit_Pixels$Select_Sunlit
  nb_Sunlit <- length(Select_Sunlit)
  # Define spatial units processed through ordination and those processed through
  # Nearest neighbor based on the first ordination batch
  print("Read Spectral Species distribution")
  RandPermKernels <- sample(seq(1, nb_Sunlit, by = 1))
  if (Nb_Units_Ordin <= nb_Sunlit) {
    Kernels_NN <- RandPermKernels[(Nb_Units_Ordin + 1):nb_Sunlit]
  } else {
    Nb_Units_Ordin <- nb_Sunlit
    Kernels_NN <- c()
  }
  # read spectral species distribution file
  SSD_All <- extract_samples_from_image(SSD_Path, Sunlit_Pixels$coordTotSort)
  # define kernels used for Ordination
  Kernels_Ordination <- RandPermKernels[1:Nb_Units_Ordin]
  Sample_Sel <- SSD_All[Kernels_Ordination, ]
  rm(SSD_All)
  gc()

  # create a Bray curtis dissimilarity matrix for each iteration
  print("compute BC dissimilarity for selected kernels")
  # create a list in with each element is an iteration
  MatBC <- matrix(0, ncol = Nb_Units_Ordin, nrow = Nb_Units_Ordin)
  SSDList <- list()
  BC.from.SSD <- list()
  for (i in 1:nb_partitions) {
    lb <- 1 + (i - 1) * nbclusters
    ub <- i * nbclusters
    SSDList[[1]] <- Sample_Sel[, lb:ub]
    SSDList[[2]] <- Sample_Sel[, lb:ub]
    BC.from.SSD <- compute_BCdiss(SSDList, pcelim)
    MatBC <- MatBC + BC.from.SSD
  }
  MatBC <- MatBC / nb_partitions

  # Perform Ordination based on BC dissimilarity matrix
  print("perform Ordination on the BC dissimilarity matrix averaged from all iterations")
  # parallel computing of Ordination can be run on 2 cores on Windows.
  # core management seems better on linux --> 4 cores possible
  MatBCdist <- as.dist(MatBC, diag = FALSE, upper = FALSE)
  BetaPCO <- NULL
  if (scaling == 'NMDS') {
    Beta_Ordination_sel <- compute_NMDS(MatBCdist)
    PCname <- 'NMDS'
  } else if (scaling == 'PCO') {
    BetaPCO <- labdsv::pco(MatBCdist, k = dimMDS)
    Beta_Ordination_sel <- BetaPCO$points
    PCname <- 'PCoA'
  }

  # Perform nearest neighbor on spatial units excluded from Ordination
  print("BC dissimilarity between samples selected for Ordination and remaining")
  coordTotSort <- Sunlit_Pixels$coordTotSort
  Ordination_est <- ordination_to_NN(Beta_Ordination_sel, SSD_Path, Sample_Sel,
                                     coordTotSort, nb_partitions, nbclusters, pcelim,
                                     nbCPU = nbCPU)

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
  list <- ls()

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

  rm(list = list[-which(list == 'BetaDiversityRGB' | list == 'Select_Sunlit' |
                          list == 'HDR_Beta' | list == 'BetaPCO')])
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

compute_BCdiss <- function(SSDList, pcelim) {
  # compute the proportion of each spectral species
  # Here, the proportion is computed with regards to the total number of sunlit pixels
  # One may want to determine if the results are similar when the proportion is computed
  # with regards to the total number of pixels (se*se)
  # however it would increase dissimilarity betwen kernels with different number of sunlit pixels
  SSD <- list()
  for (i in 1:length(SSDList)) {
    # get the total number of sunlit pixels in spatial unit
    SumSpecies <- rowSums(SSDList[[i]])
    elim <- which(SumSpecies == 0)
    if (length(elim) > 0) {
      SumSpecies[elim] <- 1
      SSDList[[i]][elim, ] <- 0
    }
    SSD[[i]] <- apply(SSDList[[i]], 2, function(x, c1) x / c1, 'c1' = SumSpecies)
    SSD[[i]][which(SSD[[i]] < pcelim)] <- 0
  }
  # matrix of bray curtis dissimilarity (size = nb kernels x nb kernels)
  # Here use the package "dissUtils" to compute dissimilarity matrices sequentially
  MatBC <- diss(SSD[[1]], SSD[[2]], method = 'braycurtis')
  # EDIT 06-Feb-2019: added this to fix problem when empty kernels occur, leading to NA BC value
  if (length(which(is.na(MatBC) == TRUE)) > 0) {
    MatBC[which(is.na(MatBC) == TRUE)] <- 1
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
  Ordin_est <- lapply(X = NN,FUN = WeightedCoordsNN,knn=knn, BetaDiversity0=BetaDiversity0)
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
