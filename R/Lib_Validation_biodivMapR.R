# ==============================================================================
# biodivMapR
# Lib_ValidationbiodivMapR.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This Library contains functions to perform validation on products from biodivMapR
# the main goal is to validate ground data associated with biodiversity metrics
# ==============================================================================


#' gets alpha diversity indicators from plot
#' @param Raster_SpectralSpecies character. path for the SpectralSpecies file computed from DiverstyMapping method
#' @param Plots list. list of paths corresponding to shapefiles defining polygons in the raster
#' @param nbclusters numeric. Number of clusters defined in k-Means.
#' @param Raster_Functional character. path for the raster file used to compute functional diversity
#' @param Selected_Features numeric. selected features for Raster_Functional. use all features if set to FALSE
#' @param Name_Plot character. Name of the plots defined in the shapefiles
#' @param Hellinger boolean. set TRUE to compute Hellinger distance matrices based on Euclidean distances
#' @param pcelim numeric. Discard spectral species when contribution is < pcelim (between 0 and 1)
#' @param FDmetric character. Functional diversity metric
#'
#' @return alpha and beta diversity metrics
#' @importFrom terra rast vect same.crs
#' @importFrom tools file_path_sans_ext
#' @export
diversity_from_plots = function(Raster_SpectralSpecies, Plots, nbclusters = 50,
                                Raster_Functional = FALSE, Selected_Features = FALSE,
                                Name_Plot = FALSE, Hellinger = FALSE, pcelim = 0.02,
                                FDmetric = c('FRic', 'FEve', 'FDiv')){

  # get hdr from Raster_SpectralSpecies
  HDR <- read_ENVI_header(paste(Raster_SpectralSpecies,'.hdr',sep=''))
  nbRepetitions <- HDR$bands
  # get the number of plots
  nbPlots <- length(Plots)
  # define variables where alpha diversity will be stored
  Richness.AllRep <- Shannon.AllRep <- Fisher.AllRep <- Simpson.AllRep <- list()
  # initialize alpha diversity for each plot
  Richness <- Shannon <- Fisher <- Simpson <- data.frame()
  Name.Vector <- list()

  ## ______________________________________________________##
  ###       Read pixel coordinates for each polygon       ###
  ## ______________________________________________________##

  # list of corodinates corresponding to each polygon to be studied
  Raster <- rast(Raster_SpectralSpecies, lyrs=1)
  XY <- list()
  for (ip in 1:nbPlots){
    File.Vector <- Plots[[ip]]
    Name.Vector[[ip]] <- tools::file_path_sans_ext(basename(File.Vector))
    print(paste('Reading pixels coordinates for polygons in ',Name.Vector[[ip]],sep=''))

    if (file.exists(paste(tools::file_path_sans_ext(File.Vector),'.shp',sep=''))){
      Plot <- terra::vect(Plots[[ip]])
      # check if vector and rasters are in the same referential
      # if not, convert vector file
      if (!same.crs(Raster, Plot)){
        stop('Raster and Plots have different projection. Plots should be reprojected to Raster CRS')
      }
    } else if (file.exists(paste(File.Vector,'kml','sep'='.'))){
      print('Please convert vector file to shpfile')
    }
    # extract data corresponding to the Raster_SpectralSpecies
    XY0 <- extract_pixels_coordinates(Raster,Plot)
    XY <- c(XY, XY0)
  }
  # total number of plots among all files
  nbPolygons <- length(XY)
  message(paste('Number of polygons: ',nbPolygons))

  ## ______________________________________________________##
  ###                 Compute alpha diversity             ###
  ## ______________________________________________________##
  # for each polygon
  Pixel_Inventory_All <- Pixel_Hellinger_All <- list()
  for (ip in 1:nbPolygons){
    # if only one polygon in the shapefile and if the polyon is not included in the Raster_SpectralSpecies
    if (all(is.na(XY[[ip]]$col))){
      # if list of individual plots provided
      if (length(Name_Plot)==nbPolygons){
        message(paste('Polygon named',Name_Plot[ip],'is out of the raster'))
        # Set name to NA
        Name_Plot[ip] <- NA
      }
      Richness <- rbind(Richness, NA, row.names = NULL, col.names = NULL)
      Fisher <- rbind(Fisher, NA, row.names = NULL, col.names = NULL)
      Shannon <- rbind(Shannon, NA, row.names = NULL, col.names = NULL)
      Simpson <- rbind(Simpson, NA, row.names = NULL, col.names = NULL)
    } else {
      ExtractIm <- extract.big_raster(Raster_SpectralSpecies, XY[[ip]])
      if (length(XY[[ip]]$col)==1){
        ExtractIm <- matrix(ExtractIm,ncol = nbRepetitions)
      }
      # compute alpha diversity for each repetition
      Pixel_Inventory <- Pixel_Hellinger <- list()
      Richness.tmp <- Shannon.tmp <- Fisher.tmp <- Simpson.tmp <- vector(length = nbRepetitions)
      # eliminate spectral species contributing to less than pcelim percent of the total valid pixels
      # pcelim <- 0.02
      for (i in 1:nbRepetitions){
        if (nbRepetitions ==1){
          Distritab <- table(ExtractIm)
        } else {
          Distritab <- table(ExtractIm[,i])
        }
        Pixel_Inventory[[i]] <- as.data.frame(Distritab)
        # fix 2022/07/20: first discard shaded pixels before application of pcelim
        if (length(which(Pixel_Inventory[[i]]$Var1==0))==1){
          Pixel_Inventory[[i]] <- Pixel_Inventory[[i]][-which(Pixel_Inventory[[i]]$Var1==0),]
        }
        SumPix <- sum(Pixel_Inventory[[i]]$Freq)
        ThreshElim <- pcelim*SumPix
        ElimZeros <- which(Pixel_Inventory[[i]]$Freq<ThreshElim)
        if (length(ElimZeros)>=1){
          Pixel_Inventory[[i]] <- Pixel_Inventory[[i]][-ElimZeros,]
        }
        # if (length(which(Pixel_Inventory[[i]]$Var1==0))==1){
        #   Pixel_Inventory[[i]] <- Pixel_Inventory[[i]][-which(Pixel_Inventory[[i]]$Var1==0),]
        # }
        Alpha <- get_alpha_metrics(Pixel_Inventory[[i]]$Freq)
        # Alpha diversity
        Richness.tmp[i] <- as.numeric(Alpha$Richness)
        Fisher.tmp[i] <- Alpha$fisher
        Shannon.tmp[i] <- Alpha$Shannon
        Simpson.tmp[i] <- Alpha$Simpson
      }
      Richness.AllRep[[ip]] <- Richness.tmp
      Shannon.AllRep[[ip]] <- Shannon.tmp
      Fisher.AllRep[[ip]] <- Fisher.tmp
      Simpson.AllRep[[ip]] <- Simpson.tmp
      Richness <- rbind(Richness, mean(Richness.tmp), row.names = NULL, col.names = NULL)
      Fisher <- rbind(Fisher, mean(Fisher.tmp), row.names = NULL, col.names = NULL)
      Shannon <- rbind(Shannon, mean(Shannon.tmp), row.names = NULL, col.names = NULL)
      Simpson <- rbind(Simpson, mean(Simpson.tmp), row.names = NULL, col.names = NULL)
      Pixel_Inventory_All[[ip]] <- Pixel_Inventory
      # compute beta with Hellinger
      if (Hellinger == TRUE){
        for (i in 1:nbRepetitions){
          # compute Hellinger distance
          Pixel_Hellinger[[i]] <- Pixel_Inventory[[i]]
          Pixel_Hellinger[[i]]$Freq <- sqrt(Pixel_Hellinger[[i]]$Freq/sum(Pixel_Hellinger[[i]]$Freq))
        }
        Pixel_Hellinger_All[[ip]] <- Pixel_Hellinger
      }
    }
  }
  Richness.AllRep <- do.call(rbind,Richness.AllRep)
  Shannon.AllRep <- do.call(rbind,Shannon.AllRep)
  Fisher.AllRep <- do.call(rbind,Fisher.AllRep)
  Simpson.AllRep <- do.call(rbind,Simpson.AllRep)

  ## ______________________________________________________##
  ###           Compute functional  diversity             ###
  ## ______________________________________________________##
  FunctionalDiversity <- NULL
  if (!Raster_Functional==FALSE){
    FunctionalDiversity <- data.frame('FEve'=vector(length = nbPolygons),
                                      'FRic'=vector(length = nbPolygons),
                                      'FDiv'=vector(length = nbPolygons))
    # for each polygon
    for (ip in 1:nbPolygons){
      # if only one polygon in the shapefile and if the polyon is not included in the raster
      if (all(is.na(XY[[ip]]$col))){
        FunctionalDiversity$FRic[ip] <- NA
        FunctionalDiversity$FEve[ip] <- NA
        FunctionalDiversity$FDiv[ip] <- NA
        # if list of individual plots provided
        if (length(Name_Plot)==nbPolygons){
          message(paste('Polygon named',Name_Plot[ip],'is out of the raster'))
          # Set name to NA
          Name_Plot[ip] <- NA
        }
      } else {
        ExtractIm <- extract.big_raster(Raster_Functional, XY[[ip]])
        if (!Selected_Features[1]==FALSE){
          ExtractIm <- ExtractIm[,Selected_Features]
        } else {
          Selected_Features <- seq(1,ncol(ExtractIm))
        }
        ij <- ExtractIm
        # keep non zero values
        ij <- matrix(ij[which(!is.na(ij[,1])),], ncol = ncol(ExtractIm))
        nbPix_Sunlit <- dim(ij)[1]
        PCsun <- nbPix_Sunlit / nrow(ExtractIm)
        spectraits <- data.frame(ij)
        if (nbPix_Sunlit>length(Selected_Features)){
          FD <- getFD(spectraits = spectraits, FDmetric = FDmetric)
          FunctionalDiversity$FRic[ip] <- FD$FRic
          FunctionalDiversity$FDiv[ip] <- FD$FDiv
          FunctionalDiversity$FEve[ip] <- FD$FEve
        } else {
          FunctionalDiversity$FRic[ip] <- 0
          FunctionalDiversity$FDiv[ip] <- 0
          FunctionalDiversity$FEve[ip] <- 0
        }
      }
    }
  }

  ## ______________________________________________________##
  ###                 Compute beta diversity              ###
  ## ______________________________________________________##
  # for each pair of plot, compute beta diversity indices
  BC <- list()
  for(i in 1:nbRepetitions){
    MergeDiversity <- matrix(0,nrow = nbclusters,ncol = nbPolygons)
    for(j in 1:nbPolygons){
      if (nbRepetitions>1){
        SelSpectralSpecies <- as.numeric(as.vector(Pixel_Inventory_All[[j]][[i]]$Var1))
        SelFrequency <- Pixel_Inventory_All[[j]][[i]]$Freq
      } else {
        SelSpectralSpecies <- as.numeric(as.vector(Pixel_Inventory_All[[j]][[i]]$ExtractIm))
        SelFrequency <- Pixel_Inventory_All[[j]][[i]]$Freq
      }
      MergeDiversity[SelSpectralSpecies,j] = SelFrequency
    }
    BC[[i]] <- vegan::vegdist(t(MergeDiversity),method="bray")
  }
  BC_mean <- 0*BC[[1]]
  for(i in 1:nbRepetitions){
    BC_mean <- BC_mean+BC[[i]]
  }

  # Hellinger
  Hellinger_mean <- Hellmat <- NULL
  if (Hellinger==TRUE){
    # for each pair of plot, compute Euclidean distance on Hellinger
    Hellmat <- list()
    for(i in 1:nbRepetitions){
      MergeDiversity <- matrix(0,nrow = nbclusters,ncol = nbPolygons)
      for(j in 1:nbPolygons){
        SelSpectralSpecies <- as.numeric(as.vector(Pixel_Hellinger_All[[j]][[i]]$Var1))
        SelFrequency <- Pixel_Hellinger_All[[j]][[i]]$Freq
        MergeDiversity[SelSpectralSpecies,j] = SelFrequency
      }
      Hellmat[[i]] <- vegan::vegdist(t(MergeDiversity),method="euclidean")
    }
    Hellinger_mean <- 0*Hellmat[[1]]
    for(i in 1:nbRepetitions){
      Hellinger_mean <- Hellinger_mean+Hellmat[[i]]
    }
    Hellinger_mean <- Hellinger_mean/nbRepetitions
  }

  # if (length(Name_Plot)>1){
  #   elim <- which(is.na(Name_Plot))
  #   if (length(elim)>0){
  #     Name_Plot <- Name_Plot[-elim]
  #   }
  # }
  BC_mean <- as.matrix(BC_mean/nbRepetitions)
  names(Richness) <- 'Richness'
  names(Fisher) <- 'Fisher'
  names(Shannon) <- 'Shannon'
  names(Simpson) <- 'Simpson'
  return(list("Richness" = Richness, "Fisher" = Fisher, "Shannon" = Shannon, "Simpson" = Simpson,
              "fisher.All" = Fisher.AllRep, "Shannon.All" = Shannon.AllRep, "Simpson.All" = Simpson.AllRep,
              "FunctionalDiversity"= FunctionalDiversity,
              'BCdiss' = BC_mean, 'BCdiss.All' = BC,
              'Hellinger' = Hellinger_mean, 'Hellinger_ALL' = Hellmat,
              'Name_Plot' = Name_Plot))
}





#' Get list of shapefiles in a directory
#'
#' @param x character or list. Directory containing shapefiles
#' @return list of shapefiles names
#' @export
list_shp <- function(x){
  if(typeof(x)=='list'){
    List.Shp <- c()
    ii <- 0
    for (shp in x){
      ii <- ii+1
      List.Shp[ii] <- dir(shp, pattern = '.shp$', full.names = TRUE, ignore.case = FALSE,include.dirs = FALSE)
    }
  } else if(typeof(x)=='character'){
    List.Shp <- dir(x, pattern = '.shp$', full.names = TRUE, ignore.case = FALSE,include.dirs = FALSE)
  }
  return(List.Shp)
}


#' Extracts pixels coordinates from raster intersecting a vector.
#' @param x terra::SpatRaster or character. Raster or raster path.
#' @param y terra::SpatVector or character. Vector or vector path.
#' @return data.frame. Pixel indices row, col of intersecting parts.
#' @import magrittr
#' @importFrom terra rast vect cells rowFromCell colFromCell
#' @importFrom dplyr mutate select group_by group_split
#' @export
extract_pixels_coordinates <- function(x, y){
  # read raster info
  if(is.character(x))
    x <- terra::rast(x, lyrs = 1)

  if(is.character(y))
    y <- terra::vect(y)

  ColRow <- cells(x, y) %>% as.data.frame() %>%
    dplyr::mutate(row = rowFromCell(x, cell),
                  col = colFromCell(x, cell)) %>%
    dplyr::select(ID, row, col) %>%
    dplyr::group_by(ID) %>% dplyr::group_split(.keep=F)

  return(ColRow)
}

#' Computes alpha diversity metrics from distribution
#' @param Distrib distribution of clusters
#'
#' @return Richness, Fisher, Shannon, Simpson
#' @importFrom vegan specnumber diversity
#' @export

get_alpha_metrics = function(Distrib){
  RichnessPlot <- vegan::specnumber(Distrib, MARGIN = 1)        # species richness
  FisherPlot <- 0
  ShannonPlot <- vegan::diversity(Distrib, index = "shannon", MARGIN = 1, base = exp(1)) # shannon's alpha
  SimpsonPlot <- vegan::diversity(Distrib, index = "simpson", MARGIN = 1, base = exp(1))
  return(list("Richness" = RichnessPlot,"fisher"=FisherPlot,"Shannon"=ShannonPlot,"Simpson"=SimpsonPlot))
}

