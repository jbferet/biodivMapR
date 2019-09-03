# ==============================================================================
# biodivMapR
# Lib_ValidationbiodivMapR.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
# ==============================================================================
# This Library contains functions to perform validation on products from biodivMapR
# the main goal is to validate ground data associated with biodiversity metrics
# ==============================================================================

#' get projection of a raster or a vector
#' @param file path for a raster or vector (shapefile)
#' @param type 'raster' or 'vector'
#' @return projection
#' @importFrom raster raster shapefile projection
#' @importFrom rgdal readOGR
#' @import tools
#' @export
get_projection <- function(file, type = 'raster'){
  if (type == 'raster'){
    obj <- raster(file)
  } else if (type == 'vector'){
    Shp.Path      = dirname(file)
    Shp.Name      = file_path_sans_ext(basename(file))
    obj  = readOGR(dsn = Shp.Path,layer = Shp.Name,verbose = FALSE)
  }
  projstr <- projection(obj)
  return(projstr)
}

#' Get list of shapefiles in a directory
#'
#' @param x character or list. Directory containing shapefiles
#' @return list of shapefiles names
#' @export
list_shp <- function(x){
  if(typeof(x)=='list'){
    List.Shp = c()
    ii = 0
    for (shp in x){
      ii = ii+1
      List.Shp[ii]  = dir(shp, pattern = '.shp$', full.names = TRUE, ignore.case = FALSE,include.dirs = FALSE)
    }
  } else if(typeof(x)=='character'){
    List.Shp  = dir(x, pattern = '.shp$', full.names = TRUE, ignore.case = FALSE,include.dirs = FALSE)
  }
  return(List.Shp)
}

# reprojects a vector file and saves it
# @param Initial.File path for a shapefile to be reprojected
# @param Projection projection to be applied to Initial.File
# @param Reprojected.File path for the reprojected shapefile
# @return
#' @importFrom rgdal readOGR writeOGR
#' @importFrom sp spTransform
#' @import tools
reproject_vector = function(Initial.File,Projection,Reprojected.File){

  Shp.Path      = dirname(Initial.File)
  Shp.Name      = file_path_sans_ext(basename(Initial.File))
  Vect.Init     = readOGR(Shp.Path,Shp.Name,verbose = FALSE)
  Proj.init     = projection(Vect.Init)

  if (!Proj.init==Projection){
    Shp.Path      = dirname(Reprojected.File)
    Shp.Name      = file_path_sans_ext(basename(Reprojected.File))
    Vect.reproj = spTransform(Vect.Init, Projection)
    writeOGR(obj = Vect.reproj, dsn = Shp.Path,layer = Shp.Name, driver="ESRI Shapefile",overwrite_layer = TRUE) #also you were missing the driver argument
  }
  return()
}

# Extracts pixels coordinates from raster corresponding to an area defined by a vector
# @param Path.Raster path for the raster file. !! BIL expected
# @param Path.Vector path for the vector file. !! SHP expected
# @return ColRow list of coordinates of pixels corresponding to each polygon in shp
#' @importFrom rgdal readOGR
#' @import tools
extract_pixels_coordinates = function(Path.Raster,Path.Vector){
  # read vector file
  Shp.Path    = dirname(Path.Vector)
  Shp.Name    = file_path_sans_ext(basename(Path.Vector))
  Shp.Crop    = readOGR(Shp.Path,Shp.Name)
  # read raster info
  Raster      = raster(Path.Raster, band = 1)
  # extract pixel coordinates from raster based on vector
  XY  = raster::cellFromPolygon (Raster,Shp.Crop)
  # for each polygon in the
  ColRow = list()
  for (i in 1:length(XY)){
    ColRow[[i]] = ind2sub(Raster,XY[[i]])
  }
  return(ColRow)
}

# Extracts pixels coordinates from raster corresponding to an area defined by a vector
# @param Path.Raster path for the raster file. !! BIL expected
# @param OGR.Vector  OGR for the vector file obtained from readOGR
# @return ColRow list of coordinates of pixels corresponding to each polygon in shp
extract_pixels_coordinates.From.OGR = function(Path.Raster,OGR.Vector){
  # read raster info
  Raster      = raster(Path.Raster, band = 1)
  # for each polygon or point in the shapefile
  ColRow = list()
  # extract pixel coordinates from raster based on vector
  if (OGR.Vector@class[1]=='SpatialPointsDataFrame'){
    XY  = raster::cellFromXY (Raster,OGR.Vector)
    ColRow[[1]] = ind2sub(Raster,XY)

  } else if  (OGR.Vector@class[1]=='SpatialPolygonsDataFrame'){
    XY  = raster::cellFromPolygon (Raster,OGR.Vector)
    for (i in 1:length(XY)){
      ColRow[[i]] = ind2sub(Raster,XY[[i]])
    }
  }
  return(ColRow)
}

# Computes alpha diversity metrics from distribution
# @param Distrib distribution of clusters
# @return Richness, Fisher, Shannon, Simpson
get_alpha_metrics = function(Distrib){
  RichnessPlot  = vegan::specnumber(Distrib, MARGIN = 1)        # species richness
  # if (length(Distrib)>1){
  #   fisherPlot    = fisher.alpha(Distrib, MARGIN = 1)      # fisher's alpha
  # } else {
  #   fisherPlot    = 0
  # }
  FisherPlot    = 0
  ShannonPlot   = vegan::diversity(Distrib, index = "shannon", MARGIN = 1, base = exp(1)) # shannon's alpha
  SimpsonPlot   = vegan::diversity(Distrib, index = "simpson", MARGIN = 1, base = exp(1))
  return(list("Richness" = RichnessPlot,"fisher"=FisherPlot,"Shannon"=ShannonPlot,"Simpson"=SimpsonPlot))
}

#' gets alpha diversity indicators from plot
#' @param Raster SpectralSpecies file computed from DiverstyMapping method
#' @param Plots list of shapefiles included in the raster
#' @param NbClusters numeric. Number of clusters defined in k-Means.
#' @param Name_Plot character. Name of teh plots defined in the shapefile
#' @return alpha and beta diversity metrics
#' @importFrom rgdal readOGR
#' @import tools
#' @export
diversity_from_plots = function(Raster, Plots,NbClusters = 50,Name_Plot = FALSE){
  # get hdr from raster
  HDR           = read_ENVI_header(paste(Raster,'.hdr',sep=''))
  nbRepetitions = HDR$bands
  # get the number of plots
  nbPlots             = length(Plots)
  # define variables where alpha diversity will be stored
  Richness.AllRep = Shannon.AllRep = Fisher.AllRep = Simpson.AllRep = list()
  # initialize alpha diversity for each plot
  Richness       = Shannon    = Fisher     = Simpson    = data.frame()
  Name.Vector       = list()
  # prepare directory to write shapefiles with correct projection
  Dir.Vector            = dirname(Plots[[1]])
  Dir.Vector.reproject  = paste(Dir.Vector,'/Reproject',sep='')
  ###                 Compute alpha diversity             ###
  # for each plot
  Pixel.Inventory.All   = list()
  PolygonID             = 0

  for (ip in 1:nbPlots){
    # prepare for possible reprojection
    File.Vector           = Plots[[ip]]
    Name.Vector[[ip]]     = file_path_sans_ext(basename(File.Vector))
    print(paste('Processing file ',Name.Vector[[ip]],sep=''))
    File.Vector.reproject = paste(Dir.Vector.reproject,'/',Name.Vector[[ip]],'.shp','sep'='')

    if (file.exists(paste(file_path_sans_ext(File.Vector),'.shp',sep=''))){
      Plot                = readOGR(Dir.Vector,Name.Vector[[ip]],verbose = FALSE)
      # check if vector and rasters are in the same referential
      Projection.Plot     = get_projection(File.Vector,'vector')
      Projection.Raster     = get_projection(Raster,'raster')
      # if not, convert vector file
      if (!Projection.Raster==Projection.Plot){
        if (!file.exists(File.Vector.reproject)){
          dir.create(Dir.Vector.reproject, showWarnings = FALSE,recursive = TRUE)
          reproject_vector(File.Vector,Projection.Raster,File.Vector.reproject)
        }
        Plot              = readOGR(Dir.Vector.reproject,Name.Vector[[ip]],verbose = FALSE)
      }
    } else if (file.exists(paste(File.Vector,'kml','sep'='.'))){
      print('Please convert vector file to shpfile')
    }

    # extract data corresponding to the raster
    XY          = extract_pixels_coordinates.From.OGR(Raster,Plot)
    # if the plot is included in the raster
    if (length(XY)==1 & length(XY[[1]]$Column)==0){
      if (length(Name_Plot)==nbPlots){
        Name_Plot[ip] = NA
      }
    }
    if (length(XY)>1 | length(XY[[1]]$Column)>0){
      ID.Poly     = list()
      idPix       = 1
      if (length(XY)==1){
        coordPix    = cbind(XY[[1]]$Row,XY[[1]]$Column)
        ID.Poly[[1]]= seq(idPix,idPix+length(XY[[1]]$Row)-1,by = 1)
      } else {
        coordPix    = list()
        ii          = 0
        for (pp in 1:length(XY)){
          ii        = ii+1
          if (length(XY[[pp]]$Row)>0){
            coordPix[[pp]]  = cbind(XY[[pp]]$Row,XY[[pp]]$Column)
            ID.Poly[[pp]]   = seq(idPix,idPix+length(XY[[pp]]$Row)-1,by = 1)
            idPix           = idPix+length(XY[[pp]]$Row)
          }
        }
        coordPix    = do.call(rbind,coordPix)
      }
      ExtractIm   = extract_pixels(coordPix,Raster,HDR)
      # compute alpha diversity for each repetition
      Pixel.Inventory = list()
      Richness.tmp = Shannon.tmp = Fisher.tmp = Simpson.tmp = vector(length = nbRepetitions)
      # eliminate spectral species contributing to less than pcelim percent of the total valid pixels
      pcelim  = 0.02
      for (ii in 1:length(ID.Poly)){
        if (length(ID.Poly[[ii]])>0){
          for (i in 1:nbRepetitions){
            Distritab     = table(ExtractIm[ID.Poly[[ii]],i])
            Pixel.Inventory[[i]]  = as.data.frame(Distritab)
            SumPix        = sum(Pixel.Inventory[[i]]$Freq)
            ThreshElim    = pcelim*SumPix
            ElimZeros = which(Pixel.Inventory[[i]]$Freq<ThreshElim)
            if (length(ElimZeros)>=1){
              Pixel.Inventory[[i]]   = Pixel.Inventory[[i]][-ElimZeros,]
            }
            if (length(which(Pixel.Inventory[[i]]$Var1==0))==1){
              Pixel.Inventory[[i]]   = Pixel.Inventory[[i]][-which(Pixel.Inventory[[i]]$Var1==0),]
            }
            Alpha = get_alpha_metrics(Pixel.Inventory[[i]]$Freq)
            # Alpha diversity
            Richness.tmp[i]   = as.numeric(Alpha$Richness)
            Fisher.tmp[i]     = Alpha$fisher
            Shannon.tmp[i]    = Alpha$Shannon
            Simpson.tmp[i]    = Alpha$Simpson
          }
          PolygonID = PolygonID+1
          Richness.AllRep[[PolygonID]] = Richness.tmp
          Shannon.AllRep[[PolygonID]]  = Shannon.tmp
          Fisher.AllRep[[PolygonID]]   = Fisher.tmp
          Simpson.AllRep[[PolygonID]]  = Simpson.tmp
          Richness   = rbind(Richness, mean(Richness.tmp), row.names = NULL, col.names = NULL)
          Fisher     = rbind(Fisher, mean(Fisher.tmp), row.names = NULL, col.names = NULL)
          Shannon    = rbind(Shannon, mean(Shannon.tmp), row.names = NULL, col.names = NULL)
          Simpson    = rbind(Simpson, mean(Simpson.tmp), row.names = NULL, col.names = NULL)
          Pixel.Inventory.All[[PolygonID]] = Pixel.Inventory
        }
      }
    }
  }
  Richness.AllRep = do.call(rbind,Richness.AllRep)
  Shannon.AllRep  = do.call(rbind,Shannon.AllRep)
  Fisher.AllRep   = do.call(rbind,Fisher.AllRep)
  Simpson.AllRep  = do.call(rbind,Simpson.AllRep)

  ###                 Compute beta diversity             ###
  # for each pair of plot, compute beta diversity indices
  BC = list()
  for(i in 1:nbRepetitions){
    MergeDiversity = matrix(0,nrow = NbClusters,ncol = PolygonID)
    for(j in 1:PolygonID){
      SelSpectralSpecies  =  as.numeric(as.vector(Pixel.Inventory.All[[j]][[i]]$Var1))
      SelFrequency        = Pixel.Inventory.All[[j]][[i]]$Freq
      MergeDiversity[SelSpectralSpecies,j] = SelFrequency
    }
    BC[[i]] = vegan::vegdist(t(MergeDiversity),method="bray")
  }
  BC_mean   = 0*BC[[1]]
  for(i in 1:nbRepetitions){
    BC_mean = BC_mean+BC[[i]]
  }
  if (length(Name_Plot)>1){
    elim = which(is.na(Name_Plot))
    if (length(elim)>0){
      Name_Plot =   Name_Plot[-elim]
    }
  }
  BC_mean = as.matrix(BC_mean/nbRepetitions)
  return(list("Richness" = Richness,"Fisher"=Fisher,"Shannon"=Shannon,"Simpson"=Simpson,'BCdiss' = BC_mean,"fisher.All"=Fisher.AllRep,"Shannon.All"=Shannon.AllRep,"Simpson.All"=Simpson.AllRep,'BCdiss.All' = BC,'Name_Plot'=Name_Plot))
}

# build a vector file from raster footprint
# borrowed from https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
# @param x path for a raster or raster object
# @param outshape path for a vector to be written
# @param gdalformat
# @param pypath
# @param readpoly
# @param quiet
# @return NULL
#' @importFrom rgdal readOGR
#' @importFrom raster writeRaster
#' @importFrom methods is
gdal_polygonizeR = function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                            pypath=NULL, readpoly=TRUE, quiet=TRUE) {

  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape  = sub('\\.shp$', '', outshape)
    f.exists  = file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists)){
      print('Footprint of raster file already defined')
    }
    # stop(sprintf('File already exists: %s',
    #              toString(paste(outshape, c('shp', 'shx', 'dbf'),
    #                             sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    writeRaster(x, {f <- tempfile(fileext='.tif')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  outshape.Shp = paste(outshape,'.shp',sep='')
  if (!file.exists(outshape.Shp)){
    system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                    pypath, rastpath, gdalformat, outshape)))
  }
  if (readpoly==TRUE) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
    return(shp)
  }
  return(NULL)
}
