## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval=FALSE
)

## ----Input / Output files------------------------------------------------
#  Input.Image.File  = system.file('extdata', 'RASTER', 'S2A_T33NUD_20180104_Subset', package = 'biodivMapR')
#  check_data(Input.Image.File)
#  
#  Input.Image.File  = raster2BIL(Raster.Path = Input.Image.File,
#                                         Sensor = 'SENTINEL_2A',
#                                         Convert.Integer = TRUE,
#                                         Output.Directory = '~/test')
#  Input.Mask.File   = FALSE
#  
#  Output.Dir        = 'RESULTS'

## ----Spatial resolution--------------------------------------------------
#  window_size = 10

## ----PCA filtering-------------------------------------------------------
#  FilterPCA = TRUE

## ----Computing options---------------------------------------------------
#  nbCPU         = 4
#  MaxRAM        = 0.5
#  nbclusters    = 50

## ----Mask non vegetated / shaded / cloudy pixels-------------------------
#  NDVI.Thresh = 0.5
#  Blue.Thresh = 500
#  NIR.Thresh  = 1500
#  print("PERFORM RADIOMETRIC FILTERING")
#  ImPathShade = perform_radiometric_filtering(Input.Image.File, Input.Mask.File, Output.Dir,
#                                              NDVI.Thresh = NDVI.Thresh, Blue.Thresh = Blue.Thresh,
#                                              NIR.Thresh = NIR.Thresh)

## ----PCA-----------------------------------------------------------------
#  print("PERFORM PCA ON RASTER")
#  PCA.Files  = perform_PCA(Input.Image.File, ImPathShade, Output.Dir,
#                                 FilterPCA = TRUE, nbCPU = nbCPU, MaxRAM = MaxRAM)
#  print("Select PCA components for diversity estimations")
#  select_PCA_components(Input.Image.File, Output.Dir, PCA.Files, File.Open = TRUE)

## ----alpha and beta diversity maps---------------------------------------
#  print("MAP SPECTRAL SPECIES")
#  map_spectral_species(Input.Image.File, Output.Dir, PCA.Files,
#                       nbCPU = nbCPU, MaxRAM = MaxRAM)
#  
#  print("MAP ALPHA DIVERSITY")
#  # Index.Alpha   = c('Shannon','Simpson')
#  Index.Alpha   = c('Shannon')
#  map_alpha_div(Input.Image.File, Output.Dir, window_size,
#                      nbCPU = nbCPU, MaxRAM = MaxRAM, Index.Alpha = Index.Alpha)
#  
#  print("MAP BETA DIVERSITY")
#  map_beta_div(Input.Image.File, Output.Dir, window_size,
#                     nbCPU = nbCPU, MaxRAM = MaxRAM)

## ----alpha and beta diversity indices from vector layer------------------
#  # location of the spectral species raster needed for validation
#  TypePCA     = 'SPCA'
#  Dir.Raster  = file.path(Output.Dir,basename(Input.Image.File),TypePCA,'SpectralSpecies')
#  Name.Raster = 'SpectralSpecies'
#  Path.Raster = file.path(Dir.Raster,Name.Raster)
#  
#  # location of the directory where shapefiles used for validation are saved
#  vect        = system.file('extdata', 'VECTOR', package = 'biodivMapR')
#  Shannon.All = list() # ??
#  
#  # list vector data
#  Path.Vector         = list.shp(vect)
#  Name.Vector         = tools::file_path_sans_ext(basename(Path.Vector))
#  
#  # read raster data including projection
#  RasterStack         = stack(Path.Raster)
#  Projection.Raster   = projection.file(Path.Raster,'raster')
#  
#  # get alpha and beta diversity indicators corresponding to shapefiles
#  Biodiv.Indicators           = diversity_from_plots(Raster = Path.Raster, Plots = Path.Vector,NbClusters = nbclusters)
#  # if no name
#  Biodiv.Indicators$Name.Plot = seq(1,length(Biodiv.Indicators$Shannon[[1]]),by = 1)
#  Shannon.RS                  = c(Biodiv.Indicators$Shannon)[[1]]

## ----Write validation----------------------------------------------------
#  # write RS indicators
#  ####################################################
#  # write indicators for alpha diversity
#  Path.Results = file.path(Output.Dir, basename(Input.Image.File), TypePCA, 'VALIDATION')
#  dir.create(Path.Results, showWarnings = FALSE, recursive = TRUE)
#  ShannonIndexFile <- file.path(Path.Results, "ShannonIndex.tab")
#  write.table(Shannon.RS, file = ShannonIndexFile, sep = "\t", dec = ".", na = " ",
#              row.names = Biodiv.Indicators$Name.Plot, col.names= F, quote=FALSE)
#  
#  Results =  data.frame(Name.Vector, Biodiv.Indicators$Richness, Biodiv.Indicators$Fisher,                                Biodiv.Indicators$Shannon, Biodiv.Indicators$Simpson)
#  names(Results)  = c("ID_Plot", "Species_Richness", "Fisher", "Shannon", "Simpson")
#  write.table(Results, file = paste(Path.Results,"AlphaDiversity.tab",sep=''), sep="\t", dec=".",               na=" ", row.names = F, col.names= T,quote=FALSE)
#  
#  # write indicators for beta diversity
#  BC_mean = Biodiv.Indicators$BCdiss
#  colnames(BC_mean) = rownames(BC_mean) = Biodiv.Indicators$Name.Plot
#  write.table(BC_mean, file = paste(Path.Results,"BrayCurtis.csv",sep=''), sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
#  

