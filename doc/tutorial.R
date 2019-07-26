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
#  FilterPCA = FALSE

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

## ----Spectral species map------------------------------------------------
#  print("MAP SPECTRAL SPECIES")
#  map_spectral_species(Input.Image.File, Output.Dir, PCA.Files,
#                       nbCPU = nbCPU, MaxRAM = MaxRAM)

## ----alpha and beta diversity maps---------------------------------------
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
#  Projection.Raster   = get_projection(Path.Raster,'raster')
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

## ----PCoA on Field Plots-------------------------------------------------
#  # apply ordination using PCoA (same as done for map_beta_div)
#  library(labdsv)
#  MatBCdist = as.dist(BC_mean, diag = FALSE, upper = FALSE)
#  BetaPCO   = pco(MatBCdist, k = 3)
#  

## ----plot PCoA & Shannon-------------------------------------------------
#  # very uglily assign vegetation type to polygons in shapefiles
#  nbSamples = c(6,4,7,7)
#  vg        = c('Forest high diversity', 'Forest low diversity', 'Forest medium diversity', 'low vegetation')
#  Type_Vegetation = c()
#  for (i in 1: length(nbSamples)){
#    for (j in 1:nbSamples[i]){
#      Type_Vegetation = c(Type_Vegetation,vg[i])
#    }
#  }
#  
#  # create data frame including alpha and beta diversity
#  library(ggplot2)
#  Results     =  data.frame('vgtype'=Type_Vegetation,'pco1'= BetaPCO$points[,1],'pco2'= BetaPCO$points[,2],'pco3' = BetaPCO$points[,3],'shannon'=Shannon.RS)
#  
#  # plot field data in the PCoA space, with size corresponding to shannon index
#  ggplot(Results, aes(x=pco1, y=pco2, color=vgtype,size=shannon)) +
#    geom_point(alpha=0.6) +
#    scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
#  filename = file.path(Path.Results,'BetaDiversity_PcoA1_vs_PcoA2.png')
#  ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
#         scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
#         dpi = 600, limitsize = TRUE)
#  
#  
#  ggplot(Results, aes(x=pco1, y=pco3, color=vgtype,size=shannon)) +
#    geom_point(alpha=0.6) +
#    scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
#  filename = file.path(Path.Results,'BetaDiversity_PcoA1_vs_PcoA3.png')
#  ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
#         scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
#         dpi = 600, limitsize = TRUE)
#  
#  ggplot(Results, aes(x=pco2, y=pco3, color=vgtype,size=shannon)) +
#    geom_point(alpha=0.6) +
#    scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
#  filename = file.path(Path.Results,'BetaDiversity_PcoA2_vs_PcoA3.png')
#  ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
#         scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
#         dpi = 600, limitsize = TRUE)
#  

