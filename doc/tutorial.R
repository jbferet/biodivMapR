## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval=FALSE
)

## ----Input / Output files------------------------------------------------
#  Input_Image_File  = system.file('extdata', 'RASTER', 'S2A_T33NUD_20180104_Subset', package = 'biodivMapR')
#  
#  # Input.Image.File  = raster2BIL(Raster.Path = Input.Image.File,
#  #                                        Sensor = 'SENTINEL_2A',
#  #                                        Convert.Integer = TRUE,
#  #                                        Output.Directory = '~/test')
#  
#  Input_Mask_File   = FALSE
#  
#  Output_Dir        = 'RESULTS'

## ----Spatial resolution--------------------------------------------------
#  window_size = 10

## ----PCA filtering-------------------------------------------------------
#  FilterPCA = FALSE

## ----Computing options---------------------------------------------------
#  nbCPU         = 2
#  MaxRAM        = 0.5
#  nbclusters    = 50

## ----Mask non vegetated / shaded / cloudy pixels-------------------------
#  NDVI_Thresh = 0.5
#  Blue_Thresh = 500
#  NIR_Thresh  = 1500
#  print("PERFORM RADIOMETRIC FILTERING")
#  Input_Mask_File = perform_radiometric_filtering(Input_Image_File, Input_Mask_File, Output_Dir,
#                                              NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
#                                              NIR_Thresh = NIR_Thresh)

## ----PCA-----------------------------------------------------------------
#  print("PERFORM PCA ON RASTER")
#  PCA_Output        = perform_PCA(Input_Image_File, Input_Mask_File, Output_Dir,
#                                 FilterPCA = TRUE, nbCPU = nbCPU,MaxRAM = MaxRAM)
#  # path for the PCA raster
#  PCA_Files         = PCA_Output$PCA_Files
#  # number of pixels used for each partition used for k-means clustering
#  Pix_Per_Partition = PCA_Output$Pix_Per_Partition
#  # number of partitions used for k-means clustering
#  nb_partitions     = PCA_Output$nb_partitions
#  # path for the updated mask
#  Input_Mask_File   = PCA_Output$MaskPath
#  # parameters of the PCA model
#  PCA_model         = PCA_Output$PCA_model
#  # definition of spectral bands to be excluded from the analysis
#  SpectralFilter    = PCA_Output$SpectralFilter
#  
#  print("Select PCA components for diversity estimations")
#  select_PCA_components(Input_Image_File, Output_Dir, PCA_Files, File_Open = TRUE)

## ----Spectral species map------------------------------------------------
#  print("MAP SPECTRAL SPECIES")
#  map_spectral_species(Input_Image_File, Output_Dir, PCA_Files,PCA_model, SpectralFilter, Input_Mask_File,
#                       Pix_Per_Partition, nb_partitions, nbCPU=nbCPU, MaxRAM=MaxRAM)

## ----alpha and beta diversity maps---------------------------------------
#  print("MAP ALPHA DIVERSITY")
#  # Index.Alpha   = c('Shannon','Simpson')
#  Index_Alpha   = c('Shannon')
#  map_alpha_div(Input_Image_File, Output_Dir, window_size,
#                nbCPU=nbCPU, MaxRAM=MaxRAM, Index_Alpha = Index_Alpha)
#  
#  print("MAP BETA DIVERSITY")
#  map_beta_div(Input_Image_File, Output_Dir, window_size, nb_partitions=nb_partitions,
#               nbCPU=nbCPU, MaxRAM=MaxRAM)

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
#  Path.Vector         = list_shp(vect)
#  Name.Vector         = tools::file_path_sans_ext(basename(Path.Vector))
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

