# ===============================================================================
# biodivMapR
# ===============================================================================
# PROGRAMMERS:
#
# Jean-Baptiste FERET <jb.feret@irstea.fr>
#
# Copyright 2019/06 Jean-Baptiste FERET
#
# biodivMapR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#
# ===============================================================================
library(raster)
library(biodivMapR)
################################################################################
##              DEFINE PARAMETERS FOR DATASET TO BE PROCESSED                 ##
################################################################################
# path (absolute or relative) for the image to process
# expected to be in ENVI HDR format, BIL interleaved
Input_Image_File  = system.file('extdata', 'RASTER', 'S2A_T33NUD_20180104_Subset', package = 'biodivMapR')

# # convert the image using Convert.Raster2BIL if not in the proper format
# Input_Image_File  = raster2BIL(Raster_Path = Input_Image_File,
#                                        Sensor = 'SENTINEL_2A',
#                                        Convert_Integer = TRUE,
#                                        Output_Dir = 'BIL_Raster')

# full path for the Mask raster corresponding to image to process
# expected to be in ENVI HDR format, 1 band, integer 8bits
# expected values in the raster: 0 = masked, 1 = selected
# set to FALSE if no mask available
Input_Mask_File   = FALSE

# relative or absolute path for the Directory where results will be stored
# For each image processed, a subdirectory will be created after its name
Output_Dir        = 'RESULTS'

# SPATIAL RESOLUTION
# resolution of spatial units for alpha and beta diversity maps (in pixels), relative to original image
# if Res.Map = 10 for images with 10 m spatial resolution, then spatial units will be 10 pixels x 10m = 100m x 100m surfaces
# rule of thumb: spatial units between 0.25 and 4 ha usually match with ground data
# too small window_size results in low number of pixels per spatial unit, hence limited range of variation of diversity in the image
window_size       = 10

# PCA FILTERING: 		Set to TRUE if you want second filtering based on PCA outliers to be processed. Slower
FilterPCA         = FALSE

# type of PCA:
# PCA: no rescaling of the data
# SPCA: rescaling of the data
TypePCA     = 'SPCA'

################################################################################
##                    DEFINE PARAMETERS FOR METHOD                            ##
################################################################################
nbCPU       = 2
MaxRAM      = 0.5
nbclusters  = 50

################################################################################
##                              PROCESS IMAGES                                ##
################################################################################
# 1- Filter data in order to discard non vegetated / shaded / cloudy pixels
NDVI_Thresh = 0.5
Blue_Thresh = 500
NIR_Thresh  = 1500
print("PERFORM RADIOMETRIC FILTERING")
Input_Mask_File = perform_radiometric_filtering(Input_Image_File,Input_Mask_File,Output_Dir,
                                                NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                NIR_Thresh = NIR_Thresh)

# 2- Compute PCA for a random selection of pixels in the raster
print("PERFORM PCA ON RASTER")
PCA_Output        = perform_PCA(Input_Image_File,Input_Mask_File,Output_Dir,FilterPCA=FilterPCA,nbCPU=nbCPU,MaxRAM = MaxRAM)
PCA_Files         = PCA_Output$PCA_Files
Pix_Per_Partition = PCA_Output$Pix_Per_Partition
nb_partitions     = PCA_Output$nb_partitions
Input_Mask_File   = PCA_Output$MaskPath
PCA_model         = PCA_Output$PCA_model
SpectralFilter    = PCA_Output$SpectralFilter

# 3- Select principal components from the PCA raster
select_PCA_components(Input_Image_File,Output_Dir,PCA_Files,File_Open = TRUE)

################################################################################
##                      MAP ALPHA AND BETA DIVERSITY                          ##
################################################################################

print("MAP SPECTRAL SPECIES")
map_spectral_species(Input_Image_File,Output_Dir,PCA_Files,PCA_model,SpectralFilter,Input_Mask_File,Pix_Per_Partition,nb_partitions,nbCPU=nbCPU,MaxRAM=MaxRAM)

print("MAP ALPHA DIVERSITY")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha   = c('Shannon')
map_alpha_div(Input_Image_File,Output_Dir,window_size,nbCPU=nbCPU,MaxRAM=MaxRAM,Index_Alpha = Index_Alpha)

print("MAP BETA DIVERSITY")
map_beta_div(Input_Image_File,Output_Dir,window_size,nb_partitions=nb_partitions,nbCPU=nbCPU,MaxRAM=MaxRAM)

################################################################################
##          COMPUTE ALPHA AND BETA DIVERSITY FROM FIELD PLOTS                 ##
################################################################################

# location of the spectral species raster needed for validation
Dir.Raster  = file.path(Output_Dir,basename(Input_Image_File),TypePCA,'SpectralSpecies')
Name.Raster = 'SpectralSpecies'
Path.Raster = file.path(Dir.Raster,Name.Raster)

# location of the directory where shapefiles used for validation are saved
vect        = system.file('extdata', 'VECTOR', package = 'biodivMapR')
Shannon.All = list()

# list vector data
Path.Vector         = list_shp(vect)
Name.Vector         = tools::file_path_sans_ext(basename(Path.Vector))

# get alpha and beta diversity indicators corresponding to shapefiles
Biodiv.Indicators           = diversity_from_plots(Raster = Path.Raster, Plots = Path.Vector,NbClusters = nbclusters)
# if no name
Biodiv.Indicators$Name_Plot = seq(1,length(Biodiv.Indicators$Shannon[[1]]),by = 1)
Shannon.RS                  = c(Biodiv.Indicators$Shannon)[[1]]

####################################################
# write RS indicators
####################################################
# write indicators for alpha diversity
Path.Results = paste(Output_Dir,'/',basename(Input_Image_File),'/',TypePCA,'/VALIDATION/',sep='')
dir.create(Path.Results, showWarnings = FALSE,recursive = TRUE)
write.table(Shannon.RS, file = paste(Path.Results,"ShannonIndex.csv",sep=''), sep="\t", dec=".", na=" ", row.names = Biodiv.Indicators$Name_Plot, col.names= F,quote=FALSE)

Results         =  data.frame(Name.Vector, Biodiv.Indicators$Richness, Biodiv.Indicators$Fisher, Biodiv.Indicators$Shannon,Biodiv.Indicators$Simpson)
names(Results)  = c("ID_Plot", "Species_Richness", "Fisher", "Shannon", "Simpson")
write.table(Results, file = paste(Path.Results,"AlphaDiversity.csv",sep=''), sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)

# write indicators for beta diversity
BC_mean = Biodiv.Indicators$BCdiss
colnames(BC_mean) = rownames(BC_mean) = Biodiv.Indicators$Name_Plot
write.table(BC_mean, file = paste(Path.Results,"BrayCurtis.csv",sep=''), sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)


####################################################
# illustrate results
####################################################
# apply ordination using PCoA (same as done for map_beta_div)
library(labdsv)
MatBCdist = as.dist(BC_mean, diag = FALSE, upper = FALSE)
BetaPCO   = pco(MatBCdist, k = 3)

# very uglily assign vegetation type to polygons in shapefiles
nbSamples = c(6,4,7,7)
vg        = c('Forest high diversity', 'Forest low diversity', 'Forest medium diversity', 'low vegetation')
Type_Vegetation = c()
for (i in 1: length(nbSamples)){
  for (j in 1:nbSamples[i]){
    Type_Vegetation = c(Type_Vegetation,vg[i])
  }
}

# create data frame including alpha and beta diversity
library(ggplot2)
Results     =  data.frame('vgtype'=Type_Vegetation,'pco1'= BetaPCO$points[,1],'pco2'= BetaPCO$points[,2],'pco3' = BetaPCO$points[,3],'shannon'=Shannon.RS)

# plot field data in the PCoA space, with size corresponding to shannon index
ggplot(Results, aes(x=pco1, y=pco2, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename = file.path(Path.Results,'BetaDiversity_PcoA1_vs_PcoA2.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 600, limitsize = TRUE)


ggplot(Results, aes(x=pco1, y=pco3, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename = file.path(Path.Results,'BetaDiversity_PcoA1_vs_PcoA3.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 600, limitsize = TRUE)

ggplot(Results, aes(x=pco2, y=pco3, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename = file.path(Path.Results,'BetaDiversity_PcoA2_vs_PcoA3.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 600, limitsize = TRUE)
