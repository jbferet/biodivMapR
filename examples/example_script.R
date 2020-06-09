# ===============================================================================
# biodivMapR
# ===============================================================================
# PROGRAMMERS:
#
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
#
# Copyright 2020/06 Jean-Baptiste FERET
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
# a ENVI format header file with information about spectral bands is expected:
# please make sure 'wavelength' and 'wavelength units' are defined
Input_Image_File <- system.file('extdata', 'RASTER', 'S2A_T33NUD_20180104_Subset', package = 'biodivMapR')

# full path for the Mask raster corresponding to image to process
# expected to be the same dimensions (columns and lines) as the image
# if possible in ENVI HDR format, 1 band, integer 8bits
# expected values in the raster: 0 = masked, 1 = selected
# set to FALSE if no mask available
Input_Mask_File <- FALSE

# Output directory: files created by script will be written there.
# For each image processed, a subdirectory will be created after its name
Output_Dir <- '~/biodiv' # fill with your own path
Output_Dir <- '../RESULTS_TEST15' # fill with your own path

# SPATIAL RESOLUTION
# resolution of spatial units for alpha and beta diversity maps (in pixels), relative to original image
# if Res.Map = 10 for images with 10 m spatial resolution, then spatial units will be 10 pixels x 10m = 100m x 100m surfaces
# rule of thumb: spatial units between 0.25 and 4 ha usually match with ground data
# too small window_size results in low number of pixels per spatial unit, hence limited range of variation of diversity in the image
window_size <- 10

# PCA FILTERING: 		Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- FALSE

# type of dimensionality reduction:
# PCA:  no rescaling of the data
# SPCA: rescaling of the data
# MNF:  minimum noise fraction
TypePCA <- 'SPCA'

# should continuum removal be performed on the image for spectral normalization?
# Continuum_Removal recommended for multi and hyperspectral data
Continuum_Removal <- TRUE

################################################################################
##                    DEFINE PARAMETERS FOR METHOD                            ##
################################################################################
nbCPU <- 4
MaxRAM <- 0.5
nbclusters <- 50

################################################################################
##                              PROCESS IMAGES                                ##
################################################################################
# 1- Filter data in order to discard non vegetated / shaded / cloudy pixels
NDVI_Thresh <- 0.5
Blue_Thresh <- 500
NIR_Thresh <- 1500
print("PERFORM RADIOMETRIC FILTERING")
Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File, Mask_Path = Input_Mask_File,
                                                 Output_Dir = Output_Dir, TypePCA = TypePCA,
                                                 NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,NIR_Thresh = NIR_Thresh)

# 2- Perform  dimensionality reduction
print("PERFORM DIMENSIONALITY REDUCTION")
PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File, Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir, TypePCA = TypePCA, FilterPCA=FilterPCA,
                          nbCPU = nbCPU, MaxRAM = MaxRAM, Continuum_Removal = Continuum_Removal)

PCA_Files <- PCA_Output$PCA_Files
Pix_Per_Partition <- PCA_Output$Pix_Per_Partition
nb_partitions <- PCA_Output$nb_partitions
Input_Mask_File <- PCA_Output$MaskPath
PCA_model <- PCA_Output$PCA_model
SpectralFilter <- PCA_Output$SpectralFilter

# 3- Select principal components from the PCA raster
# Sel_PC = path of the file where selected components are stored
Sel_PC <- select_PCA_components(Input_Image_File = Input_Image_File,
                                Output_Dir = Output_Dir, PCA_Files = PCA_Files,
                                TypePCA = TypePCA, File_Open = TRUE)

################################################################################
##    MAP ALPHA AND BETA DIVERSITY USING STANDARD SPECTRAL SPECIES APPROACH   ##
################################################################################
print("MAP SPECTRAL SPECIES")
map_spectral_species(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir,
                     PCA_Files = PCA_Files, PCA_model = PCA_model,
                     SpectralFilter = SpectralFilter, Input_Mask_File = Input_Mask_File,
                     Pix_Per_Partition = Pix_Per_Partition, nb_partitions = nb_partitions,
                     nbCPU = nbCPU, MaxRAM = MaxRAM, nbclusters = nbclusters, TypePCA = TypePCA,
                     Continuum_Removal = Continuum_Removal)

print("MAP ALPHA DIVERSITY")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')
map_alpha_div(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, TypePCA = TypePCA,
              window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha, nbclusters = nbclusters)

print("MAP BETA DIVERSITY")
map_beta_div(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, TypePCA = TypePCA,
             window_size = window_size, nb_partitions=nb_partitions, nbCPU = nbCPU, MaxRAM = MaxRAM,
             nbclusters = nbclusters)

################################################################################
##              MAP FUNCTIONAL DIVERSITY METRICS FRic, FEve, FDiv             ##
##          (Villeger et al, 2008 https://doi.org/10.1890/07-1206.1)          ##
################################################################################
## read selected features from dimensionality reduction
Selected_Features <- read.table(Sel_PC)[[1]]
## path for selected components
map_functional_div(Original_Image_File = Input_Image_File, Functional_File = PCA_Files,
                   Selected_Features = Selected_Features, Output_Dir = Output_Dir,
                   window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM,TypePCA = TypePCA)

# ################################################################################
# ## MAP PARTITIONING OF PLANT SPECTRAL DIVERSITY INTO ALPHA & BETA COMPONENTS  ##
# ##        (Laliberte et al, 2020 https://doi.org/10.1111/ele.13429)           ##
# ################################################################################
# PartitionSpectralDiversity <- map_partition_div(Original_Image_File = Input_Image_File, Partition_File = PCA_Files,
#                                                 Selected_Features = Selected_Features, Output_Dir = Output_Dir,
#                                                 window_size = window_size, TypePCA = TypePCA)

################################################################################
##                COMPUTE DIVERSITY METRICS FROM FIELD PLOTS                  ##
################################################################################
# location of the spectral species raster needed for validation
Output_SubDir <- tools::file_path_sans_ext(basename(Input_Image_File))
Dir_Raster <- file.path(Output_Dir,Output_SubDir,TypePCA,'SpectralSpecies')
Name_Raster <- 'SpectralSpecies'
Path_Raster <- file.path(Dir_Raster,Name_Raster)

# location of the directory where shapefiles used for validation are saved
vect <- system.file('extdata', 'VECTOR', package = 'biodivMapR')
Shannon_All <- list()

# list vector data
Path_Vector <- list_shp(vect)
Name_Vector <- tools::file_path_sans_ext(basename(Path_Vector))

# get alpha and beta diversity indicators corresponding to shapefiles
Biodiv_Indicators <- diversity_from_plots(Raster_SpectralSpecies = Path_Raster, Plots = Path_Vector,
                                          NbClusters = nbclusters, Raster_Functional = PCA_Files, Selected_Features = Selected_Features)
# if no name
Biodiv_Indicators$Name_Plot <- seq(1,length(Biodiv_Indicators$Shannon[[1]]),by = 1)
Shannon_RS <- c(Biodiv_Indicators$Shannon)[[1]]
FRic <- c(Biodiv_Indicators$FunctionalDiversity$FRic)
FEve <- c(Biodiv_Indicators$FunctionalDiversity$FEve)
FDiv <- c(Biodiv_Indicators$FunctionalDiversity$FDiv)

####################################################
# write RS indicators
####################################################
Path_Results <- file.path(Output_Dir,Output_SubDir,TypePCA,'VALIDATION')
dir.create(Path_Results, showWarnings = FALSE,recursive = TRUE)
write.table(Shannon_RS, file = file.path(Path_Results,"ShannonIndex.csv"),
            sep="\t", dec=".", na=" ", row.names = Biodiv_Indicators$Name_Plot, col.names= F,quote=FALSE)

Results <- data.frame(Name_Vector, Biodiv_Indicators$Richness, Biodiv_Indicators$Fisher,
                      Biodiv_Indicators$Shannon, Biodiv_Indicators$Simpson,
                      Biodiv_Indicators$FunctionalDiversity$FRic,
                      Biodiv_Indicators$FunctionalDiversity$FEve,
                      Biodiv_Indicators$FunctionalDiversity$FDiv)

names(Results)  = c("ID_Plot", "Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
write.table(Results, file = paste(Path_Results,"AlphaDiversity.csv",sep=''),
            sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)

# write indicators for beta diversity
BC_mean <- Biodiv_Indicators$BCdiss
colnames(BC_mean) <- rownames(BC_mean) <- Biodiv_Indicators$Name_Plot
write.table(BC_mean, file = paste(Path_Results,"BrayCurtis.csv",sep=''),
            sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)

####################################################
# produce plots
####################################################
# apply ordination using PCoA (same as done for map_beta_div)
library(labdsv)
MatBCdist <- as.dist(BC_mean, diag = FALSE, upper = FALSE)
BetaPCO <- pco(MatBCdist, k = 3)

# Assign vegetation type to polygons in shapefiles
nbSamples <- c(6,4,7,7)
vg <- c('Forest high diversity', 'Forest low diversity', 'Forest medium diversity', 'low vegetation')
Type_Vegetation <- c()
for (i in 1: length(nbSamples)){
  for (j in 1:nbSamples[i]){
    Type_Vegetation <- c(Type_Vegetation,vg[i])
  }
}

# create data frame including alpha and beta diversity
library(ggplot2)
Results <- data.frame('vgtype'=Type_Vegetation,'pco1'= BetaPCO$points[,1],'pco2'= BetaPCO$points[,2],'pco3' = BetaPCO$points[,3],
                      'shannon'=Shannon_RS,'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)

# plot field data in the PCoA space, with size corresponding to shannon index
ggplot(Results, aes(x=pco1, y=pco2, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'BetaDiversity_PcoA1_vs_PcoA2.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)

ggplot(Results, aes(x=pco1, y=pco3, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'BetaDiversity_PcoA1_vs_PcoA3.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)

ggplot(Results, aes(x=pco2, y=pco3, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'BetaDiversity_PcoA2_vs_PcoA3.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)

# plot field data in the PCoA space, with size corresponding to FRic
ggplot(Results, aes(x=pco1, y=pco2, color=vgtype,size=FRic)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'BetaDiversity_PcoA1_vs_PcoA2_FRic.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)


ggplot(Results, aes(x=pco1, y=pco3, color=vgtype,size=FRic)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'BetaDiversity_PcoA1_vs_PcoA3_FRic.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)

ggplot(Results, aes(x=pco2, y=pco3, color=vgtype,size=FRic)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'BetaDiversity_PcoA2_vs_PcoA3_FRic.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)

# Compare Shannon and FRic
ggplot(Results, aes(x=shannon, y=FRic, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'Shannon_Vs_FRic.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)

ggplot(Results, aes(x=shannon, y=FEve, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'Shannon_Vs_FEve.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)

ggplot(Results, aes(x=shannon, y=FDiv, color=vgtype,size=shannon)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#e6140a", "#e6d214", "#e68214", "#145ae6"))
filename <- file.path(Path_Results,'Shannon_Vs_FDiv.png')
ggsave(filename, plot = last_plot(), device = 'png', path = NULL,
       scale = 1, width = 6, height = 4, units = "in",
       dpi = 600, limitsize = TRUE)
