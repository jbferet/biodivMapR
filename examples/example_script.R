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

################################################################################
##              DEFINE PARAMETERS FOR DATASET TO BE PROCESSED                 ##
################################################################################
# path (absolute or relative) for the image to process
# expected to be in ENVI HDR format, BIL interleaved
Input.Image.File  = system.file('extdata', 'RASTER', 'S2A_T33NUD_20180104_Subset', package = 'biodivMapR')
check_data(Input.Image.File)

# convert the image using Convert.Raster2BIL if not in the proper format
Input.Image.File  = raster2BIL(Raster.Path = Input.Image.File,
                                       Sensor = 'SENTINEL_2A',
                                       Convert.Integer = TRUE,
                                       Output.Directory = '~/test')

# full path for the Mask raster corresponding to image to process
# expected to be in ENVI HDR format, 1 band, integer 8bits
# expected values in the raster: 0 = masked, 1 = selected
# set to FALSE if no mask available
Input.Mask.File   = FALSE

# relative or absolute path for the Directory where results will be stored
# For each image processed, a subdirectory will be created after its name
Output.Dir        = 'RESULTS'

# SPATIAL RESOLUTION
# resolution of spatial units for alpha and beta diversity maps (in pixels), relative to original image
# if Res.Map = 10 for images with 10 m spatial resolution, then spatial units will be 10 pixels x 10m = 100m x 100m surfaces
# rule of thumb: spatial units between 0.25 and 4 ha usually match with ground data
# too small Spatial.Res results in low number of pixels per spatial unit, hence limited range of variation of diversity in the image
Spatial.Res       = 10

# PCA FILTERING: 		Set to TRUE if you want second filtering based on PCA outliers to be processed. Slower
FilterPCA         = TRUE

################################################################################
##      Check if the image format is compatible with codes (ENVI BIL)         ##
################################################################################
check_data(Input.Image.File)

################################################################################
##                    DEFINE PARAMETERS FOR METHOD                            ##
################################################################################
nbCPU         = 4
MaxRAM        = 0.5
nbclusters    = 50

################################################################################
##                              PROCESS IMAGES                                ##
################################################################################
# 1- Filter data in order to discard non vegetated / shaded / cloudy pixels
NDVI.Thresh = 0.5
Blue.Thresh = 500
NIR.Thresh  = 1500
print("PERFORM RADIOMETRIC FILTERING")
ImPathShade         = Perform.Radiometric.Filtering(Input.Image.File,Input.Mask.File,Output.Dir,
                                                    NDVI.Thresh = NDVI.Thresh, Blue.Thresh = Blue.Thresh,
                                                    NIR.Thresh = NIR.Thresh)

# 2- Compute PCA for a random selection of pixels in the raster
print("PERFORM PCA ON RASTER")
PCA.Files           = Perform.PCA.Image(Input.Image.File,ImPathShade,Output.Dir,FilterPCA=TRUE,nbCPU=nbCPU,MaxRAM = MaxRAM)

# 3- Select principal components from the PCA raster
Select.Components(Input.Image.File,Output.Dir,PCA.Files,File.Open = TRUE)

################################################################################
##                      MAP ALPHA AND BETA DIVERSITY                          ##
################################################################################

print("MAP SPECTRAL SPECIES")
Map.Spectral.Species(Input.Image.File,Output.Dir,PCA.Files,nbCPU=nbCPU,MaxRAM=MaxRAM)


print("MAP ALPHA DIVERSITY")
# Index.Alpha   = c('Shannon','Simpson')
Index.Alpha   = c('Shannon')
map_alpha_div(Input.Image.File,Output.Dir,Spatial.Res,nbCPU=nbCPU,MaxRAM=MaxRAM,Index.Alpha = Index.Alpha)

print("MAP BETA DIVERSITY")
map_beta_div(Input.Image.File,Output.Dir,Spatial.Res,nbCPU=nbCPU,MaxRAM=MaxRAM)

################################################################################
##          COMPUTE ALPHA AND BETA DIVERSITY FROM FIELD PLOTS                 ##
################################################################################

# location of the spectral species raster needed for validation
TypePCA     = 'SPCA'
Dir.Raster  = file.path(Output.Dir,basename(Input.Image.File),TypePCA,'SpectralSpecies')
Name.Raster = 'SpectralSpecies'
Path.Raster = file.path(Dir.Raster,Name.Raster)

# location of the directory where shapefiles used for validation are saved
vect        = system.file('extdata', 'VECTOR', package = 'biodivMapR')
Shannon.All = list()

# list vector data
Path.Vector         = list.shp(vect)
Name.Vector         = tools::file_path_sans_ext(basename(Path.Vector))

# read raster data including projection
RasterStack         = stack(Path.Raster)
Projection.Raster   = projection.file(Path.Raster,'raster')

# get alpha and beta diversity indicators corresponding to shapefiles
Biodiv.Indicators           = diversity_from_plots(Raster = Path.Raster, Plots = Path.Vector,NbClusters = nbclusters)
# if no name
Biodiv.Indicators$Name.Plot = seq(1,length(Biodiv.Indicators$Shannon[[1]]),by = 1)
Shannon.RS                  = c(Biodiv.Indicators$Shannon)[[1]]

####################################################
# write RS indicators
####################################################
# write indicators for alpha diversity
Path.Results = paste(Output.Dir,'/',basename(Input.Image.File),'/',TypePCA,'/VALIDATION/',sep='')
dir.create(Path.Results, showWarnings = FALSE,recursive = TRUE)
write.table(Shannon.RS, file = paste(Path.Results,"ShannonIndex.csv",sep=''), sep="\t", dec=".", na=" ", row.names = Biodiv.Indicators$Name.Plot, col.names= F,quote=FALSE)

Results         =  data.frame(Name.Vector, Biodiv.Indicators$Richness, Biodiv.Indicators$Fisher, Biodiv.Indicators$Shannon,Biodiv.Indicators$Simpson)
names(Results)  = c("ID_Plot", "Species_Richness", "Fisher", "Shannon", "Simpson")
write.table(Results, file = paste(Path.Results,"AlphaDiversity.csv",sep=''), sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)

# write indicators for beta diversity
BC_mean = Biodiv.Indicators$BCdiss
colnames(BC_mean) = rownames(BC_mean) = Biodiv.Indicators$Name.Plot
write.table(BC_mean, file = paste(Path.Results,"BrayCurtis.csv",sep=''), sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
