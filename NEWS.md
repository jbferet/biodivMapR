# biodivMapR2 v2.3.4
## fix
- sample_from_raster adapted to handle distances when crs of input_rast in deg

# biodivMapR2 v2.3.3
## addition
- modified perform_PCA to allow PCA over any type of data

## change
- updated documentation

# biodivMapR2 v2.3.2
## change
- optimize tile processing with function get_raster_diversity_tile 

# biodivMapR2 v2.3.1
## fix
- fix tile mask management when tiles are empty

## change
- change management of parallel processing for steps where it is not optimal

# biodivMapR2 v2.3.0
## change
- added functions to process large rasters divided into tiles using get_s2_tiling function from preprocS2 package
- added functions to apply biodivMapR using moving window

# biodivMapR2 v2.2.1
## change
- clean code
- add suppressWarnings when calling progressr

# biodivMapR2 v2.2.0
## change
- merge v2 from gitlab repository
- update vignettes for use of v2 

# biodivMapR2 v2.1.4
## fix
- in 'get_diversity_from_plots' : 
- initialize Attributes with nrow = nbPlots_init
- initialize functional diversity metrics in Attributes with NA
- correct get_diversity_from_plots : assign to beta diversity values instead of mean Hill
- correct extract_vect_from_rast : update AttributeTable and discard vectors with no data

## addition
- create functions init_kmeans_samples and init_PCoA_samples to process samples 
previously extracted independently from an input raster
- possibility to provide name for output files in addition to directory

## change
- possibility to define updated mask file name as input for radiometric_filtering

# biodivMapR2 v2.1.3
## fix
- correct get_diversity_from_plots: output for FDis corrected to FDis instead of FDiv
- correct spectralspecies_per_polygon : AttributeTable not set to NULL
## addition
- computes all functional diversity when running get_diversity_from_plots 

# biodivMapR2 v2.1.2
## change
- use package fundiversity to compute diversity metrics

# biodivMapR v2.1.1
## Fix
- Fixed testthat
- fixed maps for functional diversity metrics (center reduction of trait space)

## Change
- more systematic use of 'seq_len'


# biodivMapR v2.1.0
## Addition
- added Hill numbers
- added functional diversity

## fix
- initialize input_mask to NULL in biodivMapR_Full
- handle chunk with pixels valid but no window valid

# biodivMapR v2.0.0
## Fix
- fixed SFS function when using SpatVector instead of SpatVectorCollection

# biodivMapR v1.14.3
## Fix
- fixed SFS function when using SpatVector instead of SpatVectorCollection

# biodivMapR v1.14.2

## Fix

- update fixed problem occurring when performing PCA without mask

# biodivMapR v1.14.1

## Fix

- fixed problem occurring when performing PCA without mask
- fixed issue related to plan(multisession), suggested by @rnedelec on github, 
following documented issue in https://github.com/ropensci/drake/issues/255#issuecomment-365799680

## Addition

- option added to display progress bar. Set to FALSE as default

# biodivMapR v1.14.0

## Change

- added input variables Kmeans_info_path and Kmeans_info in map_spectral_species to allow production from previously computed clusters

# biodivMapR v1.13.0

## Change

- allow use of any vector format (including gpkg, shp... ) for validation plots provided in diversity_from_plots
- allow different crs between raster and vector layers

# biodivMapR v1.12.1

## Fix

- update example_script and vignettes to fulfill removal to rgdal dependency in the package

# biodivMapR v1.12.0

## Fix

- issue #18: remove the dependency to rgdal and rgeos using terra and gdalUtilities instead

## Change

- remove extract_pixels_coordinates.From.OGR, a unique extract_pixels_coordinates now works for both character or rast/vect inputs

## Addition

- add get_gdal_info to get raster gdal info as a nested list.

# biodivMapR v1.11.1

## Change
- possibility to enable / disable progressbar in compute_spectral_species_FieldPlots and init_kmeans

# biodivMapR v1.11.0

## Fix
- fixed computation of FD metrics
- removed import progress and emstreeR

# biodivMapR v1.10.6

## Fix
- removed reference to multiprocess, multisession used exclusively from future

# biodivMapR v1.10.5

## Fix
- fixed bug when producing simpson index from map_alpha_div function
- added importFrom corresponding to progressr in continum removal

# biodivMapR v1.10.4

## Changes
- removed unused functions (VectorInRasterFootprint, get_BB, get_BB_from_fullImagem get_BB_from_Vector, get_polygonCoord_from_Shp, read_ListRasters, VectorInRasterFootprint, gdal_polygonizeR) in order to reduce dependency to rgeos
- these functions can be found in R package preprocS2

## Fix
- fixed bug when computing functional diversity metrics from 1 component only

# biodivMapR v1.10.3

## Changes
- function map_spectral_species checks if format and dimensions of Input_Mask_File are as expected
- function map_spectral_species now supports Input_Mask_File = FALSE. It writes a blank mask to ensure next processes

# biodivMapR v1.10.2

## Fix
- fixed bug when defining SelectedPCs as a vector in map_spectral_species()
- Updated vignettes

# biodivMapR v1.10.1

## Fix
- Updated alpha and beta diversity mapping when starting from a classification map

# biodivMapR v1.10.0

## Fix
- now uses red band instead of red edge band for the computation of NDVI
--> changes thresholding

## changes
- simplified inputs for functions such as map_spectral_species(), init_kmeans(), map_beta_div()
- updated vignettes

# biodivMapR v1.9.11

## changes
- optimized codes for the computation of spectral species, alpha and beta diversity
- addition of progress bars during the different steps of the computation of spectral species, alpha and beta diversity maps
- addition of functions Compute_ALPHA_SSD_per_window_list, Compute_ALPHA_SSD_per_window, Compute_ALPHA_per_window, prepare_HDR_SSD, prepare_HDR_Sunlit, RW_bytes_all, RW_bytes
- use multisession instead of multiprocess where it was forgotten

# biodivMapR v1.9.10

## Addition
- use of list for computation of beta diversity
- addition of functions getBCdiss and Normalize_SSD

# biodivMapR v1.9.9

## Addition
- added Hellinger distance as beta diversity metric

# biodivMapR v1.9.8

## Addition
- removed NMDS as possible ordination method when computing beta diversity

# biodivMapR v1.9.7

## Addition
- progress bar instead of messages
- future: multisession instead of multiprocess

# biodivMapR v1.9.6.1

## Addition
- added tutorials for the PROGYSAT workshop

# biodivMapR v1.9.6

## Fixes
- install package dissUtils directly from github ('cran/dissUtils') as it was removed from official CRAN repo

# biodivMapR v1.9.5

## Fixes
- discard marginal spectral species based on number of sunlit pixels, instead of total number of pixels considered (window size or plot used for validation)

## Changes
- harmonize default value for pcelim

# biodivMapR v1.9.4

## Fixes
- re-integrate package emstreeR (previously removed from CRAN) used to compute functional divergence

# biodivMapR v1.9.3

## Fixes
- uses future_lapply only if more than one CPU requested, otherwise use standard lapply

# biodivMapR v1.9.2

## Fixes
- fixed data type of spectral species derived from supervised classification and output directories for alpha and beta diversity

# biodivMapR v1.9.1

## Fixes
- fixed error occuring when using custom classification map instead of spectral species map

# biodivMapR v1.9.0

## Changes
- Added dimMDS defining number of dimensions to run PCoA for beta diversity
- exported functions from beta library
- updated documentation

## Fixes
- corrected weighted distance from nearest neighbors: assign exact coordinates of a sample when dissimilarity = 0
- temporary: discarded evenness from functional diversity as it uses emstreeR and requires binding to mlpack

# biodivMapR v1.8.0

## Changes
- Added functionality to allow for computation of spectral diversity maps based on spectral index stack or any raster stack
- added function to compute interquartile range (IQR) and identify outliers
- applied IQR instead of center reduction
- exported most of the functions in the Lib_ImageProcess
- Image sample was removed and placed in external repository
- Major update of vignettes for tutorial 

## Fixes
- used file.path instead of paste

# biodivMapR v1.7.0

## Changes
- capacity to produce diversity maps based on classification raster

## Fixes
- warning eliminated when testing Selected_Features in Lib_MapFunctionalDiversity.R 
- Corrected data extraction function using nbands in extract.big_raster and get_random_subset_from_image
- updated extract.big_raster and get_random_subset_from_image in order to account for 2D rasters which cannot be read with brick
- added driver definition when using read_stars (extract.big_raster)
- corrected vignettes (TypePCA used before defined in previous version) 

# biodivMapR v1.6.2

## Fixes
- added importFrom raster brick in function extract.big_raster

# biodivMapR v1.6.1

## Fixes
- Fixed bug when calling file.edit from linux terminal. initial parameter editor='internal' is not cross platforms. removed it. 

# biodivMapR v1.6.0

## Changes
- Added option to directly compute alpha and beta diversity maps from classification maps, even if it does not correspond to SpectralSpecies file produced from biodivMapR

## Fixes
- Fixed bug occuring when calling function diversity_from_plots if no functional duversity map was produced before

# biodivMapR v1.5.1

## Fixes
- Fixed bug occuring when writing image if initial raster is not a multiple of the window size. no bug occur but the raster files have wrong information
- fixed problem when performing estimation of biodiversity for plots outside of the raster: now the value is NA

# biodivMapR v1.5.0

## Changes
- Added functional diversity metrics
- updated example script
- updated vignette

# biodivMapR v1.4.0

## Changes
- implemented MNF
- discarded HDR as input variable from get_random_subset_from_image
- updated example file
- added contribution of F de Boissieu
- changed email address to teledetection.fr
- updated diversity_from_plots
- updated example script

# biodivMapR v1.3.1

## Fixes
- fixed identification and elimination of pixel samples with NA
- fixed bug occuring when input raster is BSQ interleave, and BSQ reported in PC file and following files

## Changes
- finalized preparation for MNF
- improved information in the header files

# biodivMapR v1.3.0

## Changes
- integrated stars package in order to read any file format, including TIFF format
- developed a generic function to write rasters
- prepared for MNF
- changed default red band for the computation of NDVI: closest band to 690 nm is now selected instead of closest band to 700 nm

# biodivMapR v1.2.1

## Changes
- added option to convert the mask file into proper format with raster2BIL

# biodivMapR v1.2.0

## Fixes
- fixed bug raised when processing data over large number of pixels (image products >2^31 bits)
- fixed bug by adding elimination of bands disturbed by water vapor even when ContinuumRemoval set to FALSE

## Changes
- Changed name of default directory when saving image after calling raster2BIL
- added documentation for raster data conversion using raster2BIL

# biodivMapR v1.1.0

## Fixes
- fixed bug in Lib_MapSpectralSpecies: remove constant bands from Subset$DataSubset if needed
- fix bug when continuum removal applied on pixels with constant values which were not filtered in previous stages

## Changes
- Updated continuum removal
- Named Continuum_Removal instead of CR
- corrected documetation for perform_radiometric_filtering
- Added an error report when spectral information is not conform with expectations (due to too high noise level for example, usually leading to NaN or Inf values after PCA)
- Changed name for some variables in internal functions for consistency
- will add documentation on how to analyze error reports

# biodivMapR v1.0.2

## Fixes

## Changes
- moved examples to repository root: example files are not installed with package anymore.
- removed Plots reprojection from `diversity_from_plots`: 
now, Plots and Raster must be in the same projection. Changed example Plots projection accordingly.
- removed `get_projection` (useless)

# biodivMapR v1.0.1 (Release date: 2019-09-27)

- Added NEWS.md
- Updated README.md: transfered from gitlab.irstea to github
- changed return() into return(invisible())
- updated vignettes & tutorial with latest outputs and figures from code

# biodivMapR v1.0.0 (Release date: 2019-09-26)

First release in GitHub
Submission accepted to Methods in Ecology and Evolution
