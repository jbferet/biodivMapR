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
