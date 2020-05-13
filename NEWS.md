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
