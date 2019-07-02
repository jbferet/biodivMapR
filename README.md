# DiversityMappR: an R package for α- and β-diversity mapping using remotely-sensed images

# 1 Install
After installing packages devtools and getPass, `DiversityMappR` can then be installed with the folloqing command line in R session, where `uname` is your gitlab.irstea.fr username:
```
devtools::install_git('https://gitlab.irstea.fr/jean-baptiste.feret/diversitymappr.git',
                      credentials = git2r::cred_user_pass("uname", getPass::getPass()))
```

1- link raw data included in data-raw, or put it in a directory publicly available where the package will be proposed.

2- hide all functions except:
- Check.Data.Format
- Convert.Raster2BIL
- Perform.Radiometric.Filtering
- Perform.PCA.Images
- Select.Components
- Map.Spectral.Species
- Map.Alpha.Diversity
- Map.Beta.Diversity
- Get.List.Shp
- Get.Projection
- Get.Diversity.From.Plots

3- add manual help based on 'Main_DiversityMapping.R'

4- produce a binary/compiled packages for version beta

5- put the binary package ready to be installed in a directory publicly available where the data is stored as well.

6- replace hard-coded instructions (e.g. Lib_CheckConvertData.R line 110 to get access to the HDR templates)
