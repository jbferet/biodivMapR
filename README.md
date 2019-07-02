# DiversityMappR: an R package for α- and β-diversity mapping using remotely-sensed images
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)

# 1 Install
After installing packages `devtools` and `getPass`, package `DiversityMappR` can then be installed with the folloqing command line in R session, where `uname` is your gitlab.irstea.fr username:
```
devtools::install_git('https://gitlab.irstea.fr/jean-baptiste.feret/diversitymappr.git',
                      credentials = git2r::cred_user_pass("uname", getPass::getPass()))
```

# 2 Tutorial
A tutorial script can be found in `Main_DiversityMapping.R`, showing the main steps of the processing. A tutorial vignette should soon be available.

