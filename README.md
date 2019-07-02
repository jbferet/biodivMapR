# DiversityMappR: an R package for α- and β-diversity mapping using remotely-sensed images
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)

# 1 Install
After installing packages `devtools` and `getPass`, package `DiversityMappR` can then be installed with the folloqing command line in R session, where `uname` is your gitlab.irstea.fr username:
```
devtools::install_git('https://gitlab.irstea.fr/jean-baptiste.feret/diversitymappr.git',
                      credentials = git2r::cred_user_pass("uname", getPass::getPass()))
```

# 2 Tutorial
A tutorial vignette showing the main steps of the processing can be visualised with the following command line:
```
rstudioapi::viewer(system.file('doc', 'tutorial.html', package='DiversityMappR'))
```
or for non Rstudio session:
```
vignette('tutorial', package='DiversityMappR')
```

The corresponding script is available in file `examples/tutorial.R`.

