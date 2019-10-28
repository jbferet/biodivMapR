# biodivMapR: an R package for α- and β-diversity mapping using remotely-sensed images
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)
[![Build Status](https://travis-ci.org/jbferet/biodivMapR.png?branch=master)](https://travis-ci.org/jbferet/biodivMapR)

# 1 Install
After installing packages `devtools` and `getPass`, package `biodivMapR` can then be installed with the folloqing command line in R session, where `uname` is your gitlab.irstea.fr username:
```
devtools::install_git('https://github.com/jbferet/biodivMapR.git')
```

# 2 Tutorial
A tutorial vignette showing the main steps of the processing can be visualised with the following command line:
```
rstudioapi::viewer(system.file('doc', 'tutorial.html', package='biodivMapR'))
```
or for non Rstudio session:
```
vignette('tutorial', package='biodivMapR')
```

The corresponding script is available in file `examples/tutorial.R`.

