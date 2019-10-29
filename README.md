# biodivMapR: an R package for α- and β-diversity mapping using remotely-sensed images <img src="man/figures/logo.png" align="right" alt="" width="120" />

**[https://jbferet.github.io/biodivMapR](https://jbferet.github.io/biodivMapR/index.html)**

![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)
[![Build Status](https://travis-ci.org/jbferet/biodivMapR.png?branch=master)](https://travis-ci.org/jbferet/biodivMapR)

# 1 Install
After installing packages `devtools`, package `biodivMapR` can then be installed with the following command line in R session:
```
devtools::install_github('jbferet/biodivMapR')
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

