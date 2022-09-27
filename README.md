# __biodivMapR__ <img src="man/figures/logo.png" align="right" alt="" width="200" />

# An R package for α- and β-diversity mapping using remotely-sensed images

[![build](https://img.shields.io/github/workflow/status/jbferet/biodivMapR/tic/master)](https://github.com/jbferet/biodivMapR/actions)
[![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
[![version](https://img.shields.io/github/v/release/jbferet/biodivMapR?label=version)](https://github.com/jbferet/biodivMapR)

# 1 Install

The package `remotes` first needs to be installed from the CRAN

```
install.packages("remotes")
```

Then some packages which were removed from the CRAN may need to be installed directly from authors' repository. This is the case for `dissUtils`: 

```
remotes::install_github('cran/dissUtils')
```

After installing `remotes`and `dissUtils`, `biodivMapR` should be ready for installation with the following command line in your R session:

```
remotes::install_github('jbferet/biodivMapR')
```


# 2 Tutorial

A tutorial vignette is available [here](https://jbferet.github.io/biodivMapR/articles/biodivMapR.html).

The corresponding script is available in file `examples/tutorial.R`.

# 3 Citation

If you use **biodivMapR**, please cite the following references:

Féret, J.-B., de Boissieu, F., 2019. biodivMapR: an R package for α‐ and β‐diversity mapping using remotely‐sensed images. Methods Ecol. Evol. 00:1-7. https://doi.org/10.1111/2041-210X.13310

Féret, J.-B., Asner, G.P., 2014. Mapping tropical forest canopy diversity using high-fidelity imaging spectroscopy. Ecol. Appl. 24, 1289–1296. https://doi.org/10.1890/13-1824.1


# 4 Acknowledgments / Fundings

This research was supported by the Agence Nationale de la Recherche ([ANR](https://anr.fr/en/open-calls-and-preannouncements/), France) through the young researchers project **BioCop** (ANR-17-CE32-0001)
