# DiversityMappR: an R package for α- and β-diversity mapping using remotely-sensed images

# 1 Install
After installing packages devtools and getPass, `DiversityMappR` can then be installed with the folloqing command line in R session, where `uname` is your gitlab.irstea.fr username:
```
devtools::install_git('https://gitlab.irstea.fr/jean-baptiste.feret/diversitymappr.git',
                      credentials = git2r::cred_user_pass("uname", getPass::getPass()))
```

