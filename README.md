# EpiGraphDB pQTL MR and colocalization R script
This is the github repository to provide the R scripts been used in the phenome-wide Mendelian randomization study of plasma proteome. 

<a href="http://epigraphdb.org"><img src="man/figures/logo_wide.png" alt="" height="60" style="padding:10px"/></a> <span class="pull-right"> <a href="http://www.bris.ac.uk"><img src="man/figures/ieu40.png" alt="" height="60" style="padding:10px"/></a> <a href="http://www.bris.ac.uk/ieu"><img src="man/figures/uob40.png" alt="" height="60" style="padding:10px"/></a> </span>

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/epigraphdb)](https://cran.r-project.org/package=epigraphdb)
[![Travis build status](https://travis-ci.org/MRCIEU/epigraphdb-r.svg?branch=master)](https://travis-ci.org/MRCIEU/epigraphdb-pqtl)
[![Codecov test coverage](https://codecov.io/gh/MRCIEU/epigraphdb-r/branch/master/graph/badge.svg)](https://codecov.io/gh/MRCIEU/epigraphdb-r?branch=master)

<!-- badges: end -->

[EpiGraphDB-pQTL](https://epigraphdb.org/pqtl/) is an web browser and database to show the MR and colocalization results of the pQTL MR project of EpiGraphDB. 
The browser currently contains the Mendelian Randomization and sensitivity analyses results for 989 proteins and 225 traits, i.e. diseases and risk factors. To start using this browser, simply type a protein or trait name into the "search" field, for example, ADAM19 or Lung cancer. The full list of proteins can be found by following the link on top of the "search" field.
[`epigraphdb-pqtl`](https://github.com/MRCIEU/epigraphdb-pqtl/) is an github repo to provide ease to use scripts to run MR and colocalization analysis for the pQTL MR project. 

## Installation

[`devtools`](https://devtools.r-lib.org/)
is required to install from github:

```r
###install the Two sample MR package (just need do this once) 
##source("https://bioconductor.org/biocLite.R")
#install.packages("devtools")

##to update the R package (once there is a )
#library(devtools)
#install_github("MRCIEU/TwoSampleMR")

#example of use the older version of the package
#devtools::install_github("MRCIEU/TwoSampleMR@0.3.2")
```

## Run pQTL MR analysis

`epigraphdb` provides a simple and intuitive way to run the pQTL MR analysis using "TwoSampleMR" R package

For more information on how to run the MR analysis, please check out R script `MR-pQTL-script.R`

## EpiGraphDB pQTL browser 

| link                                                | 
|-----------------------------------------------------|
| [website](http://epigraphdb.org/pqtl/)              |


## Citation

Please cite the pQTL MR analysis as

> Jie Zheng, Valeriia Haberland, Denis Baird, Venexia M Walker, Philip M Haycock, Alex X Gutteridge, Tom M Richardson, James Staley, Benjamin Elsworth, Stephen Burgess, Benjamin B Sun, John Danesh, Heiko Runz, Joseph C Maranville, Hannah Martin, James Yarmolinsky, Charles Laurin, Michael V Holmes, Jimmy Liu, Karol Estrada, Linda C McCarthy, Mark Hurle, Dawn M Waterworth, Matthew R Nelson, Adam S Butterworth, George Davey Smith, Gibran V Hemani, Robert A Scott, Tom R Gaunt. Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases.BioRXiv (2019)

```
@article {pQTL MR paper,
  author = {Jie Zheng, Valeriia Haberland, Denis Baird, Venexia M Walker, Philip M Haycock, Alex X Gutteridge, Tom M Richardson, James Staley, Benjamin Elsworth, Stephen Burgess, Benjamin B Sun, John Danesh, Heiko Runz, Joseph C Maranville, Hannah Martin, James Yarmolinsky, Charles Laurin, Michael V Holmes, Jimmy Liu, Karol Estrada, Linda C McCarthy, Mark Hurle, Dawn M Waterworth, Matthew R Nelson, Adam S Butterworth, George Davey Smith, Gibran V Hemani, Robert A Scott, Tom R Gaunt},
  title = {Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases},
  url = {https://www.biorxiv.org/content/10.1101/627398v2}
}
```

## Contact

Please get in touch with us for issues, comments, suggestions, etc. via the following methods:

- [The issue tracker on the `epigrapdb-pqtl` repo](https://github.com/MRCIEU/epigraphdb/issues)
- [The support email](mailto:feedback@epigraphdb.org)
- [The EpiGraphDB twitter](https://twitter.com/epigraphdb)