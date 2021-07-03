# EpiGraphDB pQTL MR project repository 
This is the GitHub repository to provide all information, news and updates of the EpiGraphDB pQTL MR project, which is a phenome-wide Mendelian randomization study of the plasma proteome. 

<a href="http://epigraphdb.org/pqtl/"><img src="https://epigraphdb.org/img/epigraphdb-logo.ed38e02a.svg" alt="" height="60" style="padding:10px"/></a> <span class="pull-right">&nbsp;&nbsp;&nbsp; <a href="http://www.bris.ac.uk"><img src="https://epigraphdb.org/img/uob.16744ca9.svg" alt="" height="60" style="padding:10px"/></a>&nbsp;&nbsp;&nbsp; <a href="http://www.bris.ac.uk/ieu"><img src="http://www.bristol.ac.uk/media-library/sites/integrative-epidemiology/images/mrc-ieu-logo.png" alt="" height="60" style="padding:10px"/></a> </span>

<!-- badges: start -->

<!--
[![CRAN status](https://www.r-pkg.org/badges/version/epigraphdb)](https://cran.r-project.org/package=epigraphdb)
[![Travis build status](https://travis-ci.org/MRCIEU/epigraphdb-r.svg?branch=master)](https://travis-ci.org/MRCIEU/epigraphdb-pqtl)
-->

<!-- badges: end -->

[`epigraphdb-pqtl`](https://github.com/MRCIEU/epigraphdb-pqtl/) is a GitHub repo to provide the following information: 
1. news and updates of the EpiGraphDB pQTL project, and issue reporting functionality for users. 
2. easy to use R scripts to run MR and colocalization analysis for the pQTL MR project. 

The MR and colocalization analyses results from the EpiGraphDB pQTL project have been pre-calculated and stored in the [EpiGraphDB Proteome PheWAS browser](https://epigraphdb.org/pqtl/). 

## Installation of related R packages

The scripts in this repository require the following dependencies to be installed:

[`devtools`](https://devtools.r-lib.org/)
is required to install from github:

```r
###install the Two sample MR package (just need do this once) 
source("https://bioconductor.org/biocLite.R")
install.packages("devtools")

##to install/update the R package (once there is a )
library(devtools)
install_github("MRCIEU/TwoSampleMR")

#example of use the older version of the package
devtools::install_github("MRCIEU/TwoSampleMR@0.3.2")
```

## Run pQTL MR analysis

`epigraphdb-pqtl` provides a simple and intuitive way to run pQTL MR analysis using the "TwoSampleMR" R package

For more information on how to run the MR analysis, please check out R script `MR-pQTL-script.R`

## EpiGraphDB pQTL PheWAS browser 

The [EpiGraphDB Proteome PheWAS browser](https://epigraphdb.org/pqtl/) currently contains the Mendelian randomization and sensitivity analyses results for 989 proteins and 225 traits, i.e. diseases and risk factors. To start using this browser, simply type a protein or trait name into the "search" field, for example, [ADAM19](https://epigraphdb.org/pqtl/ADAM19) or [Lung cancer](https://epigraphdb.org/pqtl/Lung%20cancer). The full list of proteins can be found by following the link on top of the "search" field.

## Citation

Please cite the pQTL MR analysis as

> Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases

> Jie Zheng, Valeriia Haberland, Denis Baird, Venexia Walker, Philip C. Haycock, Mark R. Hurle, Alex Gutteridge, Pau Erola, Yi Liu, Shan Luo, Jamie Robinson, Tom G. Richardson, James R. Staley, Benjamin Elsworth, Stephen Burgess, Benjamin B. Sun, John Danesh, Heiko Runz, Joseph C. Maranville, Hannah M. Martin, James Yarmolinsky, Charles Laurin, Michael V. Holmes, Jimmy Z. Liu, Karol Estrada, Rita Santos, Linda McCarthy, Dawn Waterworth, Matthew R. Nelson, George Davey Smith, Adam S. Butterworth, Gibran Hemani, Robert A. Scott and Tom R. Gaunt Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases. Nature Genetics. 52, pages1122–1131 (2020)

```
@article {pQTL MR paper,
  author = {Jie Zheng, Valeriia Haberland, Denis Baird, Venexia M Walker, Philip M Haycock, Alex X Gutteridge, Tom M Richardson, James Staley, Benjamin Elsworth, Stephen Burgess, Benjamin B Sun, John Danesh, Heiko Runz, Joseph C Maranville, Hannah Martin, James Yarmolinsky, Charles Laurin, Michael V Holmes, Jimmy Liu, Karol Estrada, Linda C McCarthy, Mark Hurle, Dawn M Waterworth, Matthew R Nelson, Adam S Butterworth, George Davey Smith, Gibran V Hemani, Robert A Scott, Tom R Gaunt},
  title = {Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases},
  url = {https://www.nature.com/articles/s41588-020-0682-6}
}
```

## Contact

Please get in touch with us for issues, comments, suggestions, etc. via the following methods:

- [The issue tracker on the `epigrapdb-pqtl` repo](https://github.com/MRCIEU/epigraphdb-pqtl/issues)
- [The support email](mailto:feedback@epigraphdb.org)
- [The EpiGraphDB twitter](https://twitter.com/epigraphdb)
