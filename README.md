# NitrogenUptake2016

This repository contains data and source code associated with the manuscript "Nitrogen uptake and allocation estimates for _Spartina alterniflora_ and _Distichlis spicata_" (Hill et al. 2018; DOI: https://doi.org/10.1016/j.jembe.2018.07.006).

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/NitrogenUptake2016)](https://cran.r-project.org/package=NitrogenUptake2016)


[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/NitrogenUptake2016)](https://cran.r-project.org/package=NitrogenUptake2016)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1226378.svg)](https://doi.org/10.5281/zenodo.1226378)


## Install the Package 

Install the package from CRAN:

```r
install.packages("NitrogenUptake2016")
```


Alternatively, the package can be installed from GitHub using devtools:

```r
install.packages("devtools")
devtools::install_github("troyhill/NitrogenUptake2016", build_vignettes = TRUE)
library("NitrogenUptake2016")
```

All the data used in this manuscript will then be available:

```r
?allometry
?dea
?stemHeights
?CN_mass_data
```


The manuscript vignettes can be read via:

```r
# Journal of Experimental Marine Biology and Ecology (in press)
vignette("JEMBE", package = "NitrogenUptake2016")

# Data In Brief (submitted)
vignette("DataInBrief", package = "NitrogenUptake2016")
```


## Disclaimer 

This project code is provided on an "as is" basis and the user assumes responsibility for its use. The authors have relinquished control of the information and no longer have responsibility to protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.


## Thanks 

This effort to support open science benefitted significantly from the example set by Jeff Hollister et al.'s [LakeTrophicModelling](https://github.com/USEPA/LakeTrophicModelling) manuscript.
