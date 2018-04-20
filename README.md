# NitrogenUptake2016

This repository contains data and source code associated with Hill et al. "Nitrogen uptake and allocation estimates for _Spartina alterniflora_ and _Distichlis spicata_"

[![DOI](https://zenodo.org/badge/103706037.svg)](https://zenodo.org/badge/latestdoi/103706037)


## Install the Package 

Install the package from GitHub using devtools:

```r
install.packages("devtools")
devtools::install_github("troyhill/NitrogenUptake2016", build_vignettes = TRUE)
library("NitrogenUptake2016")
```

All the data used in this manuscript will then be available using data():

```r
data("allometry")
data("dea")
data("stemHeights")
data("CN_mass_data")
```


The manuscript vignette can be read via:

```r
vignette("analysis",package="NitrogenUptake2016")
```


## Disclaimer 

This project code is provided on an "as is" basis and the user assumes responsibility for its use. The authors have relinquished control of the information and no longer have responsibility to protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.


## Thanks 

This effort to support open science benefitted significantly from the example set by Jeff Hollister et al.'s [LakeTrophicModelling](https://github.com/USEPA/LakeTrophicModelling) manuscript.
