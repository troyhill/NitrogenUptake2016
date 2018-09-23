# NitrogenUptake2016


[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/NitrogenUptake2016)](https://cran.r-project.org/package=NitrogenUptake2016) [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/NitrogenUptake2016)](https://cran.r-project.org/package=NitrogenUptake2016) [![](http://cranlogs.r-pkg.org/badges/grand-total/NitrogenUptake2016)](https://cran.r-project.org/package=NitrogenUptake2016)



This repository contains data and source code associated with the manuscript "Nitrogen uptake and allocation estimates for _Spartina alterniflora_ and _Distichlis spicata_" (Hill et al. 2018; DOI: https://doi.org/10.1016/j.jembe.2018.07.006).


## Install the package

Install our R package directly from CRAN using the following commands in the R console:

```r
install.packages("NitrogenUptake2016")
library(NitrogenUptake2016)
```


## Read our manuscripts

After loading the NitrogenUptake2016 package, the manuscript vignettes are available as pdfs:

```r
# Journal of Experimental Marine Biology and Ecology (Hill et al. 2018)
vignette("JEMBE", package = "NitrogenUptake2016")

# Data In Brief (Hill et al. submitted)
vignette("DataInBrief", package = "NitrogenUptake2016")
```



## Access our data

All the data used in our manuscripts are also available and annotated:

```r
?allometry
?dea
?stemHeights
?CN_mass_data
```


## Disclaimer 

This project code is provided on an "as is" basis and the user assumes responsibility for its use. The authors have relinquished control of the information and no longer have responsibility to protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.


## Thanks 

This effort to support open science benefitted significantly from the example set by Jeff Hollister et al.'s [LakeTrophicModelling](https://github.com/USEPA/LakeTrophicModelling) manuscript.


