#' Data: Stem masses and heights for plants collected from Colt State Park, Rhode Island, USA, during summer 2016
#'
#' A dataframe of masses and heights of stems of Spartina alterniflora and Distichlis spicata. Samples were collected from Colt State Park, Bristol, RI, USA, during May-July 2016. Column descriptions:
#'
#'
#' @format A dataframe with 170 observations of 6 variables:
#' \describe{
#' \item{site}{Study location (Colt State Park, RI)}
#'
#' \item{samplingDate}{Sampling dates}
#'
#' \item{status}{Indicates whether plant was live or dead}
#'
#' \item{height_cm}{Stem height, in centimeters, from the sediment surface to the tip of the longest leaf}
#'
#' \item{sample}{Biomass, grams}
#'
#' \item{spp}{Species (SPAL or DISP)}
#' }
#' @docType data
#' @keywords data, allometry
#' @name allometry
#' @references{
#' Hill, T.D., N.R. Sommer, C.R. Kanaskie, E.A. Santos, A.J. Oczkowski. 2018. Data and source code from: Nitrogen uptake and allocation estimates for Spartina alterniflora and Distichlis spicata. Data In Brief. 21: 466-472.
#'   \url{https://doi.org/10.1016/j.dib.2018.09.133}.
#' 
#' Hill, T.D., N.R. Sommer, C.R. Kanaskie, E.A. Santos, A.J. Oczkowski. 2018. Nitrogen uptake and allocation estimates for Spartina alterniflora and Distichlis spicata. Journal of Experimental Marine Biology and Ecology 507: 53-60.
#'   \url{https://doi.org/10.1016/j.jembe.2018.07.006}.
#' }
#' @usage allometry
#' @examples ### export to .csv:
#' write.csv(allometry, file = file.path(tempdir(), "allometry.csv"))
NULL
