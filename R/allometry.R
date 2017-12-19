#' Data: Stem masses and heights for plants collected from Colt State Park, Rhode Island, USA, during summer 2016
#'
#' A dataframe of masses and heights of stems of Spartina alterniflora and Distichlis spicata. Samples were collected from Colt State Park, Bristol, RI, USA, during May-July 2016. Column descriptions:
#'
#'
#' @format A dataframe with 170 observations of 13 variables:
#' \describe{
#' \item{site}{Study location (Colt State Park, RI)}
#'
#' \item{samplingDate}{Sampling dates}
#'
#' \item{plot}{Sampling plot (three plots per vegetation type)}
#'
#' \item{status}{Indicates whether plant was live or dead}
#'
#' \item{height_cm}{Stem height, in centimeters, from the sediment surface to the tip of the longest leaf}
#'
#' \item{height2_cm}{Stem height, in centimeters, from the sediment surface to the farthest-away leaf node}
#'
#' \item{tin_no}{Weigh boat label (administrative use)}
#'
#' \item{tin_mass}{Weigh boat mass, in grams}
#'
#' \item{tin_plant_50C}{Combined tin and plant mass after drying to 50C}
#'
#' \item{shared_base}{Does plant share a base with another stem?}
#'
#' \item{note}{Notes about plant}
#'
#' \item{sample}{Sample mass, grams (no weigh boat)}
#'
#' \item{spp}{Species (SPAL or DISP)}
#' }
#' @docType data
#' @keywords data, allometry
#' @name allometry
#' @usage allometry
#' @examples ### export to .csv:
#' write.csv(allometry, file = "allometry.csv")
#' @references Hill et al.
#' (\href{url}{text})
NULL
