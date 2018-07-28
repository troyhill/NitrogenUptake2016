#' Data: Stem heights for each mesocosm and each measurement date
#'
#' A dataframe of Spartina alterniflora and Distichlis spicata stem heights, from the mesocosms used in 15N study. Samples were collected from Colt State Park, Bristol, RI, USA, and grown in the US EPA Atlantic Ecology Division greenhouse.
#'

#' @format A dataframe with 3315 observations of 10 variables:
#'  \describe{
#' \item{date}{Measurement date}
#'
#' \item{core_num}{Mesocosm number}
#'
#' \item{species}{Species, either Spartina alterniflora (SA) or Distichlis spicata (DS)}
#'
#' \item{dead_live}{Indicates whether plant was live or dead}
#'
#' \item{plant_num}{Plants were tagged to permit growth rate calculations; this is the plant tag number}
#'
#' \item{height_cm}{Stem height, in centimeters, from the sediment surface to the tip of the longest leaf}
#'
#' \item{id}{Unique plant identifier, combining species, mesocosm number, and plant tag number}
#'
#' \item{day}{Measurement date expressed in YYYY-MM-DD format (and structured as a POSIXct column in R)}
#'
#' \item{timeSinceLast}{Days since last measurement}
#'
#' \item{new.core.id}{Mesocosm ID, including species and a unique mesocosm number (time-zero mesocosms re-numbered as mesocosms 13, 14, and 15)}
#'
#' }
#' @docType data
#' @keywords data, aboveground biomass, NAPP, stem height, stem density
#' @name stemHeights
#' @usage stemHeights
#' @examples ### export to .csv:
#' write.csv(stemHeights, file = file.path(tempdir(), "stemHeights.csv"))
NULL
