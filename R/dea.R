#' Data: Denitrification enzyme activity and in vitro N2O production rates
#'
#' @format A dataframe with six rows and five columns:
#' \describe{
#' \item{pot}{Mesocosm ID; equivalent to "new.core.id" in other datasets}
#' \item{DEA}{Denitrification enzyme activity (units = nanomoles N2O / gram dry mass / hour)}
#' \item{IV}{"In vitro" N2O production; no nutrient solution added, just filtered seawater (units = nanomoles N2O / gram dry mass / hour)}
#' \item{mcf}{Moisture correction factor (1 - gravimetric water content)}
#' \item{bd_gcm3}{Bulk density (grams per cubic centimeter)}
#'}
#' @docType data
#' @keywords data, denitrification, N2O
#' @name dea
#' @usage dea
#' @examples ### export to .csv:
#' write.csv(dea, file = file.path(tempdir(), "dea.csv"))
NULL
