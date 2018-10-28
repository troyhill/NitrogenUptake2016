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
#' @references{
#' Hill, T.D., N.R. Sommer, C.R. Kanaskie, E.A. Santos, A.J. Oczkowski. 2018. Data and source code from: Nitrogen uptake and allocation estimates for Spartina alterniflora and Distichlis spicata. Data In Brief. 21: 466-472.
#'   \url{https://doi.org/10.1016/j.dib.2018.09.133}.
#' 
#' Hill, T.D., N.R. Sommer, C.R. Kanaskie, E.A. Santos, A.J. Oczkowski. 2018. Nitrogen uptake and allocation estimates for Spartina alterniflora and Distichlis spicata. Journal of Experimental Marine Biology and Ecology 507: 53-60.
#'   \url{https://doi.org/10.1016/j.jembe.2018.07.006}.
#' }
#' @examples ### export to .csv:
#' write.csv(dea, file = file.path(tempdir(), "dea.csv"))
NULL
