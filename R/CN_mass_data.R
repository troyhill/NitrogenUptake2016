#' Data: Nutrient concentrations, stable isotope ratios, and biomass from destructive mesocosm harvests
#'
#'
#'
#' @format A dataframe with 1192 observations of 16 variables:
#' \describe{
#' \item{time}{Time point of harvest (harvested at one-week intervals)}
#' \item{new.core.id}{Unique mesocosm identifier, including species (SA or DS) and mesocosm number}
#' \item{depth_bottom}{Depth at bottom of sample (only applicable for belowground data)}
#' \item{sample.type}{Sample material; tissue type}
#' \item{interval}{Depth interval for sample; indicates the top and bottom depths (e.g., an entry of "5_10" covers the depth interval from 5-10 cm)}
#' \item{pool_label}{Label for each pool (combination of "sample.type" and "depth_bottom")}
#' \item{id}{Same as "pool_label" but with mesocosm ID included}
#' \item{species}{Spartina alterniflora (SA) or Distichlis spicata (DS)}
#' \item{d15n}{15-N isotope ratio in per mille units}
#' \item{n_pct}{Nitrogen content, decimal fraction (0.015 = 1.5 percent)}
#' \item{d13c}{13-C isotope ratio in per mille units}
#' \item{c_pct}{Carbon content, decimal fraction (0.015 = 1.5 percent)}
#' \item{total_volume_cm3}{Total volume of interval (only applicable for belowground data)}
#' \item{depth_top}{Depth at top of sample (only applicable for belowground data)}
#' \item{sample.type2}{Simplified "sample.type" column; leaf numbers dropped, and belowground stems included as "stems"}
#' \item{g_core}{Total mass (grams) in entire pool; corrects for subsampling of depth intervals}
#'}
#' @docType data
#' @keywords data, stable isotopes, 15N
#' @name CN_mass_data
#' @usage CN_mass_data
#' @references{
#' Hill, T.D., N.R. Sommer, C.R. Kanaskie, E.A. Santos, A.J. Oczkowski. 2018. Data and source code from: Nitrogen uptake and allocation estimates for Spartina alterniflora and Distichlis spicata. Data In Brief. 21: 466-472.
#'   \url{https://doi.org/10.1016/j.dib.2018.09.133}.
#' 
#' Hill, T.D., N.R. Sommer, C.R. Kanaskie, E.A. Santos, A.J. Oczkowski. 2018. Nitrogen uptake and allocation estimates for Spartina alterniflora and Distichlis spicata. Journal of Experimental Marine Biology and Ecology 507: 53-60.
#'   \url{https://doi.org/10.1016/j.jembe.2018.07.006}.
#' }
#' @examples ### export to .csv:
#' write.csv(CN_mass_data, file = file.path(tempdir(), "CN_mass_data.csv"))
NULL
