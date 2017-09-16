#' Convert per mil isotope values to atom percent
#'
#' @param perMilValues Value to be converted, in per mil notation
#' @param isotope can be 13C or 15N
#'
#' @return a numeric value or vector
#' @examples ap(10); ap(1000)
#' @export
ap <- function(perMilValues, isotope = "15N") {
  # function converts per mil enrichment values to atom percent
  if (isotope %in% "15N") {
    R_val <- 0.0036765
  } else if (isotope %in% "13C") {
    R_val <- 0.01118
  }
  (perMilValues + 1000) / ((perMilValues + 1000 + (1000 / R_val))) * 100
}
