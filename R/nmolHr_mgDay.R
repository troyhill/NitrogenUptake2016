#' Convert N2O units from nanomoles of N2O per hour to milligrams of N per day
#'
#' @param x numeric or integer value(s)
#'
#' @return numeric value
#'
#' @examples nmolHr_mgDay(dea$DEA)
#' @export
nmolHr_mgDay <- function(x) { # x = nanomoles of N20 / hr. output: mg N per day
  x / 1e6 * 14.0067 * 24 * 2
} # nmolHr_ugDay(dd.dea$flux.m2) * 365 / 1e6 * 10^4 # kg N / ha / yr (for comparison with Wigand data)
