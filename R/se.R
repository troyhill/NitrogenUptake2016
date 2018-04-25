#' Calculates standard error
#'
#' @param x numeric or integer
#'
#' @return value
#'
#' @examples
#' se(CN_mass_data$n_pct)
#' plyr::ddply(CN_mass_data, plyr::.(species, pool_label), plyr::summarise, se = se(n_pct))
#' @importFrom stats sd
#' @export
se <- function(x){
  sd(x, na.rm = T) / sqrt(sum(!is.na(x) == T))
}
