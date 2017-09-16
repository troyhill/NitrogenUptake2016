#' Calculates standard error
#'
#' @param x numeric or integer
#'
#' @return value
#'
#' @examples \dontrun{}
#' @importFrom stats sd
#' @export
se <- function(x){
  sd(x, na.rm = T) / sqrt(sum(!is.na(x) == T))
}
