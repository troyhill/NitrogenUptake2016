#' Identify lambda values from Box Cox model (based on Lu et al. 2016)
#'
#' @param df data frame
#' @param mass mass column
#' @param height height column
#' @param lam.avail possible lambda values that can be used in Box-Cox transformations
#'
#' @return lambda value
#'
#' @examples \dontrun{}
#' @importFrom MASS boxcox
#' @export
lam <- function(df, mass = "sample", height = "height_cm", lam.avail = c(-2, -1.5, -1, -2/3, -1/2, -1/3, 0, 1/3, 1/2, 2/3, 1, 1.5, 2) ) {
  dat.h <- df[!is.na(df[, height]), ]
  bc.h <- boxcox(get(mass) ~ get(height), data=dat.h, lam.avail, interp=F, plotit=F)
  lam.h <- bc.h$x[bc.h$y==max(bc.h$y)]
  lam.h
}
