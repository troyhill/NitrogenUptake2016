#' Parameterize Box-Cox model for mass-height allometry (based on Lu et al. 2016)
#'
#' @param dat dataframe with data
#' @param mass mass column
#' @param height height column
#' @param lam.avail set of possible lambda values
#' @param lam.only if TRUE, lambda is returned. If FALSE, model is returned
#'
#' @return if lam.only is FALSE, a model is returned. If lam.only is TRUE, lambda value is returned.
#' @importFrom MASS boxcox
#' @importFrom stats lm
#' @importFrom car basicPower
#' @examples \dontrun{}
#' @export
bCM <- function(dat, mass = "sample", height = "height_cm", lam.avail = c(-2, -1.5, -1, -2/3, -1/2, -1/3, 0, 1/3, 1/2, 2/3, 1, 1.5, 2),
                lam.only = FALSE) {
  mod.h <- lm(get(mass) ~ get(height), data = dat[!is.na(dat[, height]),])
  bc.h  <- boxcox(get(mass) ~ get(height), data = dat[!is.na(dat[, height]),], lam.avail, interp=F, plotit=F)
  lam.h <- bc.h$x[bc.h$y==max(bc.h$y)]
  if (lam.only == TRUE) {
    lam.h
  } else if (lam.only == FALSE) {
    mod.h.bc <- lm(basicPower(get(mass) , lam.h) ~ get(height), data = dat[!is.na(dat[, height]),])
    mod.h.bc
  }
}
