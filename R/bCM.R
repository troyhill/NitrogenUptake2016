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
#' @examples ### get allometry model for each species
#' CSP <- plyr::dlply(allometry, c("spp"), bCM)
#' CSP.coef <- plyr::ldply(CSP, stats::coef)
#' ### add lambda value
#' CSP.coef$lam <- plyr::ddply(allometry, c("spp"), function(df)  
#'                 bCM(df, lam.only = TRUE))[, "V1"]
#' @export
bCM <- function(dat, mass = "sample", height = "height_cm", lam.avail = c(-2, -1.5, -1, -2/3, -1/2, -1/3, 0, 1/3, 1/2, 2/3, 1, 1.5, 2),
                lam.only = FALSE) {
  dat.sub <- dat[(!is.na(dat[, height])) & (!is.na(dat[, mass])), ]
  height <- dat.sub[, height]
  mod.h <- stats::lm(dat.sub[, mass] ~ height)
  bc.h  <- MASS::boxcox(dat.sub[, mass] ~ height, lam.avail, interp=F, plotit=F)
  lam.h <- bc.h$x[bc.h$y==max(bc.h$y)]
  if (lam.only == TRUE) {
    lam.h
  } else if (lam.only == FALSE) {
    mod.h.bc <- lm(car::basicPower(dat.sub[, mass] , lam.h) ~ height)
    mod.h.bc
  }
}
