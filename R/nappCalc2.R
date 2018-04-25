#' Calculate net aboveground primary production
#'
#' @param dataset data
#' @param liveCol live biomass
#' @param deadCol dead biomass
#' @param yearCol year
#' @param siteCol site/plot/experimental unit identifier
#' @param timeCol time column (sequential measurements within each year)
#' @param annualReset should data be reset to zero each year
#' @param MilnerHughes If "TRUE", Milner-Hughes NAPP is calculated
#' @param EOS If "TRUE", end-of-season live biomass is reported
#' @param EOS_window window for EOSL
#' @param summarize If "TRUE", output will be a list with two elements: incremental data and summary data
#'
#' @return list
#' @importFrom zoo as.yearmon
#'
#' @export
nappCalc2 <- function(
  dataset,
  liveCol = "mass",
  deadCol = "dead",
  yearCol = "year",
  siteCol = "pot2",
  timeCol = "day",
  annualReset = "TRUE",
  MilnerHughes = "TRUE",
  EOS = "FALSE",
  EOS_window = 1,
  summarize = "TRUE"
) {
  # requires that zoo library be loaded (time is converted to yearmon for finding EOS)
  # implements Smalley (1958) and Milner and Hughes (1968)
  # runs for entire dataset, reports results by year (unique value in "yearCol") for each site (unique value in "siteCol")
  # Does not make assumptions about returns to zero. If that assumption is desired (e.g., if large gaps exist between sampling intervals),
  # it is recommended to run this function separately for each year
  #   dataset = dataframe with your data
  #   liveCol = name of the column with live biomass data
  #   deadCol = name of the column with dead biomass data
  #   yearCol = name of the column with year (4 digits)
  #   siteCol = name of the column with plot name
  #   timeCol = column with time data (in form "%b %Y"). This only matters for the summary statistics, which report peak timing
  #   MilnerHughes      = if "TRUE", also implements Millner & Hughes 1968 (sum of positive changes in standing *live* biomass)
  #   summarize = if "TRUE", summary statistics (max NAPP estimates) are reported. TODO: report peak timing
  #   EOS    = "TRUE" includes a column calculating September *live* biomass (or closest month in dataset). If there was
  #              no sampling within some number of months (+- EOS_window) of September, value is reported as NA
  #
  # Usage examples:
  # # Single site, single year
  #   test <- napp[(napp$site %in% "LUM1") & (napp$year %in% "2014"), 1:6]
  #   nappCalc(test)
  # # Single site, multiple years
  #   test2 <- napp[(napp$site %in% "LUM1"), 1:6]
  #   nappCalc(test2, summarize = "TRUE")[[2]]
  # # Multiple sites, multiple years
  #   test3 <- napp[, 1:6]
  #   nappCalc(test3)
  #   nappCalc(test3, summarize = "TRUE")[[2]]
  #   nappCalc(test3, summarize = "TRUE")[[2]]
  # # Last two columns are identical to
  #   PSC(napp[, 1:6])

  ### error checking
  countsAsTrue  <- c("T", "TRUE", "true", "True")
  countsAsFalse <- c("F", "FALSE", "false", "False")

  if (!MilnerHughes %in% c(countsAsTrue, countsAsFalse)) {
    stop ("`MilnerHughes` argument isn't recognized. Input can be either `TRUE` or `FALSE`.")
  }
  if (sum(c(liveCol, deadCol, yearCol, siteCol) %in% names(dataset)) < 4) {
    stop ("Check column names. One or more column names were not found in the dataset.")
  }
  ###

  tempData    <- dataset
  # column names as variables:
  smalley.inc <- "smalley.inc"
  smalley     <- "smalley"
  live.inc    <- "live.inc" # VTS1975's delL
  dead.inc    <- "dead.inc" # VTS1975's delD
  MH          <- "MH"
  maxMin      <- "maxMin"
  # Valiela, Teal, Sass 1975
  eV          <- "VTS1975.inc"
  VTS         <- "VTS1975"
  PSC_A       <- "psc.live"
  PSC_B       <- "psc.tot"
  EOS_col     <- "eos"

  # define acceptable window for EOS measurement
  EOS_target <- as.numeric(as.yearmon("Sep 2016", format = "%b %Y")) - 2016
  EOS_high   <- as.numeric(as.yearmon(paste0(month.abb[grep("Sep", month.abb) + EOS_window], "2016"), format = "%b %Y")) - 2016
  EOS_low   <- as.numeric(as.yearmon(paste0(month.abb[grep("Sep", month.abb) - EOS_window], "2016"), format = "%b %Y")) - 2016

  # more variables than necessary are appended to dataset
  tempData[, EOS_col] <- tempData[, PSC_B] <- tempData[, PSC_A] <- tempData[, VTS] <- tempData[, MH] <- tempData[, smalley] <- tempData[, smalley.inc] <-
    tempData[, eV] <- tempData[, dead.inc] <- tempData[, live.inc] <- as.numeric(NA)

  for (h in 1:length(unique(tempData[, siteCol]))) {
    # print(h)
    targetSite <- unique(tempData[, siteCol])[h]
    subData1 <- tempData[tempData[, siteCol] %in% targetSite, ]
    # calculates biomass increments (Bn - B(n-1)) over all available years

    if (nrow(subData1) > 1) {
      for (j in 2:(nrow(subData1))) {
        # live biomass increments
        subData1[, live.inc][j] <- subData1[, liveCol][j] - subData1[, liveCol][j - 1]
        # dead biomass increments
        subData1[, dead.inc][j] <- subData1[, deadCol][j] - subData1[, deadCol][j - 1]

        # calculate e from Valiela, Teal, Sass 1975
        if ((subData1[, live.inc][j] >= 0) & (subData1[, dead.inc][j] < 0)) {
          subData1[, eV][j] <- -subData1[, dead.inc][j]
        } else if (subData1[, live.inc][j] < 0) {
          subData1[, eV][j] <- -(subData1[, dead.inc][j] + subData1[, live.inc][j])
        }
        # change e to zero, if negative
        if (!is.na(subData1[, eV][j])) {
          if (subData1[, eV][j] < 0 ) {
            subData1[, eV][j] <- 0
          }
        }

        # apply decision rules
        # both increments positive: sum of both
        # both negative: zero
        # live + and dead -: use live
        # live - and dead +: difference (if positive, otherwise use zero)
        if ((subData1[, live.inc][j] > 0) & (subData1[, dead.inc][j] > 0)) {
          newVal <- subData1[, live.inc][j] + subData1[, dead.inc][j]
        } else if ((subData1[, live.inc][j] <= 0) & (subData1[, dead.inc][j] <= 0)) {
          newVal <- 0
        } else if ((subData1[, live.inc][j] > 0) & (subData1[, dead.inc][j] <= 0)) {
          newVal <- subData1[, live.inc][j]
        } else if ((subData1[, live.inc][j] <= 0) & (subData1[, dead.inc][j] > 0)) {
          calc <- subData1[, live.inc][j] + subData1[, dead.inc][j]
          if (calc < 0) {
            newVal <- 0
          } else {
            newVal <- calc
          }
        }
        # write values to output dataframe
        tempData[tempData[, siteCol] %in% targetSite, live.inc][j]    <- subData1[, live.inc][j]
        tempData[tempData[, siteCol] %in% targetSite, dead.inc][j]    <- subData1[, dead.inc][j]
        tempData[tempData[, siteCol] %in% targetSite, eV][j]          <- subData1[, eV][j]
        tempData[tempData[, siteCol] %in% targetSite, smalley.inc][j] <- newVal
      }
    }
    # now, calculate NAPP cumulatively for each year
    # treat NAs as zeroes in NAPP calculation
    for (k in 1:length(unique(subData1[, yearCol]))) {
      targetYear <- unique(subData1[, yearCol])[k]
      # data for a single site, single year
      subData2   <- tempData[(tempData[, yearCol] %in% targetYear) & (tempData[, siteCol] %in% targetSite), ]

      # if annualReset is TRUE, set the year's first increment to zero/NA
      # This doesn't assume NAPP goes to zero, but deals with each year separately
      if (annualReset %in% countsAsTrue) {
        subData2[1, c(live.inc, dead.inc, smalley.inc, eV)]   <- 0
      }

      # sum smalley increments
      subData2[, smalley][!is.na(subData2[, smalley.inc])] <- cumsum(subData2[, smalley.inc][!is.na(subData2[, smalley.inc])])

      # sum positive *live* biomass increments for Millner & Hughes 1968
      subData2[, "MH.inc"] <- subData2[, live.inc]
      # use only positive increments to calc MH NAPP
      subData2[, "MH.inc"][subData2[, "MH.inc"] < 0] <- 0
      subData2[, MH][!is.na(subData2[, "MH.inc"])] <- cumsum(subData2[, "MH.inc"][!is.na(subData2[, "MH.inc"])])

      # sum e from Valiela, Teal, Sass 1975
      subData2[, VTS][!is.na(subData2[, eV])] <- cumsum(subData2[, eV][!is.na(subData2[, eV])])

      # sum peak standing crop increments
      subData2[, PSC_A][!is.na(subData2[, liveCol])]                       <- cummax(subData2[, liveCol][!is.na(subData2[, liveCol])])
      # PSC_B reflects maximum summed (live + dead) biomass when live biomass is at its peak.
      subData2[, PSC_B][!is.na(subData2[, liveCol] + subData2[, deadCol])] <- max(subData2[, liveCol][!is.na(subData2[, liveCol])]) + subData2[, deadCol][subData2[, liveCol][!is.na(subData2[, liveCol])] == max(subData2[, liveCol][!is.na(subData2[, liveCol])])]

      # find "end of season" biomass
      # may need to ensure that monthly data exist? could pose a problem when assigning EOS_biomass and querying months.
      if(EOS %in% countsAsTrue) {
        # if no sampling occurred in September, widen window
        if (round(EOS_target, 2) %in% round((as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol])))), 2)) {
          EOS_biomass <- subData2[, liveCol][(as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol]))) == EOS_target)] # +
          #subData2[, deadCol][(as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol]))) == EOS_target)]
          # } else if(sum(c(round(seq(from = EOS_low, to = EOS_high, by = 1/12), 2) %in% round((as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol])))), 2))) >= 1) {
          # EOS_biomass <- max(subData2[, liveCol][(as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol]))) >= EOS_low) | (as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol]))) <= EOS_high)], na.rm = TRUE)
        } else {
          EOS_biomass <- NA
        }
        subData2[, EOS_col][subData2[, liveCol]  == EOS_biomass] <- EOS_biomass
        #subData2[, EOS_col][subData2[, liveCol] + subData2[, deadCol] == EOS_biomass] <- EOS_biomass
      }


      # add to output dataframe
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), smalley]  <- subData2[, smalley]
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), MH]       <- subData2[, MH]
      # MH should resolve to max-min if first and final biomass data = 0
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), maxMin]   <- max(subData2[, liveCol], na.rm = T) - min(subData2[, liveCol], na.rm = T)
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), VTS]      <- subData2[, VTS]
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), PSC_A]    <- subData2[, PSC_A]
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), PSC_B]    <- subData2[, PSC_B]
      if(EOS %in% countsAsTrue) {
        tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), EOS_col]   <- subData2[, EOS_col]
      }
    }
  }

  if (summarize %in% countsAsFalse) {
    output <- tempData
  } else if (summarize %in% countsAsTrue) {
    for (i in 1:length(unique(tempData[, siteCol]))) {
      # print(i)
      targetSite <- unique(tempData[, siteCol])[i]
      subData1 <- tempData[tempData[, siteCol] %in% targetSite, ]
      for (j in 1:length(unique(subData1[, yearCol]))) {
        targetYear <- unique(subData1[, yearCol])[j]
        subData2   <- subData1[subData1[, yearCol] %in% targetYear, ]

        intData <- data.frame(
          site = targetSite,
          year = targetYear,
          mean.live    = mean(subData2[, liveCol], na.rm = T),
          mean.dead    = mean(subData2[, deadCol], na.rm = T),
          napp.smalley = max(subData2[, smalley], na.rm = T),
          napp.MH      = max(subData2[, MH], na.rm = T),
          napp.maxMin  = max(subData2[, maxMin], na.rm = T),
          # napp.VTS     = max(subData2[, VTS], na.rm = T),
          napp.psc.a   = max(subData2[, PSC_A], na.rm = T),
          napp.psc.b   = max(subData2[, PSC_B], na.rm = T),
          n            = sum(!is.na(subData2[, liveCol])),
          t.smalley    = ifelse(is.finite(max(subData2[, smalley], na.rm = T)), as.character(subData2[, timeCol][which.max(subData2[, smalley])]), NA),
          t.MH         = ifelse(is.finite(max(subData2[, MH], na.rm = T)), as.character(subData2[, timeCol][which.max(subData2[, MH])]), NA),
          # t.vts        = as.character(subData2[, timeCol][which.max(subData2[, VTS])]),
          t.psc.a      = as.character(subData2[, timeCol][which.max(subData2[, PSC_A])]),
          t.psc.b      = as.character(subData2[, timeCol][which.max(subData2[, PSC_B])])
        )

        # ADD EOS to intdata if EOS == TRUE
        if (EOS %in% countsAsTrue) {
          intData$napp.EOS <- ifelse(sum(is.na(subData2[, EOS_col])) == length(subData2[, EOS_col]),
                                     NA, max(subData2[, EOS_col], na.rm = T))
        }


        if (i != 1 ) {
          finalData <- rbind(finalData, intData)
        } else if ((i == 1) & (j != 1)) {
          finalData <- rbind(finalData, intData)
        } else {
          finalData <- intData
        }

      }
    }
    finalData$napp.MH[!is.finite(finalData$napp.MH)] <- NA
    output <- list(intervalData = tempData, summary = finalData)
  }
  output
}


