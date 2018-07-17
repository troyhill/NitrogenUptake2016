### Script for analysis of GRO 2016 data published as R package

library(ggplot2)
library(gridExtra)
library(rsq)
library(reshape2)
library(plyr)
library(car) # has qqPlot() and basePower() functions
library(MuMIn)
library(MASS) # has boxcox() function
library(zoo)

library(NitrogenUptake2016)


todaysDate <- substr(as.character(Sys.time()), 1, 10)
core.area <- pot.m2 <- 0.00801185 # mesocosm surface area (m2)
top.vol   <- core.area * 0.05 * 1e6 # cm3 in core: top 5 cm only
pointSize <- 1.5 # for ggplot graphics
grayColor <- 0.55
fig2Col   <- "gray55"


# Allometry ----------------------------------------------------------------

##### Establish relationships following Lu et al.
CSP          <- dlply(allometry, c("spp"), bCM)
CSP.coef     <- ldply(CSP, coef)
CSP.coef$lam <- ddply(allometry, c("spp"), function(df)  bCM(df, lam.only = TRUE))[, "V1"]
ldply(CSP, BIC)
ldply(CSP, r.squaredGLMM)
table(allometry$spp[!is.na(allometry$sample)])




# Relative growth rates ---------------------------------------------------

### estimate plant masses from allometry
for (i in 1:nrow(stemHeights)) {
  if (stemHeights$species[i] %in% "DS") {
    stemHeights$mass[i] <- (stemHeights$height_cm[i] * CSP.coef[1, 3] + CSP.coef[1, "(Intercept)"])^(1/CSP.coef$lam[1]) # DISP
  } else if (stemHeights$species[i] %in% "SA") {
    stemHeights$mass[i] <- (stemHeights$height_cm[i] * CSP.coef[2, 3] + CSP.coef[2, "(Intercept)"])^(1/CSP.coef$lam[2]) # SPAL
  }
}


stemHeights$RGR <- stemHeights$lateStartDate <- NA

pb <- txtProgressBar(min = 0, max = length(unique(stemHeights$id[stemHeights$dead_live %in% "L"])), initial = 0, style = 3)
for(i in 1:length(unique(stemHeights$id[stemHeights$dead_live %in% "L"]))) { # 830 plants
  targPlant <- unique(stemHeights$id[stemHeights$dead_live %in% "L"])[i]
  subDat <- stemHeights[stemHeights$id %in% targPlant, ]
  if (nrow(subDat) > 1) {
    subDat$mass2 <- c(NA, subDat$mass[1:(nrow(subDat)-1)])
    subDat$day2 <- c(subDat$day[1], subDat$day[1:(nrow(subDat)-1)])
    subDat$RGR <- (log(subDat$mass) - log(subDat$mass2)) / as.numeric(difftime(subDat$day, subDat$day2, units = "days"))

    ### add RGRs to stem height dataset
    stemHeights$RGR[(which(rownames(stemHeights) %in% rownames(subDat)))] <- subDat$RGR
    if (!is.na(subDat$height[1]) & (subDat$height[1] < 6)) {
      stemHeights$lateStartDate[which(rownames(stemHeights) %in% rownames(subDat))] <- 1 # indicates if high-growth period is captured
    }
  }
  setTxtProgressBar(pb, i)
}

stemHeights$RGR[!is.finite(stemHeights$RGR)] <- NA

# mean growth rate for each plant, then each mesocosm
mean.rgr     <- ddply(stemHeights[stemHeights$dead_live %in% "L", ], .(id, new.core.id, species), summarise, meanRGR = mean(RGR, na.rm = TRUE))
mean.rgr.pot <- ddply(mean.rgr, .(new.core.id, species), summarise, RGR = mean(meanRGR, na.rm = TRUE))
mean.rgr.pot <- mean.rgr.pot[1:24, ]


# Stem density and biomass ------------------------------------------------

### summarize by pot and session
dat.live <- stemHeights[stemHeights$dead_live %in% "L", ]
ddHgt2 <- ddply(dat.live, .(date, day, core_num, species), numcolwise(sum, na.rm = TRUE))
ddHgt2$day <- as.POSIXct(ddHgt2$day, origin = "1960-01-01")
ddHgt2$cohort <- NA

for (i in 1:length(unique(ddHgt2$date))) {
  uniqueIDS <- unique(ddHgt2$core_num[ddHgt2$date == sort(unique(ddHgt2$date))[i]])
  if (length(uniqueIDS) == 3) {
    ddHgt2$cohort[ddHgt2$core_num %in% uniqueIDS] <- i
  }
}
ddHgt2$cohort <- (ddHgt2$cohort - 1) / 2
ddHgt2$cohort[is.na(ddHgt2$cohort)] <- 0


ddHgt4 <- ddply(ddHgt2, .(day, cohort, species), numcolwise(mean, na.rm = TRUE))
ddHgt4.se <- ddply(ddHgt2, .(day, cohort, species), numcolwise(se))
names(ddHgt4.se) <- paste0(names(ddHgt4.se), ".se")
ddHgt4 <- cbind(ddHgt4, ddHgt4.se[, c(5, 7, 9)])
ddHgt4$session <- as.POSIXct(ddHgt4$day, origin = "1960-01-01")
ddHgt4$species2 <- ifelse(ddHgt4$species %in% "SA", "italic(S.~alterniflora)", "italic(D.~spicata)")



####
#### Calculate Milner-Hughes NAPP
####

# get live + dead biomass for each pot at each date
dat.ld <- ddply(stemHeights, .(core_num, species, date), summarise,
                live = sum(mass[dead_live %in% "L"], na.rm = TRUE),
                dead = sum(mass[dead_live %in% "D"], na.rm = TRUE)
)

# add cohort number
dat.ld$day <- as.POSIXct(as.character(dat.ld$date), format = "%y%m%d", origin = "1960-01-01")
dat.ld$cohort <- NA

for (i in 1:length(unique(dat.ld$day))) {
  uniqueIDS <- unique(dat.ld$core_num[dat.ld$day == sort(unique(dat.ld$day))[i]])
  if (length(uniqueIDS) == 3) {
    dat.ld$cohort[dat.ld$core_num %in% uniqueIDS] <- ifelse(sum(is.finite(dat.ld$cohort)) == 0, 1, max(dat.ld$cohort, na.rm = TRUE) + 1)
  }
}
dat.ld$cohort[is.na(dat.ld$cohort)] <- 0
dat.ld$cohort[dat.ld$cohort == 5] <- 4



dat.ld$year <- 2016
dat.ld$dead <- 0
dat.ld$pot2  <- paste0(dat.ld$core_num, "-", dat.ld$species)
test <- nappCalc2(dataset = dat.ld[dat.ld$cohort > 0, ], liveCol = "live", siteCol = "pot2", timeCol = "day", summarize = "TRUE",
                  annualReset = "FALSE")
napp <- test$summary
napp$pot2 <- as.character(napp$site)
napp$core_num <- as.numeric(sapply(strsplit(napp$pot2, "-"),  "[[", 1))
napp$species  <- sapply(strsplit(napp$pot2, "-"),  "[[", 2)

napp$cohort <- dat.ld$cohort[match(napp$pot2, dat.ld$pot2)]
napp$new.core.id <- ifelse(napp$core_num > 9, paste0(napp$species, napp$core_num), paste0(napp$species, "0", napp$core_num))


### build object to merge with master: average RGR, stem characteristics
# mean growth rate for each plant, then each mesocosm
mean.rgr     <- ddply(stemHeights[stemHeights$dead_live %in% "L", ], .(id, new.core.id, species), summarise, meanRGR = mean(RGR, na.rm = TRUE))
mean.rgr.pot <- ddply(mean.rgr, .(new.core.id, species), summarise, rgr = mean(meanRGR, na.rm = TRUE))
mean.rgr.pot <- mean.rgr.pot[1:24, ]

mean.dens     <- ddply(stemHeights[stemHeights$dead_live %in% "L", ], .(new.core.id, date), summarise, stem.density = sum(!is.na(height_cm)))
mean.dens.pot <- ddply(mean.dens, .(new.core.id), summarise, dens = mean(stem.density / pot.m2, na.rm = TRUE)) # mean of all stem density msrmts for a pot
mean.dens.pot <- mean.dens.pot[1:24, ]

rgr.stems <- join_all(list(mean.rgr.pot, mean.dens.pot), by = "new.core.id")



# Denitrification enzyme activity -----------------------------------------

dea$species <- substr(dea$pot, 1, 2)
dea$species2 <- ifelse(dea$species %in% "SA", "italic(S.~alterniflora)", "italic(D.~spicata)")
dea$DEA.m2   <- dea$DEA * dea$bd_gcm3 * (top.vol / core.area) # * 5 * 1e4 # [flux per gram] * [g/cm3] * [5 cm depth] * [cm2/m2]
dea$IV.m2    <- dea$IV * dea$bd_gcm3  * (top.vol / core.area) # * 5 * 1e4
dea$gap.g    <- dea$DEA - dea$IV
dea$gap.m2   <- dea$DEA.m2 - dea$IV.m2

dd.dea <- ddply(dea, .(species, species2), summarise,
                flux.g = mean(DEA),
                flux.m2 = mean(DEA.m2),
                flux.g.se = se(DEA),
                flux.m2.se = se(DEA.m2)
)
dd.dea$type <- "DEA"

dd.iv <- ddply(dea, .(species, species2), summarise,
               flux.g = mean(IV),
               flux.m2 = mean(IV.m2),
               flux.g.se = se(IV),
               flux.m2.se = se(IV.m2)
)
dd.iv$type <- "GHG"

dd.dea <- rbind(dd.dea, dd.iv)


### stats on DEA/IV N2O flux
dea.m <- melt(dea, id.vars = c("species", "species2", "pot"), measure.vars = c(2, 3, 8:11))
head(dea.m)
dea.m$ID       <- paste0(dea.m$species, "-", dea.m$variable)

dea.g.aov <- aov(value ~ factor(ID), data = dea.m[(dea.m$variable %in% c("IV.m2", "DEA.m2")), ])
Anova(dea.g.aov, type = 3, singular = TRUE)
TukeyHSD(dea.g.aov)




# 15N data ----------------------------------------------------------------


# how much 15N was added to each pot?
KNO3 <- 2.24047           # grams of KNO3 added
soln <- 2027.5            # final volume of spike solution (mls)
KNO3_mw <- 102.1032       # g/mol
N_mw <- 15                # gN/mol KNO3
mls  <- 12                # mls added per pot
spike <- KNO3 * (N_mw / KNO3_mw) / soln * mls # total grams of 15N added to each pot


# calculate atom pct from mean values
CN_mass_data$n_ap <- ap(CN_mass_data$d15n)


# get background 15N AP
bkd <- ddply(CN_mass_data[as.character(CN_mass_data$time) %in% "t0", ], .(species, pool_label), summarise,
             n15 = mean(n_ap, na.rm = T),
             se.n15 = se(n_ap)
)



# subtract background to get 15N excess (if result is negative, 15Nxs = 0)
CN_mass_data$n15xs <- as.numeric(NA)
for(i in 1:nrow(CN_mass_data)) {
  if(!is.na(CN_mass_data$species[i])) {
    spp <- CN_mass_data$species[i]
    pool <- CN_mass_data$pool_label[i]

    a <- bkd[(bkd$species %in% spp) & (bkd$pool_label %in% pool), "n15"]

    if((is.numeric(a) & (length(a) > 0))) { # excludes all pools we don't have bkd for
      difference <- CN_mass_data$n_ap[i] - a # NOT decimal fraction
      ifelse (!is.na(difference) & (difference > 0), CN_mass_data$n15xs[i] <- difference, CN_mass_data$n15xs[i] <- 0)
    } else {
      difference <- CN_mass_data$n_ap[i] - 0.370317701 # assume bkgd of 11 d15N
      ifelse ((difference > 0) & !is.na(difference), CN_mass_data$n15xs[i] <- difference, CN_mass_data$n15xs[i] <- 0)
    }
  }
}


### Additional correction to belowground live biomass, subtracting 15Nxs of dead MOM to account for sorption
CN_mass_data$n15xs_MOM <- CN_mass_data$n15xs # n15xs still useful for total recoveries
for(i in 1:nrow(CN_mass_data)) {
  if (CN_mass_data$time[i] %in% paste0("t", 1:4)) {
    if(!is.na(CN_mass_data$depth_top[i])) {
      spp <- CN_mass_data$species[i]
      coreID <- CN_mass_data$new.core.id[i]
      pool <- CN_mass_data$pool_label[i]
      depth_bottom <- CN_mass_data$depth_bottom[i]

      a <- CN_mass_data[(CN_mass_data$new.core.id %in% coreID) & (CN_mass_data$pool_label %in% paste0("dead biomass ", depth_bottom, "cm")) &
                   (CN_mass_data$depth_bottom == depth_bottom), "n15xs"]

      if((is.numeric(a) & (length(a) > 0))) { # excludes all pools we don't have bkd for
        difference <- CN_mass_data$n15xs[i] - a # NOT decimal fraction
        ifelse (!is.na(difference) & (difference > 0), CN_mass_data$n15xs_MOM[i] <- difference, CN_mass_data$n15xs_MOM[i] <- 0)
      } else {
        difference <- CN_mass_data$n_ap[i] - 0.370317701 # assume bkgd of 11 d15N
        ifelse ((difference > 0) & !is.na(difference), CN_mass_data$n15xs_MOM[i] <- difference, CN_mass_data$n15xs_MOM[i] <- 0)
      }
    }
  }
}



CN_mass_data$n15_g_pg_recov  <- CN_mass_data$n15xs  / 100 * CN_mass_data$n_pct  # for total recovery
CN_mass_data$n15_g_recov     <- CN_mass_data$n15_g_pg_recov  * CN_mass_data$g_core  # for total recovery

CN_mass_data$n15_g_pg <- CN_mass_data$n15xs_MOM / 100 * CN_mass_data$n_pct
CN_mass_data$n15_g    <- CN_mass_data$n15_g_pg  * CN_mass_data$g_core
CN_mass_data$n_core   <- CN_mass_data$n_pct * CN_mass_data$g_core

# Q: how much sorption is there? (look just at belowground data from t1-4)
summary(CN_mass_data$n15xs[!is.na(CN_mass_data$depth_bottom) & (CN_mass_data$time %in% paste0("t", 1:4))])
summary(CN_mass_data$n15xs_MOM[!is.na(CN_mass_data$depth_bottom) & (CN_mass_data$time %in% paste0("t", 1:4))])

summary(CN_mass_data$n15xs[(CN_mass_data$depth_bottom < 11) & (CN_mass_data$sample.type %in% "dead biomass") & (CN_mass_data$time %in% paste0("t", 1:4))])
summary(CN_mass_data$n15xs_MOM[(CN_mass_data$depth_bottom < 11) & (CN_mass_data$sample.type %in% "dead biomass") & (CN_mass_data$time %in% paste0("t", 1:4))])

# look at difference between with and without MOM correction
summary(CN_mass_data$n15xs[!is.na(CN_mass_data$depth_bottom) & (CN_mass_data$time %in% paste0("t", 1:4))] - CN_mass_data$n15xs_MOM[!is.na(CN_mass_data$depth_bottom) & (CN_mass_data$time %in% paste0("t", 1:4))])
# A: actually kind of high. mean = 0.05% compared with mean 15Nxs AP = 0.25%

### does sorption vary by species?  No.
summary(CN_mass_data$n15xs[CN_mass_data$sample.type %in% "dead biomass"])
summary(aov(n15xs ~ species, data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0"), ]))
summary(aov(n15xs ~ species, data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0") & (CN_mass_data$depth_bottom < 11), ])) # doesn't change if focused only on top 10 cm

### does sorption change over time? No
summary(lm(n15xs ~ species + as.numeric(substr(time, 2, 2)), data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0"), ])) # no
summary(lm(n15xs ~ as.numeric(substr(time, 2, 2)), data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0"), ]))
summary(lm(n15xs ~ as.numeric(substr(time, 2, 2)), data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0") & (CN_mass_data$depth_bottom < 11), ])) # doesn't change if focused only on top 10 cm
plot(n15xs ~ as.numeric(substr(time, 2, 2)), data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0") & (CN_mass_data$depth_bottom < 6), ])





# get 15Nxs mass in each pool in each core
# total recovery using belowground pools (rather than bulk sediment)
mat2 <- ddply(CN_mass_data[!CN_mass_data$sample.type %in% c("bulk sediment"), ], .(time, species, new.core.id), summarise,
              n15_2 = sum(n15_g_recov, na.rm  = T)
)
mat2$recovery2 <- mat2$n15_2 / (spike) # "recovery" uses bulk belowground data, "recovery2" uses root pools
summary(mat2$recovery2[mat2$recovery2 > 0.05])



### relps between napp and 15n recoveries (excludes dead biomass)
bgd <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("coarse roots", "fine roots", "rhizomes"), ],
             .(time, species, new.core.id), summarise,
             g        = sum(g_core, na.rm = T),
             n15    = sum(n15_g, na.rm  = T),
             n_core   = sum(n_core, na.rm = T)
)
bgd$g_pg <- bgd$n15 / bgd$g


# species-level belowground inventory over time
abv <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("decomp layer", "leaf 1", "leaf 2", "leaf 3", "leaf 4",
                                            "leaf 5", "leaf 6", "leaf 7", "leaf 8", "leaf 9", "leaf 10", "leaf 11",
                                            # "dead leaf", "standing dead", "microbe mat",
                                            "standing dead", "stems"), ], .(time, species, new.core.id), summarise,
             g        = sum(g_core, na.rm = T),
             n15      = sum(n15_g, na.rm  = T),
             n_core   = sum(n_core, na.rm = T)
)
abv$g_pg <- abv$n15 / abv$g

names(abv)[4:7] <- paste0(names(abv)[4:7], ".ag")
names(bgd)[4:7] <- paste0(names(bgd)[4:7], ".bg")

comb <- join_all(list(napp, abv, bgd, mat2), by = "new.core.id")




###
### compare predicted (allometry) and observed final masses
###

### allometry estimates
for (i in 1:length(unique(dat.ld$pot2))) {
  targPot <- unique(dat.ld$pot2)[i]
  maxDate <- max(dat.ld$day[dat.ld$pot2 %in% targPot])
  tempDat <- dat.ld[(dat.ld$pot2 %in% targPot) & (dat.ld$day == maxDate), ]
  if (i == 1) {
    dat.ld.sub  <- tempDat
  }
  if (i > 1) {
    dat.ld.sub <- rbind(dat.ld.sub, tempDat)
  }
}
dat.ld.sub$id <- paste0(dat.ld.sub$species, "-", dat.ld.sub$core_num)

### observed values
obs <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("leaf 1", "leaf 2", "leaf 3", "leaf 4", # "dead leaf",
                                            "leaf 5", "leaf 6", "leaf 7", "leaf 8", "leaf 9", "leaf 10", "leaf 11",
                                            "stems"), ], .(time, species, new.core.id), summarise,
             g        = sum(g_core, na.rm = T)
)
obs$core_num <- as.character(as.integer(substr(obs$new.core.id, 3, 4)))
obs$core_num[as.integer(obs$core_num) > 12] <- paste0(obs$species[as.integer(obs$core_num) > 12], as.integer(obs$core_num[as.integer(obs$core_num) > 12]) - 12, "_T0")
obs$id <- paste0(obs$species, "-", obs$core_num)

obs <- join_all(list(obs, dat.ld.sub[, c(4,9:10)]), by = "id")
names(obs)[which(names(obs) %in% "live")] <- "allom.est"
obs$diff <- (obs$allom.est - obs$g)             # magnitude accuracy of allometry prediction (g)
obs$diff.pct <- (obs$allom.est - obs$g) / obs$g # percent accuracy of allometry prediction
summary(obs$diff)
summary(obs[obs$species %in% "SA", ])
summary(obs[obs$species %in% "DS", ])

obs$species2 <- ifelse(obs$species %in% "SA", "italic(S.~alterniflora)", "italic(D.~spicata)")



###
### aboveground data
###

### mass over time in aboveground compartments. mean +- se by species (sum across depth intervals)

ag2 <- ddply(CN_mass_data[CN_mass_data$sample.type2 %in% c("stems", "leaf"), ],
             .(time, species, new.core.id, sample.type2), summarise,
             g   = sum(g_core / pot.m2, na.rm  = T),
             n_pct = mean(n_pct, na.rm = TRUE),
             n   = sum(n_core / pot.m2, na.rm  = T),
             n15 = sum(n15_g, na.rm  = T),
             n15_pg = n15 / (g * pot.m2),
             c13 = mean(d13c, na.rm = TRUE),
             cn  = sum(c_pct * g_core, na.rm = TRUE) / n
)

ag2.sp <- ddply(ag2, .(time, species, sample.type2), numcolwise(mean))
ag2.se <- ddply(ag2, .(time, species, sample.type2), numcolwise(se))
names(ag2.se)[4:10] <- paste0(names(ag2.sp)[4:10], ".se")
ag2.sp <- cbind(ag2.sp, ag2.se[, 4:10])
ag2.sp$t <- as.integer(substr(ag2.sp$time, 2, 2))
ag2.sp$session <- unique(ddHgt4$session)[c(1, 3, 5, 7, 9)][as.numeric(as.factor(ag2.sp$time))]
ag2.sp$species2 <- ifelse(ag2.sp$species %in% "SA", "italic(S.~alterniflora)", "italic(D.~spicata)")





###
### belowground data
###

### mass over time in belowground compartments. mean +- se by species (sum across depth intervals)
bg <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("coarse roots", "fine roots", "rhizomes"), ],
            .(time, species, new.core.id, sample.type), summarise,
            g   = sum(g_core / pot.m2, na.rm  = T),
            n_pct = mean(n_pct, na.rm = TRUE),
            n   = sum(n_core / pot.m2, na.rm  = T),
            n15 = sum(n15_g, na.rm  = T),
            n15_pg = n15 / (g * pot.m2)
)

bg.sp <- ddply(bg, .(time, species, sample.type), numcolwise(mean))
bg.se <- ddply(bg, .(time, species, sample.type), numcolwise(se))
names(bg.se)[4:8] <- paste0(names(bg.se)[4:8], ".se")
bg.sp <- cbind(bg.sp, bg.se[, 4:8])
bg.sp$t <- as.integer(substr(bg.sp$time, 2, 2))
bg.sp$session <- unique(ddHgt4$session)[c(1, 3, 5, 7, 9)][as.numeric(as.factor(bg.sp$time))]
bg.sp$species2 <- ifelse(bg.sp$species %in% "SA", "italic(S.~alterniflora)", "italic(D.~spicata)")





# 15N inventories ---------------------------------------------------------
# all vegetative compartments
all_pools <- ddply(CN_mass_data[CN_mass_data$sample.type2 %in% 
                  c("stems", "leaf", "coarse roots", "fine roots", "rhizomes"), ],
            .(time, species, sample.type2, new.core.id), summarise,
            g   = sum(g_core / pot.m2, na.rm  = T),
            n_pct = mean(n_pct, na.rm = TRUE),
            n   = sum(n_core / pot.m2, na.rm  = T),
            n15 = sum(n15_g, na.rm  = T) # g per pot
)

all_pools1 <- ddply(all_pools, .(time, species, sample.type2), summarise,
              n15.mean    = mean(n15, na.rm  = T),
              n15.se      = se(n15),
              recovery    = mean(n15 / spike, na.rm = TRUE),
              recovery.se = se(n15 / spike)
)
all_pools2 <- ddply(all_pools[!all_pools$time %in% "t0", ], .(species, sample.type2), summarise,
                    n15.mean    = mean(n15, na.rm  = T),
                    n15.se      = se(n15),
                    recovery    = mean(n15 / spike, na.rm = TRUE),
                    recovery.se = se(n15 / spike)
)

# aboveground
ag <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("belowground stems", "stems", "leaf", paste0("leaf ", 1:10)), ],
            .(time, species, new.core.id), summarise,
            g   = sum(g_core / pot.m2, na.rm  = T),
            n_pct = mean(n_pct, na.rm = TRUE),
            n   = sum(n_core / pot.m2, na.rm  = T),
            n15 = sum(n15_g, na.rm  = T) # g per pot
)

agd1 <- ddply(ag, .(time, species), summarise,
              n15.mean    = mean(n15, na.rm  = T),
              n15.se      = se(n15),
              recovery    = mean(n15 / spike, na.rm = TRUE),
              recovery.se = se(n15 / spike)
)
# species-level belowground inventory over time
bgd <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("coarse roots", "fine roots", "rhizomes"), ], #, "dead biomass"
             .(time, species, new.core.id), summarise,
             n15    = sum(n15_g, na.rm  = T)
)
bgd1 <- ddply(bgd, .(time, species), summarise,
              n15.mean    = mean(n15, na.rm  = T),
              n15.se      = se(n15),
              recovery    = mean(n15 / spike, na.rm = TRUE),
              recovery.se = se(n15 / spike)
)

bgd1$type <- "Belowground"
agd1$type <- "Aboveground"
mgd <- rbind(agd1, bgd1)

tot <- ddply(mat2, .(time, species), summarise,
             n15.mean    = mean(n15_2, na.rm  = T),
             n15.se      = se(n15_2),
             recovery    = mean(n15_2 / spike, na.rm = TRUE),
             recovery.se = se(n15_2 / spike)
)
tot$type <- "Total"

# total recovery combined above + belowground (includes dead aboveground but only live belowground)
tot.live <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("coarse roots", "fine roots", "rhizomes",
                                                 "leaf 1", "leaf 2", "leaf 3", "leaf 4",
                                                 "leaf 5", "leaf 6", "leaf 7", "leaf 8", "leaf 9", "leaf 10", "leaf 11",
                                                 "stems"), ], #, "dead biomass"
                  .(time, species, new.core.id), summarise,
                  n15    = sum(n15_g, na.rm  = T)
)
tot.live1 <- ddply(tot.live, .(time, species), summarise,
                   n15.mean    = mean(n15, na.rm  = T),
                   n15.se      = se(n15),
                   recovery    = mean(n15 / spike, na.rm = TRUE),
                   recovery.se = se(n15 / spike)
)


mgd <- rbind(mgd, tot)
# replace unique times with actual sampling dates using indexing
mgd$session <- unique(ddHgt4$session)[c(1, 3, 5, 7, 9)][as.numeric(as.factor(mgd$time))]


# N uptake rates ----------------------------------------------------------
# Aboveground: estimate from biomass and weighted average biomass N concentration
# get weighted average of aboveground N concentrations for each core
n_mean <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("belowground stems", "stems", "leaf", paste0("leaf ", 1:10)), ],
                .(species, time, new.core.id), summarise,
                tot_mass = sum(g_core / pot.m2, na.rm = TRUE),
                tot_n    = sum(n_core / pot.m2, na.rm = TRUE),
                tot_c    = sum(c_pct * g_core / pot.m2, na.rm = TRUE),
                n_wa = tot_n / tot_mass,
                c_wa = tot_c / tot_mass
)
n_bg_mean <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("fine roots", "coarse roots", "rhizomes"), ],
                   .(species, time, new.core.id), summarise,
                   tot_mass_bg = sum(g_core / pot.m2, na.rm = TRUE),
                   tot_n_bg    = sum(n_core / pot.m2, na.rm = TRUE),
                   tot_c_bg    = sum(c_pct * g_core  / pot.m2, na.rm = TRUE),
                   n_wa_bg = tot_n_bg / tot_mass_bg,
                   c_wa_bg = tot_c_bg / tot_mass_bg
)
bg_del <- ddply(bg, .(new.core.id), summarise, # I think all these values are per m2
                g_bg     = sum(g, na.rm = TRUE),
                g_froots = sum(g[sample.type %in% "fine roots"], na.rm = TRUE),
                g_roots  = sum(g[sample.type %in% c("fine roots", "coarse roots")], na.rm = TRUE),
                n_bg     = sum(n, na.rm = TRUE),
                n_pct_bg = mean(n_pct * g, na.rm = TRUE) / g_bg,
                n_pct_bg2 = n_bg / g_bg,
                n15_bg   = sum(n15, na.rm = TRUE)
)

# all relevant values are per m2
master <- join_all(list(napp[, c("napp.MH", "new.core.id")], # primary production
                        rgr.stems,
                        n_mean, n_bg_mean[, -c(1:2)], # N inventories above and belowground
                        ag[, 3:7], bg_del # 15N recoveries and unweighted average pct N
), by = "new.core.id")
master$napp.MH  <- master$napp.MH / pot.m2 # g/m2
master$n15_core <- master$n15 + master$n15_bg # total live biomass inventory
master$recovery <- master$n15_core / spike
master$session  <- unique(ddHgt4$session)[c(3, 5, 7, 9)][as.numeric(as.factor(master$time))]

#  15N applied on 24 June; first stem height msrmts on 22 June
# use master$timeDiff for 15N rate calculations
# use master$timeDiff + 2 for NAPP rate calculations
master$timeDiff <- as.numeric(difftime(master$session, "2016-06-24", units = "days"))

master$n_uptake_15n     <- master$n15 * 1e3 / master$timeDiff / pot.m2 # mg N/day/m2 accumulating in aboveground tissue, using 15N
master$n_uptake_15n_bg  <- master$n15_bg * 1e3 / master$timeDiff /pot.m2 # mg N/day/m2
master$prodn_rate       <- master$napp.MH / (master$timeDiff + 2) # g/m2/day
master$n_uptake_biomass <- master$prodn_rate * 1e3 * master$n_wa # mg N/m2/day accumulating in aboveground tissue, using primary production rate and weighted average pct N
master$cn_ag    <- master$c_wa / master$n_wa # C:N ratio in aboveground tissue
master$cn_bg    <- master$c_wa_bg / master$n_wa_bg
master$cluster  <- ifelse(master$time %in% c("t1", "t2"), "1-2 weeks", "3-4 weeks")
master$n15_pgBG <- (master$n_uptake_15n + master$n_uptake_15n_bg) / master$g_bg # 15N uptake per gram belowground biomass: mg 15N/day/g TOTAL BG biomass
master$n15_pgFR <- (master$n_uptake_15n + master$n_uptake_15n_bg) / master$g_froots # 15N uptake per gram coarse roots: mg 15N/day/g fine roots
master$n15_pgCR <- (master$n_uptake_15n + master$n_uptake_15n_bg) / (master$g_roots - master$g_froots) # 15N uptake per gram coarse roots: mg 15N/day/g coarse roots
master$n15_pgR  <- (master$n_uptake_15n + master$n_uptake_15n_bg) / master$g_roots # 15N uptake per gram fine+coarse root biomass: mg 15N/day/g root biomass
master$n15_pgFR2 <- master$n15_core * 1e3 / master$g_froots           # mg 15N / g fine roots
master$tot_15n_uptake <- master$n_uptake_15n + master$n_uptake_15n_bg # mg N/day/m2

### apply relationship between total N uptake and 15N uptake to estimate belowground production
summary(lm3_4 <- lm(n_uptake_biomass ~ I(n_uptake_15n ) , data = master[master$time %in% c("t3", "t4"), ]))
master$bg_n_est <- as.numeric(master$n_uptake_15n_bg) * coefficients(lm3_4)[2] + coefficients(lm3_4)[1] # total N uptake belowground: mg N/m2/day
master$bg_n_est[master$time %in% c("t1", "t2")] <- NA
master$bgp_est <- master$bg_n_est / master$n_pct_bg2 / 1e3 # [estimated total N uptake (mg N/day)] / [mg N / mg biomass] / [1000 mg/g] = [g belowground biomass / day (/m2)]
master$bgp_biomass_est <- master$bgp_est * master$timeDiff # [g belowground biomass / day] * [timeDiff (days)] = total belowground production during period [g/m2]

summary(master$bgp_est / master$prodn_rate)

### differences in total nitrogen interception between species
master$tot_N_int <- master$bg_n_est + master$n_uptake_biomass + ifelse(master$species %in% "SA", nmolHr_mgDay(dd.dea$flux.m2[2]), nmolHr_mgDay(dd.dea$flux.m2[1]))
master$tot_N_int[master$time %in% c("t1", "t2")] <- NA


dd.master <- ddply(master, .(time, species), numcolwise(mean))
dd.master.se <- ddply(master, .(time, species), numcolwise(se))
dd.master$n_uptake_biomass.se <- dd.master.se$n_uptake_biomass
dd.master$prodn_rate.se <- dd.master.se$prodn_rate
dd.master$session <- unique(ddHgt4$session)[c(3, 5, 7, 9)][as.numeric(as.factor(dd.master$time))]





bgp.out <- ddply(master[master$time %in% c("t3", "t4"), ], .(species), summarise,
                 bg_n_uptake         = mean(bg_n_est), # mg/m2/d
                 bg_n_uptake.se      = se(bg_n_est),
                 ag_n_uptake         = mean(n_uptake_biomass),
                 ag_n_uptake.se      = se(n_uptake_biomass),

                 bg_production       = mean(bgp_est), # g/m2/day
                 bg_production.se    = se(bgp_est),
                 ag_production       = mean(prodn_rate), # g/m2/day
                 ag_production.se    = se(prodn_rate),

                 new_bg_biomass      = mean(bgp_biomass_est), # g/m2
                 new_bg_biomass.se   = se(bgp_biomass_est),
                 root_shoot_ratio    = mean(bgp_est / prodn_rate),
                 root_shoot_ratio.se = se(bgp_est / prodn_rate)
)
# write.csv(bgp.out, file = "C:/RDATA/greenhouse/output/GRO/bgp_estimates.csv", row.names = FALSE)









# Documentation of statements made in text --------------------------------

# Abstract ----------------------------------------------------------------
ddply(master, .(species), summarise,
      Nint = mean(tot_N_int, na.rm = TRUE),
      Nint.se = se(tot_N_int),
      production    = mean(bgp_est + prodn_rate, na.rm = TRUE),
      production.se = se(bgp_est + prodn_rate)
      )


# Results: a. Stem allometry, aboveground biomass and production ----------
table(allometry$spp[!is.na(allometry$sample)])
ldply(CSP, r.squaredGLMM)

# change over time in production?
summary(lm.sa.prod <- lm(napp.MH ~ session, data = master[master$species %in% "SA",]))
summary(lm.ds.prod <- lm(napp.MH ~ session, data = master[master$species %in% "DS",]))

ddply(master, .(species), summarise,
      production    = mean(prodn_rate, na.rm = TRUE),
      production.se = se(prodn_rate)
)
t.test(prodn_rate ~ species, data = master)

for (i in 1:length(unique(master$session))) {
  print(unique(master$session)[i])
  print(t.test(prodn_rate ~ species, data = master[(master$session %in% c(unique(master$session)[i])), ]))
}




# Results: b.	Leaf and stem biomass, N pools, and uptake ------------------
# leaf, stem biomass differences
leaf <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("leaf", paste0("leaf ", 1:10)), ],
              .(species, time, new.core.id), summarise,
              tot_mass = sum(g_core / pot.m2, na.rm = TRUE),
              tot_n    = sum(n_core / pot.m2, na.rm = TRUE),
              tot_c    = sum(c_pct * g_core / pot.m2, na.rm = TRUE),
              n_wa = tot_n / tot_mass,
              c_wa = tot_c / tot_mass)
leaf$t <- as.numeric(substr(leaf$time, 2, 2)) * 7

stem <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("belowground stems", "stems"), ],
              .(species, time, new.core.id), summarise,
              tot_mass = sum(g_core / pot.m2, na.rm = TRUE),
              tot_n    = sum(n_core / pot.m2, na.rm = TRUE),
              tot_c    = sum(c_pct * g_core / pot.m2, na.rm = TRUE),
              n_wa = tot_n / tot_mass,
              c_wa = tot_c / tot_mass
)
stem$t <- as.numeric(substr(stem$time, 2, 2)) * 7

leaf2 <- leaf
names(leaf2)[4:8] <- paste0(names(leaf)[4:8], ".leaf")
leaf2 <- cbind(leaf2, stem[4:8])
leaf2$n_wa_ag <- (leaf2$tot_n.leaf + leaf2$tot_n) / (leaf2$tot_mass.leaf + leaf2$tot_mass)

ddply(leaf[, -c(2, 3)], .(species), numcolwise(mean))
ddply(leaf[, -c(2, 3)], .(species), numcolwise(se))

# leaf biomass
summary(lm.sa.leaf <- lm(tot_mass ~ t, data = leaf[leaf$species %in% "SA",]))
summary(lm.ds.leaf <- lm(tot_mass ~ t, data = leaf[leaf$species %in% "DS",]))
t.test(tot_mass ~ species, data = leaf)

# leaf N content
summary(lm.sa.leafN <- lm(n_wa ~ t, data = leaf[leaf$species %in% "SA",]))
summary(lm.ds.leafN <- lm(n_wa ~ t, data = leaf[leaf$species %in% "DS",]))
t.test(n_wa ~ species, data = leaf)


# stems 
ddply(stem[, -c(2, 3)], .(species), numcolwise(mean))
ddply(stem[, -c(2, 3)], .(species), numcolwise(se))

# stem biomass
summary(lm.sa.stem <- lm(tot_mass ~ t, data = stem[stem$species %in% "SA",]))
summary(lm.ds.stem <- lm(tot_mass ~ t, data = stem[stem$species %in% "DS",]))
t.test(tot_mass ~ species, data = stem)

# stem N content
summary(lm.sa.leafN <- lm(n_wa ~ t, data = leaf[leaf$species %in% "SA",]))
summary(lm.ds.leafN <- lm(n_wa ~ t, data = leaf[leaf$species %in% "DS",]))
t.test(n_wa ~ species, data = leaf)

# mean N content for all
ddply(stem[, -c(2, 3)], .(), numcolwise(mean))
ddply(stem[, -c(2, 3)], .(), numcolwise(se))

# stem N stock
t.test(tot_n ~ species, data = stem)


# aboveground weighted average N content
summary(lm.sa.N.ag <- lm(n_wa_ag ~ t, data = leaf2[leaf2$species %in% "SA",]))
summary(lm.ds.N.ag <- lm(n_wa_ag ~ t, data = leaf2[leaf2$species %in% "DS",]))
t.test(n_wa_ag ~ species, data = leaf2)

ddply(leaf2, .(species, time), summarise,
      n_pct = mean(n_wa_ag),
      n_pct.se = se(n_wa_ag))
ddply(leaf2, .(species), summarise,
      n_pct = mean(n_wa_ag),
      n_pct.se = se(n_wa_ag))
ddply(leaf2, .(), summarise,
      n_pct = mean(n_wa_ag),
      n_pct.se = se(n_wa_ag))

### total N uptake (NAPP * %N)
summary(lm.sa.Nupt <- lm(n_uptake_biomass ~ session, data = master[master$species %in% "SA",]))
summary(lm.ds.Nupt <- lm(n_uptake_biomass ~ session, data = master[master$species %in% "DS",]))
t.test(n_uptake_biomass ~ species, data = master) # no difference bt species as a whole

for (i in 1:length(unique(master$session))) {
  print(unique(master$session)[i])
  print(t.test(n_uptake_biomass ~ species, data = master[(master$session %in% c(unique(master$session)[i])), ]))
}
ddply(master, .(species), summarise,
      n_uptake = mean(n_uptake_biomass),
      n_uptake.se = se(n_uptake_biomass))




# Results: c.	Belowground biomass and N pools -----------------------------
# belowground biomass and tissue N
rts <- ddply(CN_mass_data[CN_mass_data$sample.type %in% c("fine roots", "coarse roots", "rhizomes"), ],
             .(species, time, new.core.id, sample.type), summarise,
             tot_mass = sum(g_core / pot.m2, na.rm = TRUE),
             tot_n    = sum(n_core / pot.m2, na.rm = TRUE),
             tot_c    = sum(c_pct * g_core / pot.m2, na.rm = TRUE),
             n_wa = tot_n / tot_mass,
             c_wa = tot_c / tot_mass
)
rts$t <- as.numeric(substr(rts$time, 2, 2)) * 7

ddply(rts[, -c(2, 3)], .(species, sample.type), numcolwise(mean))
ddply(rts[, -c(2, 3)], .(species, sample.type), numcolwise(se))

# fine roots
organ <- "fine roots"
summary(lm(tot_mass ~ t, data = rts[(rts$sample.type %in% organ) & (rts$species %in% "SA"), ]))
summary(lm(tot_mass ~ t, data = rts[(rts$sample.type %in% organ) & (rts$species %in% "DS"), ]))
t.test(tot_mass ~ species, data = rts[(rts$sample.type %in% organ), ])
t.test(n_wa ~ species, data = rts[(rts$sample.type %in% organ), ])
t.test(tot_n ~ species, data = rts[(rts$sample.type %in% organ), ])

# coarse roots
organ <- "coarse roots"
t.test(tot_mass ~ species, data = rts[(rts$sample.type %in% organ), ])
summary(lm(tot_mass ~ t, data = rts[(rts$sample.type %in% organ) & (rts$species %in% "SA"), ]))
summary(lm(tot_mass ~ t, data = rts[(rts$sample.type %in% organ) & (rts$species %in% "DS"), ]))


t.test(n_wa ~ species, data = rts[(rts$sample.type %in% organ), ])
t.test(tot_n ~ species, data = rts[(rts$sample.type %in% organ), ])



# rhizomes
organ <- "rhizomes"
summary(lm(tot_mass ~ t, data = rts[(rts$sample.type %in% organ) & (rts$species %in% "SA"), ]))
summary(lm(tot_mass ~ t, data = rts[(rts$sample.type %in% organ) & (rts$species %in% "DS"), ]))
t.test(tot_mass ~ species, data = rts[(rts$sample.type %in% organ), ])
t.test(n_wa ~ species, data = rts[(rts$sample.type %in% organ), ])
t.test(tot_n ~ species, data = rts[(rts$sample.type %in% organ), ])




# Results d.	15N recoveries and N uptake estimates -----------------------
# t.test(recovery2 ~ species, data = mat2[mat2$time %in% "t4", ])
summary(mat2$recovery2[mat2$recovery2 > 0.05])

# comparison of total 15N recovery by species
t.test(recovery2 ~ species, data = mat2)

# variation over time in total 15N recovery
mat2$days <- as.numeric(substr(mat2$time, 2, 2)) * 7
Anova(aov.tot <- aov(recovery2 ~ factor(species) * days, data = mat2), type = 2)
summary(lm(recovery2 ~ days, data = mat2[mat2$species %in% "DS",]))
summary(lm(recovery2 ~ days, data = mat2[mat2$species %in% "SA",]))

# comparison of belowground 15N recoveries between species
t.test(n15_bg ~ species, data = master)
Anova(aov.bg <- aov(n15_bg ~ species * days, data = master), type = 2)


# comparison of aboveground 15N recoveries between species
t.test(n15 ~ species, data = master)
master$days <- as.numeric(substr(master$time, 2, 2)) * 7
Anova(aov.ag.n15 <- aov(n15 ~ species * days, data = master), type = 2)
summary(lm(n15 ~ days, data = master[master$species %in% "DS",]))
summary(lm(n15 ~ days, data = master[master$species %in% "SA",]))


t.test(n15 ~ species, data = master)


# Q: how much sorption is there? (look just at belowground data from t1-4)
summary(CN_mass_data$n15xs[!is.na(CN_mass_data$depth_bottom) & (CN_mass_data$time %in% paste0("t", 1:4))])
summary(CN_mass_data$n15xs_MOM[!is.na(CN_mass_data$depth_bottom) & (CN_mass_data$time %in% paste0("t", 1:4))])

# Q: does sorption vary by species?  No.
summary(CN_mass_data$n15xs[CN_mass_data$sample.type %in% "dead biomass"])
summary(stats::aov(n15xs ~ species, data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0"), ]))
summary(stats::aov(n15xs ~ species, data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0") & (CN_mass_data$depth_bottom < 11), ])) # doesn't change if focused only on top 10 cm

### does sorption change over time? No
summary(stats::lm(n15xs ~ species + as.numeric(substr(time, 2, 2)), data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0"), ])) # no
summary(stats::lm(n15xs ~ as.numeric(substr(time, 2, 2)), data = CN_mass_data[(CN_mass_data$sample.type %in% "dead biomass") & (!CN_mass_data$time %in% "t0"), ]))


### differences in 15N-nitrogen interception between species (weeks 3 & 4)
t.test(n_uptake_15n ~ species, data = master[!is.na(master$bgp_biomass_est), ])
t.test(n_uptake_15n_bg ~ species, data = master[!is.na(master$bgp_biomass_est), ])

ddply(master[!is.na(master$bgp_biomass_est), ], .(species), summarise,
      bg.15n.uptake = mean(n_uptake_15n_bg),
      bg.15n.uptake.se = se(n_uptake_15n_bg),
      ag.15n.uptake = mean(n_uptake_15n),
      ag.15n.uptake.se = se(n_uptake_15n))

# 15N uptake per gram fine-roots 
master$spp2 <- factor(master$species)
contrasts(master$spp2) <- contr.poly(2) 
Anova(lm(n15_pgFR ~ timeDiff * spp2, data = master), type = "III")
summary(lm.sa <- lm(n15_pgFR ~ timeDiff, data = master[master$species %in% "SA", ]))
summary(lm.ds <- lm(n15_pgFR ~ timeDiff, data = master[master$species %in% "DS", ]))



# Calculate intersection of lines of best fit
rightSide <- lm.ds$coefficients[2] - lm.sa$coefficients[2]
leftSide  <- lm.sa$coefficients[1] - lm.ds$coefficients[1]
leftSide / rightSide # days until no difference detectable


# Analysis of N uptake drivers --------------------------------------------
summary(lm.mv <- lm(tot_15n_uptake ~ (napp.MH + rgr + g_roots + timeDiff) * factor(species), data = master))
rsq::rsq.partial(lm.mv) # partial r2 in a pooled model doesn't make sense

summary(lm.mv.sa <- lm(tot_15n_uptake ~ napp.MH + rgr + g_roots + timeDiff, data = master[master$species %in% "SA", ]))
summary(lm.mv.ds <- lm(tot_15n_uptake ~ napp.MH + rgr + g_roots + timeDiff, data = master[master$species %in% "DS", ]))
rsq::rsq.partial(lm.mv.sa)
rsq::rsq.partial(lm.mv.ds)

# belowground production
summary(lm3_4)


### differences in total nitrogen interception between species
t.test(tot_N_int ~ species, data = master)







# Results e.	Denitrification enzyme activity and total N interception ----
t.test(DEA.m2 ~ species, data = dea)
t.test(IV.m2 ~ species, data = dea)

dea.x <- ddply(dea, .(species), numcolwise(mean))
dea.se <- ddply(dea, .(species), numcolwise(se))
dea.x
dea.se
dea.se[, 6:7] / dea.x[, 6:7]

tableCols <- c("species", "bg_n_est", "n_uptake_biomass", "bgp_est", "prodn_rate")
ddply(master[master$time %in% c("t3", "t4"), tableCols], .(species), numcolwise(mean))
ddply(master[master$time %in% c("t3", "t4"), tableCols], .(species), numcolwise(se))

nmolHr_mgDay(dea.x[, c("DEA.m2", "IV.m2")]) 
nmolHr_mgDay(dea.se[, c("DEA.m2", "IV.m2")])

t.test(tot_N_int ~ species, data = master) # weeks three and four

### assemble Table 1: 15N recovery and uptake by compartment
CN_mass_data$day.no <- as.numeric(substr(CN_mass_data$time, 2, 2)) * 7
tbl.prep <- ddply(CN_mass_data[(!CN_mass_data$time %in% "t0") & 
                                 (CN_mass_data$sample.type2 %in% c("coarse roots", "fine roots", "rhizomes",
                                 "stems", "leaf")), ],
      .(day.no, species, new.core.id, sample.type2), summarise,
      n15xs          = sum(n15xs_MOM * n_core, na.rm = T) / sum(n_core, na.rm = TRUE),
      g              = sum(g_core / pot.m2, na.rm  = T),
      n15.mg         = sum(n15_g, na.rm  = T) * 1000 / pot.m2, # mg per m2
      n15.perg       = n15.mg * pot.m2 / sum(g_core, na.rm  = T)) # mg 15N per g dw
      # uptake_rate    = mean(ap(d15n) * n_pct * g_core, na.rm  = T) * 1000 / day.no,
      # sp_uptake_rate = sum(ap(d15n) * n_pct * g_core, na.rm  = T) * 1000 / sum(g_core, na.rm  = T) / day.no
tbl.prep$uptake_rate <- tbl.prep$n15.mg / tbl.prep$day.no
tbl.prep$sp_uptake_rate <- tbl.prep$n15.perg / tbl.prep$day.no

tbl.prep.for.stats <- ddply(tbl.prep, .(species, sample.type2), summarise,
                            n15 = mean(n15.mg),
                            n15.se   = se(n15.mg),
                            n15.pg = mean(n15.perg),
                            n15.pg.se = se(n15.perg),
                            uptake = mean(uptake_rate),
                            uptake.se = se(uptake_rate),
                            uptake.pg = mean(sp_uptake_rate),
                            uptake.pg.se = se(sp_uptake_rate)
                            )

tbl.prep.all <- ddply(CN_mass_data[(!CN_mass_data$time %in% "t0") & 
                                    (CN_mass_data$sample.type2 %in% c("coarse roots", "fine roots", "rhizomes",
                                                                      "stems", "leaf")), ],
                     .(day.no, species, new.core.id), summarise,
                     g              = sum(g_core / pot.m2, na.rm  = T),
                     n15.mg         = sum(n15_g, na.rm  = T) * 1000 / pot.m2, # mg per m2
                     n15.perg       = n15.mg * pot.m2 / sum(g_core, na.rm  = T)) # mg 15N per g dw
tbl.prep.all$uptake_rate <- tbl.prep.all$n15.mg / tbl.prep.all$day.no
tbl.prep.all$sp_uptake_rate <- tbl.prep.all$n15.perg / tbl.prep.all$day.no

tbl.prep.all2 <- ddply(tbl.prep.all, .(species), summarise,
                      n15 = mean(n15.mg),
                      n15.se   = se(n15.mg),
                      n15.pg = mean(n15.perg),
                      n15.pg.se = se(n15.perg),
                      uptake = mean(uptake_rate),
                      uptake.se = se(uptake_rate),
                      uptake.pg = mean(sp_uptake_rate),
                      uptake.pg.se = se(sp_uptake_rate)
)
tbl.prep.all2$sample.type2 <- "Total"

tbl.prep.ag <- ddply(CN_mass_data[(!CN_mass_data$time %in% "t0") & 
                                 (CN_mass_data$sample.type2 %in% c("stems", "leaf")), ],
                  .(day.no, species, new.core.id), summarise,
                  g              = sum(g_core / pot.m2, na.rm  = T),
                  n15.mg         = sum(n15_g, na.rm  = T) * 1000  / pot.m2, # mg per m2
                  n15.perg       = n15.mg  * pot.m2 / sum(g_core, na.rm  = T)) # mg 15N per g dw
# uptake_rate    = mean(ap(d15n) * n_pct * g_core, na.rm  = T) * 1000 / day.no,
# sp_uptake_rate = sum(ap(d15n) * n_pct * g_core, na.rm  = T) * 1000 / sum(g_core, na.rm  = T) / day.no
tbl.prep.ag$uptake_rate <- tbl.prep.ag$n15.mg / tbl.prep.ag$day.no
tbl.prep.ag$sp_uptake_rate <- tbl.prep.ag$n15.perg / tbl.prep.ag$day.no

tbl.prep.ag2 <- ddply(tbl.prep.ag, .(species), summarise,
                            n15 = mean(n15.mg),
                            n15.se   = se(n15.mg),
                            n15.pg = mean(n15.perg),
                            n15.pg.se = se(n15.perg),
                            uptake = mean(uptake_rate),
                            uptake.se = se(uptake_rate),
                            uptake.pg = mean(sp_uptake_rate),
                            uptake.pg.se = se(sp_uptake_rate)
)
tbl.prep.ag2$sample.type2 <- "Total aboveground"



tbl.prep.bg <- ddply(CN_mass_data[(!CN_mass_data$time %in% "t0") & 
                                    (CN_mass_data$sample.type2 %in% c("fine roots", "coarse roots", "rhizomes")), ],
                     .(day.no, species, new.core.id), summarise,
                     g              = sum(g_core / pot.m2, na.rm  = T),
                     n15.mg         = sum(n15_g, na.rm  = T) * 1000  / pot.m2, # mg per m2
                     n15.perg       = n15.mg  * pot.m2 / sum(g_core, na.rm  = T)) # mg 15N per g dw
# uptake_rate    = mean(ap(d15n) * n_pct * g_core, na.rm  = T) * 1000 / day.no,
# sp_uptake_rate = sum(ap(d15n) * n_pct * g_core, na.rm  = T) * 1000 / sum(g_core, na.rm  = T) / day.no
tbl.prep.bg$uptake_rate <- tbl.prep.bg$n15.mg / tbl.prep.bg$day.no
tbl.prep.bg$sp_uptake_rate <- tbl.prep.bg$n15.perg / tbl.prep.bg$day.no

tbl.prep.bg2 <- ddply(tbl.prep.bg, .(species), summarise,
                      n15 = mean(n15.mg),
                      n15.se   = se(n15.mg),
                      n15.pg = mean(n15.perg),
                      n15.pg.se = se(n15.perg),
                      uptake = mean(uptake_rate),
                      uptake.se = se(uptake_rate),
                      uptake.pg = mean(sp_uptake_rate),
                      uptake.pg.se = se(sp_uptake_rate)
)
tbl.prep.bg2$sample.type2 <- "Total belowground"

tbl.final <- plyr::rbind.fill(list(tbl.prep.for.stats, tbl.prep.ag2, tbl.prep.bg2, tbl.prep.all2))
# write.csv(tbl.final, "C:/RDATA/greenhouse/output/GRO/table_15n.csv")

for (i in 1:length(unique(tbl.prep$sample.type2))) {
  print(unique(tbl.prep$sample.type2)[i])
  print(paste("n15.mg P-value: ", t.test(n15.mg ~ species, data = tbl.prep[tbl.prep$sample.type2 %in% unique(tbl.prep$sample.type2)[i], ])$p.value))
  print(paste("uptake_rate P-value: ", t.test(uptake_rate ~ species, data = tbl.prep[tbl.prep$sample.type2 %in% unique(tbl.prep$sample.type2)[i], ])$p.value))
  }

setDat <- tbl.prep.bg # tbl.prep.all, tbl.prep.ag, tbl.prep.bg
paste("n15.mg P-value: ", t.test(n15.mg ~ species, data = setDat)$p.value)
paste("uptake_rate P-value: ", t.test(uptake_rate ~ species, data = setDat)$p.value)

# average recovery across species
ddply(tbl.prep.ag, .(), summarise,
      n15 = mean(n15.mg) / (spike * 1e3 / pot.m2),
      n15.se   = se(n15.mg) / (spike * 1e3 / pot.m2)
)
ddply(tbl.prep.bg, .(), summarise,
      n15 = mean(n15.mg) / (spike * 1e3 / pot.m2),
      n15.se   = se(n15.mg) / (spike * 1e3 / pot.m2)
)

# difference between above and belowground recovery by species
ag$type <- "aboveground"
bgd$type <- "belowground"
bgag <- rbind(ag[(!ag$time %in% "t0"), c(1:3, 7:8)], bgd[(!bgd$time %in% "t0"), ])
Anova(aov2 <- aov(n15 ~ type * species, data = bgag), type = 2)
TukeyHSD(aov2)

# stats for table 2
t.test(bg_n_est ~ species, data = master[master$time %in% c("t3", "t4"), ])
t.test(n_uptake_biomass ~ species, data = master[master$time %in% c("t3", "t4"), ]) # DS higher than SA
t.test(n_uptake_biomass ~ species, data = master[master$time %in% c("t3", "t4"), ])

t.test(bgp_est ~ species, data = master[master$time %in% c("t3", "t4"), ])
t.test(prodn_rate ~ species, data = master[master$time %in% c("t3", "t4"), ])  # DISP is more productive than SPAL
t.test(prodn_rate ~ species, data = master[master$time %in% c("t3", "t4"), ])

t.test(DEA.m2 ~ species, data = dea)
t.test(IV.m2 ~ species, data = dea)


# Table 3: literature estimates
ddply(tbl.prep.all, .(species, day.no), summarise,
      n15 = mean(n15.mg) / 243,
      n15.se   = se(n15.mg) / 243
)
t.test(n15.mg ~ species, data = tbl.prep.all)

# Discussion ---------------------------------------------------------------

ddply(master[, c("time", "species", "n15_pgFR")], .(species, time), numcolwise(mean))







# Figure 1 - Allometry ----------------------------------------------------

# png(filename = paste0("C:/RDATA/greenhouse/output/GRO/Figure1_", todaysDate, ".png"), width = 80, height = 80, units = "mm", res = 400)
par(fig = c(0, 1, 0, 1), mar = c(4, 4, 0.5, 0.5))
plot(sample ~ height_cm, data = allometry[(allometry$spp %in% "DISP"), ],
     cex = pointSize / 2, pch = 19, col = fig2Col,
     ylab = "Total mass (g)", xlab = "Height (cm)", xlim = c(0, 80), ylim = c(0, 0.65), xaxt = "n",
     las = 1, tcl = 0.25, tck = 0.01, bty = "n", yaxs = "i", xaxs = "i")
abline(h = 0)
abline(v = 0)
axis(side = 1, tcl = 0.25, tck = 0.01, at = axTicks(1), labels = axTicks(1))
points(sample ~ height_cm, data = allometry[(allometry$spp %in% "SPAL"), ],  cex = pointSize / 2, pch = 17)

x <- allometry[(allometry$spp %in% "DISP"), "height_cm"]
x.spal <- allometry[(allometry$spp %in% "SPAL"), "height_cm"]

y.pred2 <- (x * CSP.coef[1, 3] + CSP.coef[1, "(Intercept)"])^(1/CSP.coef$lam[1]) # DISP
y.pred3 <- (x.spal * CSP.coef[2, 3] + CSP.coef[2, "(Intercept)"])^(1/CSP.coef$lam[2]) # SPAL
lines(x = x[order(y.pred2)], y = y.pred2[order(y.pred2)], lty = 1, lwd = 2, col = fig2Col)
lines(x = x.spal[order(y.pred3)], y = y.pred3[order(y.pred3)], lty = 1, lwd = 2)

text(x = 20, y = 0.52, cex = 0.95, expression(italic("Spartina")))
text(x = 65, y = 0.3, cex = 0.95, expression(italic("Distichlis")))
# dev.off()




# Figure 2 - Aboveground biomass over time ----------------------------------
pd <- position_dodge(1.2)

ggplot(ddHgt4[ddHgt4$cohort > 0, ], aes(y = mass / pot.m2, x = as.Date(session), colour = species, shape = species)) +
  geom_point(size = pointSize, position = pd) + theme_classic() + 
  geom_errorbar(aes(ymin = (mass - mass.se) / pot.m2, ymax = (mass + mass.se) / pot.m2), width = 0, position = pd) +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = c(unique(ddHgt4$species[ddHgt4$cohort > 0])[1], unique(ddHgt4$species[ddHgt4$cohort > 0])[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(ddHgt4$species[ddHgt4$cohort > 0])[1], unique(ddHgt4$species[ddHgt4$cohort > 0])[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  ylim(0, 600) + facet_grid(. ~ cohort, labeller = label_parsed) +
  labs(y = expression("Biomass (g "%.%m^-2~")"), x = "") +
  scale_x_date(breaks = as.Date(unique(ddHgt4$session))[c(1, 2, 4, 6, 8)], labels = date_format("%b-%d")) +
  theme(legend.position = c(0.125, 0.87), legend.text.align = 0,
        legend.background = element_rect(fill = NA, colour = NA), axis.text.x=element_text(angle=45,hjust=1))
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/Figure2_", todaysDate, ".png"), width = 120, height= 70, units = "mm", dpi = 400)

unique(ddHgt4$day[ddHgt4$cohort == 1])

# Figure 3A - NAPP ------------------------------------------
# set dodge width for points & error bars
pd <- position_dodge(0.6)
pointSize2 <- 2

fig3A <- ggplot(dd.master, aes(y = prodn_rate, x = as.Date(session), shape = species, colour = species)) + 
  geom_point(size = pointSize2, position = pd) + theme_classic() +
  # facet_wrap(~ species) +
  geom_errorbar(aes(ymin = prodn_rate - prodn_rate.se, ymax = prodn_rate + prodn_rate.se), width = 0, position = pd) +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = c(unique(dd.master$species)[1], unique(dd.master$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(dd.master$species)[1], unique(dd.master$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  ylim(0, 10.2) +
  ylab(expression("Primary production (g"%.%m^{-2}%.%"d"^{-1}~")")) + xlab("") +
  scale_x_date(breaks = as.Date(unique(dd.master$session)), labels = date_format("%b-%d")) +
  theme(legend.text.align = 0, legend.position = c(5,5),
        legend.background = element_rect(fill = NA, colour = NA)) +
  annotate("text", x = as.Date(unique(dd.master$session)[3]), y = 10, label = "*") +
  annotate("text", x = as.Date(unique(dd.master$session)[4]), y = 8, label = "*") +
  annotate("text", x = as.Date(unique(dd.master$session)[1]), y = 10.2, label = "A")
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/Figure3A_", todaysDate, ".png"), width = 80, height= 80, units = "mm", dpi = 400)



# Figure 3B - Aboveground N uptake from production  ----------------------------------------

fig3B <- ggplot(dd.master, aes(y = n_uptake_biomass, x = as.Date(session), shape = species, colour = species)) + geom_point(size = pointSize2, position = pd) + theme_classic() +
  # facet_wrap(~ species) +
  geom_errorbar(aes(ymin = n_uptake_biomass - n_uptake_biomass.se, ymax = n_uptake_biomass + n_uptake_biomass.se), width = 0, position = pd) +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = c(unique(dd.master$species)[1], unique(dd.master$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(dd.master$species)[1], unique(dd.master$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  ylim(0, 150) +
  ylab(expression("N uptake (mg N"%.%d^{-1}%.%"m"^{-2}*")")) + xlab("") +
  scale_x_date(breaks = as.Date(unique(dd.master$session)), labels = date_format("%b-%d")) +
  theme(legend.text.align = 0, legend.position = c(5,5),
        legend.background = element_rect(fill = NA, colour = NA)) +
  annotate("text", x = as.Date(unique(dd.master$session)[3]), y = 147, label = "*") +
  annotate("text", x = as.Date(unique(dd.master$session)[4]), y = 125, label = "*") +
  annotate("text", x = as.Date(unique(dd.master$session)[1]), y = 150, label = "B")
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/Figure3B_.png"), width = 80, height= 80, units = "mm", dpi = 400)

ggsave(paste0("C:/RDATA/greenhouse/output/GRO/Figure3_", todaysDate,".png"), width = 80, height= 160, units = "mm", dpi = 400, 
       gridExtra::arrangeGrob(fig3A, fig3B))


# Figure 4 - 15N recoveries ------------------------------------------------
pd2 <- position_dodge(1.2)

ggplot(mgd[!mgd$time %in% "t0", ], aes(x = as.Date(session), y = recovery, colour = species, shape = species)) + 
  geom_point(size = pointSize2, position = pd2) + ylab(expression(" "^15~"N recovery")) + xlab("") +
  geom_errorbar(aes(ymin = recovery - recovery.se, ymax = recovery + recovery.se), width = 0, position = pd2) + theme_classic()  +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = unique(mgd$species), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = unique(mgd$species), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  facet_grid(. ~ type)  +
  scale_y_continuous(labels = scales::percent, lim = c(0, 0.85)) + 
  scale_x_date(breaks = as.Date(unique(mgd[!mgd$time %in% "t0", "session"])), labels = date_format("%b-%d")) +
  theme(legend.position = c(0.15, 0.9), legend.text.align = 0, axis.title.x=element_blank(),
        legend.background = element_rect(fill = NA, colour = NA), axis.text.x=element_text(angle=45,hjust=1))
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/Figure4_", todaysDate, ".png"), width = 120, height= 90, units = "mm", dpi = 400)



# Figure 5 - 15N uptake per gram fine roots ---------------------------------
rt2 <- ddply(master, .(species, timeDiff), summarise,
             xbar = mean(n15_pgFR, na.rm = TRUE),
             se = se(n15_pgFR),
             session = mean(session))
ggplot(rt2, aes(x = timeDiff, y = xbar, colour = species, shape = species)) + geom_point(size = pointSize2) + theme_classic() +
  ylim(0, 0.15) +
  geom_errorbar(aes(ymin = xbar - se, ymax = xbar + se), width = 0) +
  xlab("Experiment duration (days)") +
  ylab(expression(" "^15~"N uptake (mg "^{15}~N%.%d^{-1}%.%"gdw"^{-1}~")")) +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = unique(rt2$species), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(rt2$species)), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_x_continuous(breaks = unique(rt2$timeDiff), lim = c(7, 28)) +
  theme(legend.text.align = 0, legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = NA, colour = NA)) +
  geom_smooth(method = "lm", se = FALSE)
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/Figure5_", todaysDate, ".png"), width = 80, height= 80, units = "mm", dpi = 400)




# Figure 6 - N uptake estimated by 2 methods --------------------------------
ggplot(master[master$cluster %in% "3-4 weeks", ], aes(x = n_uptake_15n, y = n_uptake_biomass, colour = species, shape = species)) + geom_point(size = pointSize2) + theme_classic() +
  ylim(0, 150)  + # facet_wrap(~ cluster, scales = "free_x") +
  ylab(expression("N uptake (mg N "%.%d^{-1}%.%"m"^{-2}*")")) +
  xlab(expression(""^15*"N uptake (mg "^15*"N"%.%d^{-1}%.%"m"^{-2}*")")) +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = c(unique(master$species)[1], unique(master$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(master$species)[1], unique(master$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  theme(legend.text.align = 0, legend.position = c(0.2, 0.9),
        legend.background = element_rect(fill = NA, colour = NA)) +
  geom_smooth(data = subset(master, cluster %in% "3-4 weeks"), aes(group = 1), method = "lm", se = FALSE, colour = "black") +
  geom_text(data = data.frame(cluster = "3-4 weeks", n_uptake_biomass = 50, n_uptake_15n = 3, species = "SA"), label =
              "paste(italic(R) ^ 2, \" = 0.45\")", parse = TRUE, colour = "black")
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/Figure6_", todaysDate, ".png"), width = 80, height= 80, units = "mm", dpi = 400)
summary(lm3_4)

# Figure S1 - predicted vs obs biomass ---------------------------------------
ggplot(obs, aes(x = allom.est / pot.m2, y = g / pot.m2, colour = species2, shape = species2)) +
  geom_point(size = pointSize) + theme_classic() %+replace% theme(legend.title = element_blank()) +
  labs(x = expression("Predicted biomass (allometry; g"%.%m^-2~")"), y = expression("Measured biomass (g "%.%m^-2~")")) +
  ylim(0, 650) + xlim(0, 650) +
  scale_colour_grey(start = 0.5, end = 0.1, name = "", breaks = c(unique(obs$species2)[1], unique(obs$species2)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(obs$species2)[1], unique(obs$species2)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  theme(legend.position = c(0.3, 0.9), legend.text.align = 0,
        legend.background = element_rect(fill = NA, colour = NA)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/FigureS1_", todaysDate, ".png"), width = 80, height= 80, units = "mm", dpi = 400)




# Figure S2 - aboveground pools ---------------------------------------
ggplot(ag2.sp, aes(y = g, x = as.Date(session), colour = species2, shape = species2)) + geom_point(size = pointSize2, position = pd2) + theme_classic() %+replace% theme(legend.title = element_blank()) +
  facet_grid(sample.type2 ~ .) + geom_errorbar(aes(ymin = g - g.se, ymax = g + g.se), width = 0, position = pd2) +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = c(unique(ag2.sp$species2)[1], unique(ag2.sp$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(ag2.sp$species2)[1], unique(ag2.sp$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  ylab(expression("Biomass (g "%.%m^-2~")")) +  xlab("") + ylim(0, 325) +
  scale_x_date(breaks = as.Date(unique(ag2.sp$session)), labels = date_format("%b-%d")) +
  theme(legend.position = c(2.2, 0.9), legend.text.align = 0,
        legend.background = element_rect(fill = NA, colour = NA))
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/FigureS2_part1_", todaysDate, ".png"), width = 75, height= 90, units = "mm", dpi = 400)


ggplot(ag2.sp, aes(y = n_pct, x = as.Date(session), colour = species2, shape = species2)) + geom_point(size = pointSize2, position = pd2) + theme_classic() %+replace% theme(legend.title = element_blank()) +
  facet_grid(sample.type2 ~ .) + geom_errorbar(aes(ymin = n_pct - n_pct.se, ymax = n_pct + n_pct.se), width = 0, position = pd2) +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = c(unique(ag2.sp$species2)[1], unique(ag2.sp$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(ag2.sp$species2)[1], unique(ag2.sp$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  ylab("Tissue N") +  xlab("") + scale_y_continuous(labels = scales::percent) +
  scale_x_date(breaks = as.Date(unique(ag2.sp$session)), labels = date_format("%b-%d")) +
  theme(legend.position = c(2.2, 0.9), legend.text.align = 0,
        legend.background = element_rect(fill = NA, colour = NA))
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/FigureS2_part2_", todaysDate, ".png"), width = 75, height= 90, units = "mm", dpi = 400)



# Figure S3 - belowground biomass pools --------------------------------------
ggplot(bg.sp, aes(y = g, x = as.Date(session), colour = species2, shape = species2)) + geom_point(size = pointSize2, position = pd2) + theme_classic() +
  facet_grid(sample.type ~ .) + geom_errorbar(aes(ymin = g - g.se, ymax = g + g.se), width = 0, position = pd2) +
  scale_colour_grey(start = grayColor, end = 0.1, name = "", breaks = c(unique(bg.sp$species2)[1], unique(bg.sp$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  scale_shape_manual(values = c(16, 17), name = "", breaks = c(unique(bg.sp$species2)[1], unique(bg.sp$species)[2]), labels = c(expression(italic(Distichlis)), expression(italic(Spartina)))) +
  ylab(expression("Biomass (g "%.%m^-2~")")) +  xlab("") +
  scale_x_date(breaks = as.Date(unique(bg.sp$session)), labels = date_format("%b-%d")) +
  theme(legend.text.align = 0, legend.position = c(2.2, 1),
        legend.background = element_rect(fill = NA, colour = NA), axis.text.x=element_text(angle=45,hjust=1))
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/FigureS3_part1_", todaysDate, ".png"), width = 75, height= 90, units = "mm", dpi = 400)


ggplot(bg.sp, aes(y = n_pct, x = as.Date(session), colour = species2, shape = species2)) + geom_point(size = pointSize2, position = pd2) + theme_classic() +
  facet_grid(sample.type ~ .) + geom_errorbar(aes(ymin = n_pct - n_pct.se, ymax = n_pct + n_pct.se), width = 0, position = pd2) +
  scale_colour_grey(start = 0.5, end = 0.1, name = "", breaks = c(unique(bg.sp$species2)[1], unique(bg.sp$species)[2]), labels = c(expression(italic(D.~spicata)), expression(italic(S.~alterniflora)))) +
  scale_x_date(breaks = as.Date(unique(bg.sp$session)), labels = date_format("%b-%d")) +
  ylab("Tissue N content") +  xlab("") + scale_y_continuous(labels = scales::percent) +
  theme(legend.text.align = 0, legend.position = c(2.2, 1),
        legend.background = element_rect(fill = NA, colour = NA), axis.text.x=element_text(angle=45,hjust=1))
# ggsave(paste0("C:/RDATA/greenhouse/output/GRO/FigureS3_part2_", todaysDate, ".png"), width = 75, height= 90, units = "mm", dpi = 400)




