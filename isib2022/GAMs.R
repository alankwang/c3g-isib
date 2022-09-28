# GAMs - NewData
devtools::load_all("../isib2022")
library("data.table")
library("ggplot2")
library("magrittr")
library("mgcv")
library("rstan")
library("splines")
library("visreg")
c3g <- readRDS('data/c3g_imputed_0711.rds')

# clean up the data
c3g <- c3g[GFR>=0 ,] # remove observations with no GFR/negative GFR
c3g <- c3g[GFR_12m6w >= 0, ] # remove observations with negative GFR_12m6w
c3g <- c3g[(c3g$Etiology == 'C3G: DDD' | c3g$Etiology == 'C3G: C3GN' | c3g$Etiology == 'Immune Complex GN'), ]
c3g <- c3g[!is.na(c3g$soluble_level_12m6w), ]

# GFR as a function of GFR_12m6w
fit.GFR <- gam(GFR ~ s(GFR_12m6w), data = c3g)
visreg(fit.GFR, 'GFR_12m6w')
abline(coef=c(0,1), lwd = 1.5)
summary(fit.GFR) #R-sq = 0.83


# see if any of the other variables are significant
# UPCR_12m6w
fit <- gam(GFR ~ s(GFR_12m6w) + s(UPCR_12m6w), data = c3g)
visreg(fit, 'UPCR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit)

fit <- gam(GFR ~ s(UPCR_12m6w, GFR_12m6w), data=c3g)
visreg(fit, 'UPCR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit)


# logUPCR_12m6w
fit <- gam(GFR ~ s(GFR_12m6w) + s(logUPCR_12m6w), data=c3g)
visreg(fit, 'logUPCR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25), ylim=c(0,120))
summary(fit)

fit <- gam(GFR ~ s(logUPCR_12m6w, GFR_12m6w), data=c3g)
visreg(fit, 'logUPCR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit) # R-sq = 0.824
# looks like they change at different rates


# C3G_12m6w
fit <- gam(GFR ~ s(GFR_12m6w) + s(C3_12m6w), data = c3g)
visreg(fit, 'C3_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit)

fit <- gam(GFR ~ s(C3_12m6w, GFR_12m6w), data=c3g)
visreg(fit, 'C3_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit)
# no additional effect


# UACR_12m6w
fit <- gam(GFR ~ s(GFR_12m6w) + s(UACR_12m6w), data = c3g)
visreg(fit, 'UACR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit)

fit <- gam(GFR ~ s(UACR_12m6w, GFR_12m6w), data=c3g)
visreg(fit, 'UACR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit)


# soluble_level_12m6w
c3g <- c3g[c3g$soluble_level_12m6w<=20, ]
fit <- gam(GFR ~ s(soluble_level_12m6w) + s(GFR_12m6w), data=c3g[soluble_level_12m6w<=20])
visreg(fit, 'soluble_level_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, 
       breaks=seq(25,100,25), xlim=c(0, 20), data=c3g[soluble_level_12m6w<=20])
summary(fit) # 0.829

fit <- gam(GFR ~ s(soluble_level_12m6w, GFR_12m6w), data=c3g)
visreg(fit, 'soluble_level_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, 
       breaks=seq(25,100,25), xlim=c(0, 20))
summary(fit)
# no additional effect

# Etiology
c3g$Etiology <- unclass(as.factor(c3g$Etiology))
fit <- gam(GFR ~ s(GFR_12m6w) + Etiology, data=c3g)
fit <- gam(GFR ~ s(GFR_12m6w, Etiology), data=c3g)
visreg(fit, 'GFR_12m6w', by='Etiology', gg=TRUE, overlay=TRUE)
summary(fit)
# no difference


# AgeBiopsy
fit <- gam(GFR ~ s(GFR_12m6w) + s(AgeBiopsy), data=c3g)
visreg(fit, 'GFR_12m6w', by='AgeBiopsy', gg=TRUE, overlay=TRUE, breaks=seq(5,40,5))
summary(fit)
# no additional effect


# YrsDisease
fit <- gam(GFR ~ s(GFR_12m6w) + s(YrsDisease), data=c3g)
visreg(fit, 'GFR_12m6w', by='YrsDisease', gg=TRUE, overlay=TRUE, breaks=seq(1,10,1))
summary(fit)
# no additional effect


# RaceEthn
#c3g$RaceEthn <- unclass(c3g$RaceEthn)
fit <- gam(GFR ~ s(GFR_12m6w) + RaceEthn, data=c3g)
visreg(fit, 'GFR_12m6w', by='RaceEthn', gg=TRUE, overlay=TRUE)
summary(fit)
# no additional effect


# logUPCR_12m6w and UACR_12m6w
fit <- gam(GFR ~ s(GFR_12m6w) + s(logUPCR_12m6w) + s(UACR_12m6w), data=c3g)
visreg(fit)
summary(fit) # R-sq = 0.798

fit <- gam(GFR ~ s(GFR_12m6w, logUPCR_12m6w) + s(UACR_12m6w), data=c3g)
visreg(fit, ylim=c(40, 125))
# now logUPCR doesn't seem to make a difference



# different smooth functions
fit <- gam(GFR ~ s(GFR_12m6w) + s(UACR_12m6w), data=c3g)
visreg(fit, 'UACR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))

fit <- gam(GFR ~ te(GFR_12m6w,UACR_12m6w), data=c3g)
visreg(fit, 'UACR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))

fit <- gam(GFR ~ ti(GFR_12m6w, UACR_12m6w), data=c3g)
visreg(fit, 'UACR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
# this one looks flat but the rest look like it makes a difference

fit <- gam(GFR ~ t2(GFR_12m6w, UACR_12m6w), data=c3g)
visreg(fit, 'UACR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))


fit <- gam(GFR ~ t2(GFR_12m6w, logUPCR_12m6w), data=c3g)
visreg(fit, 'logUPCR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
# doesn't look like it makes a difference with the other ones
