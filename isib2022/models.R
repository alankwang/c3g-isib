c3g <- readRDS('data/c3g_imputed_0713.rds')

# clean up the data
c3g <- c3g[GFR>=0, ] # remove observations with no GFR/negative GFR or GFR >120
c3g <- c3g[GFR_12m6w >= 0, ] # remove observations with negative GFR_12m6w
c3g <- c3g[(c3g$Etiology == 'C3G: DDD' | c3g$Etiology == 'C3G: C3GN' | c3g$Etiology == 'Immune Complex GN'), ]
c3g <- c3g[!is.na(c3g$soluble_level_12m6w), ]


# GFR_12m6w
fit.GFR <- gam(GFR ~ s(GFR_12m6w), data = c3g)
visreg(fit.GFR, 'GFR_12m6w', xlab='GFR 12months ago', ylab='GFR')
abline(coef=c(0,1), lwd = 1.5, lty = 2)
summary(fit.GFR) #R-sq = 0.84


# logUPCR_12m6w
fit <- gam(GFR ~ s(GFR_12m6w) + s(logUPCR_12m6w), data=c3g)
visreg(fit, 'logUPCR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit)

fit <- gam(GFR ~ s(logUPCR_12m6w, GFR_12m6w), data=c3g)
visreg(fit, 'logUPCR_12m6w', by='GFR_12m6w', gg=TRUE, overlay=TRUE, breaks=seq(25,100,25))
summary(fit) # R-sq = 0.824


# YrsDisease
fit <- gam(GFR ~ s(GFR_12m6w) + s(YrsDisease), data=c3g)
visreg(fit, 'GFR_12m6w', by='YrsDisease', gg=TRUE, overlay=TRUE, breaks=seq(1,10,1))
summary(fit)



# trying models to find the best one
fit <- gam(GFR ~ s(GFR_12m6w) + 
             s(logUPCR_12m6w) + 
             s(UACR_12m6w) + 
             s(soluble_level_12m6w) +
             s(YrsDisease),
             data=c3g[soluble_level_12m6w<20])
visreg(fit, data=c3g[soluble_level_12m6w<20])
summary(fit)



fit <- gam(GFR ~ s(GFR_12m6w, logUPCR_12m6w) + 
             s(UACR_12m6w) + 
             s(soluble_level_12m6w) +
             s(YrsDisease),
           data=c3g[soluble_level_12m6w<20])
visreg(fit, data=c3g[soluble_level_12m6w<20])
summary(fit)



# models with eGFR
# all variables
fit <- gam(eGFR2 ~ s(eGFR2_12m6w) +
             s(logUPCR_12m6w) +
             s(UACR_12m6w) +
             s(soluble_level_12m6w) +
             s(YrsDisease), 
             data=c3g)
visreg(fit)
summary(fit)
# logUPCR_12m6w not significant here (p-value=0.109)


# everything but YrsDisease
fit <- gam(eGFR2 ~ s(eGFR2_12m6w) +
             s(logUPCR_12m6w) +
             s(UACR_12m6w) +
             s(soluble_level_12m6w),
             data=c3g)
visreg(fit)
summary(fit) 
# now logUPCR_12m6w p-value=0.22


# everything but UACR_12m6w
fit <- gam(eGFR2 ~ s(eGFR2_12m6w) +
             s(logUPCR_12m6w) +
             s(soluble_level_12m6w) +
             s(YrsDisease),
             data=c3g)
visreg(fit)
summary(fit) 
# logUPCR_12m6w p-value=0.00159


# everything but soluble_level_12m6w
fit <- gam(eGFR2 ~ s(eGFR2_12m6w) +
             s(logUPCR_12m6w) +
             s(UACR_12m6w) +
             s(YrsDisease),
             data=c3g)
summary(fit)
# logUPCR_12m6w p-value=0.0845
# UACR_12m6w p-value=0.12


# everything but logUPCR_12m6w
fit <- gam(eGFR2 ~ s(eGFR2_12m6w) + 
             s(UACR_12m6w) +
             s(soluble_level_12m6w) +
             s(YrsDisease),
             data=c3g)
visreg(fit)
summary(fit)
# UACR_12m6w p-value=7.85e-05
#### USE THIS MODEL ####


# UACR and UPCR are correlated
cor.test(c3g$UACR_12m6w, c3g$UPCR_12m6w) # p-value=0.2391
#so we don't need both UACR and UPCR



# eGFR uncapped
fit <- gam(eGFR2 ~ s(eGFR2_12m6w) + 
             s(UACR_12m6w) +
             s(soluble_level_12m6w) +
             s(YrsDisease),
           data=c3g[soluble_level_12m6w <=20])
visreg(fit, data=c3g[soluble_level_12m6w <=20])
summary(fit)


# GFR uncapped 
D <- c3g[GFR>=0 & GFR_12m6w >= 0 & !is.na(soluble_level_12m6w) & soluble_level_12m6w < 20]
fit <- gam(GFR ~ s(GFR_12m6w) +
             s(UACR_12m6w) +
             s(soluble_level_12m6w),
           data=D)
visreg(fit, data=D)
summary(fit)


# GFR capped
c3g.cap <- c3g[GFR <=120, ]
fit <- gam(GFR ~ s(GFR_12m6w) +
             s(logUPCR_12m6w) +
             s(soluble_level_12m6w) +
             s(YrsDisease), 
           data=c3g.cap[soluble_level_12m6w <= 20])
visreg(fit, data=c3g.cap[soluble_level_12m6w <=20])
summary(fit)

# eGFR2 capped
c3g[, egfr_cap := pmin(eGFR2, 120)]
c3g[, egfr_cap_12m6w := pmin(eGFR2_12m6w, 120)]
fit <- gam(egfr_cap ~ s(egfr_cap_12m6w) +
             s(logUPCR_12m6w) +
             s(log(soluble_level_12m6w)) +
             s(YrsDisease), 
           data=c3g[soluble_level_12m6w <= 20])
visreg(fit, data=c3g[soluble_level_12m6w <=20])
summary(fit)
