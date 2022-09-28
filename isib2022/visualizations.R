## Slides 
c3g <- readRDS('data/c3g_imputed_0713.rds')
c3g <- c3g[(c3g$Etiology == 'C3G: DDD' | c3g$Etiology == 'C3G: C3GN' | c3g$Etiology == 'Immune Complex GN'), ]
c3g <- c3g[!is.na(c3g$soluble_level_12m6w), ]
c3g[, soluble_level := pmin(soluble_level_12m6w, 3)]
c3g <- c3g[!is.na(c3g$eGFR2), ]
c3g <- c3g[!is.na(c3g$eGFR2_12m6w), ]
c3g[ , eGFR2_cap := pmin(eGFR2, 120)]
c3g[ , eGFR2_12m6w_cap := pmin(eGFR2_12m6w, 120)]




# eGFR2(eGFR2_12m6w)
plot(c3g$eGFR2_12m6w_cap, c3g$eGFR2_cap, pch=20,col='darkgrey',
     main='Baseline Prediction', 
     xlab='eGFR', 
     ylab='eGFR prediction (12 months from now)')
abline(coef=c(0,1), lty=1, lwd=2, col='dodgerblue')


fit <- gam(eGFR2_cap ~ eGFR2_12m6w_cap, data=c3g)
visreg(fit, xlab='eGFR', ylab='eGFR prediction (12 months from now)',
       main='baseline eGFR prediction')
summary(fit)
abline(coef=c(0,1), lty=2, lwd=1.5)

predict(fit, newdata=c3g[66])

lm_fit <- lm(eGFR2_cap ~ eGFR2_12m6w_cap, data=c3g)
visreg(lm_fit)
summary(lm_fit)


# simple model (only based on eGFR_12m6w_cap and )
lm_fit <- lm(eGFR2_cap ~ 0 + eGFR2_12m6w_cap:I(YrsDisease<4), data=c3g)
plot(c3g$YrsDisease, c3g$eGFR2_cap, xlim=c(0,10), pch=20, col='darkgrey',
     xlab='Years with Disease', ylab='predicted eGFR', main='Years Prediction')
# DF <- data.frame(eGFR2_12m6w_cap = 75, YrsDisease = seq(0, 10, len=100))
# lines(DF$YrsDisease, predict(lm_fit, newdata=DF), col='dodgerblue')
lines(predict(lm_fit, c3g), col='dodgerblue', lwd=2)
abline(v=4, lty=2, lwd=1.5)
summary(lm_fit)

visreg(lm_fit, 'eGFR2_12m6w_cap', by='YrsDisease', gg=TRUE, overlay=TRUE, breaks=c(0, 4)) + geom_abline(intercept=0, slope=1, col='darkgrey', lty=2)
visreg(lm_fit, 'YrsDisease', cond=list(eGFR2_12m6w_cap=80), xlim=c(0,10), ylim=c(60,100))
abline(h=80, lty=2, col='darkgrey', lwd=1.5)

predict(lm_fit, c3g[10,])
c3g$eGFR2_12m6w_cap[10]
c3g$YrsDisease[10]

c3g$YrsDisease[403]
c3g$eGFR2_12m6w_cap[403]
predict(lm_fit, c3g[403, ])


# UACR
fit <- lm(eGFR2_cap ~ UACR_12m6w, data=c3g)
fit <- gam(eGFR2_cap ~ s(UACR_12m6w), data=c3g)
visreg(fit, ylim=c(0, 120))



# Model 2 (years with disease)
c3g[, DiseaseTime := cut(YrsDisease, breaks=c(0, 4, Inf))]
gam(eGFR2_cap ~ 0 + eGFR2_12m6w_cap:DiseaseTime, data=c3g) %>% summary
gam(eGFR2_cap ~ 0 + eGFR2_12m6w_cap + I(eGFR2_12m6w_cap * (DiseaseTime=='(0,4]')), data=c3g) %>% summary
fit <- gam(eGFR2_cap ~ 0 + eGFR2_12m6w_cap:DiseaseTime, data=c3g)

visreg(fit, 'eGFR2_12m6w_cap', by='DiseaseTime', gg=TRUE, overlay=TRUE,
       xlab='eGFR', ylab='predicted eGFR') +
  geom_abline(intercept=0, slope=1, lty=2, col='gray60', lwd=0.75)

visreg(fit, 'DiseaseTime', xlab='Years with disease',
       ylab='predicted eGFR', ylim=c(60,100), 
       cond=list(eGFR2_12m6w_cap=80))
abline(h=80, lty=2, col='darkgrey', lwd=1.5)




# Our Model

fit <- gam(eGFR2_cap ~ eGFR2_12m6w_cap +
             logUPCR_12m6w +
             s(soluble_level, sp=0.07) + 
             s(YrsDisease), 
           data=c3g)

visreg(fit, ylab='predicted eGFR')
summary(fit)

visreg(fit, 'eGFR2_12m6w_cap',
       xlab='eGFR',
       ylab='predicted eGFR')
abline(coef=c(0,1), lty=2, lwd=1.5, col='darkgrey')

visreg(fit, 'logUPCR_12m6w', cond=list(eGFR2_12m6w_cap=80), 
       xlab='logUPCR',
       ylab='predicted eGFR',
       ylim=c(50,100),
       xlim=c(-3,3))
abline(h=80, lty=2, lwd=1.5)

visreg(fit, 'soluble_level', cond=list(eGFR2_12m6w_cap=80),
       xlab='soluble C5b-9 level',
       ylab='predicted eGFR',
       ylim=c(50, 100))
abline(h=80, lty=2, lwd=1.5)
abline(v=0.3, lty=2, lwd=1.5)

visreg(fit, 'YrsDisease', cond=list(eGFR2_12m6w_cap=80), 
       xlab='Years with disease',
       ylab='predicted eGFR',
       xlim=c(0,10),
       ylim=c(50,100))
abline(h=80, lty=2, lwd=1.5)
abline(v=4, lty=2, lwd=1.5)


# prediction residuals
c3g$eGFR_pred <- predict(fit, newdata=c3g)

old_diff <- sqrt(mean((c3g$eGFR2_cap - c3g$eGFR2_12m6w_cap)^2, na.rm=TRUE))
old_diff
new_diff <- sqrt(mean((c3g$eGFR2_cap - c3g$eGFR_pred)^2, na.rm=TRUE))
new_diff



## Split data
library(dplyr)
set.seed(123)
train <- c3g[sample(1:nrow(c3g), round(nrow(c3g)*.8)),]
test <- anti_join(c3g, train)

# prediction residuals
test$eGFR_pred <- predict(fit, newdata=test)

(old_rMSE <- sqrt(mean((test$eGFR2_cap - test$eGFR2_12m6w_cap)^2, na.rm=TRUE)))
(new_rMSE <- sqrt(mean((test$eGFR2_cap - test$eGFR_pred)^2, na.rm=TRUE)))

# try multiple training/testing splits
nsim <- 1000
improve <- data.frame(old_rMSE=rep(0, nsim), 
                      new_rMSE=rep(0,nsim), 
                      improve=rep(0,nsim))
for (i in 1:nsim) {
  train <- c3g[sample(1:nrow(c3g), round(nrow(c3g)*.8)),]
  test <- anti_join(c3g, train)
  test$eGFR_pred <- predict(fit, newdata=test)
  improve$old_rMSE[i] <- sqrt(mean((test$eGFR2_cap - test$eGFR2_12m6w_cap)^2, na.rm=TRUE))
  improve$new_rMSE[i] <- sqrt(mean((test$eGFR2_cap - test$eGFR_pred)^2, na.rm=TRUE))
  improve$improve[i] <- ifelse(improve$old_rMSE[i] > improve$new_rMSE[i], 1, 0)
}
mean(improve$improve)
