##############################################################
# 
#         SENSITIVITY ANALYSIS FOR FINAL REGRESSION MODELS
#                     Code by: Andrea Lund
#                 Last updated: January 22, 2020
#

# packages
library(lme4)
#library(nnet) # for multinomial/polytamous logistic regression
#library(mlogit) # for multinomial logistic regression with random effects
#library(MCMCglmm) # for multinomial logistic regresion with random effects without data manipulation
#library(ordinal) # for ordinal logistic regression with random effects
library(multgee) # for multinomial logistic regression using ordLORgee()

# set working directory
setwd("C:/Users/andre/Box Sync/Schisto (alund2@stanford.edu)/Data")

# import merged data
allData <- read.csv("allData.csv")
allData$headSchoolYrsCat <- factor(allData$headSchoolYrsCat)
allData$nWivesCat <- factor(allData$nWivesCat)
allData$wealthF <- factor(allData$wealthF)

shData <- read.csv("shData.csv")

smData <- read.csv("smData.csv")

allData$Sh_low <- as.factor(ifelse(allData$Sh_cat == 1, 1,
                         ifelse(allData$Sh_cat == 0, 0, "NA")))
allData$Sh_hi <- as.factor(ifelse(allData$Sh_cat == 2, 1,
                        ifelse(allData$Sh_cat == 0, 0, "NA")))

# functions
ci_calc_mixed <- function(mod) {
  se <- sqrt(diag(vcov(mod)))
  estCI <- cbind(OR = fixef(mod), 
                 LL = fixef(mod)-1.96*se, 
                 UL = fixef(mod)+1.96*se)
  orCI <- exp(estCI)
  return(orCI)
}

###
#
#           SECONDARY OUTCOMES
#       SCHISTO INTENSITY CATEGORIES
#


###### USING MULTINOMIAL GEE 
##### S. HAEMATOBIUM
sh_multgee <- ordLORgee(formula = Sh_cat ~ irrigArea + 
                          locF + 
                          DEM01F + 
                          DEM04c +
                          headSchoolYrsCat + 
                          nWivesCat + 
                          pumpYN + 
                          wealthF + 
                          irrigAreaV,
                          # factor(village), 
                     data = allData, id = C_ConcessionId,
                     link = "logit", LORstr = "category.exch")

summary(sh_multgee)

# generate ORs and confidence intervals
# calculate standard errors
se <- diag(sh_multgee[["robust.variance"]]); se

# low intensity infection compared to no infection
lowOR <- exp(coefficients(sh_multgee)[1] + coefficients(sh_multgee)[3]); lowOR
lowOR_LL <- exp((coefficients(sh_multgee)[1] - 1.96*se[1]) + 
                 coefficients(sh_multgee)[3] - 1.96*se[3]); lowOR_LL
lowOR_UL <- exp((coefficients(sh_multgee)[1] + 1.96*se[1]) + 
                 coefficients(sh_multgee)[3] + 1.96*se[3]); lowOR_UL
low_ORCI <- cbind(lowOR, lowOR_LL, lowOR_UL); low_ORCI

# high intensity infection compared to no infection
hiOR <- exp(coefficients(sh_multgee)[2] + coefficients(sh_multgee)[3]); hiOR
hiOR_LL <- exp((coefficients(sh_multgee)[2] - 1.96*se[2]) + 
                 coefficients(sh_multgee)[3] - 1.96*se[3]); hiOR_LL
hiOR_UL <- exp((coefficients(sh_multgee)[2] + 1.96*se[2]) + 
                 coefficients(sh_multgee)[3] + 1.96*se[3]); hiOR_UL
hi_ORCI <- cbind(hiOR, hiOR_LL, hiOR_UL); hiOR

multgee_ORCI <- rbind(low_ORCI, hi_ORCI)
colnames(multgee_ORCI) <- c("OR", "LL", "UL"); multgee_ORCI
rownames(multgee_ORCI) <- c("Low vs no infection", "High vs no infection"); multgee_ORCI
multgee_ORCI <- as.data.frame(multgee_ORCI)
multgee_ORCI$species <- "sh"
multgee_ORCI$model <- "multgee"
multgee_ORCI$level <- c("low", "high")


##### S. MANSONI
allData$Sm_catR <- ifelse(allData$Sm_cat == 2, 1,
                          ifelse(allData$Sm_cat %in% c(3:4), 2, 0))
n <- dim(allData)[1]
temp_rows <- sample(1:n,size=n,replace=TRUE)
model_temp <-  ordLORgee(formula = Sm_catR ~ irrigArea + 
                             locF + 
                             DEM01F + 
                             DEM04c + 
                             headSchoolYrsCat + 
                             nWivesCat + 
                             pumpYN +  
                             wealthF + 
                             irrigAreaV,
                             #factor(village), 
                           data = allData[temp_rows,],
                           id = C_ConcessionId, repeated = NULL,
                           link = "logit", LORstr = "category.exch", 
                           bstart = c(1.8, 2.9, -0.25, -0.25, 0,
                                      -16, -0.25, -0.25, 0, -0.2,
                                      -0.35, -0.5, -0.2, -0.3, 0, 0))


sm_multgee <- ordLORgee(formula = Sm_catR ~ irrigArea + 
                          #locF + # doesnt run by itself
                          #DEM01F + # doesnt run by itself
                          #DEM04c + # doesnt run by itself
                          #headSchoolYrsCat + # robust cov matrix not positive definite
                          #nWivesCat + # doesnt run by itself
                          #pumpYN +  
                          #wealthF + 
                          #irrigAreaV,
                          factor(village), 
                        data = allData, id = C_ConcessionId,
                        link = "logit", LORstr = "category.exch",
                        bstart = c(1.8, 2.9, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0, 
                                   0, 0, 0))

# logit starting points
initShp <- c(0, 0, 0, -0.25, -0.25, 0,
          -16, -0.25, -0.25, 0, -0.2,
          -0.35, -0.5, -0.2, -0.3, 0)
initBiv <- c(1.8, 2.9, 0, 0, 0, 0,
             0, 0 , 0, 0,
             0.1, -0.3, 0, 0)
initBiv <- c(1.8, 2.9, 0.01, -0.07, -0.03, 0.01,
             -0.2, -0.07, -0.06, -0.01, 0.02,
             -0.01, -0.03, 0.05, 0, 0)


summary(sm_multgee)


###### USING SEPARATE LOGITS
# low intensity infection vs no infection
shLow <- dplyr::filter(allData, Sh_low != "NA")
shLow_mod <- glmer(Sh_hi ~ irrigArea + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = shLow, family = binomial(link = "logit"))
# warning: max grad = 0.0055

# high intensity infection vs no infection
shHi <- dplyr::filter(allData, Sh_hi != "NA")
shHi_mod <- glmer(Sh_hi ~ irrigArea + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = shHi, family = binomial(link = "logit"))
# max|grad = 0.0362319
summary(shHi_mod)

# generate point estimate and confidence interval
shLow_ORCI <- ci_calc(shLow_mod)[2,]; shLow_ORCI
shHi_ORCI <- ci_calc(shHi_mod)[2,]; shHi_ORCI
shCat_ORCI <- rbind(shLow_ORCI, shHi_ORCI)
rownames(shCat_ORCI) <- c("low", "high"); shCat_ORCI
shCat_ORCI <- as.data.frame(shCat_ORCI)
shCat_ORCI$species <- "sh"
shCat_ORCI$model <- "seplogit"


##### S. MANSONI
table(allData$Sm_cat) 
# 0 -> 1014; 2 -> 131; 3 -> 54; 4 -> 23 
# ^ these categories ^ aren't quite right; recode to 0, 1, 2, 3
# and combine medium and high categories to improve chances of convergence
allData$Sm_cat <- ifelse(allData$Sm_cat == 2, 1, 
                         ifelse(allData$Sm_cat == 3, 2,
                                ifelse(allData$Sm_cat == 4, 3, allData$Sm_cat)))
allData$Sm_catR <- ifelse(allData$Sm_cat %in% c(2,3), 2, allData$Sm_cat)

# create variable and data set for low vs no infection comparison
allData$Sm_low <- as.factor(ifelse(allData$Sm_cat %in% 2:3, "NA",
                                   ifelse(allData$Sm_cat == 0, 0, 1)))
table(allData$Sm_cat, allData$Sm_low)
smLow <- dplyr::filter(allData, Sm_low != "NA") # n = 1145

# create variable and data set for med/hi vs no infection comparison
allData$Sm_medhi <- as.factor(ifelse(allData$Sm_cat %in% c(2, 3), 1,
                                   ifelse(allData$Sm_cat == 0, 0, "NA")))
table(allData$Sm_cat, allData$Sm_medhi)
smMedhi <- dplyr::filter(allData, Sm_medhi != "NA") # n = 1068

# fit separate logit models
# low intensity vs no infection
smLow_mod <- glmer(Sm_low ~ irrigArea + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = smLow, family = binomial(link = "logit"))
# unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue
summary(smLow_mod)

smLow_ORCI <- ci_calc_mixed(smLow_mod)[2,]

# medium/high intensity vs no infection
smMedhi_mod <- glmer(Sm_medhi ~ irrigArea + locF + DEM01F + DEM04c +
                     factor(headSchoolYrsCat) + factor(nWivesCat) + 
                     pumpYN + wealthF + irrigAreaV +
                     (1|village/C_ConcessionId), 
                   data = smMedhi, family = binomial(link = "logit"))
# unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue
summary(smMedhi_mod)
smMedhi_ORCI <- ci_calc_mixed(smMedhi_mode)[2,]

# point estimates and confidence intervals for separate sm models
smLow_ORCI <- ci_calc(smLow_mod)[2,]; smLow_ORCI
smMedhi_ORCI <- ci_calc(smMedhi_mod)[2,]; smMedhi_ORCI
smCat_ORCI <- rbind(smLow_ORCI, smMedhi_ORCI); smCat_ORCI
rownames(smCat_ORCI) <- c("low", "medium/high"); smCat_ORCI
smCat_ORCI <- as.data.frame(smCat_ORCI)
smCat_ORCI$species <- "sm"
smCat_ORCI$model <- "seplogit"

seplogit_ORCI <- rbind(shCat_ORCI, smCat_ORCI); seplogit_ORCI
seplogit_ORCI$level = c("low", "high", "low", "medium/high"); seplogit_ORCI

# multinomial model with fixed effects for village
# Sh
ShCat <- multinom(Sh_cat ~ irrigArea + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV + factor(village), 
                  data = allData)

# generate point estimates and confidence intervals
shcat_coef <- summary(ShCat)$coefficients
shcat_coef <- shcat_coef[,2]
shcat_se <- summary(ShCat)$standard.errors
shcat_se <- shcat_se[,2]
shcat_or <- exp(shcat_coef)
shcat_ll <- exp(shcat_coef - 1.96*shcat_se)
shcat_ul <- exp(shcat_coef + 1.96*shcat_se)
shcat_orci <- cbind(shcat_or, shcat_ll, shcat_ul); shcat_orci
shcat_orci <- as.data.frame(shcat_orci)
colnames(shcat_orci) <- c("OR", "LL", "UL")
shcat_orci$species <- "sh"
shcat_orci$model <- "fixed effects"
shcat_orci <- as.data.frame(shcat_orci)
shcat_orci$level <- c("low", "high")


# Sm
SmCat <- multinom(Sm_cat ~ irrigArea + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV + factor(village), 
                  data = allData)
summary(SmCat)

# generate point estimates and confidence intervals
smcat_coef <- summary(SmCat)$coefficients
smcat_coef <- smcat_coef[,2]
smcat_se <- summary(SmCat)$standard.errors
smcat_se <- smcat_se[,2]
smcat_or <- exp(smcat_coef)
smcat_ll <- exp(smcat_coef - 1.96*smcat_se)
smcat_ul <- exp(smcat_coef + 1.96*smcat_se)
smcat_orci <- cbind(smcat_or, smcat_ll, smcat_ul); smcat_orci
colnames(smcat_orci) <- c("OR", "LL", "UL")
smcat_orci <- as.data.frame(smcat_orci)
smcat_orci$species <- "sm"
smcat_orci$model <- "fixed effects"
smcat_orci$level <- c("low", "medium", "high")

multinom_ORCI <- rbind(shcat_orci, smcat_orci); multinom_ORCI


# combine estimates from all multinomial models
allMultinom_ORCI <- rbind(multgee_ORCI, seplogit_ORCI, multinom_ORCI)
write.csv(allMultinom_ORCI, "C:/Users/andre/Box Sync/Papers/Ch2IRisk-IrrigatedAg/Analysis/multinom_output.csv")


# ZERO-INFLATED MODELS
library(pscl)
# S. HAEMATOBIUM
# poisson
shc_zip <- zeroinfl(round(Sh_median, 0) ~ irrigArea + locF + DEM01F + DEM04c +
                       factor(headSchoolYrsCat) + factor(nWivesCat) + 
                       pumpYN + wealthF + irrigAreaV + factor(village),
                     data = allData, dist = "poisson")
summary(shc_zip)

shc_zip_IRRCI_all <- cbind(exp(coef(shc_zip)), exp(confint(shc_zip)))
shc_zip_IRRCI <- shc_zip_IRRCI_all[2,]

# negative binomial
shc_zinb <- zeroinfl(round(Sh_median, 0) ~ irrigArea + locF + DEM01F + DEM04c +
                       factor(headSchoolYrsCat) + factor(nWivesCat) + 
                       pumpYN + wealthF + irrigAreaV + factor(village),
                     data = shData, dist = "negbin")

summary(shc_zinb)

# generate point estimates and confidence intervals
shc_zinb_IRRCI_all <- cbind(exp(coef(shc_zinb)), exp(confint(shc_zinb)))
shc_zinb_IRRCI <- shc_zinb_IRRCI_all[2,]

# Vuong test to compare ZIP and ZINB models (Hilbe 2011, p 378)
#   V gt 1.96 indicates ZIP is preferred; V lt -1.96 indicates ZINB is preferred; 
#   values of V between 1.96 and -1.96 indicate neither model is preferred over the other

pp <- predict(shc_zip)
pnb <- predict(shc_zinb)
u <- log(pp/pnb)
umean <- mean(u)
nobs <- dim(shData)[1]
stdu <- sqrt(var(u))
v <- (sqrt(nobs)*umean)/stdu
dnorm(v) # 0.006056562

# LR test to compare ZINB and ZIP models
LLp <-logLik(shc_zip)
LLnb <- logLik(shc_zinb)
LRtest <- -2*(LLp -- LLnb)
dchisq(1, LRtest)/2
# retain point estimates and confidence intervals for best-fitting model


smc_zinb <- zeroinfl(round(Sm_median, 0) ~ irrigArea + locF + DEM01F + DEM04c +
                       factor(headSchoolYrsCat) + factor(nWivesCat) + 
                       pumpYN + wealthF + irrigAreaV ,
                     data = allData, dist = "negbin")

summary(smc_zinb) # NaNs in sqrt(daig(vcov(mod))), for headSchoolYrsCat1
smc_zinb_IRRCI_all <- cbind(exp(coef(smc_zinb)), exp(confint(smc_zinb)))
smc_zinb_IRRCI <- smc_zinb_IRRCI_all[c(2, 17),]
###
#
#           SECONDARY OUTCOMES
#          COINFECTION PRESENCE
#
table(allData$coinf) # 1 --> n = 167, 0 --> n = 1060
coinf <- glmer(coinf ~ irrigArea + locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) + 
                        pumpYN + wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))
# singular fit
summary(coinf)

###
#
#           SECONDARY OUTCOMES
#          NUMBER OF INFECTIONS
#
table(allData$nInf)
nInf <- multinom(nInf ~ irrigArea + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV, 
                  data = allData)

ninf_coef <- summary(nInf)$coefficients
ninf_coef <- ninf_coef[,2]
ninf_se <- summary(nInf)$standard.errors
ninf_se <- ninf_se[,2]
ninf_or <- exp(ninf_coef)
ninf_ll <- exp(ninf_coef - 1.96*ninf_se)[,2]
ninf_ul <- exp(ninf_coef + 1.96*ninf_se)[,2]
ninf_orci <- cbind(ninf_or, ninf_ll, ninf_ul); ninf_orci


###
#
#           SECONDARY EXPOSURE
#       IRRIGATED AREA CATEGORIES
#

##### S. HAEMATOBIUM PRESENCE
# 
table(allData$irrigAreaCat)
table(allData$irrigAreaCat, allData$locF)

shp.irrigcat <- glmer(Sh_presence ~ factor(irrigAreaCat) + locF + DEM01F + DEM04c +
                     factor(headSchoolYrsCat) + factor(nWivesCat) + 
                     pumpYN + wealthF + irrigAreaV +
                     (1|village/C_ConcessionId), 
                   data = allData, family = binomial(link = "logit"))
# max|grad = 0.0203
summary(shp.irrigcat)

shp.irrigcat_loc <- glmer(Sh_presence ~ factor(irrigAreaCat)*locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) + 
                        pumpYN + wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))

anova(shp.irrigcat, shp.irrigcat_loc, test = "Chisq") # Chisq = 5.1, df = 3, p = 0.16
# no signficant interaction

shp.irrigcat_pump <- glmer(Sh_presence ~ factor(irrigAreaCat)*pumpYN + locF + DEM01F + DEM04c +
                            factor(headSchoolYrsCat) + factor(nWivesCat) + 
                            wealthF + irrigAreaV +
                            (1|village/C_ConcessionId), 
                          data = allData, family = binomial(link = "logit"))

anova(shp.irrigcat, shp.irrigcat_pump, test = "Chisq") # Chisq = 6.3179, df = 3, p = 0.09713
# no significant interaction

# point estimates + confidence intervals of interest
shp.irrigcat_ci <- ci_calc(shp.irrigcat)
shp.irrigcat_ci <- shp.irrigcat_ci[2:4,]

# model diagnostics
binned_residuals(shp.irrigcat) # 94% of residuals inside error bounds
plot_model(shp.irrigcat, type = "resid")

#####  S. MANSONI PRESENCE
#
smp.irrigcat <- glmer(Sm_presence ~ factor(irrigAreaCat) + locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) + 
                        pumpYN + wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))

smp.irrigcat_loc <- glmer(Sm_presence ~ factor(irrigAreaCat)*locF + DEM01F + DEM04c +
                     factor(headSchoolYrsCat) + factor(nWivesCat) + 
                     pumpYN + wealthF + irrigAreaV +
                     (1|village/C_ConcessionId), 
                   data = allData, family = binomial(link = "logit"))
# large eigenvalue ratio

anova(smp.irrigcat, smp.irrigcat_loc, test = "Chisq") # Chisq = 1.82, df = 3, p = 0.61
# no interaction

smp.irrigcat_pump <- glmer(Sm_presence ~ factor(irrigAreaCat)*pumpYN + locF + DEM01F + DEM04c +
                            factor(headSchoolYrsCat) + factor(nWivesCat) + 
                            wealthF + irrigAreaV +
                            (1|village/C_ConcessionId), 
                          data = allData, family = binomial(link = "logit"))

anova(smp.irrigcat, smp.irrigcat_pump, test = "Chisq") # Chisq = 4.8819, df = 3, p = 0.1807
# no interaction 

# point estimates + confidence intervals of interest
smp.irrigcat_ci <- ci_calc(smp.irrigcat)
smp.irrigcat_ci <- smp.irrigcat_ci[2:4,]

# model diagnostics
binned_residuals(smp.irrigcat) # 89% of residuals inside error bounds
plot_model(smp.irrigcat, type = "resid")

###
#
#           SECONDARY EXPOSURE
#       BINARY (Y/N) IRRIGATED AREA 
#

# create binary variable for irrigated area presence 
allData$irrigAreaYN <- ifelse(allData$irrigArea > 0, 1, 0)
table(allData$irrigAreaYN, allData$irrigArea)
table(allData$irrigAreaYN)
table(allData$irrigAreaYN, allData$locF)

##### 
#     S. HAEMATOBIUM PRESENCE
# 

shp.irrigyn <- glmer(Sh_presence ~ irrigAreaYN + locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) + 
                        pumpYN + wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))

shp.irrigyn_loc <- glmer(Sh_presence ~ irrigAreaYN*locF + DEM01F + DEM04c +
                       factor(headSchoolYrsCat) + factor(nWivesCat) + 
                       pumpYN + wealthF + irrigAreaV +
                       (1|village/C_ConcessionId), 
                     data = allData, family = binomial(link = "logit"))

anova(shp.irrigyn, shp.irrigyn_loc, test = "Chisq") # Chisq = 2.8, df = 1, p = 0.09
# no interaction

shp.irrigyn_pump <- glmer(Sh_presence ~ irrigAreaYN*pumpYN + locF + DEM01F + DEM04c +
                            factor(headSchoolYrsCat) + factor(nWivesCat) + 
                            wealthF + irrigAreaV +
                            (1|village/C_ConcessionId), 
                          data = allData, family = binomial(link = "logit"))
# max grad warning 0.00117799

anova(shp.irrigyn, shp.irrigyn_pump, test = "Chisq") # Chisq = 1.3278, df = 1, p = 0.2492
# no interaction

ci_calc(shp.irrigyn)[2,] # irrigAreaYN OR = 1.04, 95% CI (0.73, 1.46)

binned_residuals(shp.irrigyn) # 89%
plot_models(shp.irrigyn, type = "resid")

##### 
#     S. MANSONI PRESENCE
# 

smp.irrigyn <- glmer(Sm_presence ~ irrigAreaYN + locF + DEM01F + DEM04c +
                       factor(headSchoolYrsCat) + factor(nWivesCat) + 
                       pumpYN + wealthF + irrigAreaV +
                       (1|village/C_ConcessionId), 
                     data = allData, family = binomial(link = "logit"))

smp.irrigyn_loc <- glmer(Sm_presence ~ irrigAreaYN*locF + DEM01F + DEM04c +
                       factor(headSchoolYrsCat) + factor(nWivesCat) + 
                       pumpYN + wealthF + irrigAreaV +
                       (1|village/C_ConcessionId), 
                     data = allData, family = binomial(link = "logit"))
# max grad = 0.00731978; large eigenvalue ratio

anova(smp.irrigyn, smp.irrigyn_loc, test = "Chisq") # chisq = 0.27, df = 1, p = 0.60

smp.irrigyn_pump <- glmer(Sm_presence ~ irrigAreaYN*pumpYN + locF + DEM01F + DEM04c +
                           factor(headSchoolYrsCat) + factor(nWivesCat) + 
                           wealthF + irrigAreaV +
                           (1|village/C_ConcessionId), 
                         data = allData, family = binomial(link = "logit"))
# unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue

anova(smp.irrigyn, smp.irrigyn_pump, test = "Chisq") # Chisq = 4.3338, df = 1, p = 0.03736

summary(smp.irrigyn_pump)

ci_calc(smp.irrigyn)[2,]
#        OR        LL        UL 
# 0.9678376 0.6297084 1.4875291 

ci_calc_int_mixed(smp.irrigyn_pump)

binned_residuals(smp.irrigyn) # 80%
plot_models(smp.irrigyn, type = "resid")

###
#
#           SECONDARY EXPOSURE
#       AREA OF MARKET GARDEN PLOTS
#

table(allData$gardenArea)
dplyr::summarize(allData, mean = mean(gardenArea), sd = sd(gardenArea))

allData %>% 
  group_by(locF) %>% 
  summarize(mean = mean(gardenArea), sd = sd(gardenArea))

###### S. haematobium
shp.garden <- glmer(Sh_presence ~ gardenArea + locF + DEM01F + DEM04c +
                      factor(headSchoolYrsCat) + factor(nWivesCat) + 
                      pumpYN + wealthF + irrigAreaV +
                      (1|village/C_ConcessionId), 
                    data = allData, family = binomial(link = "logit"))

shp.garden_loc <- glmer(Sh_presence ~ gardenArea*locF + DEM01F + DEM04c +
                      factor(headSchoolYrsCat) + factor(nWivesCat) + 
                      pumpYN + wealthF + irrigAreaV +
                      (1|village/C_ConcessionId), 
                    data = allData, family = binomial(link = "logit"))
# max grad = 0.00109824

anova(shp.garden, shp.garden_loc, test = "Chisq") # Chisq = 0.0768, df = 1, p = 0.7817
# no interaction

shp.garden_pump <- glmer(Sh_presence ~ gardenArea*pumpYN + locF + DEM01F + DEM04c +
                          factor(headSchoolYrsCat) + factor(nWivesCat) + 
                          wealthF + irrigAreaV +
                          (1|village/C_ConcessionId), 
                        data = allData, family = binomial(link = "logit"))

anova(shp.garden, shp.garden_pump, test = "Chisq") # Chisq = 0.0157, df = 1, p = 0.9004

ci_calc(shp.garden)[2,]
#        OR        LL        UL 
# 1.1118795 0.9068996 1.3631897 

###### S. mansoni

smp.garden <- glmer(Sm_presence ~ gardenArea + locF + DEM01F + DEM04c +
                          factor(headSchoolYrsCat) + factor(nWivesCat) + 
                          pumpYN + wealthF + irrigAreaV +
                          (1|village/C_ConcessionId), 
                        data = allData, family = binomial(link = "logit"))
# unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue

smp.garden_loc <- glmer(Sm_presence ~ gardenArea*locF + DEM01F + DEM04c +
                      factor(headSchoolYrsCat) + factor(nWivesCat) + 
                      pumpYN + wealthF + irrigAreaV +
                      (1|village/C_ConcessionId), 
                    data = allData, family = binomial(link = "logit"))
# large eigenvalue ratio

anova(smp.garden, smp.garden_loc, test = "Chisq") # Chisq = 0.6417,  df = 1, p = 0.4231
# no loc interaction

smp.garden_pump <- glmer(Sm_presence ~ gardenArea*pumpYN + locF + DEM01F + DEM04c +
                          factor(headSchoolYrsCat) + factor(nWivesCat) + 
                          wealthF + irrigAreaV +
                          (1|village/C_ConcessionId), 
                        data = allData, family = binomial(link = "logit"))
# unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue

anova(smp.garden, smp.garden_pump, test = "Chisq") # Chisq = 5.6429, df = 1, p = 0.01753

ci_calc_int_mixed(smp.garden_pump)
#           strata     estR       llR       ulR
#gardenArea      1 1.317707 0.8934236 1.9434802
#gardenArea      0 0.809466 0.6806814 0.9626166

ci_calc(smp.garden)[2,]
#        OR        LL        UL 
# 1.1212285 0.9363429 1.3426207 


###
#
#           SECONDARY EXPOSURE
#          AREA OF RICE FIELDS
#

table(allData$riz, allData$locF)
dplyr::summarize(allData, mean = mean(riz), sd = sd(riz))

allData %>% 
  group_by(locF) %>% 
  summarize(mean = mean(riz), sd = sd(riz))

# S. haematobium
shp.rice <- glmer(Sh_presence ~ riz + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = allData, family = binomial(link = "logit"))
# max grad = 0.00766617

riverData <- filter(allData, locF == "river")
shp.rice_r <- glmer(Sh_presence ~ riz + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = riverData, family = binomial(link = "logit"))
# max grad = 0.00121786

ci_calc(shp.rice_r)
shp.rice_loc <- glmer(Sh_presence ~ riz*locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) + 
                        pumpYN + wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))
# max grad = 0.00922713

anova(shp.rice, shp.rice_loc, test = "Chisq") # Chisq = 1.2238, df = 1, p = 0.2686
# no loc interaction

shp.rice_pump <- glmer(Sh_presence ~ riz*pumpYN + locF + DEM01F + DEM04c +
                         factor(headSchoolYrsCat) + factor(nWivesCat) + 
                         wealthF + irrigAreaV +
                         (1|village/C_ConcessionId), 
                       data = allData, family = binomial(link = "logit"))

anova(shp.rice, shp.rice_pump, test = "Chisq") # Chisq = 0, df = 1, p = 0.996
# no pump interaction

ci_calc(shp.rice)[2,]
#        OR        LL        UL 
# 1.0739802 0.9480313 1.2166618 

# S. mansoni
smp.rice <- glmer(Sm_presence ~ riz + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = allData, family = binomial(link = "logit"))
# large eigenvalue ratio

smp.rice_loc <- glmer(Sm_presence ~ riz*locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = allData, family = binomial(link = "logit"))
# max grad = 0.016; large eigenvalue ratio

anova(smp.rice, smp.rice_loc, test = "Chisq") # Chisq = 1.1972, df = 1, p = 0.2739

smp.rice_pump <- glmer(Sm_presence ~ riz*pumpYN + locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) + 
                        wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))
# unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue

anova(smp.rice, smp.rice_pump, test = "Chisq") # Chisq = 5.3539, df = 1, p = 0.02068

ci_calc_int_mixed(smp.rice_pump)
#    strata         estR       llR     ulR
#riz      1 3.062909e-11 0.0000000     Inf <-- weird numbers result of fit warnings?
#riz      0 9.871017e-01 0.8964177 1.08696


ci_calc(smp.rice)[2,]
# OR        LL        UL 
# 0.9599319 0.7927017 1.1624413 

###
#
#           SECONDARY EXPOSURE
#       AREA OF IRRIGATED MONOCROPS 
#

dplyr::summarize(allData, mean = mean(irrigmonoArea), sd = sd(irrigmonoArea))

allData %>% 
  group_by(locF) %>% 
  summarize(mean = mean(irrigmonoArea), sd = sd(irrigmonoArea))


# S. haematobium
shp.mono <- glmer(Sh_presence ~ irrigmonoArea + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = allData, family = binomial(link = "logit"))
# unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalues

shp.mono_loc <- glmer(Sh_presence ~ irrigmonoArea*locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = allData, family = binomial(link = "logit"))
# max grad = 0.0175813

anova(shp.mono, shp.mono_loc, test = "Chisq") # Chisq = 1.0696, df = 1, p = 0.301
# no location interaction

shp.mono_pump <- glmer(Sh_presence ~ irrigmonoArea*pumpYN + locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) + 
                        wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))
# unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue

anova(shp.mono, shp.mono_pump, test = "Chisq") # Chisq = 0.5476, df = 1, p = 0.4593

ci_calc(shp.mono)[2,]
#       OR       LL       UL 
# 1.137880 1.007616 1.284983 

# S. mansoni
smp.mono <- glmer(Sm_presence ~ irrigmonoArea + locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = allData, family = binomial(link = "logit"))
# max grad = 0.00142217, large eigenvalue ratio

smp.mono_loc <- glmer(Sm_presence ~ irrigmonoArea*locF + DEM01F + DEM04c +
                    factor(headSchoolYrsCat) + factor(nWivesCat) + 
                    pumpYN + wealthF + irrigAreaV +
                    (1|village/C_ConcessionId), 
                  data = allData, family = binomial(link = "logit"))
# max grad = 0.0026222, large eigenvalue ratio

anova(smp.mono, smp.mono_loc, test = "Chisq") # Chisq = 0.2473, df = 1, p = 0.619
# no location interaciton

smp.mono_pump <- glmer(Sm_presence ~ irrigmonoArea*pumpYN + locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) + 
                        wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))
# max grad = 0.00155025, large eigenvalue ratio

anova(smp.mono, smp.mono_pump, test = "Chisq") # Chisq = 0.058, df = 1, p = 0.8096
# no pump interaction

ci_calc(smp.mono)[2,]
#        OR        LL        UL 
# 0.9803055 0.8662416 1.1093891 

