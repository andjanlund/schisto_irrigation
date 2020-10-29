###############################################################
# 
#    regression analysis of irrigated agriculture as
#    a risk factor for schistosomiasis reinfection
#    coded by: Andrea Lund
#    last updated: October 29, 2019
#
###############################################################

# install packages
library(stats) # for a variety of statistical functions
library(lme4) # for mixed effects regression
library(car) # for vif/collinearity assessment
library(stargazer) # for publication-quality regression output tables
library(sjPlot) # custom correlation matrices 
library(msm) # multi-state markov in continuous time required for additive interactions
library(ggplot2) # for aesthetically-pleasing plots
library(MASS) # for negative binomial regression 
library(performance) # for generating binned residuals, other regression diagnostics
library(lmtest) # for lrtest()
library(visreg) # for visualizing regression output
library(broom)
library(dotwhisker)
library(lattice) # for dotplots of random effects
library(DescTools) # for Hosmer-Lemeshow goodness-of-fit test
library(ResourceSelection) # for Hosmer-Lemeshow test
library(MKmisc) # for Hosmer-Lemeshow test
library(ggpubr) # for combining plots 

# set working directory
setwd("C:/Users/andre/Box Sync/Schisto (alund2@stanford.edu)/Data")

# import merged data
allData <- read.csv("allData.csv")
allData$headSchoolYrsCat <- factor(allData$headSchoolYrsCat)
allData$nWivesCat <- factor(allData$nWivesCat)
allData$wealthF <- factor(allData$wealthF)
shData <- read.csv("shData.csv")
smData <- read.csv("smData.csv")


# subset data by location
shRiver <- shData[which(shData$locF == "river"),]
shLake <- shData[which(shData$locF == "lake"),]
smRiver <- smData[which(smData$locF == "river"),]
smLake <- smData[which(smData$locF == "lake"),]

# add vRice variable to river data sets
shRiver$vRice <- ifelse(shRiver$village %in% c("LR", "MD", "NM"), 1, 0)
smRiver$vRice <- ifelse(smRiver$village %in% c("LR", "MD", "NM"), 1, 0)

# designate functions
# source additive interaction script 
source("C:/Users/andre/Box Sync/Papers/Ch2IRisk-IrrigatedAg/Analysis/Additive interactions function.R")

####  FUNCTIONS   ####
# explained deviance
explainedD <- function(mod) {
  (((mod$null.deviance - mod$deviance)/mod$null.deviance) * 100)
}

# dispersion parameter, based on Zuur et al 2009 p 224
dispersion <- function(mod) {
  n = nrow(mod$data)
  p = length(mod$coefficients)
  phi = mod$deviance/(n-p)
  phi
}

# function to generate variables needed for confounding assessment
confound_output <- function(mod, omit) {
  
  coefE <- coef(mod)
  coefE <- coefE[2]
  ciE <- confint(mod)
  ciE <- ciE[2,]
  expE <- exp(coef(mod))
  expE <- expE[2]
  expCI <- exp(confint(mod))
  expCI <- expCI[2,]
  
  out <- c(omit, as.numeric(coefE), as.numeric(ciE), as.numeric(expE), as.numeric(expCI))
  out
}


# assess presence of separation for fitted models
sep <- function(mod) {
  vals <- predict(mod, mod$data, type = 'response')
  min <- min(vals) 
  max <- max(vals)
  return(c(min, max))
}

## calculate confidence intervals for fitted models
# https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
ci_calc <- function(mod) {
  se <- sqrt(diag(vcov(mod)))
  estCI <- cbind(est = fixef(mod), 
                 LL = fixef(mod)-1.96*se, 
                 UL = fixef(mod)+1.96*se)
  orCI <- exp(estCI)
  return(orCI)
  }


#############################################################
#
#                   PRESENCE/ABSENCE OUTCOMES

############################################################
#
#                   SCHISTOSOMA HAEMATOBIUM

## bivariate linear model
shp1<- glm(Sh_presence ~ irrigArea, data = allData, family = binomial(link = "logit"))
summary(shp1)
shp_cOR <- cbind(est = exp(coef(shp1)),
                 ci = exp(confint(shp1)))
shp_cOR <- shp_cOR[2,] # irrigarea only

## TEST INTERACTION TERMS ON FULL FIXED EFFECTS MODEL
shp2 <- glm(Sh_presence ~ irrigArea + pumpYN + locF + DEM01F + DEM04c +
              headSchoolYrsCat + nWivesCat + eth_mode3 + nFishAnyC + 
              kidsAnyTask + nWatptsC + nWatPtsAgKidsC + 
              wealthF +
              distanceWatPtC + marketDistance_kmC + irrigAreaV , 
              data = allData, family = binomial(link = "logit")) 


shp_aOR <- cbind(est = exp(coef(shp2)),
                 ci = exp(confint(shp2)))
shp_aOR <- shp_aOR[2,]

extractAIC(shp2) 
summary(shp2) 
exp(coef(shp2))
exp(confint(shp2))
vif(shp2)

## ADD/TEST RANDOM EFFECTS
# household random intercepts
shp3.h <- glmer(Sh_presence ~ irrigArea + locF + DEM01F + DEM04c +
                headSchoolYrsCat + 
                nWivesCat + eth_mode3 + 
                nFishAnyC +
                nWatptsC + nWatPtsAgKidsC + kidsAnyTask + 
                pumpYN + 
                wealthF +
                distanceWatPtC + 
                marketDistance_kmC + irrigAreaV +
                (1|C_ConcessionId),
              data = allData, family = binomial(link = "logit")) 

anova(shp3.h, shp2, test = "Chisq") # Chisq = 0.1364, df = 1, p = 0.7118

# village random intercepts
shp3.v <- glmer(Sh_presence ~ irrigArea + locF + DEM01F + DEM04c +
                  headSchoolYrsCat + 
                  nWivesCat + eth_mode3 + 
                  nFishAnyC +
                  nWatptsC + nWatPtsAgKidsC + kidsAnyTask + 
                  pumpYN + 
                  wealthF +
                  distanceWatPtC + 
                  marketDistance_kmC + irrigAreaV +
                  (1|village),
                data = allData, family = binomial(link = "logit")) 

# warnings: max|grad| = 0.00260283; very large eigenvalue - rescale variables?
anova(shp3.v, shp2, test = "Chisq") # Chisq = 18.556, df = 1, p =  1.65e-05 

# nested random intercepts
shp3.hv <- glmer(Sh_presence ~ irrigArea + locF + DEM01F + DEM04c +
                   headSchoolYrsCat + 
                   nWivesCat + eth_mode3 + 
                   nFishAnyC +
                   nWatptsC + nWatPtsAgKidsC + 
                   kidsAnyTask + 
                   pumpYN + 
                   wealthF +
                   distanceWatPtC + 
                   marketDistance_kmC + irrigAreaV +
                  (1|village/C_ConcessionId), 
                data = allData, family = binomial(link = "logit"))

anova(shp3.hv, shp3.h, test = "Chisq") # Chisq = 18.419, df = 1, p = 1.773e-05

shp3.hv_loc <- glmer(Sh_presence ~ irrigArea*locF + DEM01F + DEM04c +
                       headSchoolYrsCat + 
                       nWivesCat + eth_mode3 + 
                       nFishAnyC +
                       nWatptsC + nWatPtsAgKidsC + 
                       kidsAnyTask + 
                       pumpYN + 
                       wealthF +
                       distanceWatPtC + 
                       marketDistance_kmC + irrigAreaV +
                       (1|village/C_ConcessionId), 
                     data = allData, family = binomial(link = "logit"))

anova(shp3.hv, shp3.hv_loc, test = "Chisq") # Chisq = 0.2367, df = 1, p = 0.6266

shp3.hv_pump <- glmer(Sh_presence ~ irrigArea*pumpYN + locF + DEM01F + DEM04c +
                       headSchoolYrsCat + 
                       nWivesCat + eth_mode3 + 
                       nFishAnyC +
                       nWatptsC + nWatPtsAgKidsC + 
                       kidsAnyTask + 
                       wealthF +
                       distanceWatPtC + 
                       marketDistance_kmC + irrigAreaV +
                       (1|village/C_ConcessionId), 
                     data = allData, family = binomial(link = "logit"))

anova(shp3.hv, shp3.hv_pump, test = "Chisq") # Chisq = 0.2198, df = 1, p = 0.6392

# report non-stratified estimates
shp_mOR <- ci_calc(shp3.hv)[2,]

## MODEL FIT DIAGNOSTICS
# plot binned residuals (from Gelman & Hill p 97)
binned_residuals(shp3.hv) #91%

# plot residuals against numeric predictors
plot_model(shp3.hv, type = "resid" )
binned_residuals(model = shp3.hv, term = "irrigArea") # 86%
binned_residuals(model = shp3.hv, term = "DEM04c") # 100%
binned_residuals(model = shp3.hv, term = "nFishAnyC") # 100%
binned_residuals(model = shp3.hv, term = "nWatptsC") # 83%
binned_residuals(model = shp3.hv, term = "nWatPtsAgKidsC") # 67% - probably bad fit
binned_residuals(model = shp3.hv, term = "kidsAnyTask") # 100%
binned_residuals(model = shp3.hv, term = "pumpYN") # 100%
binned_residuals(model = shp3.hv, term = "distanceWatPtC") # 94%
binned_residuals(model = shp3.hv, term = "marketDistance_kmC") # 94%
binned_residuals(model = shp3.hv, term = "irrigAreaV") # 93%

# collinearity assessment
vif(shp3.hv) # highest value is 1.63 for eth_mode3

# influential observations
# calculate Cooks distance for model 
cooksd <- cooks.distance(shp3.hv)
# Plot the Cook's Distance 
sample_size <- nrow(allData)
plot(cooksd, pch="*", cex=1.5, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 2, col = "red") # add cutoff line
abline(h = 4/sample_size, col="red")  # add cutoff line - 4/n = 0.00325, seems small
#text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd >= 2, names(cooksd),""), col="red")  # add labels
length(which(cooksd >= 2)) # n = 0
length(which(cooksd >= 4/sample_size)) # n = 639 

# plot random intercepts
dotplot(ranef(shp3.hv, condVar = T))

# check normality of random intercepts
qqnorm(ranef(shp3.hv)$village[[1]], main = "Village Effects")
qqnorm(ranef(shp3.hv)$'village/C_ConcessionId'[[1]], main = "Household Effects") # error y is empty or has only NAs

# plot random effects against predictors
plot_model(shp3.hv, type = "pred", pred.type = "re")
plot_model(shp3.hv, type = "diag", grid) %>% plot_grid()

## GENERATE FINAL ESTIMATES
# generate confidence intervals 
#exp(confint.merMod(shp3.hv, method = "boot")) # this takes forever to run
shp_mOR <- ci_calc(shp3.hv)[2,] # but the "by-hand" calculation doesnt...

# export crude, adjusted and mixed effects estimates for plotting
shp_OR <- as.data.frame(rbind(shp_cOR, shp_aOR, shp_mOR))
colnames(shp_OR) <- c("OR", "LL", "UL")
shp_OR <- shp_OR %>% 
  dplyr::mutate(Estimate = c("crude", "adjusted", "mixed"))


############################################################
#
#                   SCHISTOSOMA MANSONI

## bivariate linear model
smp1 <- glm(Sm_presence ~ irrigArea, data = allData, family = binomial)

# crude odds ratio for presence of S mansoni
smp_cOR <- cbind(est = exp(coef(smp1)),
                 ci = exp(confint(smp1)))
smp_cOR <- smp_cOR[2,] # irrigarea only

# full model
smp2 <- glm(Sm_presence ~ irrigArea + pumpYN + locF + DEM01F + DEM04c +
              headSchoolYrsCat + nWivesCat + eth_mode3 + nFishAnyC + 
              kidsAnyTask + nWatptsC + nWatPtsAgKidsC + 
              wealthF +
              distanceWatPtC + marketDistance_kmC + irrigAreaV , 
            data = allData, family = binomial(link = "logit"))
summary(smp2) 

# rewrite code to calculate stratum-specific aOR from interaction term
smp_aOR <- cbind(exp(coef(smp2)), exp(confint(smp2)))
smp_aOR <- smp_aOR[2,]


# additive interaction - Excel file saved to Box > Papers > Ch2IrrigatedAg > Analysis > additive-interaction-smp2.1
d =  data.frame(allData$irrigArea, allData$loc, allData$Sm_presence)
additive_interactions(smp2.1, dat = d, recode = T) # preventive exposure - need to recode


## FIT RANDOM EFFECTS
# household-level intercepts
smp3.h <- glmer(Sm_presence ~ irrigArea + pumpYN + locF + DEM01F + DEM04c +
                   headSchoolYrsCat + nWivesCat + eth_mode3 + nFishAnyC + 
                   kidsAnyTask + nWatptsC + nWatPtsAgKidsC + 
                   wealthF +
                   distanceWatPtC + marketDistance_kmC + irrigAreaV + 
               (1|C_ConcessionId), 
             data = allData, family = binomial(link = "logit"))
# warnings: max|grad| = 0.246688, very large eigenvalue - rescale variables?
anova(smp3.h, smp2, test = "Chisq") # chisq = 2.1205, df = 1, p = 0.1453

# village-level intercepts
smp3.v <-  glmer(Sm_presence ~ irrigArea + pumpYN + locF + DEM01F + DEM04c +
                   headSchoolYrsCat + nWivesCat + eth_mode3 + nFishAnyC + 
                   kidsAnyTask + nWatptsC + nWatPtsAgKidsC + 
                   wealthF +
                   distanceWatPtC + marketDistance_kmC + irrigAreaV + 
                   (1|village), 
                 data = allData, family = binomial(link = "logit"))

anova(smp3.v, smp2, test = "Chisq") # chisq = 32.571, df = 1, p = 1.149e-08

# households nested in village-level intercepts
smp3.hv <-  glmer(Sm_presence ~ irrigArea + pumpYN + locF + DEM01F + DEM04c +
                    headSchoolYrsCat + nWivesCat + eth_mode3 + nFishAnyC + 
                    kidsAnyTask + nWatptsC + nWatPtsAgKidsC + 
                    wealthF +
                    distanceWatPtC + marketDistance_kmC + irrigAreaV + 
                   (1|village/C_ConcessionId), 
                 data = allData, family = binomial(link = "logit"))

isSingular(smp3.hv) # FALSE

anova(smp3.hv, smp3.h, test = "Chisq") # Chisq = 30.45, df = 1, p = 3.426e-08

smp3.hv_loc <-  glmer(Sm_presence ~ irrigArea*locF + pumpYN + DEM01F + DEM04c +
                    headSchoolYrsCat + nWivesCat + eth_mode3 + nFishAnyC + 
                    kidsAnyTask + nWatptsC + nWatPtsAgKidsC + 
                    wealthF +
                    distanceWatPtC + marketDistance_kmC + irrigAreaV + 
                    (1|village/C_ConcessionId), 
                  data = allData, family = binomial(link = "logit"))

anova(smp3.hv_loc, smp3.hv, test = "Chisq") # Chisq = 2.9742, df = 1, p = 0.0846

smp3.hv_pump <-  glmer(Sm_presence ~ irrigArea*pumpYN + locF + DEM01F + DEM04c +
                        headSchoolYrsCat + nWivesCat + eth_mode3 + nFishAnyC + 
                        kidsAnyTask + nWatptsC + nWatPtsAgKidsC + 
                        wealthF +
                        distanceWatPtC + marketDistance_kmC + irrigAreaV + 
                        (1|village/C_ConcessionId), 
                      data = allData, family = binomial(link = "logit"))

anova(smp3.hv_pump, smp3.hv, test = "Chisq") # Chisq = 4.0692, df = 1, p = 0.04367


## MODEL FIT DIAGNOSTICS
# plot binned residuals (from Gelman & Hill p 97)
binned_residuals(smp3.hv) # 83%

# plot random intercepts
dotplot(ranef(smp3.hv, condVar = T))

# check normality of random intercepts
qqnorm(ranef(smp3.hv)$village[[1]], main = "Village Effects")

# plot residuals against numeric predictors
plot_model(smp3.hv, type = "resid" )
binned_residuals(model = smp3.hv, term = "irrigArea") # 92%
binned_residuals(model = smp3.hv, term = "DEM04c") # 100%
binned_residuals(model = smp3.hv, term = "headSchoolYrsC") # 71% - probably bad fit
binned_residuals(model = smp3.hv, term = "nWivesC") # 80%
binned_residuals(model = smp3.hv, term = "nFishAnyC") # 100%
binned_residuals(model = smp3.hv, term = "nWatptsC") # 100%
binned_residuals(model = smp3.hv, term = "nWatPtsAgKidsC") # 67% - probably bad fit
binned_residuals(model = smp3.hv, term = "kidsAnyTask") # 100%
binned_residuals(model = smp3.hv, term = "pumpYN") # 100%
binned_residuals(model = smp3.hv, term = "wealthF") #100%
binned_residuals(model = smp3.hv, term = "distanceWatPtC") # 91%
binned_residuals(model = smp3.hv, term = "marketDistance_kmC") # 94%
binned_residuals(model = smp3.hv, term = "irrigAreaV") # 93%

# collinearity assessment
vif(smp3.hv) # highest values are for irrigArea (1.537845), locF (1.434780) and interact term (1.467746)

# influential observations
# calculate Cooks distance for model 
cooksd <- cooks.distance(smp3.hv)
# Plot the Cook's Distance 
sample_size <- nrow(allData)
plot(cooksd, pch="*", cex=1.5, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 2, col = "red") # add cutoff line
#abline(h = 4/sample_size, col="red")  # add cutoff line - 4/n = 0.00325, seems small
#text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd >= 2, names(cooksd),""), col="red")  # add labels
length(which(cooksd >= 2)) # n = 0 -- obs 1055, 543, 22, 23, 846
#length(which(cooksd >= 4/sample_size)) # n = 598 

# subset data to exclude 12 most influential observations
smDataR <- allData[which(cooksd < 2),]

# refit model with subsetted data
smp3.hv_R <-  glmer(Sm_presence ~ irrigArea + pumpYN + locF + DEM01F + DEM04c +
                      headSchoolYrsCat + nWivesCat + eth_mode3 + nFishAnyC + 
                      kidsAnyTask + nWatptsC + nWatPtsAgKidsC + 
                      wealthF +
                      distanceWatPtC + marketDistance_kmC + irrigAreaV +  
                    (1|village/C_ConcessionId), 
                  data = smDataR, family = binomial(link = "logit"))

# plot binned residuals
binned_residuals(smp3.hv_R) # 82% of residuals inside error bounds - probably bad fit

# residuals by predictor 
plot_model(smp3.hv_R, type = "resid")
binned_residuals(smp3.hv_R, term = "irrigArea") # 92%
binned_residuals(smp3.hv_R, term = "DEM01F") # 100%
binned_residuals(smp3.hv_R, term = "DEM04c") # 100%
binned_residuals(smp3.hv_R, term = "headSchoolYrsC") # 71%
binned_residuals(smp3.hv_R, term = "nWivesC") # 80%
binned_residuals(smp3.hv_R, term = "nFishAnyC") # 100%
binned_residuals(smp3.hv_R, term = "kidsAnyTask") # 100%
binned_residuals(smp3.hv_R, term = "nWatptsC")# 100%
binned_residuals(smp3.hv_R, term = "nWatPtsAgKidsC") # 67%
binned_residuals(smp3.hv_R, term = "pumpYN") # 100%
binned_residuals(smp3.hv_R, term = "wealthF") #100%
binned_residuals(smp3.hv_R, term = "distanceWatPtC") # 89%
binned_residuals(smp3.hv_R, term = "marketDistance_kmC") # 88%
binned_residuals(smp3.hv_R, term = "irrigAreaV") # 93%

# collinearity assessment
vif(smp3.hv_R)

# influential observations
# calculate Cooks distance for model 
cooksd <- cooks.distance(smp3.hv_R)
# Plot the Cook's Distance 
sample_size <- nrow(smData)
plot(cooksd, pch="*", cex=1.5, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 2, col = "red") # add cutoff line
#abline(h = 4/sample_size, col="red")  # add cutoff line - 4/n = 0.00325, seems small
#text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd >= 2, names(cooksd),""), col="red")  # add labels
length(which(cooksd >= 2)) # n = 3 
#length(which(cooksd >= 4/sample_size)) # n = 598 

# plot random intercepts
dotplot(ranef(smp3.hv_R, condVar = T))

# check normality of random intercepts
qqnorm(ranef(smp3.hv_R)$village[[1]], main = "Village Effects")
qqnorm(ranef(shp3.hv_R)$'village/C_ConcessionId'[[1]], main = "Household Effects") # error y is empty or has only NAs

# plot random effects against predictors
plot_model(smp3.hv_R, type = "pred", pred.type = "re") 
plot_model(smp3.hv_R, type = "diag", grid) %>% plot_grid()
plot_model(smp3.hv_R, type = "int")


# calculate confidence intervals
smp_mOR <- ci_calc(smp3.hv)[2,]

# full regression output for all sh and sm presence models 
stargazer(shp1, shp2, shp3.hv,
          smp1, smp2, smp3.hv,
          type="html",out="presence_regression_output_nonDAG.doc",
          intercept.bottom = F,
          intercept.top = T,
          digits=2)

################################################
# combine Sh_presence and Sm_pressence results into figures/tables
# 1. combine all shp model estimates 
shp_OR <- as.data.frame(rbind(shp_cOR, shp_aOR, shp_mOR))
colnames(shp_OR) <- c("OR", "LL", "UL")
shp_OR <- shp_OR %>% 
  dplyr::mutate(Estimate = c("crude", "adjusted", "mixed")) 
# 2. combine all smp estimates into single data frame
smp_OR <- rbind(smp_cOR, smp_aOR, smp_mOR); smp_OR
smp_OR <- as.data.frame(smp_OR) 
colnames(smp_OR) <- c("OR", "LL", "UL")
smp_OR$Estimate <- c("crude", "adjusted", "mixed"); smp_OR
smp_OR$Estimate <- factor(smp_OR$Estimate, levels = c("crude", "adjusted", "mixed"))
# 3. add column identifying parasite species
shp_OR <- dplyr::mutate(shp_OR, species = "Sh")
smp_OR <- dplyr::mutate(smp_OR, species = "Sm")
# 4. combine sh and sm presence estimates
presence_OR <- rbind(shp_OR, smp_OR)
# 5. make sure categories displayed in desired order
presence_OR$Estimate <- factor(presence_OR$Estimate, levels = c("crude", "adjusted", "mixed"))
presence_OR$species <- factor(presence_OR$species, levels = c("Sm", "Sh")) # puts Sh on top
presence_OR <- dplyr::mutate(presence_OR, outcome = "presence")

# plot presence results
presence_plot <- ggplot(presence_OR,
                        aes(x = species, y = OR, ymin = LL, ymax = UL)) +
  geom_pointrange(aes(lty = Estimate),
                  position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = LL, ymax = UL, lty = Estimate), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = 2) +
  xlab("Schistosome species") +
  ylab("Odds Ratio") +
  scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0.5, 1.5, 0.25)) +
  coord_flip() +
  #scale_colour_manual(values = c("grey41", "darkslategrey", "darkred")) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  guides(col = guide_legend(reverse = T)) +
  theme_bw() 

# prediction plots
# plot predicted probabilities of shp
shp.pred <- plot_model(shp3.hv, type = "pred", pred.type = "fe",
                       term = "irrigArea [all]", title = " " , 
                       colors = 1,
                       axis.title = c("Area of irrigated land (hectares)",
                                      "Probability of S. haematobium infection"))

smp.pred <- plot_model(smp3.hv, type = "pred", pred.type = "fe",
                       term = "irrigArea [all]", title = " " , 
                       #colors = c("darkred", "slategrey"),
                       axis.title = c("Area of irrigated land (hectares)",
                                      "Probability of S. mansoni infection"))

# combine prediction plots
ggarrange(shp.pred, smp.pred, labels = c("A", "B"))


#############################################################
#
#                   EGG COUNTS OUTCOMES

#############################################################

#                
#               SCHISTOSOMA HAEMATOBIUM
#


shc.mixedTMB <- glmmTMB(round(Sh_median, 0) ~ irrigArea + DEM01F + DEM04c + locF +
                          factor(headSchoolYrsCat) + factor(eth_mode3) + nFishAnyC + 
                          factor(nWivesCat) + nWatptsC + nWatPtsAgKidsC + pumpYN + wealthF +
                          kidsAnyTask + distanceWatPtC + marketDistance_kmC + irrigAreaV +
                          (1|village/C_concession),
                        data = allData,
                        family = nbinom2(link = "log"))

summary(shc.mixedTMB)  
shc.mixed_ORCI <- exp(confint(shc.mixedTMB, method = "Wald")) # 1.05 (0.97, 1.13)

# crude and adjusted models
shc.adjTMB <- glmmTMB(round(Sh_median, 0) ~ irrigArea + DEM01F + DEM04c + locF +
                        factor(headSchoolYrsCat) + factor(eth_mode3) + nFishAnyC + 
                        factor(nWivesCat) + nWatptsC + nWatPtsAgKidsC + pumpYN + wealthF +
                        kidsAnyTask + distanceWatPtC + marketDistance_kmC + irrigAreaV,
                      data = allData,
                      family = nbinom2(link = "log"))
summary(shc.adjTMB)
shc.adj_ORCI <- exp(confint(smc.adjTMB, method = "Wald")) # 1.04 (0.85, 1.28)


shc.crudeTMB <- glmmTMB(round(Sh_median, 0) ~ irrigArea,
                        data = allData,
                        family = nbinom2(link = "log"))
summary(shc.crudeTMB)
shc.crude_ORCI <- exp(confint(smc.crudeTMB, method = "Wald")) # 1.02 (0.85, 1.22)


# prediction plot 
shc.pred <- plot_model(shc.mixedTMB, type = "pred", pred.type = "fe",
                       term = "irrigArea [all]", title = " " , 
                       #colors = c("darkred", "slategrey"),
                       axis.title = c("Area of irrigated land (hectares)",
                                      "Intensity of S. mansoni infection"))

# residual plots
shc_res = simulateResiduals(shc.mixedTMB)
shcResid <- plot(shc_res)
testResiduals(shc_res)

#                
#               SCHISTOSOMA MANSONI
#

smc.mixedTMB <- glmmTMB(round(Sm_median, 0) ~ irrigArea + DEM01F + DEM04c + locF +
                          factor(headSchoolYrsCat) + factor(eth_mode3) + nFishAnyC + 
                          factor(nWivesCat) + nWatptsC + nWatPtsAgKidsC + pumpYN + wealthF +
                          kidsAnyTask + distanceWatPtC + marketDistance_kmC + irrigAreaV +
                          (1|village/C_concession),
                        data = allData,
                        family = nbinom2(link = "log"))

summary(smc.mixedTMB)  
smc.mixed_ORCI <- exp(confint(smc.mixedTMB, method = "Wald")) # 1.07 (0.82, 1.39)

# crude and adjusted
smc.adjTMB <- glmmTMB(round(Sm_median, 0) ~ irrigArea + DEM01F + DEM04c + locF +
                          factor(headSchoolYrsCat) + factor(eth_mode3) + nFishAnyC + 
                          factor(nWivesCat) + nWatptsC + nWatPtsAgKidsC + pumpYN + wealthF +
                          kidsAnyTask + distanceWatPtC + marketDistance_kmC + irrigAreaV,
                        data = allData,
                        family = nbinom2(link = "log"))

summary(smc.adjTMB)  
smc.adj_ORCI <- exp(confint(smc.adjTMB, method = "Wald")) # 1.04 (0.85, 1.28)

smc.crudeTMB <- glmmTMB(round(Sm_median, 0) ~ irrigArea,
                      data = allData,
                      family = nbinom2(link = "log"))

summary(smc.crudeTMB)  
smc.crude_ORCI <- exp(confint(smc.crudeTMB, method = "Wald")) # 1.02 (0.85, 1.22)



# combine crude, adjusted and mixed estimates
shc_ORCI <- rbind(shc.crude_ORCI[2,], shc.adj_ORCI[2,], shc.mixed_ORCI[2,])
rownames(shc_ORCI) <- c("crude", "adjusted", "mixed")
colnames(shc_ORCI) <- c("LL", "UL", "IRR")
shc_ORCI <- as.data.frame(shc_ORCI)
shc_ORCI[ "Estimate" ] <- rownames(shc_ORCI)
shc_ORCI$species <- "Sh"

smc_ORCI <- rbind(smc.crude_ORCI[2,], smc.adj_ORCI[2,], smc.mixed_ORCI[2,])
rownames(smc_ORCI) <- c("crude", "adjusted", "mixed")
colnames(smc_ORCI) <- c("LL", "UL", "IRR")
smc_ORCI <- as.data.frame(smc_ORCI)
smc_ORCI[ "Estimate" ] <- rownames(smc_ORCI)
smc_ORCI$species <- "Sm"

intensity_IRRCI <- rbind(shc_ORCI, smc_ORCI)
intensity_IRRCI$species <- factor(intensity_IRRCI$species, levels = c("Sm", "Sh")) # puts Sh on top

# plot intensity results
intensity_plot <- ggplot(intensity_IRRCI,
                         aes(x = species, y = IRR, ymin = LL, ymax = UL)) +
                  geom_pointrange(aes(lty = Estimate),
                                  position = position_dodge(width = 0.5)) + 
                  geom_errorbar(aes(ymin = LL, ymax = UL, lty = Estimate), 
                                width = 0.2, position = position_dodge(width = 0.5)) +
                  geom_hline(yintercept = 1, linetype = 2) +
                  xlab(" ") + 
                  ylab("Rate Ratio") +
                  scale_y_continuous(trans = 'log2', limits = c(0.5, 1.5), breaks = seq(0.5, 1.5, 0.25)) +
                  coord_flip() +
                  #scale_colour_manual(values = c("grey41", "darkslategrey", "darkred")) +
                  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
                  guides(col = guide_legend(reverse = T)) +
                  theme_bw()

# prediction plot 
smc.pred <- plot_model(smc.mixedTMB, type = "pred", pred.type = "fe",
                       term = "irrigArea [all]", title = " " , 
                       #colors = c("darkred", "slategrey"),
                       axis.title = c("Area of irrigated land (hectares)",
                                      "Intensity of S. mansoni infection"))

# residual plots
smc_res = simulateResiduals(smc.mixedTMB)
smcResid <- plot(smc_res)
testResiduals(smc_res)

#                
#               COMBINED PLOTS
#

# import intensity results generated in State
#intensity_irr <- read.csv("C:/Users/andre/Box Sync/Papers/Ch2IRisk-IrrigatedAg/Analysis/intensity_irr_nonDAG.csv")
# make sure categories displayed in desired order
#intensity_irr$Estimate <- factor(intensity_irr$Estimate, levels = c("crude", "adjusted", "mixed"))
#intensity_irr$species <- factor(intensity_irr$species, levels = c("Sm", "Sh")) # puts Sh on top
#intensity_irr <- dplyr::mutate(intensity_irr, outcome = "presence")


# combine presence and intensity plots
ggarrange(presence_plot, intensity_plot, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")

# prediction plots
ggarrange(shc.pred, smc.pred, labels = c("A", "B"))

# regression output for all count models
tab_model(shc.crudeTMB,
          transform = NULL, show.intercept = T,
          show.est = T, show.se = T, show.p = T,
          show.aic = T, show.loglik = T, show.obs = T,
          show.ci = F,
          file = "shcCrude_tab.doc")
tab_model(shc.adjTMB,
          transform = NULL, show.intercept = T,
          show.est = T, show.se = T, show.p = T,
          show.aic = T, show.loglik = T, show.obs = T,
          show.ci = F, 
          file = "shcAdj_tab.doc")
tab_model(shc.mixedTMB,
          transform = NULL, show.intercept = T,
          show.est = T, show.se = T, show.p = T,
          show.aic = T, show.loglik = T, show.obs = T,
          show.ci = F, 
          file = "shcMix_tab.doc")
tab_model(smc.crudeTMB,
          transform = NULL, show.intercept = T,
          show.est = T, show.se = T, show.p = T,
          show.aic = T, show.loglik = T, show.obs = T,
          show.ci = F, 
          file = "smcCrude_tab.doc")
tab_model(smc.adjTMB,
          transform = NULL, show.intercept = T,
          show.est = T, show.se = T, show.p = T,
          show.aic = T, show.loglik = T, show.obs = T,
          show.ci = F, 
          file = "smcAdj_tab.doc")
tab_model(smc.mixedTMB,
          transform = NULL, show.intercept = T,
          show.est = T, show.se = T, show.p = T,
          show.aic = T, show.loglik = T, show.obs = T,
          show.ci = F, 
          file = "smcMix_tab.doc") 

###################################################################
# designate outcome variables
allData$Sh_medianR = round(allData$Sh_median, 0)
shData$Sh_medianR = round(shData$Sh_median, 0)

## bivariate linear model - discern Poisson versus negative binomial
shcp1 <- glm(shOutcome ~ irrigArea, data = shData, family = poisson)
summary(shcp1)
# irrigArea   0.012372   0.002029   6.098 1.07e-09 ***
confint(shcp1)
# irrigArea   0.008343243 0.01629627
exp(confint(shcp1))
# irrigArea    1.008378  1.01643
explainedD(shcp1) # 0.02802147
dispersion(shcp1) # 99.85911
# compare to dispersion parameter to manual calculation
manPhi = shcp1$deviance / (1232-2); manPhi # 99.85911!
# apparent overdispersion is definitely present

# plot fitted values
MyData <- data.frame(irrigArea = seq(from = 0, to = 39, by = 0.5))
Pred <- predict(shcp1, newdata = MyData, type = "response")
plot(x = allData$irrigArea, y = allData$Sh_meanR)
lines(MyData$irrigArea, Pred) 

#
#     IS THE APPARENT OVERDISPERSION REAL?    per Hilbe 2011 page 142
#

# 1. add predictors - fit model with all pre-screened covariates
shcp2.0 <- glm(Sh_median ~ irrigArea + locF + DEM01F + DEM04 +
                 headSchoolYrs + eth_mode3 + nFishAny +
                 nWatPtsAgKids + pump +
                 distanceWatPt + marketDistance_km, data = shData, family = poisson)
summary(shcp2.0)
explainedD(shcp2.0) # 20.55695 (with median = 20.57768)
dispersion(shcp2.0) # 80.13516 (with median = 79.86144)
# take-away: dispersion parameter decreases with inclusion of additiona predictors but remains well-above 1


# 2. test for interaction
# 2.1 river-lake geography
shcp2.1 <- glm(Sh_median ~ irrigArea*locF + DEM01F + DEM04 +
                 headSchoolYrs + eth_mode3 + nFishAny +
                 nWatPtsAgKids + pump +
                 distanceWatPt + marketDistance_km, data = shData, family = poisson)
summary(shcp2.1) # interaction term is highly significant
explainedD(shcp2.1) # 21.92029 (with median = 21.94346)
dispersion(shcp2.1) # 78.82466 (with median = 78.55261)
anova(shcp2.0, shcp2.1, test = "Chisq") # chisq = 1675, df = 1, p < 2.2e-16
# with median, chisq = 1672.7, df = 1, p < 2.2e-16


# 2.2 pump ownership
shcp2.2 <- glm(Sh_median ~ irrigArea*pump + locF + DEM01F + DEM04 +
                headSchoolYrs + eth_mode3 + nFishAny +
                nWatPtsAgKids +
                distanceWatPt + marketDistance_km, data = shData, family = poisson)
summary(shcp2.2) # interaction is highly significant
explainedD(shcp2.2) # 20.71099 (with median = 20.73202)
dispersion(shcp2.2) # 80.04549 (with median = 79.77174)
anova(shcp2.0, shcp2.2, test = "Chisq") # Chisq = 189.26, df = 1, p < 2.2e-16
# median --> chisq = 189.02, same conclusion

# 2.3 three-way loc*pump*irrigArea interaction
shcp2.3 <- glm(Sh_median ~ irrigArea*pump*locF + DEM01F + DEM04 +
                 headSchoolYrs + eth_mode3 + nFishAny +
                 nWatPtsAgKids +
                 distanceWatPt + marketDistance_km, data = shData, family = poisson)
summary(shcp2.3) # three-way interaction is not significant (with median p = 0.047)
explainedD(shcp2.3) # 22.14103 (median --> 22.16189)
dispersion(shcp2.3) # 78.79605 (median --> 78.52636)
# take-away: none of the interaction terms reduced the dispersion parameter substantially

# 3. check for influential observations
plot(shcp2.1, which = 5) # in residuals vs leverage plot, there are several observations outside the Cooks D line

# investigate influential observations further, following guidance provided at link below: 
# https://stats.stackexchange.com/questions/164099/removing-outliers-based-on-cooks-distance-in-r-language
# calculate Cooks distance for model 
cooksd <- cooks.distance(shcp2.1)
# Plot the Cook's Distance using the traditional 4/n criterion 
# note: I do not know anything about this 4/n criterion. I recall Kleinbaum's threshold was much higher, ~ 2
sample_size <- nrow(shData)
plot(cooksd, pch="*", cex=1.5, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 2, col = "red") # add cutoff line
#abline(h = 4/sample_size, col="red")  # add cutoff line
#text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd >= 2, names(cooksd),""), col="red")  # add labels
length(which(cooksd >= 2)) # n = 16
length(which(cooksd >= 4/sample_size)) # n = 922 - over HALF the observations. this doesn't seem reasonable

# subset data to exclude 12 most influential observations
shData_r <- shData[which(cooksd < 2),] # n = 1220

# fit same model with reduced data set
shcp2.1r <- glm(Sh_median ~ irrigArea*locF + DEM01F + DEM04 +
                 headSchoolYrs + eth_mode3 + nFishAny +
                 nWatPtsAgKids + pump  +
                 distanceWatPt + marketDistance_km, data = shData_r, family = poisson)
summary(shcp2.1r)
explainedD(shcp2.1r) # 26.31191 (median --> 26.3444)
dispersion(shcp2.1r) # 59.03235 (median --> 58.77585)
# take-away: removing influential observations brings dispersion parameter down, but it still does not approach 1

# 4. center/scale continuous variables with no logical zero value
shcp4 <- glm(Sh_median ~ irrigArea*locF + DEM01F + DEM04c +
               headSchoolYrs + eth_mode3 + nFishAny +
               nWatPtsAgKids + pump +
               distanceWatPtC + marketDistance_kmC, data = shData_r, family = poisson)
summary(shcp4)
explainedD(shcp4) # 26.31191 (median --> 26.3444)
dispersion(shcp4) # 59.03235 (median --> 58.77585)
# take-away: dispersion parameter remains high even after centering/scaling relevant variables

# 5. check alternative link functions by running Tukey-Pregibon link test (Hilbe 2011 p 156)
# if hat2 term is statistically significant, link to original model is not well-specified
hat <- hatvalues(shcp4)
hat2 <- hat*hat
dat <- as.data.frame(cbind(shData_r$Sh_meanR, hat, hat2))
tplt <- glm(V1 ~ hat + hat2, data = dat, family = poisson)
summary(tplt)
confint(tplt) # hat2 is significant. 
# conclusion: link function is not appropriate. the apparent overdispersion is real.



#
#               FIT NEGATIVE BINOMIAL MODEL
#

# bivariate model
shcnb1.0 <- glm.nb(outcome ~ irrigArea, data = shData, link = log)
summary(shcnb1.0)
plot(shcnb1.0, which = 5) # all observations are within reasonable thresholds of Cook's distance

# plot fitted values
MyData <- data.frame(irrigArea = seq(from = 0, to = 39, by = 0.5))
Pred <- predict(shcnb1.0, newdata = MyData, type = "response")
plot(x = shData$irrigArea, y = outcome, 
     xlab = "Area of irrigated land (hectares)", 
     ylab = "Intensity of urogenital schistosome (S. haematobium) infection ")
lines(MyData$irrigArea, Pred)

# fit saturated model with all pre-screened covariates
  shcnb2.0 <- glm.nb(outcome ~ irrigArea + locF + DEM01F + DEM04c +
                       headSchoolYrs + nWives + eth_mode3 + nFishAny +
                       nWatpts + nWatPtsAgKids + kidsAnyTask + pumpYN +
                       distanceWatPtC + marketDistance_kmC + wealthF, 
                    data = shData)

# test negative binomial model assumptions via LR test comparing nb and poisson models
as.numeric(pchisq(2*(logLik(shcnb2.0)-logLik(shcp2.0)), df = 1, lower.tail = F)) # 0 
lrtest(shcp2.0, shcnb2.0) # chisq = 92586, df = 1, p < 2.2e-16
# negative binomial confirmed superior to poisson regression

### PROCEED WITH PRE-SPECIFIED MODEL SELECTION STEPS
## test pre-specified interaction terms
# 1. location
shcnb2.1 <- glm.nb(outcome ~ irrigArea*locF + DEM01F + DEM04c +
                  headSchoolYrs + nWives + eth_mode3 + nFishAny +
                  nWatpts + nWatPtsAgKids + kidsAnyTask + pumpYN +
                  distanceWatPtC + marketDistance_kmC + wealthF, data = shData)
# use likelihood ratio test to decide if interactions should be retained? 
lrtest(shcnb2.0, shcnb2.1) # Chisq = 7.8025, df = 1, p = 0.005217; yes, retain interaction

# 2. pump
shcnb2.2 <- glm.nb(outcome ~ irrigArea*pumpYN + locF + DEM01F + DEM04c +
                     headSchoolYrs + nWives + eth_mode3 + nFishAny +
                     nWatpts + nWatPtsAgKids + kidsAnyTask +
                     distanceWatPtC + marketDistance_kmC + wealthF, data = shData)
lrtest(shcnb2.0, shcnb2.2) # chisq = 0.6992, df = 1, p = 0.4031; don't retain interaction
# take away: proceed with models stratified on location 

### STRATIFIED MODELS
## river
shcnb3r <- glm.nb(outcomeR ~ irrigArea + pumpYN + DEM01F + DEM04c +
                    headSchoolYrs + nWives + eth_mode3 + nFishAny +
                    nWatpts + nWatPtsAgKids + kidsAnyTask +
                    distanceWatPtC + marketDistance_kmC + wealthF, 
                 data = shRiver) 
# Error: NA/NaN/Inf in 'x'
# https://stats.stackexchange.com/questions/52527/unable-to-fit-negative-binomial-regression-in-r-attempting-to-replicate-publish
# may be a convergence problem that arises from specifying too many varaibles for the number of observations

## lake
shcnb3l <- glm.nb(outcomeL ~ irrigArea + pumpYN + DEM01F + DEM04c +
                    headSchoolYrs + nWives + eth_mode3 + nFishAny +
                    nWatpts + nWatPtsAgKids + kidsAnyTask +
                    distanceWatPtC + marketDistance_kmC + wealthF, 
                 data = shLake)

step(shcnb3l)

# fit random effects
# household level random intercepts
shcnb3l.h <- glmer.nb(outcomeL ~ irrigArea + pumpYN + DEM01F + DEM04c +
                       headSchoolYrsC + nWivesC + eth_mode3 + nFishAnyC +
                       nWatptsC + nWatPtsAgKidsC + kidsAnyTask +
                       distanceWatPtC + marketDistance_kmC + wealthF + 
                       (1|C_concession), 
                     data = shLake,
                    control= glmerControl(optimizer="bobyqa", 
                                          tolPwrss=1e-3,
                                          optCtrl = list(maxfun = 100000)))

# took a very long time to run (with or without bobyqa optimizer)
# boundary (singular) fit)
# generated 10 warnings: (1, 3, 5, 7, 9) max|grad| convergence warning, 
# (2, 4, 6, 8, 10) model nearly unidentifiable - rescale variables?,
isSingular(shcnb3.h) # FALSE

lrtest(shcnb3l, shcnb3l.h) # Chisq = 14.502, df = 1, p = 0.00014

# households nested in village level random intercepts
shcnb3l.hv <- glmer.nb(outcomeL ~ irrigArea + pumpYN + DEM01F + DEM04c +
                         headSchoolYrs + nWives + eth_mode3 + nFishAny +
                         nWatpts + nWatPtsAgKids + kidsAnyTask +
                         distanceWatPtC + marketDistance_kmC + wealthF +
                         (1|village/C_concession), 
                     data = shLake,
                     control= glmerControl(optimizer="bobyqa", 
                                           tolPwrss=1e-3,
                                           optCtrl = list(maxfun = 100000)))

# boundary (singular) fit
isSingular(shcnb3l.hv) # FALSE
# compare mixed effects model to fixed effects only
lrtest(shcnb3l, shcnb3l.hv) # Chisq = 75.096, df = 2, p < 2.2e-16
# model with random effects (households nested in villages) fits better than model without them
lrtest(shcnb3.h, shcnb3l.hv) # Chisq = 60.594, df = 1, p < 7.015e-15
# model with nested random effects fits better than household-level random effects only

### PROCEED WITH NESTED RANDOM EFFECTS MODEL - shcnb3l.hv
exp(confint.merMod(shcnb3l.hv, devtol = 1e-01))
# Error in zeta(shiftpar, start = opt[seqpar1][-w]) : 
# profiling detected new, lower deviance
# suggested solution: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q3/022394.html
# tweaking devtol all the way up to 1e-02 generates the same warning, which is probably an untenable tolerance

ci_calc(shcnb3l.hv) # IRR = 0.97481913, 95% CI (0.87781651, 1.0825410)

## CONTINUE THIS WORK IN STATA.
## BUT ALSO EXAMINE DISTRIBUTION OF OUTCOME 
summary(allData$Sh_medianR) # max 1745
histogram(allData$Sh_medianR, xlim = c(0,600), breaks = 600)
which(allData$Sh_medianR > 600) # single observation (1319) gt 600

allData[1319,] # ST/1/030 in household CC/ST/19 
# this individual has urine egg counts of 310 and 3180 - could this be a typo?

# recode the outcome for this observation to 345 (average of 310 and 380)
allData$Sh_medianR <- ifelse(allData$Sh_medianR > 600, 340, allData$Sh_medianR)

# rerun model with recoded variable
shc_hv_recode <- glmer.nb(Sh_medianR ~ irrigArea*locF + pumpYN + DEM01F + DEM04c +
                         factor(headSchoolYrsCat) + factor(nWivesCat) + wealthF + irrigAreaV +
                         (1|village/C_concession), 
                       data = allData) 
# still takes several minutes to run
# 28 warnings
confint(shc_hv_recode)



#############################################################
#
#                 SCHISTOSOMA MANSONI
#

# designate outcome variables
allData$Sm_medianR = round(allData$Sm_median, 0)
smData$Sm_medianR = round(smData$Sm_median, 0)

## bivariate linear model - discern Poisson versus negative binomial
smcp1 <- glm(smOutcome ~ irrigArea, data = smData, family = poisson)
summary(smcp1)
explainedD(smcp1) # 0.02241996
dispersion(smcp1) # 199.6639 - holy cow! 

###       IS THE APPARENT OVERDISPERSION REAL - per Hilbe 2011 page 142
# 1. adding predictors - fit model with all pre-screened covariates
smcp2.0 <- glm(smOutcome ~ irrigArea + locF + DEM01F + DEM04 +
                 headSchoolYrs + nWives + eth_mode3 + nFishAny + nWatpts +
                 nWatPtsAgKids + pumpYN + wealthF + 
                 distanceWatPt + marketDistance_km, data = smData, family = poisson)
summary(smcp2.0)
explainedD(smcp2.0) # 32.68599
dispersion(smcp2.0) # 136.1054
# take-away: dispersion parameter decreases with inclusion of additiona predictors but remains well-above 1

# 2. test for interaction
# 2.1 river-lake geography
smcp2.1 <- glm(smOutcome ~ irrigArea*locF + DEM01F + DEM04 +
                 headSchoolYrs + nWives + eth_mode3 + nFishAny + nWatpts +
                 nWatPtsAgKids + pumpYN + wealthF + 
                 distanceWatPt + marketDistance_km, data = smData, family = poisson)
summary(smcp2.1) # interaction term is highly significant
explainedD(smcp2.1) # 33.03943
dispersion(smcp2.1) # 135.5032

# 2.2 pump ownership
smcp2.2 <- glm(smOutcome ~ irrigArea*pump + locF + DEM01F + DEM04 +
                 headSchoolYrs + nWives + eth_mode3 + nFishAny + nWatpts +
                 nWatPtsAgKids + pumpYN + wealthF + 
                 distanceWatPt + marketDistance_km, data = smData, family = poisson)
summary(smcp2.2) # interaction is highly significant
explainedD(smcp2.2) # 36.00523
dispersion(smcp2.2) # 129.6092

# 2.3 three-way loc*pump*irrigArea interaction
smcp2.3 <- glm(smOutcome ~ irrigArea*pump*locF + DEM01F + DEM04 +
                 headSchoolYrs + nWives + eth_mode3 + nFishAny + nWatpts +
                 nWatPtsAgKids + pumpYN + wealthF + 
                 distanceWatPt + marketDistance_km, data = smData, family = poisson)
# warning: fitted rates numerically 0 occurred
summary(smcp2.3) # three-way interaction is not significant
explainedD(smcp2.3) # 36.28965
dispersion(smcp2.3) # 129.3557
# take-away: none of the interaction terms reduced the dispersion parameter substantially

# 3. check for influential observations
plot(smcp2.1, which = 5) # in residuals vs leverage plot, there are several observations outside the Cooks D line

# investigate influential observations further, following guidance provided at link below: 
# https://stats.stackexchange.com/questions/164099/removing-outliers-based-on-cooks-distance-in-r-language
# calculate Cooks distance for model 
cooksd <- cooks.distance(smcp2.1)

# Plot the Cook's Distance using the traditional 4/n criterion 
# note: I do not know anything about this 4/n criterion. I recall Kleinbaum's threshold was much higher, ~ 2
sample_size <- nrow(smData)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 2, col = "red") # add cutoff line
#abline(h = 4/sample_size, col="red")  # add cutoff line
#text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd >= 2, names(cooksd),""), col="red")  # add labels
length(which(cooksd >= 2)) # n = 23
length(which(cooksd >= 4/sample_size)) # n = 611 - over HALF the observations. this doesn't seem reasonable

# subset data to exclude 12 most influential observations
smData_r <- smData[which(cooksd < 2),] # n = 1198
smData_r$smOutcome <- round(smData_r$Sm_median, 0)
# fit same model with reduced data set
smcp2.1r <- glm(smOutcome ~ irrigArea*locF + DEM01F + DEM04 +
                  headSchoolYrs + nWives + eth_mode3 + nFishAny + nWatpts +
                  nWatPtsAgKids + pumpYN + wealthF + 
                  distanceWatPt + marketDistance_km, data = smData_r, family = poisson)
summary(smcp2.1r)
explainedD(smcp2.1r) # 36.26953
dispersion(smcp2.1r) # 43.84749

# 4. center/scale continuous variables with no logical zero value
smcp4 <- glm(smOutcome ~ irrigArea*locF + DEM01F + DEM04 +
               headSchoolYrs + nWives + eth_mode3 + nFishAny + nWatpts +
               nWatPtsAgKids + pumpYN + wealthF + 
               distanceWatPt + marketDistance_km, data = smData_r, family = poisson)
summary(smcp4)
explainedD(smcp4) # 36.26953
dispersion(smcp4) # 43.84749

# 5. check alternative link functions by running Tukey-Pregibon link test (Hilbe 2011 p 156)
# if hat2 term is statistically significant, link to original model is not well-specified
hat <- hatvalues(smcp4)
hat2 <- hat*hat
dat <- as.data.frame(cbind(smData_r$Sh_meanR, hat, hat2))
tplt <- glm(V1 ~ hat + hat2, data = dat, family = poisson)
summary(tplt)
confint(tplt) # hat2 is significant. 
# conclusion: link function is not appropriate. the apparent overdispersion is real.

### FIT NEGATIVE BINOMIAL MODEL
# bivariate model
smcnb1 <- glm.nb(smOutcome ~ irrigArea, data = smData)
summary(smcnb1)
plot(smcnb1, which = 5) # all observations are within reasonable thresholds of Cook's distance


# fit saturated model with all pre-screened covariates
smcnb2.0 <- glm.nb(smOutcome ~ irrigArea + locF + DEM01F + DEM04c +
                     headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                     nWatPtsAgKids + pumpYN + wealthF + 
                     distanceWatPtC + marketDistance_kmC, 
                 data = smData)
# warnings: algorithm did not converge; alternation limit reached (default number of iterations is 25)
summary(smcnb2.0)
plot(smcnb2.0, which = 5) # no observations with Cooks D ge 1

# test pre-specified interaction terms
# location
smcnb2.1 <- glm.nb(smOutcome ~ irrigArea*locF + DEM01F + DEM04c +
                     headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                     nWatPtsAgKids + pumpYN + wealthF + 
                     distanceWatPtC + marketDistance_kmC, 
                   data = smData)
# error: NA/NaN/Inf in 'x'; also algorithm did not converge, step size truncated due to divergence
anova(smcnb2.1, smcnb2.0, test = "Chisq") 
# Chisq = 49.09462, df = 2, p = 2.183931e-11 <-- weird. why df = 2??? 
# significance of LR test suggests improved fit with inclusion of irrigArea*locF interaction

# pump
smcnb2.2 <- glm.nb(smOutcome ~ irrigArea*pumpYN + locF + DEM01F + DEM04c +
                     headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                     nWatPtsAgKids + wealthF + 
                     distanceWatPtC + marketDistance_kmC, 
                   data = smData)
# error: NA/NaN/Inf in 'x'; also step size truncated due to divergence
anova(smcnb2.0, smcnb2.2) 

## PROCEED WITH STRATIFIED MODELS
# river
smcnb_r1 <- glm.nb(smOutcomeR ~ irrigArea + DEM01F + DEM04c +
                     headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                     nWatPtsAgKids + pumpYN + wealthF + 
                     distanceWatPtC + marketDistance_kmC, 
                   data = smRiver)
# error: NA/NaN/Inf in 'x'
summary(smcnb_r1)

# lake
smcnb_l1 <- glm.nb(smOutcomeL ~ irrigArea + DEM01F + DEM04c +
                     headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                     nWatPtsAgKids + pumpYN + wealthF + 
                     distanceWatPtC + marketDistance_kmC, 
                   data = smLake)
# error: NA/NaN/Inf in 'x'


## test random effects
# full fixed effects only model
smcnb3 <- glm.nb(smOutcome ~ irrigArea*locF + DEM01F + DEM04c +
                   headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                   nWatPtsAgKids + kidsAnyTask + pump + wealthF +
                   distanceWatPtC + marketDistance_kmC, 
                 data = smData, 
                 control = glm.control(maxit = 50))


visreg(smcnb3, "irrigArea", by = "locF", scale = "response")

smcnb3r <- glm.nb(Sm_meanR ~ irrigArea + DEM01F + DEM04c +
                    headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                    nWatPtsAgKids + kidsAnyTask + pump + wealthF +
                    distanceWatPtC + marketDistance_kmC, 
                  data = smRiver, 
                  control = glm.control(maxit = 50))
# Error: NA/NaN/Inf in 'x'
# smRiver has 527 observations 

smcnb3l <- glm.nb(Sm_meanR ~ irrigArea + DEM01F + DEM04c +
                    headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                    nWatPtsAgKids + kidsAnyTask + pump + wealthF +
                    distanceWatPtC + marketDistance_kmC, 
                  data = smLake, 
                  control = glm.control(maxit = 50))
# Error: NA/NaN/Inf in 'x'
# smLake has 695 observations

smcnb3r <- glm.nb(Sm_meanR ~ irrigArea + DEM01F + DEM04c +
                    headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                    nWatPtsAgKids + kidsAnyTask + pump + wealthF +
                    distanceWatPtC + marketDistance_kmC, 
                  data = smRiver, 
                  control = glm.control(maxit = 50))

# household level
smcnb3.h <- glmer.nb(Sm_meanR ~ irrigArea*locF + DEM01F + DEM04c +
                       headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                       nWatPtsAgKids + kidsAnyTask + pump + wealthF +
                       distanceWatPtC + marketDistance_kmC +
                       (1|C_concession), 
                     data = smData, 
                     control=glmerControl(optimizer="bobyqa", 
                                          tolPwrss=1e-3,
                                          optCtrl = list(maxfun = 100000)))

# 14 warnings: (1) iteration limit reached, (2) NaNs produced, (3-14) max|grad
lrtest(smcnb3, smcnb3.h) # Chisq = 51508 (?), df = 1, p < 2.2e-16

# village level
smcnb3.v <- glmer.nb(Sm_meanR ~ irrigArea*locF + DEM01F + DEM04c +
                       headSchoolYrs + nWivesC + eth_mode3 + nFishAny + nWatpts +
                       nWatPtsAgKids + kidsAnyTask + pump + wealthF +
                       distanceWatPtC + marketDistance_kmC +
                       (1|village), 
                     data = smData, 
                     control=glmerControl(optimizer="bobyqa", 
                                          tolPwrss=1e-3,
                                          optCtrl = list(maxfun = 100000)))
# boundary (singular) fit 
# plus 7 max|grad| warnings, one warning about model being nearly identifiable - rescale variables, 
# one 'unable to evaluate scaled gradient' warning and one 'degenerate Hessian with 1 negative eigenvalues' warning
lrtest(smcnb3, smcnb3.v) # Chisq = 6.274, df = 1, p = 0.01225

# CONTINUE IN STATA
summary(allData$Sm_medianR) # max 4636
histogram(allData$Sm_medianR, xlim = c(0,5000), breaks = 5000)
which(allData$Sm_medianR > 1000) # n = 7 observations with intensities gt 1000
which(allData$Sm_medianR > 2000) # n = 3 observations with intensities gt 2000

allData[187,] # GG/2/002 in CC/GG/02 with Sm egg counts of 1584 and 1452
allData[197,] # GG/4/006 in CC/GG/06 with Sm egg counts of 4950 and 4323
allData[201,] # GG/1/010 in CC/GG/06 with Sm egg counts of 2233 and 1815
allData[209,] # GG/2/002 in CC/GG/14 with Sm egg counts of 1584 and 1452 - same kid as 187 :-/
allData[821,] # MG/2/022 in CC/MG/03 with Sm egg counts of 1848 and 1023
allData[831,] # MG/2/018 in CC/MG/04 with Sm egg counts of 3993 and 3663
allData[1156,] # NE/2/021 in CC/NE/51 with Sm egg counts of 396, 561, 2013 and 1716

# recode values exceeding 2000 to 2000
allData$Sm_medianR = ifelse(allData$Sm_medianR > 2000, 2000, allData$Sm_medianR)

# refit model with recoded data
smc_nb_hv <- glmer.nb(Sm_medianR ~ irrigArea*locF + DEM01F + DEM04c +
                       factor(headSchoolYrsCat) + factor(nWivesCat) +
                       pumpYN + wealthF + irrigAreaV +
                       (1|village/C_concession), 
                     data = allData, 
                     control=glmerControl(optimizer="bobyqa"))
# 50 warnings
summary(smc_nb_hv) # very small SE values and convergence code = 0
