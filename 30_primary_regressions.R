###############################################################
# 
#    DAG-based model to estimate the relationship between
#    irrigated agriculture and schistosomiasis reinfection
#    coded by: Andrea Lund
#    last updated: January 15, 2020
#
###############################################################

# install packages
library(stats) # for a variety of statistical functions
library(lme4) # for mixed effects regression
library(MASS) # for negative binomial regression 
library(glmmTMB) # for negative binomial regression
library(car) # for vif/collinearity assessment
library(stargazer) # for publication-quality regression output tables
library(sjPlot) # custom correlation matrices 
library(msm) # multi-state markov in continuous time required for additive interactions
library(ggplot2) # for aesthetically-pleasing plots
library(performance) # for generating binned residuals, other regression diagnostics
library(DHARMa) # for diagnostics on glmmTMB models
library(lmtest) # for lrtest()
library(visreg) # for visualizing regression output
library(lattice) # for dotplots of random effects
library(ggpubr) # for combining plots 

# set working directory
setwd("C:/Users/andre/Box Sync/Papers/Ch2IRisk-IrrigatedAg")

# import merged data
allData <- read.csv("allData.csv")

# convert categorical variables to factors
allData$headSchoolYrsCat <- factor(allData$headSchoolYrsCat)
allData$nWivesCat <- factor(allData$nWivesCat)
allData$wealthF <- factor(allData$wealthF)


## FUNCTIONS
'%notin%' <- function(x,y) {!(x %in% y)}

## calculate confidence intervals for fitted models
# https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
ci_calc <- function(mod) {
  se <- sqrt(diag(vcov(mod)))
  estCI <- cbind(OR = fixef(mod), 
                 LL = fixef(mod)-1.96*se, 
                 UL = fixef(mod)+1.96*se)
  orCI <- exp(estCI)
  return(orCI)
}

# calculate confidence intervals for fitted models with irrigArea X loc interactions
ci_calc_int <- function(mod) {
  # all coefficients and ses
  se <- sqrt(diag(vcov(mod))); coefs <- coef(mod)
  # subset coefficients and standard errors
  beta1 <- coefs[2]; se1 <- se[2]
  beta3 <- coefs[13]; se3 <- se[13]
  # estimates for river stratum
  estR <- beta1 + beta3
  llR <- (beta1 - se1) + (beta3 - se3)
  ulR <- (beta1 + se1) + (beta3 + se3)
  estCI_R <- cbind(estR, llR, ulR); ORCI_R <- exp(estCI_R)
  # estimates for lake stratum
  estL <- beta1; llL <- beta1 - se1; ulL <- beta1 + se1
  estCI_L <- cbind(estL, llL, ulL); ORCI_L <- exp(estCI_L)
  # package output into single matrix
  strata <- c(1, 0) # 1 = river, 0 = lake
  ORCI <- rbind(ORCI_R, ORCI_L);  ORCI <- cbind(strata, ORCI)
  return(ORCI)
}

# calculate confidence intervals for fitted mixed models with irrigArea x loc interaction
ci_calc_int_mixed <- function(mod) {
  # all coefficients and ses
  se <- sqrt(diag(vcov(mod))); coef <- fixef(mod)
  # subset coefficients and standard errors
  beta1 <- coef[2]; se1 <- se[2]
  beta3 <- coef[13]; se3 <- se[13]
  # estimates for river stratum
  estR <- beta1 + beta3
  llR <- (beta1 - se1) + (beta3 - se3)
  ulR <- (beta1 + se1) + (beta3 + se3)
  estCI_R <- cbind(estR, llR, ulR); ORCI_R <- exp(estCI_R)
  # estimates for lake stratum
  estL <- beta1; llL <- beta1 - se1; ulL <- beta1 + se1
  estCI_L <- cbind(estL, llL, ulL); ORCI_L <- exp(estCI_L)
  # package output into single matrix
  strata <- c(1, 0) # 1 = river, 0 = lake
  ORCI <- rbind(ORCI_R, ORCI_L);  ORCI <- cbind(strata, ORCI)
  return(ORCI)
}

###############################################
#
#                OUTCOME #1
#               Sh presence
#
   
# crude model           
shp.crude <- glm(Sh_presence ~ irrigArea, data = allData, family = binomial(link = "logit"))

# estimate point estimate and confidence intervals for crude model
est_shp.crude_all <- exp(coef(shp.crude))
est_shp.crude_exp <- est_shp.crude_all[2]
lim_shp.crude_all <- exp(confint(shp.crude))
lim_shp.crude_exp <- lim_shp.crude_all[2,]
shp.crude_ORCI <- c(est_shp.crude_exp, lim_shp.crude_exp)

# adjusted model
shp.adj <- glm(Sh_presence ~ irrigArea + locF + DEM01F + DEM04c +
                   headSchoolYrsCat + nWivesCat + 
                   pumpYN + wealthF + irrigAreaV, 
               data = allData, family = binomial(link = "logit"))


summary(shp.adj)

est_shp.adj_all <- exp(coef(shp.adj))
est_shp.adj_exp <- est_shp.adj_all[2]
lim_shp.adj_all <- exp(confint(shp.adj))
lim_shp.adj_exp <- lim_shp.adj_all[2,]
shp.adj_ORCI <- c(est_shp.adj_exp, lim_shp.adj_exp)

shp.adj_loc <- glm(Sh_presence ~ irrigArea*locF + DEM01F + DEM04c +
                     headSchoolYrsCat + nWivesCat + 
                     pumpYN + wealthF + irrigAreaV, 
                   data = allData, family = binomial(link = "logit"))

anova(shp.adj, shp.adj_loc, test = "Chisq") # Chisq = 2.1162, df = 1, p = 0.1457

shp.adj_pump <- glm(Sh_presence ~ irrigArea*pumpYN + locF + DEM01F + DEM04c +
                     headSchoolYrsCat + nWivesCat + 
                     wealthF + irrigAreaV, 
                   data = allData, family = binomial(link = "logit"))

anova(shp.adj, shp.adj_pump, test = "Chisq") # Chisq = 1.0481, df = 1, p = 0.3059


# mixed model
shp.mixed <- glmer(Sh_presence ~ irrigArea + locF + DEM01F + DEM04c +
                   headSchoolYrsCat + nWivesCat + 
                   pumpYN + wealthF + irrigAreaV +
                   (1|village/C_ConcessionId), 
                 data = allData, family = binomial(link = "logit"))
# max grad = 0.00904828

summary(shp.mixed)

anova(shp.mixed, shp.adj, test = "Chisq") # Chisq = 94.689, df = 2, p < 2.2e-16
anova(shp.adj, shp.crude, test = "Chisq") # Chisq = 223.25, df = 13, p < 2.2e-16
anova(shp.mixed, shp.crude, test = "Chisq") # Chisq = 321.18, df = 15, p < 2.2e-16

shp.mixed_loc <- glmer(Sh_presence ~ irrigArea*locF + DEM01F + DEM04c +
                         headSchoolYrsCat + nWivesCat + 
                         pumpYN + wealthF + irrigAreaV +
                         (1|village/C_ConcessionId), 
                       data = allData, family = binomial(link = "logit"))

anova(shp.mixed_loc, shp.mixed, test = "Chisq") # Chisq = 0.6009, df = 1, p = 0.4382
# no location interaction in mixed model

shp.mixed_pump <- glmer(Sh_presence ~ irrigArea*pumpYN + locF + DEM01F + DEM04c +
                         headSchoolYrsCat + nWivesCat + 
                         wealthF + irrigAreaV +
                         (1|village/C_ConcessionId), 
                       data = allData, family = binomial(link = "logit"))
anova(shp.mixed_pump, shp.mixed, test = "Chisq") # Chisq = 0.44, df = 1, p = 0.44
# no pump interaction in mixed model

shp.mixed_r <- glmer(Sh_presence ~ locF + DEM01F + DEM04c +
                       factor(headSchoolYrsCat) + factor(nWivesCat) + 
                       pumpYN + wealthF + irrigAreaV +
                       (1|village/C_ConcessionId), 
                     data = allData, family = binomial(link = "logit"))

anova(shp.mixed, shp.mixed_r, test = "Chisq") # Chisq = 7.3339, df = 1, p = 0.006767
# irrigated area is a significant predicted of sh infection presence by LRT

### model diagnostics
# residuals
shp.resid <- binned_residuals(shp.mixed) # 91%
# residuals by predictor
plot_model(shp.mixed, type = "resid")
binned_residuals(model = shp.mixed, term = "irrigArea") # 93%
binned_residuals(model = shp.mixed, term = "DEM04c") # 100%
binned_residuals(model = shp.mixed, term = "pumpYN") # 100%
binned_residuals(model = shp.mixed, term = "wealthF") #100%
binned_residuals(model = shp.mixed, term = "irrigAreaV") # 93%
# collinearity
vif(shp.mixed) # highest value 1.27 for nWivesCat and 1.39 for wealthF
# outliers/influence
cooksd <- cooks.distance(shp.mixed) # no obs exceeds 2; one obs exceeds 1.5
plot(cooksd, pch="*", cex=1.5, main="Influential Obs by Cooks distance")

# random effects
dotplot(ranef(shp.mixed, condVar = T, transf = exp))
qqnorm(ranef(shp.mixed)$village[[1]], main = "Village Effects")

# model estimates
#exp(confint.merMod(shp.mixed, devtol = 1e-05)) # doesnt immediately produce zeta error when devtol = 1e-06 but runs for several minutes
shp.mixed_aOR_all <- ci_calc(shp.mixed)
shp.mixed_ORCI <- shp.mixed_aOR_all[2,]

# combine all shp model estimates 
shp_OR <- as.data.frame(rbind(shp.crude_ORCI, shp.adj_ORCI, shp.mixed_ORCI))
colnames(shp_OR) <- c("OR", "LL", "UL")
shp_OR <- shp_OR %>% 
  dplyr::mutate(Estimate = c("crude", "adjusted", "mixed")) 

# generate forest plot for all Sh presence estimates
ggplot(shp_OR, 
       aes(x = Estimate, y = OR, ymin = LL, ymax = UL)) +
  geom_pointrange(aes(col = Estimate)) + 
  geom_errorbar(aes(ymin = LL, ymax = UL, col = Estimate), width = 0.2) +
  geom_hline(yintercept =1, linetype=2)+
  ylab("Odds ratio") +
  scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0.5, 1.5, 0.25)) +
  #facet_grid(strata ~ ., scales = "free_y" ) +
  coord_flip() +
  scale_colour_manual(values = c("grey41", "darkslategrey", "darkred")) +
  guides(col = guide_legend(reverse = T)) +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20)) +
  theme_classic()

# plot predicted probabilities
shp.pred <- plot_model(shp.mixed, type = "pred", pred.type = "fe",
           term = "irrigArea [all]", title = " " , 
           colors = 1,
           axis.title = c("Area of irrigated land (hectares)",
                           "Probability of S. haematobium infection"))


###############################################
#
#                OUTCOME #2
#               Sm presence
#

# crude model           
smp.crude <- glm(Sm_presence ~ irrigArea, data = allData, family = binomial(link = "logit"))
summary(smp.crude)

smp.crude_OR <- confint(smp.crude)

# adjusted model
smp.adj <- glm(Sm_presence ~ irrigArea + locF + DEM01F + DEM04c +
                     headSchoolYrsCat + nWivesCat + pumpYN + wealthF +
                     irrigAreaV, data = allData, family = binomial(link = "logit"))
summary(smp.adj)

smp.adj_int <- glm(Sm_presence ~ irrigArea*locF + DEM01F + DEM04c +
                 headSchoolYrsCat + nWivesCat + pumpYN + wealthF +
                 irrigAreaV, data = allData, family = binomial(link = "logit"))

anova(smp.adj, smp.adj_int, test = "Chisq") # Chisq = 12.948, df = 1, p = 0.0003202
# retain interaction model

smp.adj_ORCI <- c(exp(coef(smp.adj)[2]), exp(confint(smp.adj)[2,]))

# mixed model 
smp.mixed <- glmer(Sm_presence ~ irrigArea + locF + DEM01F + DEM04c +
                     headSchoolYrsCat + nWivesCat + pumpYN + wealthF +
                     irrigAreaV +
                     (1|village/C_ConcessionId), 
                   data = allData, family = binomial(link = "logit"),
                   control= glmerControl(optimizer="bobyqa"))
# with Nelder Mead optimizer: max grad = 0.00660773; large eigenvalue ratio
# with bobyqa optimizer: max number of function evaluations exceeded; max grad = 0.0365815; large eigenvalue ratio 

summary(smp.mixed)

anova(smp.mixed, smp.adj, test = "Chisq") # Chisq = 99.804, df = 2, p = <2.2e-16

smp.mixed_loc <- glmer(Sm_presence ~ irrigArea*locF + DEM01F + DEM04c +
                         headSchoolYrsCat + nWivesCat + pumpYN + wealthF +
                         irrigAreaV +
                         (1|village/C_ConcessionId), 
                       data = allData, family = binomial(link = "logit"),
                       control= glmerControl(optimizer="bobyqa"))
# convergence warning max grad = 0.0164299; large eigenvalue ratio

anova(smp.mixed, smp.mixed_loc, test = "Chisq") # Chisq = 1.7556, df = 1, p = 0.1852
# interaction not significant when accounting for data structure via random effects
# proceed without it

smp.mixed_pump <- glmer(Sm_presence ~ irrigArea*pumpYN + locF + DEM01F + DEM04c +
                         headSchoolYrsCat + nWivesCat + wealthF +
                         irrigAreaV +
                         (1|village/C_ConcessionId), 
                       data = allData, family = binomial(link = "logit"),
                       control= glmerControl(optimizer="bobyqa"))
# warnings: maximum number of function evaluations exceeded; unable to evaluate scaled gradient; degenerate Hessian

anova(smp.mixed, smp.mixed_pump, test = "Chisq") # Chisq = 2.2355, df = 1, p = 0.1349
# interaction not significant; proceed without it

# compare crude, adjusted and mixed models
anova(smp.mixed, smp.adj, test = "Chisq")
anova(smp.crude, smp.adj, test = "Chisq")

### model diagnostics
# residuals
smp.resid <- binned_residuals(smp.mixed) # 74% - yikes. porbably bad model fit
# residuals by predictor
plot_model(smp.mixed, type = "resid") # poor fit likely driven by headSchoolYrsCat
binned_residuals(model = smp.mixed, term = "irrigArea") # 92%
binned_residuals(model = smp.mixed, term = "DEM04c") # 100%
binned_residuals(model = smp.mixed, term = "pumpYN") # 100%
binned_residuals(model = smp.mixed, term = "irrigAreaV") # 100%
# collinearity
vif(smp.mixed) # highest value 1.25 for nWivesCat and 1.31 for wealth
# outliers/influence
cooksd <- cooks.distance(smp.mixed) # largest Cooks D is ~1 
plot(cooksd, pch="*", cex=1.5, main="Influential Obs by Cooks distance")
#text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd >= 2, names(cooksd),""), col="red") 
#length(which(cooksd > 10)) # n = 1
#allData[634,] # CC/MD/07 - MD/1/006 - household with irrigArea = 10 ha

# random effects
dotplot(ranef(smp.mixed, condVar = T, transf = exp))
qqnorm(ranef(smp.mixed)$village[[1]], main = "Village Effects")

# crude
smp.crude_OR <-c(exp(coef(smp.crude)[2]), exp(confint(smp.crude)[2,]))

# adjusted
smp.adj_OR <- c(exp(coef(smp.adj)[2]), exp(confint(smp.adj)[2,]))

# mixed
smp.mixed_OR <- ci_calc(smp.mixed)[2,]

# combine all estimates into single data frame
smp_OR <- rbind(smp.crude_OR, smp.adj_OR, smp.mixed_OR); smp_OR
smp_OR <- as.data.frame(smp_OR) 
colnames(smp_OR) <- c("OR", "LL", "UL")
smp_OR$Estimate <- c("crude", "adjusted", "mixed"); smp_OR
smp_OR$Estimate <- factor(smp_OR$Estimate, levels = c("crude", "adjusted", "mixed"))

# plot smp results
ggplot(smp_OR, 
       aes(x = Estimate, y = OR, ymin = LL, ymax = UL)) +
  geom_pointrange(aes(col = Estimate)) + 
  geom_errorbar(aes(ymin = LL, ymax = UL, col = Estimate), width = 0.2) +
  geom_hline(yintercept =1, linetype=2)+
  ylab("Odds ratio") +
  scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0.5, 1.5, 0.25)) +
  coord_flip() +
  scale_colour_manual(values = c("grey41", "darkslategrey", "darkred")) +
  guides(col = guide_legend(reverse = T)) +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20)) +
  theme_classic()

smp.pred <- plot_model(smp.mixed, type = "pred", pred.type = "fe",
           term = "irrigArea [all]", title = " " , 
           #colors = c("darkred", "slategrey"),
           axis.title = c("Area of irrigated land (hectares)",
                          "Probability of S. mansoni infection"))

################################################
# combine Sh_presence and Sm_pressence results into figures/tables
shp_OR <- dplyr::mutate(shp_OR, species = "Sh")
smp_OR <- dplyr::mutate(smp_OR, species = "Sm")
presence_OR <- rbind(shp_OR, smp_OR)
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
  scale_y_continuous(trans = 'log2', limits = c(0.5, 1.5), breaks = seq(0.5, 1.5, 0.25)) +
  coord_flip() +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  #scale_colour_manual(values = c("grey41", "darkslategrey", "darkred")) +
  theme_classic() +
  theme(legend.position = "bottom", 
        text = element_text(size=12), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12))

# save Figure 3a
ggsave("Figure 3a.pdf", plot = presence_plot, device = "pdf",
       path = "./Submission 3 - Lancet Planet Health/Figures",
       width = 107, units = "mm", dpi = 300)
  

# full regression output for all sh and sm presence models 
stargazer(shp.crude, shp.adj, shp.mixed,
          smp.crude, smp.adj, smp.mixed,
          type="html",out="presence_regression_output_7FEB.doc",
          intercept.bottom = F,
          intercept.top = T,
          digits=2)


# diagnostics in performance package
check_collinearity(shp.mixed)
check_convergence(shp.mixed) # FALSE, gradient = 0.01949683
icc(shp.mixed) # adjusted ICC = 0.204; conditional ICC = 0.159
performance_hosmer(shp.mixed) # Chisq = 6.201, df = 8, p = 0.625 - good fit
model_performance(shp.mixed, metrics = "all")

check_collinearity(smp.mixed)
check_convergence(smp.mixed) # FALSE, gradient = 0.01949683
icc(smp.mixed) # adjusted ICC = 0.256; conditional ICC = 0.084
performance_hosmer(smp.mixed) # Chisq = 7.982, df = 8, p = 0.435 - good fit
model_performance(smp.mixed, metrics = "all")

###############################################
#
#                OUTCOME #3
#               Sh count
#

shc.mixed <- glmer.nb(Sh_median ~ irrigArea*locF + DEM01F + DEM04c +
                        factor(headSchoolYrsCat) + factor(nWivesCat) +
                        pumpYN + wealthF + irrigAreaV +
                        (1|village/C_ConcessionId), 
                      data = allData)
control= glmerControl(optimizer="bobyqa", 
                      tolPwrss=1e-3,
                      optCtrl = list(maxfun = 100000))
# doesn't run

# try glmmTMB
shc.mixedTMB <- glmmTMB(round(Sh_median, 0) ~ irrigArea + locF + DEM01F + DEM04c +
                          factor(headSchoolYrsCat) + factor(nWivesCat) +
                          pumpYN + wealthF + irrigAreaV +
                          (1|village/C_ConcessionId),
                        data = allData, se = T, 
                        family = nbinom2(link = "log"))

summary(shc.mixedTMB)
confint(shc.mixedTMB, method = "profile") # runs for several minutes without finishing
waldORCI <- exp(confint(shc.mixedTMB, method = "Wald")) # 1.05 (0.98, 1.13)

# crude and adjusted models
shc.crudeTMB <- glmmTMB(round(Sh_median, 0) ~ irrigArea,
                      data = allData, se = T, 
                      family = nbinom2(link = "log"))

shc.adjTMB <- glmmTMB(round(Sh_median, 0) ~ irrigArea + locF + DEM01F + DEM04c +
                          factor(headSchoolYrsCat) + factor(nWivesCat) +
                          pumpYN + wealthF + irrigAreaV,
                        data = allData, se = T, 
                        family = nbinom2(link = "log"))

# compare crude, adjusted and mixed models for shc
anova(shc.adjTMB, shc.crudeTMB, test = "Chisq")
anova(shc.adjTMB, shc.mixedTMB, test = "Chisq")

# prediction plots
shc.pred <- plot_model(shc.mixedTMB, type = "pred", pred.type = "fe",
                       term = "irrigArea [all]", title = " " , 
                       #colors = c("darkred", "slategrey"),
                       axis.title = c("Area of irrigated land (hectares)",
                                      "Intensity of S. haematobium infection"))


# diagnostics in performance package
check_collinearity(shc.mixedTMB)
icc(shc.mixedTMB) # adjusted ICC = 1.00; conditional ICC = 0.567
model_performance(shc.mixedTMB, metrics = "all")

# residuals from DHARMa package
shc_res = simulateResiduals(shc.mixedTMB)
shcResid <- plot(shc_res)
testResiduals(shc_res)

###############################################)
#
#                OUTCOME #4
#               Sm count
#

smc.mixedTMB <- glmmTMB(round(Sm_median, 0) ~ irrigArea + locF + DEM01F + DEM04c +
                          factor(headSchoolYrsCat) + factor(nWivesCat) +
                          pumpYN + wealthF + irrigAreaV +
                          (1|village),
                        data = allData,
                        family = nbinom2(link = "log"))

# warning: model convergence problem; extreme or very small eigenvalues detected
# model output contains coefficients but no standard errors
# model runs with random intercepts for village only

summary(smc.mixedTMB)  
waldORCI <- exp(confint(smc.mixedTMB, method = "Wald")) # 1.09 (0.89, 1.32)

# crude and adjusted models
smc.adjTMB <- glmmTMB(round(Sm_median, 0) ~ irrigArea + locF + DEM01F + DEM04c +
                          factor(headSchoolYrsCat) + factor(nWivesCat) +
                          pumpYN + wealthF + irrigAreaV,
                        data = allData,
                        family = nbinom2(link = "log"))

smc.crudeTMB <- glmmTMB(round(Sm_median, 0) ~ irrigArea,
                        data = allData,
                        family = nbinom2(link = "log"))

# compare crude, adjusted and mixed models
anova(smc.adjTMB, smc.crudeTMB, test = "Chisq")
anova(smc.adjTMB, smc.mixedTMB, test = "Chisq")

# plot predicted values
smc.pred <- plot_model(smc.mixedTMB, type = "pred", pred.type = "fe",
                       term = "irrigArea [all]", title = " " , 
                       #colors = c("darkred", "slategrey"),
                       axis.title = c("Area of irrigated land (hectares)",
                                      "Intensity of S. mansoni infection"))

# residuals from DHARMa package
smc_res = simulateResiduals(smc.mixedTMB)
smcResid <- plot(smc_res)
testResiduals(smc_res)


# combined plots and tables for both count outcomes
# residuals --> shcResid, smcResid
ggarrange(shcResid, smcResid, labels = c("A", "B"))
# prediction
ggarrange(shc.pred, smc.pred, labels = c("A", "B"))
# regression output for all count models
stargazer(shc.crudeTMB, shc.adjTMB, shc.mixedTMB,
          smc.crudeTMB, smc.adjTMB, smc.mixedTMB,
          type="html", out="count_regression_output_2APR.doc",
          intercept.bottom = F,
          intercept.top = T,
          digits=2) 
# error: unrecognized object type -- stargazer doesn't like tmb objects

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

###############################################
#
#             PLOT INTENSITY RESULTS
#
# import point estimates and confidence intervals generated in stata
intensity_irr <- read.csv("C:/Users/andre/Box Sync/Papers/Ch2IRisk-IrrigatedAg/Analysis/intensity_irr.csv")
intensity_irr$species <- factor(intensity_irr$species, levels = c("Sm", "Sh"))
intensity_irr$Estimate <- factor(intensity_irr$Estimate, levels = c("crude", "adjusted", "mixed"))
intensity_irr <- dplyr::mutate(intensity_irr, outcome = "intensity")

# plot shc results
shc_irr <- dplyr::filter(intensity_irr, species == "Sh")
ggplot(shc_irr, 
       aes(x = Estimate, y = IRR, ymin = LL, ymax = UL)) +
  geom_pointrange(aes(col = Estimate)) + 
  geom_errorbar(aes(ymin = LL, ymax = UL, col = Estimate), width = 0.2) +
  geom_hline(yintercept =1, linetype=2)+
  ylab("Rate ratio") +
  scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0.5, 1.5, 0.25)) +
  #facet_grid(stratum ~ . ) +
  coord_flip() +
  scale_colour_manual(values = c("grey41", "darkslategrey", "darkred")) +
  guides(col = guide_legend(reverse = T)) +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20)) +
  theme_bw()


# plot smc results
smc_irr <- dplyr::filter(intensity_irr, species == "Sm")
ggplot(smc_irr, 
       aes(x = Estimate, y = IRR, ymin = LL, ymax = UL)) +
  geom_pointrange(aes(col = Estimate)) + 
  geom_errorbar(aes(ymin = LL, ymax = UL, col = Estimate), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 2)+
  ylab("Rate ratio") +
  scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0.5, 1.5, 0.25)) +
  #facet_grid(stratum ~ . ) +
  coord_flip() +
  scale_colour_manual(values = c("grey41", "darkslategrey", "darkred")) +
  guides(col = guide_legend(reverse = T)) +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20)) +
  theme_bw()



################################################
# combine Sh_presence and Sm_pressence results into figures/tables

intensity_plot <- ggplot(intensity_irr,
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
  theme_classic() + 
  theme(legend.position = "bottom",
        text = element_text(size=12), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12))


# save Figure 3b
ggsave("Figure 3b.pdf", plot = intensity_plot, device = "pdf",
       path = "./Submission 3 - Lancet Planet Health/Figures",
       width = 107, units = "mm", dpi = 300)


# combine presence and intensity plots
presence_intensity <- ggarrange(presence_plot, intensity_plot, labels = c("A", "B"),
                      common.legend = TRUE, legend = "bottom")

# save Figure 3
ggsave("Figure 3.jpeg", plot = presence_intensity, device = "jpeg",
       path = "./Submission 4 - Infect Dis Poverty/Figures",
       width = 214, units = "mm", dpi = 300)

# combine prediction plots
ggarrange(shp.pred, smp.pred, labels = c("A", "B"))

# plot binned residual plots side-by-side
source("C:/Users/andre/Documents/GitHub/Schisto-Irrigation/multiplot.R")
multiplot(binned_residuals(shp.mixed), 
          binned_residuals(smp.mixed), cols = 2)

# plot random effects dot plot for both shp and smp
ggarrange(dotplot(ranef(shp.mixed, condVar = T, transf = exp)),
          dotplot(ranef(smp.mixed, condVar = T, transf = exp)), 
          labels = c("A", "B"))



