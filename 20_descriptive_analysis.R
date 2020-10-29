###############################################################
# 
#    descriptive analysis preparing for regression
#
#    (1) univariate analyses of outcome variables and covariates
#    (2) bivariate analysis for screening pre-specified covariates
#    (3) mean/SD and n/% table for screened covariates
#    (4) summary of missingness for all screened covariates
#
#    coded by: Andrea Lund
#    last updated: September 18, 2019
#
###############################################################

# install packages
library(tidyverse)
library(plyr)
library(PerformanceAnalytics) # for correlation visualization
library(sjPlot) # custom correlation matrices
library(Hmisc) # for correlation matrices
library(stargazer) # for generating descriptive tables
library(ggplot2)
library(ggpubr) # for combining plots 
library(GGally)

# set working directory
setwd("~/GitHub/Schisto-Irrigation")

# source script that generates merged data sets
source("10 merging datasets.R")

# import data
allData <- read.csv("allData.csv")
shData <- read.csv("shData.csv")
smData <- read.csv("smData.csv")

# function
neat.table <- function(x, name){
  xx <- data.frame(x)
  names(xx) <- c("Value", "Count")
  xx$Percent <- with(xx, round(Count/sum(Count)*100, digits = 1))
  data.frame(Variable = name, xx)
}
#############################################################
#
#               univariate analyses
#         primary and secondary outcome data
#

### primary outcomes
# Sh presence 
table(allData$Sh_presence) # 0: n = 427; 1: n = 805
# Sm presence
table(allData$Sm_presence) # 0: n = 1014; 1: n = 208

# Sh mean
summary(allData$Sh_mean) # min = 0, mean = 35.29, med = 2.50, max = 1745, n = 140 missing
hist(allData$Sh_mean, breaks = 1000)

# Sm mean
summary(allData$Sm_mean) # min = 0, mean = 34.4, med = 0.0, max = 4636, n = 150 missing
hist(allData$Sm_mean, breaks = 1000)

# coinfection
summary(allData$coinf) # n = 145 missing
table(allData$coinf) # 0: n = 1060; 1: n = 167

# secondary outcomes
summary(allData$Sh_cat)
#    0    1    2 NA's 
#  427  580  225  140 

summary(allData$Sm_cat)
#    0    2    3    4 NA's 
# 1014  131   54   23  150

### subset primary and secondary outcome data by data type for descriptive tables

## categorical outcome variables
outcomeCat <- allData %>%
  # transform categorical outcomes into factors
  mutate(Sh_presenceF = as.factor(Sh_presence), 
         Sm_presenceF = as.factor(Sm_presence),
         Sh_catF = as.factor(Sh_cat), Sm_catF = as.factor(Sm_cat),
         coinfF = as.factor(coinf), nInfF = as.factor(nInf)) %>% 
  # subset data to categorical outcome variables
  select(locF, Sh_presenceF:Sm_presenceF, Sh_catF:Sm_catF, coinfF, nInfF) 

# code for categorical descriptive table adapted from: https://stackoverflow.com/questions/15111629/create-summary-table-of-categorical-variables-of-different-lengths
# overall 
x <- lapply(outcomeCat[,2:7], table) 
catTabO <- do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))
# river
outcomeCatR <- outcomeCat[which(outcomeCat$locF == "river"),]
y <- lapply(outcomeCatR[,2:7], table)
catTabR <- do.call(rbind, lapply(seq_along(y), function(i)neat.table(y[i], names(y[i]))))
# lake
outcomeCatL <- outcomeCat[which(outcomeCat$locF == "lake"),]
z <- lapply(outcomeCatL[,2:7], table)
catTabL <- do.call(rbind, lapply(seq_along(z), function(i)neat.table(z[i], names(z[i]))))

# create and export final table with overall and stratified summaries
catTab <- bind_cols(catTabO, catTabR[,3:4], catTabL[,3:4])
write.csv(catTab, "catTabORL.csv")

## numeric outcome variables
# overall
numTabO <- allData %>% 
  filter(!is.na(Sh_mean) & !is.na(Sm_mean)) %>% 
  select(Sh_mean, Sm_mean) %>%
  rename("shc" = "Sh_mean", "smc" = "Sm_mean") %>% 
  summarise_each(funs(mean, sd, 
                      min, median, 
                      q1 = quantile(., 0.25),
                      q3 = quantile(., 0.75), max)) %>% 
  gather(stat, val) %>% 
  separate(stat, into = c("var", "stat"), sep = "_") %>% 
  spread(var, val) %>% 
  arrange(match(stat, c("mean", "sd", "min", "q1", "median", "q3", "max")))

# river
numTabR <- allData %>% 
  filter(!is.na(Sh_mean) & !is.na(Sm_mean) & locF == "river") %>% 
  select(Sh_mean, Sm_mean) %>%
  rename("shc" = "Sh_mean", "smc" = "Sm_mean") %>% 
  summarise_each(funs(mean, sd, 
                      min, median, 
                      q1 = quantile(., 0.25),
                      q3 = quantile(., 0.75), max)) %>% 
  gather(stat, val) %>% 
  separate(stat, into = c("var", "stat"), sep = "_") %>% 
  spread(var, val) %>% 
  arrange(match(stat, c("mean", "sd", "min", "q1", "median", "q3", "max")))

# lake
numTabL <- allData %>% 
  filter(!is.na(Sh_mean) & !is.na(Sm_mean) & locF == "lake") %>% 
  select(Sh_mean, Sm_mean) %>%
  rename("shc" = "Sh_mean", "smc" = "Sm_mean") %>% 
  summarise_each(funs(mean, sd, 
                      min, median, 
                      q1 = quantile(., 0.25),
                      q3 = quantile(., 0.75), max)) %>% 
  gather(stat, val) %>% 
  separate(stat, into = c("var", "stat"), sep = "_") %>% 
  spread(var, val) %>% 
  arrange(match(stat, c("mean", "sd", "min", "q1", "median", "q3", "max")))

# create and export final table with overall and stratified summaries
numTab <- bind_cols(numTabO, numTabR[,2:3], numTabL[,2:3])
write.csv(numTab, "numTabORL.csv")

### TO DO ### 
#   - two-panel ggplot histogram of Sh and Sm count distribution
#       - faceted by location

# Sh count distribution

mu <- plyr::ddply(shData, "locF", summarise, grp.mean = mean(Sh_mean))
Sh_hist <- ggplot(shData, aes(x = Sh_median, col = locF)) +
  geom_histogram(binwidth = 1) +
  facet_grid(locF ~ .) +
  geom_vline(data = mu, aes(xintercept = grp.mean),
             linetype = 2) +
  labs(x = "Median intensity of S. haematobium infection (eggs/10 mL)",
       y = "Number of school-aged children") + 
  scale_x_continuous(limits = c(0, 1800), breaks = seq(0, 1800, 200)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  scale_color_manual(values = c("darkred", "darkslategrey")) +
  theme_bw()

mu <- plyr::ddply(smData, "locF", summarise, grp.mean = mean(Sm_mean))
Sm_hist <- ggplot(smData, aes(x = Sm_median, col = locF)) +
  geom_histogram(binwidth = 1) +
  facet_grid(locF ~ .) +
  geom_vline(data = mu, aes(xintercept = grp.mean),
             linetype = 2) +
  labs(x = "Median intensity of S. mansoni infection (epg)",
       y = "Number of school-aged children") + 
  #scale_x_continuous(limits = c(-1, 1000)) + # n = 8 observations beyond this
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, 100)) +
  scale_color_manual(values = c("darkred", "darkslategrey")) +
  theme_bw()

# geom_errorbar(aes(ymin = LL, ymax = UL, col = Estimate), width = 0.25) +
# geom_hline(yintercept =1, linetype=2)+
# ylab("Odds ratio") +
# facet_grid(strata ~ ., scales = "free_y" ) +
#  scale_y_continuous(limits = c(0.5, 1.5), breaks = seq(0.5, 1.5, 0.25)) +
#  coord_flip() +
#  scale_colour_manual(values = c("grey41", "darkslategrey", "darkred")) +
#  guides(col = guide_legend(reverse = T)) +
#  theme(axis.text.y = element_blank()) +
#  theme_bw()
  
ggarrange(Sh_hist, Sm_hist, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")


#############################################################
#
#               univariate analyses
#         exposure and pre-specified covariates
#

# individual-level
# age
allData %>% 
  group_by(village) %>% 
  summarise(mean_age = mean(DEM04), sd_age = sd(DEM04))

# sex 
 count(allData, DEM01)

#     DEM01     n
#      <dbl> <int>
# 1     1   707
# 2     2   665

sexByVillage <- allData %>% 
  group_by(village) %>% 
  count(DEM01) %>% 
  spread(DEM01, n) %>% 
  rename("male" = "1", "female" = "2")

#### generate descriptive tables for sharing
## numeric covariates
# overall
numCovarsO <- allData %>%
  select(irrigArea, DEM04, headSchoolYrs, nWives, nFishAny, nWatpts,
         nWatPtsAgKids, pump, marketDistance_km, distanceWatPt) %>%
  rename("age" = "DEM04", "marketDistanceK" = "marketDistance_km") %>% 
  summarise_each(funs(mean, sd,
                      min, median, 
                      q1 = quantile(., 0.25),
                      q3 = quantile(., 0.75), max)) %>%
  gather(stat, val) %>% 
  separate(stat, into = c("var", "stat"), sep = "_") %>% 
  spread(stat, val) %>% 
  select(var, mean, sd) %>% 
  arrange(match(var, c("age", "irrigArea", "pump", "nFishAny", "headSchoolYrs",
                       "nWives", "nWatpts", "nWatPtsAgKids", "distanceWatPt", "marketDistanceK")))

# river
numCovarsR <- allData %>%
  filter(locF == "river") %>% 
  select(irrigArea, DEM04, headSchoolYrs, nWives, nFishAny, nWatpts,
         nWatPtsAgKids, pump, marketDistance_km, distanceWatPt) %>%
  rename("age" = "DEM04", "marketDistanceK" = "marketDistance_km") %>% 
  summarise_each(funs(mean, sd,
                      min, median, 
                      q1 = quantile(., 0.25),
                      q3 = quantile(., 0.75), max)) %>%
  gather(stat, val) %>% 
  separate(stat, into = c("var", "stat"), sep = "_") %>% 
  spread(stat, val) %>% 
  select(var, mean, sd) %>% 
  arrange(match(var, c("age", "irrigArea", "pump", "nFishAny", "headSchoolYrs",
                       "nWives", "nWatpts", "nWatPtsAgKids", "distanceWatPt", "marketDistanceK")))

# lake
numCovarsL <- allData %>%
  filter(locF == "lake") %>% 
  select(irrigArea, DEM04, headSchoolYrs, nWives, nFishAny, nWatpts,
         nWatPtsAgKids, pump, marketDistance_km, distanceWatPt) %>%
  rename("age" = "DEM04", "marketDistanceK" = "marketDistance_km") %>% 
  summarise_each(funs(mean, sd,
                      min, median, 
                      q1 = quantile(., 0.25),
                      q3 = quantile(., 0.75), max)) %>%
  gather(stat, val) %>% 
  separate(stat, into = c("var", "stat"), sep = "_") %>% 
  spread(stat, val) %>% 
  select(var, mean, sd) %>% 
  arrange(match(var, c("age", "irrigArea", "pump", "nFishAny", "headSchoolYrs",
                       "nWives", "nWatpts", "nWatPtsAgKids", "distanceWatPt", "marketDistanceK")))
                   
# create and export final table with overall and stratified summaries
numCovars <- bind_cols(numCovarsO, numCovarsR[,2:3], numCovarsL[,2:3])
write.csv(numCovars, "numCovars.csv")

## categorical covariates
catCovars <- allData %>% 
  select(locF, DEM01F, eth_modeF, kidsAnyTaskF, wealthF, locF)

# overall
x <- lapply(catCovars[,2:5], table) 
covarsTabO <- do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

# river
catCovarsR <- catCovars[which(catCovars$locF == "river"),]
y <- lapply(catCovarsR[,2:5], table)
covarsTabR <- do.call(rbind, lapply(seq_along(y), function(i)neat.table(y[i], names(y[i]))))

# lake
catCovarsL <- catCovars[which(catCovars$locF == "lake"),]
z <- lapply(catCovarsL[,2:5], table)
covarsTabL <- do.call(rbind, lapply(seq_along(z), function(i)neat.table(z[i], names(z[i]))))

# create and export final table with overall and stratified summaries
catCovarsTab <- bind_cols(covarsTabO, covarsTabR[,3:4], covarsTabL[,3:4])
write.csv(catCovarsTab, "catCovarsTab.csv")

#############################################################
#
#                     bivariate analyses
#          primary outcomes and pre-specified covariates 
#

### visualize correlations between primary outcomes and pre-specified covariates
## outcomes: Sh_presence, Sh_mean, Sm_presence, Sm_mean
## individual-level covariates: DEM01 (sex), DEM04 (age)
ShP_ind <- allData %>% 
  select(Sh_presence, DEM01, DEM04) %>% 
  rename("Age" = "DEM04", "Sex" = "DEM01")
chart.Correlation(ShP_ind) # what are the default significance thresholds for the asterisks?
# stars are derived from cor.test
#cor.test(allData$Sh_presence, allData$DEM01) # one star: p-value = 0.01112
#cor.test(allData$Sh_presence, allData$DEM04) # three stars: p-value = 0.0002539
#cor.test(allData$DEM01, allData$DEM04) # no stars: p-value = 0.1256

ShC_ind <- allData %>% 
  select(Sh_mean, DEM01, DEM04) %>% 
  rename("Age" = "DEM04", "Sex" = "DEM01")
chart.Correlation(ShC_ind)

cor.test(allData$Sh_mean, allData$DEM04) # one dot: p-value = 0.08454

SmP_ind <- allData %>% 
  select(Sm_presence, DEM01, DEM04) %>% 
  rename("Age" = "DEM04", "Sex" = "DEM01")
chart.Correlation(SmP_ind)

SmC_ind <- allData %>% 
  select(Sm_mean, DEM01, DEM04) %>% 
  rename("Age" = "DEM04", "Sex" = "DEM01")
chart.Correlation(SmC_ind)

## primary exposure variable: irrigArea 
ShP_exp <- allData %>% 
  select(Sh_presence, irrigArea, pumpYN)
chart.Correlation(ShP_exp)

cor.test(allData$Sh_presence, allData$irrigArea) # two stars: p-value = 0.001208

ShC_exp <- allData %>% 
  select(Sh_mean, irrigArea, pumpYN)
chart.Correlation(ShC_exp)

SmP_exp <- allData %>% 
  select(Sm_presence, irrigArea, pumpYN)
chart.Correlation(SmP_exp)

SmC_exp <- allData %>% 
  select(Sm_mean, irrigArea, pumpYN)
chart.Correlation(SmC_exp)

## household-level socio-demographic covariates: head_age, headSchoolYrs, eth_mode (non-num()), nWives, wealthQuintile1, pump
ShP_dem <- allData %>% 
  select(Sh_presence, head_age, hhsize, headSchoolYrs, nWives, wealthQuintile1, pumpYN, marketDistance_km)
chart.Correlation(ShP_dem)

ShC_dem <- allData %>% 
  select(Sh_mean, head_age, hhsize, headSchoolYrs, nWives, wealthQuintile1, pumpYN, marketDistance_km)
chart.Correlation(ShC_dem)

cor.test(allData$Sh_mean, allData$nWives) # t = 0.65577, df = 1230, p-value = 0.5121

SmP_dem <- allData %>% 
  select(Sm_presence, head_age, hhsize, headSchoolYrs, nWives, wealthQuintile1, pumpYN, marketDistance_km)
chart.Correlation(SmP_dem)

cor.test(allData$Sm_presence, allData$nWives) # t = -0.011149, df = 1220, p-value = 0.9911

SmC_dem <- allData %>% 
  select(Sm_mean, head_age, hhsize, headSchoolYrs, nWives, wealthQuintile1, pumpYN, marketDistance_km) 
chart.Correlation(SmC_dem)

## household-level water-related covariates: nFishAny, nWatpts, surface, distanceWatPt, nWatPtsAgKids, kidsAnyTask 
ShP_wat <- allData %>% 
  select(Sh_presence, nFishAny, nWatpts, surface, distanceWatPt, nWatPtsAgKids, kidsAnyTask, pumpYN)
chart.Correlation(ShP_wat)

ShC_wat <- allData %>% 
  select(Sh_mean, nFishAny, nWatpts, surface, distanceWatPt, nWatPtsAgKids, kidsAnyTask, pumpYN)
chart.Correlation(ShC_wat)

SmP_wat <- allData %>% 
  select(Sm_presence, nFishAny, nWatpts, surface, distanceWatPt, nWatPtsAgKids, kidsAnyTask, pumpYN)
chart.Correlation(SmP_wat)

SmC_wat <- allData %>% 
  select(Sm_mean, nFishAny, nWatpts, surface, distanceWatPt, nWatPtsAgKids, kidsAnyTask, pumpYN)
chart.Correlation(SmC_wat)

# bivariate measures of association between primary outcomes and categorical covariates
# non-ordered categorical variables: sex (DEM01), ethnicity (eth_mode), river-lake geography (loc)

# chi-square test for presence/absence outcomes
# Sh_presence
chisq.test(x = allData$Sh_presence, y = allData$DEM01F) # X-squared = 6.3404, df = 1, p-value = 0.0118
chisq.test(x = allData$Sh_presence, y = allData$eth_mode3, # sparse counts generate warning of incorrect approximation
           simulate.p.value = T) 
# chi-squared result: X-squared = 169.22, df = 8, p-value < 2.2e-16
# chi-squared with simulated p-values: X-squared = 169.22, df = NA, p-value = 0.0004998
# fishers exact test
fisher.test(allData$Sh_presence, allData$eth_mode, simulate.p.value = T) # p-value = 0.0004998

chisq.test(x = allData$Sh_presence, y = allData$locF) # X-squared = 204.02, df = 1, p-value < 2.2e-16

# Sm_presence
chisq.test(x = allData$Sm_presence, y = allData$DEM01) # X-squared = 1.5092, df = 1, p-value = 0.2193
chisq.test(x = allData$Sm_presence, y = allData$eth_mode, # sparse counts generate warning of incorrect approximation
           simulate.p.value = T)
# chi-squared result: X-squared = 31.6, df = 8, p-value = 0.0001098
# chi-squared result with simulated p-values: X-squared = 31.6, df = NA, p-value = 0.0009995
# fishers exact test
fisher.test(allData$Sm_presence, y = allData$eth_mode, simulate.p.value = T) # p-value = 0.0004998
chisq.test(x = allData$Sm_presence, y = allData$locF) # X-squared = 10.614, df = 1, p-value = 0.001122

# one-way ANOVA/bivariate linear regression for count/intensity outcomes
# Sh_mean
ShSex_aov <- aov(Sh_mean ~ DEM01F, data = allData)
summary(ShSex_aov)
#               Df  Sum Sq Mean Sq F value Pr(>F)
# DEM01F         1     443     443   0.059  0.808
# Residuals   1246 9401263    7545               
# 148 observations deleted due to missingness

ShEth_aov <- aov(Sh_mean ~ eth_mode3, data = allData)
summary(ShEth_aov)
#              Df  Sum Sq Mean Sq F value  Pr(>F)    
# eth_mode3      5  221458   44292   6.002 1.7e-05 ***
# Residuals   1226 9047336    7380                    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#  164 observations deleted due to missingness                  

ShLoc_aov <- aov(Sh_mean ~ locF, data = allData)
summary(ShLoc_aov)
#               Df  Sum Sq Mean Sq F value Pr(>F)    
# locF           1  529366  529366   74.34 <2e-16 ***
# Residuals   1246 8872340    7121                   
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#  148 observations deleted due to missingness

# Sm_mean
SmSex_aov <- aov(Sm_mean ~ DEM01F, data = allData)
summary(SmSex_aov)
#               Df   Sum Sq Mean Sq F value Pr(>F)
# DEM01F         1    32566   32566   0.708    0.4
# Residuals   1236 56886900   46025               
# 158 observations deleted due to missingness
# aov(Sh_mean ~ eth_mode3, data = allData)

SmEth_aov <- aov(Sm_mean ~ eth_mode3, data = allData)
summary(SmEth_aov)
#              Df   Sum Sq Mean Sq F value Pr(>F)
# eth_mode3      5   366221   73244   1.579  0.163
# Residuals   1216 56417093   46396               
# 174 observations deleted due to missingness
SmEth_lm <- lm(Sm_mean ~ eth_mode3, data = allData)
summary(SmEth_lm)

SmLoc_aov <- aov(Sm_mean ~ locF, data = allData)
summary(SmLoc_aov)
#               Df   Sum Sq Mean Sq F value Pr(>F)   
# locF           1   337734  337734   7.378 0.0067 **
#  Residuals   1236 56581732   45778                  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#  158 observations deleted due to missingness

### correlation matrix with p-values for selecting covariates
## outcomes: Sh_presence, Sh_mean, Sm_presence, Sm_mean
## individual-level covariates: DEM01 (sex), DEM04 (age)
## primary exposure variable: irrigArea 
## household-level socio-demographic covariates: 
#     head_age, headSchoolYrs, eth_mode (non-num()), nWives, wealthQuintile1, pump
## household-level water-related covariates: 
#     nFishAny, nWatpts, surface, distanceWatPt, nWatPtsAgKids, kidsAnyTask 

## UPSTREAM DEMOGRAPHIC VARIABLES
# Sh presence and upstream demographic variables
# headSchoolYrs, nWives, wealth, pump, marketDistance, kidsAnyTask
GGally::ggpairs(shData, columns = c(25, 58, 129, 46, 49, 108, 183, 164),
                ggplot2::aes(col = locF, position = "dodge", alpha = 0.5),
                columnLabels = c("Infect", "irrig", "pump", "headEd", "nWives",
                                 "kidsAg", "market", "wealth"))

# Sm presence and upstream demographic variables
# headSchoolYrs, nWives, wealth, pump, marketDistance, kidsAnyTask
GGally::ggpairs(smData, columns = c(26, 58, 129, 46, 49, 108, 183, 164),
                ggplot2::aes(col = locF, position = "dodge", alpha = 0.5),
                columnLabels = c("Infect", "irrig", "pump", "headEd", "nWives",
                                 "kidsAg", "market", "wealth"))

# Sh intensity and upstream demographic variables
# headSchoolYrs, nWives, wealth, pump, marketDistance, kidsAnyTask
GGally::ggpairs(shData, columns = c(28, 58, 129, 46, 49, 108, 183, 164),
                ggplot2::aes(col = locF, position = "dodge", alpha = 0.5),
                columnLabels = c("Infect", "irrig", "pump", "headEd", "nWives",
                                 "kidsAg", "market", "wealth"))


# Sm intensity and upstream demographic variables
# headSchoolYrs, nWives, wealth, pump, marketDistance, kidsAnyTask
GGally::ggpairs(smData, columns = c(29, 58, 129, 46, 49, 108, 183, 164),
                ggplot2::aes(col = locF, position = "dodge", alpha = 0.5),
                columnLabels = c("Infect", "irrig", "pump", "headEd", "nWives",
                                 "kidsAg", "market", "wealth"))

# WATER-RELATED VARIABLES
# Sh presence and water-related variables
# nFishAny, nWatpts, nWatPtsAgKids, kidsAnyTask, distanceWatPt 
GGally::ggpairs(shData, columns = c(25, 58, 129, 55, 102, 105, 108, 179),
                ggplot2::aes(col = locF, position = "dodge", alpha = 0.5),
                columnLabels = c("Infect", "irrig", "pump", "fishers", "nWatPts", 
                                 "nWatPtsAgkids", "kidsAg", "distanceWatpt"))

# Sm presence and water-related variables
# nFishAny, nWatpts, nWatPtsAgKids, kidsAnyTask, distanceWatPt 
GGally::ggpairs(smData, columns = c(26, 58, 129, 55, 102, 105, 108, 179),
                ggplot2::aes(col = locF, position = "dodge", alpha = 0.5),
                columnLabels = c("Infect", "irrig", "pump", "fishers", "nWatPts", 
                                 "nWatPtsAgkids", "kidsAg", "distanceWatpt"))


# all regression variables saved to a data frame
regressionVars = allData[,c(2,6,15,29:31,34,47:49,57,60,101,104,
                            107:108,128,131,163,178,180)]

# generate correlation matrix wtih p-values for all regression variables , saved to a word file
sjt.corr(regressionVars, show.p = T, p.numeric = T,
         file = "bivariate correlation - regression vars.doc")

# save post-screening regression variables to a data frame
screenedVars <- allData[,c(60, 2, 6, 48, 51, 57, 101, 104, 107:108, 128, 163, 178, 180 )]
# generate correlation plot for all regression variables
# any variables not meeting significance threshold are omitted
rvCor <- rcorr(as.matrix(screenedVars))
M <- cor(screenedVars)
pv <- cor.mtest(screenedVars, conf.level = 0.8)
corrplot(rvCor$r, is.corr = T, p.mat = rvCor$P, insig = "label_sig",  sig.level = 0.2,
         #         sig.level = c(0.01, 0.05, 0.1, 0.2), pch = c("***", "**", "*", "."), pch.cex = 4,
         method = "shade", type = "upper", diag = F, tl.col = "black")

corrplot(rvCor$r, is.corr = T, p.mat = rvCor$P, insig = "label_sig", sig.level = 0.05,
         #         sig.level = c(0.01, 0.05, 0.1, 0.2), pch = c("***", "**", "*", "."), pch.cex = 4,
         method = "shade", type = "upper", diag = F, tl.col = "black")

corrplot(rvCor$r, is.corr = T, p.mat = rvCor$P, insig = "label_sig", sig.level = 0.01,
         #         sig.level = c(0.01, 0.05, 0.1, 0.2), pch = c("***", "**", "*", "."), pch.cex = 4,
         method = "shade", type = "upper", diag = F, tl.col = "black")


sjt.corr(screenedVars, show.p = T, p.numeric = T,
         file = "bivariate correlation - screened vars.doc")
#############################################################
#                assessment of missing data
#
# counts of missing observations for all columns
missingData <- colSums((is.na(allData)))

# subsets of data missing for each outcome
missingShp <- allData[which(is.na(allData$Sh_presence)),]
missingShm <- allData[which(is.na(allData$Sh_mean)),]
missingSmp <- allData[which(is.na(allData$Sm_presence)),]
missingSmm <- allData[which(is.na(allData$Sm_mean)),]

# which villages are missing outcome data?
naShpV <- table(missingShp$village)
# DF DT FS GG GK LR MA MB MD ME MG MO MT NE NM ST 
#  7 37  6  8 12  6  5  9  5  2 10 13  4  8  2 12
naShmV <- table(missingShm$village)
# DF DT FS GG GK LR MA MB MD ME MG MO MT NE NM ST 
#  7 37  6  8 12  6  5  9  5  2 10 13  4  8  2 12 
naSmpV <- table(missingSmp$village)
# DF DT FS GG GK LR MA MB MD ME MG MO MT NE NM ST 
#  7 40  6  7 12  6  5  9  5  2 14 15  5  8  2 13
naSmmV <- table(missingSmm$village)
# DF DT FS GG GK LR MA MB MD ME MG MO MT NE NM ST 
#  7 40  6  7 12  6  5  9  5  2 14 15  5  8  2 13

# what is the age distribution of children missing outcome data?
hist(missingShp$DEM04)
hist(missingShm$DEM04)
hist(missingSmp$DEM04)
hist(missingSmm$DEM04)

# extract character variables plus numeric concession ID key
charData <- allData %>% 
  select_if(is.character) 

# calculate geometric means of infection intensity overall and by location strata
allData$Sh_gMean <- ifelse(is.na(allData$Sh_count1) & is.na(allData$Sh_count2), NA,
                           ifelse(!is.na(allData$Sh_count1) & !is.na(allData$Sh_count2), 
                                  exp(((log(allData$Sh_count1 + 1) + log(allData$Sh_count2+ 1)))/2) - 1, 
                                  ifelse(is.na(allData$Sh_count1), 
                                         exp(log(allData$Sh_count2 + 1)/1)-1, 
                                         exp(log(allData$Sh_count1 + 1)/1)-1)))

gMean = exp(mean(log(allData$Sh_mean+1), na.rm = T) - 1)

gMean_all <- allData %>% 
  mutate(Sh_trans = log(Sh_mean + 1),
         Sm_trans = log(Sm_mean + 1)) %>% 
  summarise(Sh_gMean = (exp(mean(Sh_trans, na.rm = T))-1),
            Sm_gMean = (exp(mean(Sm_trans, na.rm = T))-1),
            Sh_sd  = sd(Sh_mean, na.rm = T), 
            Sm_sd = sd(Sm_mean, na.rm = T)) %>% 
  mutate(locF = "overall")

gMean_loc <- allData %>%
  group_by(locF) %>% 
  mutate(Sh_trans = log(Sh_mean + 1),
         Sm_trans = log(Sm_mean + 1)) %>% 
  summarise(Sh_gMean = (exp(mean(Sh_trans, na.rm = T))-1),
            Sm_gMean = (exp(mean(Sm_trans, na.rm = T))-1),
            Sh_sd  = sd(Sh_mean, na.rm = T), 
            Sm_sd = sd(Sm_mean, na.rm = T))

gMeans <- rbind(gMean_all, gMean_loc)

