############################0 SETUP################################################################
# 0 Setup

library(tidyverse)
library(gt)
library(gtsummary)
library(mice)
library(tableone)
# library(devtools)
# devtools::install_github("Exeter-Diabetes/EHRBiomarkr")
library(EHRBiomarkr)
library(broom)
library(flextable)
library(survival)
library(survminer)
library(rms)
library(tidyverse)
library(PSweight)
library(patchwork)
library(pROC)
library(riskRegression)
library(grid)
library(gridExtra)
library(dcurves)
library(mgcv)
library(ggpattern)
library(ggsignif)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())

# set random seed
set.seed(123)

# set number of imputations
n.imp <- 10

# set number of quantiles
n.quantiles <- 10

# set number of bootstrap samples
n.bootstrap <- 500

#today <- as.character(Sys.Date(), format="%Y%m%d")
today <- "2024-12-21"


#write function to pool HRs from multiple imputations later on
pool.rubin.HR <- function(COEFS,SE,n.imp){
  mean.coef <- mean(COEFS)
  W <- mean(SE^2)
  B <- var(COEFS)
  T.var <- W + (1+1/n.imp)*B
  se.coef <- sqrt(T.var)
  rm <- (1+1/n.imp)*B/W
  df <- (n.imp - 1)*(1+1/rm)^2
  LB.CI <- mean.coef - (se.coef*1.96)
  UB.CI <- mean.coef + (se.coef*1.96)
  F <- (-mean.coef)^2/T.var
  P <- pf(q=F, df1=1, df2=df, lower.tail = FALSE)
  HR <- exp(mean.coef)
  HR.lower <- exp(LB.CI)
  HR.upper <- exp(UB.CI)
  output <- c(HR, HR.lower, HR.upper, df, F, P)
  names(output) <- c('HR', 'lower bound', 'upper bound', 'df', 'F', 'P')
  return(output)}

# function to pool estimates from multiple imputations further down
pool.rubin.KM <- function(EST,SE,n.imp){
  mean.est <- mean(EST)
  W <- mean(SE^2)
  B <- var(EST)
  T.var <- W + (1+1/n.imp)*B
  se.est <- sqrt(T.var)
  rm <- (1+1/n.imp)*B/W
  df <- (n.imp - 1)*(1+1/rm)^2
  LB.CI <- mean.est - (se.est*1.96)
  UB.CI <- mean.est + (se.est*1.96)
  F <- (-mean.est)^2/T.var
  P <- pf(q=F, df1=1, df2=df, lower.tail = FALSE)
  observed_survival <- 1-mean.est
  lower_ci <- 1-UB.CI
  higher_ci <- 1-LB.CI
  output <- c(observed_survival, lower_ci, higher_ci, df, F, P)
  names(output) <- c('observed survival', 'lower bound', 'upper bound', 'df', 'F', 'P')
  return(output)}

# set vector of covariates used for multivariable adjustments and weighting
covariates <- c("dstartdate_age", "malesex", "imd2015_10", "ethnicity_5cat", "initiation_year", "prebmi", "prehba1c",
                "pretotalcholesterol", "preegfr", "uacr", "presbp", "ncurrtx", "statin", "INS", 
                "ACEi_or_ARB", "smoking_status", "dstartdate_dm_dur_all", "predrug_hypertension", "predrug_af", "hosp_admission_prev_year")

#outcomes to be studied:
outcomes <- c("ckd_egfr40", "ckd_egfr50", "ckd_egfr50_5y", "death", "macroalb", "dka", "amputation", "side_effect")
outcomes_per_drugclass <- c("ckd_egfr50", "ckd_egfr50_pp")

# set default colour-blind accessible colours for figures later on
cols <- c("SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00")
#in further analyses, the dpp4/su group will be combined, and we will use the dpp4 colour for this (strongest contrast)
cols <- c(cols, "DPP4i/SU" = "#0072B2")

# variables to be shown in tables
vars <- c("dstartdate_age", "malesex", "ethnicity_5cat", "imd2015_10",             # sociodemographic variables
          "prebmi", "presbp", "predbp", "pretotalcholesterol", "prehdl", "preldl", # vital signs and laboratory measurements
          "pretriglyceride", "prehba1c",  "preegfr",
          "uacr", "albuminuria_unconfirmed", "albuminuria",
          "dstartdate_dm_dur_all", "smoking_status", "predrug_hypertension",   # comorbidities
          "predrug_af", "predrug_dka", "genital_infection", "hosp_admission_prev_year",
          "initiation_year",
          "ncurrtx", "statin", "INS", "ACEi_or_ARB"
          
)

#categorical variables
factors <- c("malesex", "ethnicity_5cat", "imd2015_10", "albuminuria_unconfirmed", "albuminuria", 
             "smoking_status", "predrug_hypertension",
             "predrug_af", "predrug_dka", "genital_infection", "hosp_admission_prev_year",
             "initiation_year",
             "ncurrtx", "statin", "INS", "ACEi_or_ARB")

nonnormal <- c("uacr", "dstartdate_dm_dur_all")