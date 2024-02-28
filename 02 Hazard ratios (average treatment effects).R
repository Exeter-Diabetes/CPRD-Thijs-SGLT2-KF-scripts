## In this script our aim is to see whether the average treatment for SGLT2i vs SU is similar
## in the literature and in this CPRD cohort
## we will compare hazard ratios for SGLT2i vs SU in CPRD to those from trials.

## for reference, the trial meta-analysis HR reported in Lancet. 2022 Nov 19; 400(10365): 1788–1801 is:
# HR 0.62 (0.56 - 0.68)

## Here we will be calculating HRs for SGLT2i and DPP4i vs SU for:
## CKD: onset of CKD stage 3a-5
## CKD: fall in eGFR of <=40% from baseline or onset of CKD stage 5
## all-cause mortality

## originally, there were plans to add other outcomes, e.g. 
## prescription of medications for advanced kidney disease (e.g. injections for renal anaemia, phosphate binders).
## however, due to the low number of subjects progressing to advanced CKD it is unfeasible to look into this.

## Contents:
# 0 setup
# 1 calculate weights (IPTW, overlap weights)
# 2 calculate hazard ratios
# 3 store and display hazard ratios
# 4 sensitivity analyses for subgroups

############################0 SETUP################################################################
# 0 Setup

library(broom)
library(flextable)
library(survival)
library(survminer)
library(rms)
library(tidyverse)
library(PSweight)
library(patchwork)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-02-26_t2d_ckdpc_imputed_data_incl_egfr_below_60.Rda")

n.imp <- 10
set.seed(123)
today <- as.character(Sys.Date(), format="%Y%m%d")

#create variable studydrug2 where we compare SGLT2i vs combined group of DPP4i/SU
temp <- temp %>% mutate(
  studydrug = ifelse(studydrug == "SGLT2", "SGLT2i", 
                     ifelse(studydrug == "DPP4", "DPP4i",
                            "SU")),
  studydrug2 = ifelse(!studydrug == "SGLT2i", "DPP4i/SU", "SGLT2i"))

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

############################1 CALCULATE WEIGHTS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in them, then stack them again at the end
temp$IPTW <- temp$overlap <- 
temp$IPTW2 <- temp$overlap2 <- NA

# take non-imputed dataset to append the imputed datasets to later on
temp2 <- temp[temp$.imp == 0,]

# propensity score formula
ps.formula <- formula(paste("studydrug ~ 

                               #sociodemographic characteristics:
                               dstartdate_age + malesex + imd2015_10 + 
                               ethnicity_5cat + 
                               
                               #laboratory and vital signs measurements:
                               preweight + prebmi + prehba1c + pretriglyceride + prehdl + preldl +
                               pretotalcholesterol + prealt + preegfr + uacr + presbp + predbp +
                               ckdpc_40egfr_score +
                               
                               #medications:
                               ncurrtx + MFN + statin + INS +
                               ACEi_or_ARB + BB + CCB + ThZD + loopD + MRA + 
                               steroids + immunosuppr +
                               
                               #Comorbidities
                               qrisk2_smoking_cat + dstartdate_dm_dur_all + 
                               predrug_hypertension + osteoporosis + predrug_dka +
                               predrug_falls + predrug_urinary_frequency + predrug_volume_depletion +
                               predrug_acutepancreatitis + predrug_micturition_control + 
                               predrug_dementia + hosp_admission_prev_year"))

ps.formula2 <- formula(paste("studydrug2 ~ 

                               #sociodemographic characteristics:
                               dstartdate_age + malesex + imd2015_10 + 
                               ethnicity_5cat + 
                               
                               #laboratory and vital signs measurements:
                               preweight + prebmi + prehba1c + pretriglyceride + prehdl + preldl +
                               pretotalcholesterol + prealt + preegfr + uacr + presbp + predbp + 
                               ckdpc_40egfr_score +
                               
                               #medications:
                               ncurrtx + MFN + statin + INS +
                               ACEi_or_ARB + BB + CCB + ThZD + loopD + MRA + 
                               steroids + immunosuppr +
                               
                               #Comorbidities
                               qrisk2_smoking_cat + dstartdate_dm_dur_all + 
                               predrug_hypertension + osteoporosis + predrug_dka +
                               predrug_falls + predrug_urinary_frequency + predrug_volume_depletion +
                               predrug_acutepancreatitis + predrug_micturition_control + 
                               predrug_dementia + hosp_admission_prev_year"))


#calculate weights in each imputed dataset
for (i in 1:n.imp) {
  
  print(paste0("Calculating weights for imputed dataset number ", i))
  
  imp.x <- temp[temp$.imp == i,]
  
  w.overlap <- SumStat(ps.formula=ps.formula, 
                       data = imp.x, 
                       weight = c("IPW", "overlap"))
  
  w.overlap2 <- SumStat(ps.formula=ps.formula2, 
                       data = imp.x, 
                       weight = c("IPW", "overlap"))
  
  ## Add weights to data frame
  weights <- w.overlap$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
  imp.x$IPTW <- weights$IPW
  imp.x$overlap <- weights$overlap
  
  weights2 <- w.overlap2$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
  imp.x$IPTW2 <- weights2$IPW
  imp.x$overlap2 <- weights2$overlap
  
  # append imputed dataset with weights to dataframe with the original dataset + appended imputed datasets
  temp2 <- rbind(temp2, imp.x)
  
  rm(imp.x)
  rm(weights)
  
}

temp <- temp2
rm(temp2)

#plot SMD in weighted / unweighted populations
plot(w.overlap)
#density plot of propensity scores
plot(w.overlap, type = "density")

summary(w.overlap, weighted.var = TRUE, metric = "ASD")

# save dataset with weights so this can be used in subsequent scripts
cohort <- temp
rm(temp)
save(cohort, file=paste0(today, "_t2d_ckdpc_imputed_data_withweights_incl_egfr_below_60.Rda"))
#load("2024-02-26_t2d_ckdpc_imputed_data_withweights_incl_egfr_below_60.Rda")
############################2 CALCULATE HAZARD RATIOS################################################################

## 2 calculate hazard ratios (unadjusted, adjusted, weighted) and n events per study drug

#outcomes to be studied:
kf_outcomes <- c(#"ckd_345", 
                 "ckd_egfr40", 
                 #"death", 
                 #"ckd_345_pp", 
                 "ckd_egfr40_pp"
                 #, "death_pp"
                 )

kf_key_outcomes <- c(#"ckd_345", 
  "ckd_egfr40")

#create empty data frame to which we can append the hazard ratios once calculated
all_sglt2i_hrs <- 
  all_dpp4i_hrs <-
  all_SGLT2ivsDPP4i_hrs <- 
  all_hrs <- 
  data.frame()

for (k in kf_outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  # calculate number of subjects in each group
  count <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug, events) %>%
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug"))
  
  f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug + 
                               
                               #sociodemographic characteristics:
                               dstartdate_age + malesex + imd2015_10 + 
                               ethnicity_5cat + initiation_year +
                               
                               #laboratory and vital signs measurements:
                               preweight + prebmi + prehba1c + pretriglyceride + prehdl + preldl +
                               pretotalcholesterol + prealt + preegfr + uacr + presbp + predbp + 
                               ckdpc_40egfr_score +
                               
                               #medications:
                               ncurrtx + MFN + statin + INS +
                               ACEi_or_ARB + BB + CCB + ThZD + loopD + MRA + 
                               steroids + immunosuppr +
                               
                               #Comorbidities
                               qrisk2_smoking_cat + dstartdate_dm_dur_all + 
                               predrug_hypertension + osteoporosis + predrug_dka +
                               predrug_falls + predrug_urinary_frequency + predrug_volume_depletion +
                               predrug_acutepancreatitis + predrug_micturition_control + 
                               predrug_dementia + hosp_admission_prev_year"))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.SGLT2i.unadj <- SE.SGLT2i.unadj <-
    COEFS.DPP4i.unadj <- SE.DPP4i.unadj <-
    COEFS.SGLT2ivsDPP4i.unadj <- SE.SGLT2ivsDPP4i.unadj <-
  # for the adjusted survival models
    COEFS.SGLT2i.adj <- SE.SGLT2i.adj <-
    COEFS.DPP4i.adj <- SE.DPP4i.adj <-
    COEFS.SGLT2ivsDPP4i.adj <- SE.SGLT2ivsDPP4i.adj <-
  # for the overlap weighted survival models
    COEFS.SGLT2i.ow <- SE.SGLT2i.ow <-
    COEFS.DPP4i.ow <- SE.DPP4i.ow <-
    COEFS.SGLT2ivsDPP4i.ow <- SE.SGLT2ivsDPP4i.ow <-
  # for the inverse probability of treatment weighted survival models
    COEFS.SGLT2i.iptw <- SE.SGLT2i.iptw <-
    COEFS.DPP4i.iptw <- SE.DPP4i.iptw <-
    COEFS.SGLT2ivsDPP4i.iptw <- SE.SGLT2ivsDPP4i.iptw <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.DPP4i.unadj[i] <- fit.unadj$coefficients[1]
    COEFS.SGLT2i.unadj[i] <- fit.unadj$coefficients[2]
    SE.DPP4i.unadj[i] <- sqrt(fit.unadj$var[1,1])
    SE.SGLT2i.unadj[i] <- sqrt(fit.unadj$var[2,2])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted, cohort[cohort$.imp == i,])
    
    COEFS.DPP4i.adj[i] <- fit.adj$coefficients[1]
    COEFS.SGLT2i.adj[i] <- fit.adj$coefficients[2]
    SE.DPP4i.adj[i] <- sqrt(fit.adj$var[1,1])
    SE.SGLT2i.adj[i] <- sqrt(fit.adj$var[2,2])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = overlap)
    
    COEFS.DPP4i.ow[i] <- fit.ow$coefficients[1]
    COEFS.SGLT2i.ow[i] <- fit.ow$coefficients[2]
    SE.DPP4i.ow[i] <- sqrt(fit.ow$var[1,1])
    SE.SGLT2i.ow[i] <- sqrt(fit.ow$var[2,2])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = IPTW)
    
    COEFS.DPP4i.iptw[i] <- fit.iptw$coefficients[1]
    COEFS.SGLT2i.iptw[i] <- fit.iptw$coefficients[2]
    SE.DPP4i.iptw[i] <- sqrt(fit.iptw$var[1,1])
    SE.SGLT2i.iptw[i] <- sqrt(fit.iptw$var[2,2])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
    
    #calculate HRs for SGLT2i vs DPP4
    cohort$studydrug <- relevel(factor(cohort$studydrug), ref = "DPP4i")
    
    #unadjusted analyses first
    fit.unadj <- coxph(f, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.SGLT2ivsDPP4i.unadj[i] <- fit.unadj$coefficients[2]
    SE.SGLT2ivsDPP4i.unadj[i] <- sqrt(fit.unadj$var[2,2])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted, cohort[cohort$.imp == i,])
    
    COEFS.SGLT2ivsDPP4i.adj[i] <- fit.adj$coefficients[2]
    SE.SGLT2ivsDPP4i.adj[i] <- sqrt(fit.adj$var[2,2])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = overlap)
    
    COEFS.SGLT2ivsDPP4i.ow[i] <- fit.ow$coefficients[2]
    SE.SGLT2ivsDPP4i.ow[i] <- sqrt(fit.ow$var[2,2])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = IPTW)
    
    COEFS.SGLT2ivsDPP4i.iptw[i] <- fit.iptw$coefficients[2]
    SE.SGLT2ivsDPP4i.iptw[i] <- sqrt(fit.iptw$var[2,2])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
    
    cohort$studydrug <- relevel(factor(cohort$studydrug), ref = "SU")
    
  }
  
  # pool hazard ratios
  unadjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.unadj, SE.SGLT2i.unadj, n.imp)
  unadjusted_DPP4i <- pool.rubin.HR(COEFS.DPP4i.unadj, SE.DPP4i.unadj, n.imp)
  unadjusted_SGLT2ivsDPP4i <- pool.rubin.HR(COEFS.SGLT2ivsDPP4i.unadj, SE.SGLT2ivsDPP4i.unadj, n.imp)
  
  adjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.adj, SE.SGLT2i.adj, n.imp)
  adjusted_DPP4i <- pool.rubin.HR(COEFS.DPP4i.adj, SE.DPP4i.adj, n.imp)  
  adjusted_SGLT2ivsDPP4i <- pool.rubin.HR(COEFS.SGLT2ivsDPP4i.adj, SE.SGLT2ivsDPP4i.adj, n.imp)
  
  overlapweighted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.ow, SE.SGLT2i.ow, n.imp)
  overlapweighted_DPP4i <- pool.rubin.HR(COEFS.DPP4i.ow, SE.DPP4i.ow, n.imp)
  overlapweighted_SGLT2ivsDPP4i <- pool.rubin.HR(COEFS.SGLT2ivsDPP4i.ow, SE.SGLT2ivsDPP4i.ow, n.imp)
  
  iptw_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.iptw, SE.SGLT2i.iptw, n.imp)
  iptw_DPP4i <- pool.rubin.HR(COEFS.DPP4i.iptw, SE.DPP4i.iptw, n.imp)
  iptw_SGLT2ivsDPP4i <- pool.rubin.HR(COEFS.SGLT2ivsDPP4i.iptw, SE.SGLT2ivsDPP4i.iptw, n.imp)
  
  
  # save pooled HR and 95% confidence interval
  unadjusted_SGLT2i_string <- paste0(round(unadjusted_SGLT2i[1], 2), " (", round(unadjusted_SGLT2i[2], 2), ", ", round(unadjusted_SGLT2i[3], 2), ")")
  unadjusted_DPP4i_string <- paste0(round(unadjusted_DPP4i[1], 2), " (", round(unadjusted_DPP4i[2], 2), ", ", round(unadjusted_DPP4i[3], 2), ")")
  unadjusted_SGLT2ivsDPP4i_string <- paste0(round(unadjusted_SGLT2ivsDPP4i[1], 2), " (", round(unadjusted_SGLT2ivsDPP4i[2], 2), ", ", round(unadjusted_SGLT2ivsDPP4i[3], 2), ")")

  adjusted_SGLT2i_string <- paste0(round(adjusted_SGLT2i[1], 2), " (", round(adjusted_SGLT2i[2], 2), ", ", round(adjusted_SGLT2i[3], 2), ")")
  adjusted_DPP4i_string <- paste0(round(adjusted_DPP4i[1], 2), " (", round(adjusted_DPP4i[2], 2), ", ", round(adjusted_DPP4i[3], 2), ")")
  adjusted_SGLT2ivsDPP4i_string <- paste0(round(adjusted_SGLT2ivsDPP4i[1], 2), " (", round(adjusted_SGLT2ivsDPP4i[2], 2), ", ", round(adjusted_SGLT2ivsDPP4i[3], 2), ")")
  
  overlapweighted_SGLT2i_string <- paste0(round(overlapweighted_SGLT2i[1], 2), " (", round(overlapweighted_SGLT2i[2], 2), ", ", round(overlapweighted_SGLT2i[3], 2), ")")
  overlapweighted_DPP4i_string <- paste0(round(overlapweighted_DPP4i[1], 2), " (", round(overlapweighted_DPP4i[2], 2), ", ", round(overlapweighted_DPP4i[3], 2), ")")
  overlapweighted_SGLT2ivsDPP4i_string <- paste0(round(overlapweighted_SGLT2ivsDPP4i[1], 2), " (", round(overlapweighted_SGLT2ivsDPP4i[2], 2), ", ", round(overlapweighted_SGLT2ivsDPP4i[3], 2), ")")

  iptw_SGLT2i_string <- paste0(round(iptw_SGLT2i[1], 2), " (", round(iptw_SGLT2i[2], 2), ", ", round(iptw_SGLT2i[3], 2), ")")
  iptw_DPP4i_string <- paste0(round(iptw_DPP4i[1], 2), " (", round(iptw_DPP4i[2], 2), ", ", round(iptw_DPP4i[3], 2), ")")
  iptw_SGLT2ivsDPP4i_string <- paste0(round(iptw_SGLT2ivsDPP4i[1], 2), " (", round(iptw_SGLT2ivsDPP4i[2], 2), ", ", round(iptw_SGLT2ivsDPP4i[3], 2), ")")
  
  
  # combine in dataframe that we can tabulate
  outcome_sglt2i_hr <- cbind(outcome=k, 
                             count[c("SU_count", "SGLT2i_count")], 
                             followup[c("SU_followup", "SGLT2i_followup")], 
                             events[c("SU_events", "SGLT2i_events")], 
                            unadjusted_SGLT2i_string, 
                            adjusted_SGLT2i_string, 
                            overlapweighted_SGLT2i_string, 
                            iptw_SGLT2i_string)
  
  outcome_dpp4i_hr <- cbind(outcome=k, 
                            count[c("SU_count", "DPP4i_count")], 
                            followup[c("SU_followup", "DPP4i_followup")], 
                            events[c("SU_events", "DPP4i_events")], 
                            unadjusted_DPP4i_string, 
                            adjusted_DPP4i_string, 
                            overlapweighted_DPP4i_string, 
                            iptw_DPP4i_string)
  
  outcome_SGLT2ivsDPP4i_hr <- cbind(outcome=k, 
                                    count[c("DPP4i_count", "SGLT2i_count")], 
                                    followup[c("DPP4i_followup", "SGLT2i_followup")], 
                                    events[c("DPP4i_events", "SGLT2i_events")],  
                                    unadjusted_SGLT2ivsDPP4i_string, 
                                    adjusted_SGLT2ivsDPP4i_string, 
                                    overlapweighted_SGLT2ivsDPP4i_string, 
                                    iptw_SGLT2ivsDPP4i_string)
  
  all_sglt2i_hrs <- rbind(all_sglt2i_hrs, outcome_sglt2i_hr)
  all_dpp4i_hrs <- rbind(all_dpp4i_hrs, outcome_dpp4i_hr)
  all_SGLT2ivsDPP4i_hrs <- rbind(all_SGLT2ivsDPP4i_hrs, outcome_SGLT2ivsDPP4i_hr)
  
  outcome_hr <- rbind(
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "unadjusted", 
          HR = unadjusted_SGLT2i[1], LB = unadjusted_SGLT2i[2], UB = unadjusted_SGLT2i[3], 
          string = unadjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "adjusted", 
          HR = adjusted_SGLT2i[1], LB = adjusted_SGLT2i[2], UB = adjusted_SGLT2i[3], 
          string = adjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "overlap-weighted", 
          HR = overlapweighted_SGLT2i[1], LB = overlapweighted_SGLT2i[2], UB = overlapweighted_SGLT2i[3], 
          string = overlapweighted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "IPTW", 
          HR = iptw_SGLT2i[1], LB = iptw_SGLT2i[2], UB = iptw_SGLT2i[3], 
          string = iptw_SGLT2i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "unadjusted", 
          HR = unadjusted_DPP4i[1], LB = unadjusted_DPP4i[2], UB = unadjusted_DPP4i[3], 
          string = unadjusted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "adjusted", 
          HR = adjusted_DPP4i[1], LB = adjusted_DPP4i[2], UB = adjusted_DPP4i[3], 
          string = adjusted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "overlap-weighted", 
          HR = overlapweighted_DPP4i[1], LB = overlapweighted_DPP4i[2], UB = overlapweighted_DPP4i[3], 
          string = overlapweighted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "IPTW", 
          HR = iptw_DPP4i[1], LB = iptw_DPP4i[2], UB = iptw_DPP4i[3], 
          string = iptw_DPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "unadjusted", 
          HR = unadjusted_SGLT2ivsDPP4i[1], LB = unadjusted_SGLT2ivsDPP4i[2], UB = unadjusted_SGLT2ivsDPP4i[3], 
          string = unadjusted_SGLT2ivsDPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "adjusted", 
          HR = adjusted_SGLT2ivsDPP4i[1], LB = adjusted_SGLT2ivsDPP4i[2], UB = adjusted_SGLT2ivsDPP4i[3], 
          string = adjusted_SGLT2ivsDPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "overlap-weighted", 
          HR = overlapweighted_SGLT2ivsDPP4i[1], LB = overlapweighted_SGLT2ivsDPP4i[2], UB = overlapweighted_SGLT2ivsDPP4i[3], 
          string = overlapweighted_SGLT2ivsDPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "IPTW", 
          HR = iptw_SGLT2ivsDPP4i[1], LB = iptw_SGLT2ivsDPP4i[2], UB = iptw_SGLT2ivsDPP4i[3], 
          string = iptw_SGLT2ivsDPP4i_string)
  )
  
  all_hrs <- rbind(all_hrs, outcome_hr)
  
}

# review HRs for different treatment contrasts
flextable(all_sglt2i_hrs)
flextable(all_dpp4i_hrs)
flextable(all_SGLT2ivsDPP4i_hrs) # need to compare pp outcomes as these censor DPP4i when starting SU and vice versa
# we can see that the DPP4i vs SU HRs are non-significant, and we can therefore pool these 2 groups together.

## compare SGLT2i vs combined group of DPP4i/SU

#create empty data frame to which we can append the hazard ratios once calculated
all_SGLT2ivsDPP4iSU_hrs <- data.frame()

## remove double overlapping entries for DPP4i and SU that overlap (take one only)
cohort <- cohort %>% group_by(.imp, patid) %>% filter(
  !duplicated(studydrug2)
) %>% ungroup()


for (k in kf_key_outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  # calculate number of subjects in each group
  count <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug2, events) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2"))
  
  f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + 
                               
                               #sociodemographic characteristics:
                               dstartdate_age + malesex + imd2015_10 + 
                               ethnicity_5cat + initiation_year +
                               
                               #laboratory and vital signs measurements:
                               preweight + prebmi + prehba1c + pretriglyceride + prehdl + preldl +
                               pretotalcholesterol + prealt + preegfr + uacr + presbp + predbp + 
                               ckdpc_40egfr_score +
                               
                               #medications:
                               ncurrtx + MFN + statin + INS +
                               ACEi_or_ARB + BB + CCB + ThZD + loopD + MRA + 
                               steroids + immunosuppr +
                               
                               #Comorbidities
                               qrisk2_smoking_cat + dstartdate_dm_dur_all + 
                               predrug_hypertension + osteoporosis + predrug_dka +
                               predrug_falls + predrug_urinary_frequency + predrug_volume_depletion +
                               predrug_acutepancreatitis + predrug_micturition_control + 
                               predrug_dementia + hosp_admission_prev_year"))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.SGLT2i.unadj <- SE.SGLT2i.unadj <-
    # for the adjusted survival models
    COEFS.SGLT2i.adj <- SE.SGLT2i.adj <-
    # for the overlap weighted survival models
    COEFS.SGLT2i.ow <- SE.SGLT2i.ow <-
    # for the inverse probability of treatment weighted survival models
    COEFS.SGLT2i.iptw <- SE.SGLT2i.iptw <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f2, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.SGLT2i.unadj[i] <- fit.unadj$coefficients[1]
    SE.SGLT2i.unadj[i] <- sqrt(fit.unadj$var[1,1])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted2, cohort[cohort$.imp == i,])
    
    COEFS.SGLT2i.adj[i] <- fit.adj$coefficients[1]
    SE.SGLT2i.adj[i] <- sqrt(fit.adj$var[1,1])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted2, cohort[cohort$.imp ==i,], weights = overlap2)
    
    COEFS.SGLT2i.ow[i] <- fit.ow$coefficients[1]
    SE.SGLT2i.ow[i] <- sqrt(fit.ow$var[1,1])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted2, cohort[cohort$.imp ==i,], weights = IPTW2)
    
    COEFS.SGLT2i.iptw[i] <- fit.iptw$coefficients[1]
    SE.SGLT2i.iptw[i] <- sqrt(fit.iptw$var[1,1])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
  }
  
  # pool hazard ratios
  unadjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.unadj, SE.SGLT2i.unadj, n.imp)
  adjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.adj, SE.SGLT2i.adj, n.imp)
  overlapweighted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.ow, SE.SGLT2i.ow, n.imp)
  iptw_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.iptw, SE.SGLT2i.iptw, n.imp)
  
  # save pooled HR and 95% confidence interval
  unadjusted_SGLT2i_string <- paste0(round(unadjusted_SGLT2i[1], 2), " (", round(unadjusted_SGLT2i[2], 2), ", ", round(unadjusted_SGLT2i[3], 2), ")")
  adjusted_SGLT2i_string <- paste0(round(adjusted_SGLT2i[1], 2), " (", round(adjusted_SGLT2i[2], 2), ", ", round(adjusted_SGLT2i[3], 2), ")")
  overlapweighted_SGLT2i_string <- paste0(round(overlapweighted_SGLT2i[1], 2), " (", round(overlapweighted_SGLT2i[2], 2), ", ", round(overlapweighted_SGLT2i[3], 2), ")")
  iptw_SGLT2i_string <- paste0(round(iptw_SGLT2i[1], 2), " (", round(iptw_SGLT2i[2], 2), ", ", round(iptw_SGLT2i[3], 2), ")")
 
  # combine in dataframe that we can tabulate
#  outcome_SGLT2ivsDPP4iSU_hr <- cbind(outcome=k, count[c(1:2)], followup[c(1:2)], events[c(1:2)], 
#                            unadjusted_SGLT2i_string, adjusted_SGLT2i_string, overlapweighted_SGLT2i_string, iptw_SGLT2i_string)
  
#  all_SGLT2ivsDPP4iSU_hrs <- rbind(all_SGLT2ivsDPP4iSU_hrs, outcome_SGLT2ivsDPP4iSU_hr)
  
  outcome_SGLT2ivsDPP4iSU_hr <- rbind(
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "unadjusted", 
          HR = unadjusted_SGLT2i[1], LB = unadjusted_SGLT2i[2], UB = unadjusted_SGLT2i[3], string = unadjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "adjusted", 
          HR = adjusted_SGLT2i[1], LB = adjusted_SGLT2i[2], UB = adjusted_SGLT2i[3], string = adjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "overlap-weighted", 
          HR = overlapweighted_SGLT2i[1], LB = overlapweighted_SGLT2i[2], UB = overlapweighted_SGLT2i[3], string = overlapweighted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "IPTW", 
          HR = iptw_SGLT2i[1], LB = iptw_SGLT2i[2], UB = iptw_SGLT2i[3], string = iptw_SGLT2i_string))
    
    all_SGLT2ivsDPP4iSU_hrs <- rbind(all_SGLT2ivsDPP4iSU_hrs, outcome_SGLT2ivsDPP4iSU_hr)
  
}
############################3 STORE AND DISPLAY HAZARD RATIOS################################################################

# save all_hrs table and SGLT2i vs DPP4i/su table
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")

save(all_hrs, file=paste0(today, "_all_hrs_incl_egfr_below_60.Rda"))
save(all_SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_all_SGLT2ivsDPP4iSU_hrs_incl_egfr_below_60.Rda"))
#load("2024-02-26_all_hrs_incl_egfr_below_60.Rda")
#load("2024-02-26_all_SGLT2ivsDPP4iSU_hrs_incl_egfr_below_60.Rda")

# show table with events, follow up time, and hazard ratios
flextable(all_sglt2i_hrs)
flextable(all_dpp4i_hrs)
flextable(all_SGLT2ivsDPP4i_hrs)
flextable(all_SGLT2ivsDPP4iSU_hrs)

# prep all_hrs dataframe for forest plot to show hazard ratios and add literature-reported HR
trial_hr <- cbind(outcome = "ckd_egfr40", contrast = "SGLT2i vs SU", analysis = "meta-analysis of RCTs", 
                  HR = 0.62, LB = 0.56, UB = 0.68, string = "0.62 (0.56, 0.68)")
# all_hrs <- rbind(trial_hr, all_hrs) # we do not want to have the individual 
all_SGLT2ivsDPP4iSU_hrs <- rbind(trial_hr, all_SGLT2ivsDPP4iSU_hrs)

all_hrs$model <- paste0(all_hrs$string, " [", all_hrs$analysis, "]")
all_hrs$model <- factor(all_hrs$model, levels = unique(all_hrs$model))

all_SGLT2ivsDPP4iSU_hrs$model <- paste0(all_SGLT2ivsDPP4iSU_hrs$string, " [", all_SGLT2ivsDPP4iSU_hrs$analysis, "]")
all_SGLT2ivsDPP4iSU_hrs$model <- factor(all_SGLT2ivsDPP4iSU_hrs$model, levels = unique(all_SGLT2ivsDPP4iSU_hrs$model))

# have to coerce HR and CI to class numeric as they sometimes default to character
class(all_hrs$HR) <- class(all_hrs$LB) <- class(all_hrs$UB) <-
  class(all_SGLT2ivsDPP4iSU_hrs$HR) <- class(all_SGLT2ivsDPP4iSU_hrs$LB) <- class(all_SGLT2ivsDPP4iSU_hrs$UB) <- "numeric"

# remove unadjusted HR as we do not want to plot these
all_hrs <- all_hrs[!all_hrs$analysis == "unadjusted",]
all_SGLT2ivsDPP4iSU_hrs <- all_SGLT2ivsDPP4iSU_hrs[!all_SGLT2ivsDPP4iSU_hrs$analysis == "unadjusted",]

# create labels
labels <- cbind(outcome = "", contrast = "", analysis = "", HR = "", LB = "", UB = "", string = "", model = "Hazard Ratio (95% CI)")
labels_plot <- all_hrs

for (i in unique(all_hrs$contrast)) {
  for (m in unique(all_hrs$outcome)) {
    labels_temp <- labels %>% data.frame()
    labels_temp$outcome <- m
    labels_temp$contrast <- i
    labels_plot <- rbind(labels_temp, labels_plot)
  }
}

labels_plot$model <- factor(labels_plot$model, levels = rev(unique(labels_plot$model)))

# create labels for forest plot with studydrug2 (SGLT2i vs DPP4i/SU combined)
labels_plot2 <- all_SGLT2ivsDPP4iSU_hrs

for (i in unique(all_SGLT2ivsDPP4iSU_hrs$contrast)) {
  for (m in unique(all_SGLT2ivsDPP4iSU_hrs$outcome)) {
    labels_temp <- labels %>% data.frame()
    labels_temp$outcome <- m
    labels_temp$contrast <- i
    labels_plot2 <- rbind(labels_temp, labels_plot2)
  }
}

labels_plot2$model <- factor(labels_plot2$model, levels = rev(unique(labels_plot2$model)))

# layout for plots below
layout <- c(
  area(t = 0, l = 0, b = 30, r = 4), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 1, l = 5, b = 30, r = 9) # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
)
# #############3A FOREST PLOT FOR HRs FOR NEW CKD#############
# p_ckd345_1 <- 
#   all_hrs %>%
#   filter(outcome == "ckd_345") %>%
#   filter(contrast == "SGLT2i vs SU") %>%
#   ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
#   scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
#   coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_345",]$analysis)) + 1), 
#                   xlim=c(0.5, 2)) +
#   theme_classic() +
#   geom_point(aes(x=HR), shape=15, size=3) +
#   geom_linerange(aes(xmin=LB, xmax=UB)) +
#   geom_vline(xintercept = 1, linetype="dashed") +
#   annotate("text", x = .65, 
#            y = length(unique(all_hrs[all_hrs$outcome == "ckd_345",]$analysis)) + 1, 
#            label = "SGLT2i protective") +
#   annotate("text", x = 1.5, 
#            y = length(unique(all_hrs[all_hrs$outcome == "ckd_345",]$analysis)) + 1, 
#            label = "SGLT2i harmful") +
#   labs(x="Hazard Ratio", y="") +
#   ggtitle("New CKD (SGLT2i vs SU)") +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y= element_blank(),
#         axis.text.y= element_blank(),
#         axis.title.y= element_blank()) 
# 
# 
# p_left_1 <-
#   labels_plot %>%
#   filter(outcome == "ckd_345") %>%
#   filter(contrast == "SGLT2i vs SU") %>%
#   ggplot(aes(y = (model))) + 
#   # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
#   geom_text(
#     aes(x = 1, label = model),
#     hjust = 0,
#     fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_345" & labels_plot$contrast == "SGLT2i vs SU",]$
#                         model == "Hazard Ratio (95% CI)", "bold", "plain")
#   ) +
#   theme_void() +
#   coord_cartesian(xlim = c(0, 4))
# 
# # final plot arrangement
# p_left_1 + p_ckd345_1 + plot_layout(design = layout)
# 
# 
# ## same for SGLT2i vs DPP4
# 
# p_ckd345_2 <- 
#   all_hrs %>%
#   filter(outcome == "ckd_345") %>%
#   filter(contrast == "SGLT2i vs DPP4i") %>%
#   ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
#   scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
#   coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_345",]$analysis)) + 1), 
#                   xlim=c(0.5, 2)) +
#   theme_classic() +
#   geom_point(aes(x=HR), shape=15, size=3) +
#   geom_linerange(aes(xmin=LB, xmax=UB)) +
#   geom_vline(xintercept = 1, linetype="dashed") +
#   annotate("text", x = .65, 
#            y = length(unique(all_hrs[all_hrs$outcome == "ckd_345",]$analysis)) + 1, 
#            label = "SGLT2i protective") +
#   annotate("text", x = 1.5, 
#            y = length(unique(all_hrs[all_hrs$outcome == "ckd_345",]$analysis)) + 1, 
#            label = "SGLT2i harmful") +
#   labs(x="Hazard Ratio", y="") +
#   ggtitle("New CKD (SGLT2i vs DPP4i)") +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y= element_blank(),
#         axis.text.y= element_blank(),
#         axis.title.y= element_blank()) 
# 
# p_left_2 <-
#   labels_plot %>%
#   filter(outcome == "ckd_345") %>%
#   filter(contrast == "SGLT2i vs DPP4i") %>%
#   ggplot(aes(y = (model))) + 
#   # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
#   geom_text(
#     aes(x = 1, label = model),
#     hjust = 0,
#     fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_345" & labels_plot$contrast == "SGLT2i vs DPP4i",]$
#                         model == "Hazard Ratio (95% CI)", "bold", "plain")
#   ) +
#   theme_void() +
#   coord_cartesian(xlim = c(0, 4))
# 
# # final plot arrangement
# p_left_2 + p_ckd345_2 + plot_layout(design = layout)
# 
# 
# ## same for DPP4i vs SU
# # for this outcome we will use the _pp outcome as patients were then censored if they switched from DPP4i to SU
# # (no overlapping intervals)
# p_ckd345_3 <- 
#   all_hrs %>%
#   filter(outcome == "ckd_345_pp") %>%
#   filter(contrast == "DPP4i vs SU") %>%
#   ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
#   scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
#   coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_345_pp",]$analysis)) + 1), 
#                   xlim=c(0.5, 2)) +
#   theme_classic() +
#   geom_point(aes(x=HR), shape=15, size=3) +
#   geom_linerange(aes(xmin=LB, xmax=UB)) +
#   geom_vline(xintercept = 1, linetype="dashed") +
#   annotate("text", x = .65, 
#            y = length(unique(all_hrs[all_hrs$outcome == "ckd_345_pp",]$analysis)) + 1, 
#            label = "DPP4i protective") +
#   annotate("text", x = 1.5, 
#            y = length(unique(all_hrs[all_hrs$outcome == "ckd_345_pp",]$analysis)) + 1, 
#            label = "DPP4i harmful") +
#   labs(x="Hazard Ratio", y="") +
#   ggtitle("New CKD (DPP4i vs SU)") +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y= element_blank(),
#         axis.text.y= element_blank(),
#         axis.title.y= element_blank()) 
# 
# p_left_3 <-
#   labels_plot %>%
#   filter(outcome == "ckd_345_pp") %>%
#   filter(contrast == "DPP4i vs SU") %>%
#   ggplot(aes(y = (model))) + 
#   # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
#   geom_text(
#     aes(x = 1, label = model),
#     hjust = 0,
#     fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_345_pp" & labels_plot$contrast == "DPP4i vs SU",]$
#                         model == "Hazard Ratio (95% CI)", "bold", "plain")
#   ) +
#   theme_void() +
#   coord_cartesian(xlim = c(0, 4))
# 
# # final plot arrangement
# p_left_3 + p_ckd345_3 + plot_layout(design = layout)
# 
# 
# ## and for SGLT2i vs DPP4i/SU combined
# p_ckd345_4 <- 
#   all_SGLT2ivsDPP4iSU_hrs %>%
#   filter(outcome == "ckd_345") %>%
#   ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
#   scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
#   coord_cartesian(ylim=c(1,length(unique(all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$outcome == "ckd_345",]$analysis)) + 1), 
#                   xlim=c(0.5, 2)) +
#   theme_classic() +
#   geom_point(aes(x=HR), shape=15, size=3) +
#   geom_linerange(aes(xmin=LB, xmax=UB)) +
#   geom_vline(xintercept = 1, linetype="dashed") +
#   annotate("text", x = .65, 
#            y = length(unique(all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$outcome == "ckd_345",]$analysis)) + 1, 
#            label = "SGLT2i protective") +
#   annotate("text", x = 1.5, 
#            y = length(unique(all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$outcome == "ckd_345",]$analysis)) + 1, 
#            label = "SGLT2i harmful") +
#   labs(x="Hazard Ratio", y="") +
#   ggtitle("New CKD (SGLT2i vs DPP4i/SU)") +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y= element_blank(),
#         axis.text.y= element_blank(),
#         axis.title.y= element_blank()) 
# 
# 
# 
# p_left_4 <-
#   labels_plot2 %>%
#   filter(outcome == "ckd_345") %>%
#   ggplot(aes(y = (model))) + 
#   # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
#   geom_text(
#     aes(x = 1, label = model),
#     hjust = 0,
#     fontface = ifelse(labels_plot2[labels_plot2$outcome == "ckd_345",]$
#                         model == "Hazard Ratio (95% CI)", "bold", "plain")
#   ) +
#   theme_void() +
#   coord_cartesian(xlim = c(0, 4))
# 
# # final plot arrangement
# p_left_4 + p_ckd345_4 + plot_layout(design = layout)

#############3B FOREST PLOT FOR HRs FOR 40% DECLINE IN EGFR#############

# plot
p_egfr40_1 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr40") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40",]$analysis)) + 1), 
                  xlim=c(0.5, 2)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.5, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40",]$analysis)) + 1, 
           label = "Favours SU") +
  labs(x="", y="") +
  ggtitle("40% decline in eGFR or ESRD", subtitle = "SGLT2i vs SU") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_1 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr40") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = (model))) + 
  geom_text(
    aes(x = 1, label = model),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr40" & labels_plot$contrast == "SGLT2i vs SU",]$
                        string == "", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))


## same for SGLT2i vs DPP4

p_egfr40_2 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr40") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40" & !all_hrs$analysis == "meta-analysis of RCTs",]$analysis)) + 1), 
                  xlim=c(0.5, 2)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40" & !all_hrs$analysis == "meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.5, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40" & !all_hrs$analysis == "meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours DPP4i") +
  labs(x="", y="") +
  ggtitle(label = "", subtitle="SGLT2i vs DPP4i") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 

p_left_2 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr40") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = (model))) + 
  geom_text(
    aes(x = 1, label = model),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr40" & labels_plot$contrast == "SGLT2i vs DPP4i",]$
                        string == "", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))




## same for DPP4i vs SU

p_egfr40_3 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr40_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40_pp" & !all_hrs$analysis == "meta-analysis of RCTs",]$analysis)) + 1), 
                  xlim=c(0.5, 2)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40_pp" & !all_hrs$analysis == "meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours DPP4i") +
  annotate("text", x = 1.5, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr40_pp" & !all_hrs$analysis == "meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours SU") +
  labs(x="Hazard Ratio", y="") +
  ggtitle(label = "", subtitle="DPP4i vs SU") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 

p_left_3 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr40_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  ggplot(aes(y = (model))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 1, label = model),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr40_pp" & labels_plot$contrast == "DPP4i vs SU",]$
                        string == "", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))


layout2 <- c(
  area(t = 2, l = 0, b = 7, r = 4), 
  area(t = 2, l = 5, b = 7, r = 9),
  area(t = 8, l = 0, b = 13, r = 4), 
  area(t = 8, l = 5, b = 13, r = 9),
  area(t = 14, l = 0, b = 19, r = 4), 
  area(t = 14, l = 5, b = 19, r = 9)
)

p_left_1 + p_egfr40_1 + p_left_2 + p_egfr40_2 + p_left_3 + p_egfr40_3 + plot_layout(design = layout2)

## now for SGLT2i vs DPP4i/SU combined

# plot
p_egfr40_4 <- 
  all_SGLT2ivsDPP4iSU_hrs %>%
  filter(outcome == "ckd_egfr40") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
  coord_cartesian(ylim=c(1,length(unique(all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$outcome == "ckd_egfr40",]$analysis)) + 1), 
                  xlim=c(0.5, 2)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$outcome == "ckd_egfr40",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.5, 
           y = length(unique(all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$outcome == "ckd_egfr40",]$analysis)) + 1, 
           label = "Favours DPP4i/SU") +
  labs(x="Hazard Ratio", y="") +
  ggtitle("40% decline in eGFR or ESRD") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5)) 


p_left_4 <-
  labels_plot2 %>%
  filter(outcome == "ckd_egfr40") %>%
  ggplot(aes(y = (model))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 1, label = model),
    hjust = 0,
    fontface = ifelse(labels_plot2[labels_plot2$outcome == "ckd_egfr40",]$
                        model == "Hazard Ratio (95% CI)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

# final plot arrangement
p_left_4 + p_egfr40_4 + plot_layout(design = layout)

############################4 SENSITIVITY ANALYSES################################################################

cohort <- cohort %>% mutate(
  egfr_below_60 = ifelse(preegfr < 60, T, F),
  macroalbuminuria = ifelse(uacr < 30, F, T),
  microalbuminuria = ifelse(uacr <3, F, ifelse(macroalbuminuria == T, F, T))
)

kf_key_outcomes <- c(#"ckd_345", 
  "ckd_egfr40")

# compare HRs between different groups (albuminuria, egfr below 60)

#create empty data frame to which we can append the hazard ratios once calculated
alb_SGLT2ivsDPP4iSU_hrs <- alb_hrs <-
  egfr_SGLT2ivsDPP4iSU_hrs <- egfr_hrs <- data.frame()

cohort <- cohort %>% mutate(albuminuria_status = ifelse(macroalbuminuria == T, "Macroalbuminuria", ifelse(
  microalbuminuria == T, "Microalbuminuria", "No albuminuria"
)),
albuminuria_status = relevel(as.factor(albuminuria_status), ref = "No albuminuria")
)

## analyses stratified by albuminuria status
for (k in kf_key_outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  # calculate number of subjects in each group
  count <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2,albuminuria_status) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2,albuminuria_status) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2,albuminuria_status) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug2, albuminuria_status, events) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*albuminuria_status"))
  
  f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*albuminuria_status +

                               #sociodemographic characteristics:
                               dstartdate_age + malesex + imd2015_10 +
                               ethnicity_5cat + initiation_year +

                               #laboratory and vital signs measurements:
                               preweight + prebmi + prehba1c + pretriglyceride + prehdl + preldl +
                               pretotalcholesterol + prealt + preegfr + uacr + presbp + predbp +
                               ckdpc_40egfr_score +

                               #medications:
                               ncurrtx + MFN + statin + INS +
                               ACEi_or_ARB + BB + CCB + ThZD + loopD + MRA +
                               steroids + immunosuppr +

                               #Comorbidities
                               qrisk2_smoking_cat + dstartdate_dm_dur_all +
                               predrug_hypertension + osteoporosis + predrug_dka +
                               predrug_falls + predrug_urinary_frequency + predrug_volume_depletion +
                               predrug_acutepancreatitis + predrug_micturition_control +
                               predrug_dementia + hosp_admission_prev_year"))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.noalb.unadj <- SE.noalb.unadj <-
    COEFS.microalb.unadj <- SE.microalb.unadj <-
    COEFS.macroalb.unadj <- SE.macroalb.unadj <-
    # for the adjusted survival models
    COEFS.noalb.adj <- SE.noalb.adj <-
    COEFS.microalb.adj <- SE.microalb.adj <-
    COEFS.macroalb.adj <- SE.macroalb.adj <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f2, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.noalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"]
    SE.noalb.unadj[i] <- sqrt(fit.unadj$var[1,1])
    
    COEFS.microalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:albuminuria_statusMicroalbuminuria"]
    SE.microalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[5,5]) + 2 * vcov(fit.unadj)[1,5])
    
    COEFS.macroalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:albuminuria_statusMacroalbuminuria"]
    SE.macroalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[4,4]) + 2 * vcov(fit.unadj)[1,4])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted2, cohort[cohort$.imp == i,])
    
    COEFS.noalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"]
    SE.noalb.adj[i] <- sqrt(fit.adj$var[1,1])
    
    COEFS.microalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:albuminuria_statusMicroalbuminuria"]
    SE.microalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[60,60]) + 2 * vcov(fit.adj)[1,60])
    
    COEFS.macroalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:albuminuria_statusMacroalbuminuria"]
    SE.macroalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[59,59]) + 2 * vcov(fit.adj)[1,59])
    
    p_value_interaction_alb <- summary(fit.adj)$coefficients[59,5]
    
  }
  
  # pool hazard ratios
  unadjusted_noalb <- pool.rubin.HR(COEFS.noalb.unadj, SE.noalb.unadj, n.imp)
  adjusted_noalb <- pool.rubin.HR(COEFS.noalb.adj, SE.noalb.adj, n.imp)

  unadjusted_microalb <- pool.rubin.HR(COEFS.microalb.unadj, SE.microalb.unadj, n.imp)
  adjusted_microalb <- pool.rubin.HR(COEFS.microalb.adj, SE.microalb.adj, n.imp)
  
  unadjusted_macroalb <- pool.rubin.HR(COEFS.macroalb.unadj, SE.macroalb.unadj, n.imp)
  adjusted_macroalb <- pool.rubin.HR(COEFS.macroalb.adj, SE.macroalb.adj, n.imp)

  # save pooled HR and 95% confidence interval
  unadjusted_noalb_string <- paste0(round(unadjusted_noalb[1], 2), " (", round(unadjusted_noalb[2], 2), ", ", round(unadjusted_noalb[3], 2), ")")
  adjusted_noalb_string <- paste0(round(adjusted_noalb[1], 2), " (", round(adjusted_noalb[2], 2), ", ", round(adjusted_noalb[3], 2), ")")

  unadjusted_microalb_string <- paste0(round(unadjusted_microalb[1], 2), " (", round(unadjusted_microalb[2], 2), ", ", round(unadjusted_microalb[3], 2), ")")
  adjusted_microalb_string <- paste0(round(adjusted_microalb[1], 2), " (", round(adjusted_microalb[2], 2), ", ", round(adjusted_microalb[3], 2), ")")

  unadjusted_macroalb_string <- paste0(round(unadjusted_macroalb[1], 2), " (", round(unadjusted_macroalb[2], 2), ", ", round(unadjusted_macroalb[3], 2), ")")
  adjusted_macroalb_string <- paste0(round(adjusted_macroalb[1], 2), " (", round(adjusted_macroalb[2], 2), ", ", round(adjusted_macroalb[3], 2), ")")
 
  # combine in dataframe that we can tabulate
  noalb_hr <- cbind(outcome=k, count[1,c(2:3)], followup[1,c(2:3)], events[1,c(2:3)],
                    unadjusted=unadjusted_noalb_string, adjusted=adjusted_noalb_string#, overlapweighted_noalb_string, iptw_noalb_string
  )
  microalb_hr <- cbind(outcome=k, count[2,c(2:3)], followup[2,c(2:3)], events[2,c(2:3)],
                       unadjusted=unadjusted_microalb_string, adjusted=adjusted_microalb_string#, overlapweighted_microalb_string, iptw_microalb_string
  )
  macroalb_hr <- cbind(outcome=k, count[3,c(2:3)], followup[3,c(2:3)], events[3,c(2:3)],
                       unadjusted=unadjusted_macroalb_string, adjusted=adjusted_macroalb_string#, overlapweighted_macroalb_string, iptw_macroalb_string
  )
  
  alb_SGLT2ivsDPP4iSU_hrs <- rbind(noalb_hr, microalb_hr, macroalb_hr)
  
  temp <- rbind(
    cbind(outcome = k, contrast = "No albuminuria", analysis = "adjusted",
          HR = adjusted_noalb[1], LB = adjusted_noalb[2], UB = adjusted_noalb[3], string = adjusted_noalb_string),
    cbind(outcome = k, contrast = "Microalbuminuria", analysis = "adjusted",
          HR = adjusted_microalb[1], LB = adjusted_microalb[2], UB = adjusted_microalb[3], string = adjusted_microalb_string),
    cbind(outcome = k, contrast = "Macroalbuminuria", analysis = "adjusted",
          HR = adjusted_macroalb[1], LB = adjusted_macroalb[2], UB = adjusted_macroalb[3], string = adjusted_macroalb_string)
  )
  
  alb_hrs <- rbind(alb_hrs, temp)
  
}

## analyses stratified by egfr status
for (k in kf_key_outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  # calculate number of subjects in each group
  count <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2,egfr_below_60) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2,egfr_below_60) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- cohort[cohort$.imp > 0,] %>%
    group_by(studydrug2,egfr_below_60) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug2, egfr_below_60, events) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*egfr_below_60"))
  
  f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*egfr_below_60 +

                               #sociodemographic characteristics:
                               dstartdate_age + malesex + imd2015_10 +
                               ethnicity_5cat + initiation_year +

                               #laboratory and vital signs measurements:
                               preweight + prebmi + prehba1c + pretriglyceride + prehdl + preldl +
                               pretotalcholesterol + prealt + preegfr + uacr + presbp + predbp +
                               ckdpc_40egfr_score +

                               #medications:
                               ncurrtx + MFN + statin + INS + 
                               ACEi_or_ARB + BB + CCB + ThZD + loopD + MRA +
                               steroids + immunosuppr +

                               #Comorbidities
                               qrisk2_smoking_cat + dstartdate_dm_dur_all +
                               predrug_hypertension + osteoporosis + predrug_dka +
                               predrug_falls + predrug_urinary_frequency + predrug_volume_depletion +
                               predrug_acutepancreatitis + predrug_micturition_control +
                               predrug_dementia + hosp_admission_prev_year"))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.above60.unadj <- SE.above60.unadj <-
    COEFS.below60.unadj <- SE.below60.unadj <-
    # for the adjusted survival models
    COEFS.above60.adj <- SE.above60.adj <-
    COEFS.below60.adj <- SE.below60.adj <-
    
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f2, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.above60.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"]
    SE.above60.unadj[i] <- sqrt(fit.unadj$var[1,1])
    
    COEFS.below60.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:egfr_below_60TRUE"]
    SE.below60.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[3,3]) + 2 * vcov(fit.unadj)[1,3])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted2, cohort[cohort$.imp == i,])
    
    COEFS.above60.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"]
    SE.above60.adj[i] <- sqrt(fit.adj$var[1,1])
    
    COEFS.below60.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:egfr_below_60TRUE"]
    SE.below60.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[58,58]) + 2 * vcov(fit.adj)[1,58])
    
    p_value_interaction_egfr <- summary(fit.adj)$coefficients[58,5]
    # rm(fit.unadj)
    # rm(fit.adj)
    
  }
  
  # pool hazard ratios
  unadjusted_above60 <- pool.rubin.HR(COEFS.above60.unadj, SE.above60.unadj, n.imp)
  adjusted_above60 <- pool.rubin.HR(COEFS.above60.adj, SE.above60.adj, n.imp)
  
  unadjusted_below60 <- pool.rubin.HR(COEFS.below60.unadj, SE.below60.unadj, n.imp)
  adjusted_below60 <- pool.rubin.HR(COEFS.below60.adj, SE.below60.adj, n.imp)
  
  
  # save pooled HR and 95% confidence interval
  unadjusted_above60_string <- paste0(round(unadjusted_above60[1], 2), " (", round(unadjusted_above60[2], 2), ", ", round(unadjusted_above60[3], 2), ")")
  adjusted_above60_string <- paste0(round(adjusted_above60[1], 2), " (", round(adjusted_above60[2], 2), ", ", round(adjusted_above60[3], 2), ")")
  
  unadjusted_below60_string <- paste0(round(unadjusted_below60[1], 2), " (", round(unadjusted_below60[2], 2), ", ", round(unadjusted_below60[3], 2), ")")
  adjusted_below60_string <- paste0(round(adjusted_below60[1], 2), " (", round(adjusted_below60[2], 2), ", ", round(adjusted_below60[3], 2), ")")
  
  
  # combine in dataframe that we can tabulate
  above60_hr <- cbind(outcome=k, count[1,c(2:3)], followup[1,c(2:3)], events[1,c(2:3)],
                      unadjusted=unadjusted_above60_string, adjusted=adjusted_above60_string
  )
  below60_hr <- cbind(outcome=k, count[2,c(2:3)], followup[2,c(2:3)], events[2,c(2:3)],
                      unadjusted=unadjusted_below60_string, adjusted=adjusted_below60_string
  )
  
  
  egfr_SGLT2ivsDPP4iSU_hrs <- rbind(above60_hr, below60_hr)
  
  temp <- rbind(
    cbind(outcome = k, contrast = "eGFR above 60", analysis = "adjusted",
          HR = adjusted_above60[1], LB = adjusted_above60[2], UB = adjusted_above60[3], 
          string = adjusted_above60_string),
    cbind(outcome = k, contrast = "eGFR below 60", analysis = "adjusted",
          HR = adjusted_below60[1], LB = adjusted_below60[2], UB = adjusted_below60[3], 
          string = adjusted_below60_string)
  )
  
  egfr_hrs <- rbind(egfr_hrs, temp)
  
}

## save HRs and display
flextable(alb_SGLT2ivsDPP4iSU_hrs)
flextable(egfr_SGLT2ivsDPP4iSU_hrs)
save(alb_hrs, file=paste0(today, "_alb_hrs_incl_egfr_below_60.Rda"))
save(egfr_hrs, file=paste0(today, "_egfr_hrs_incl_egfr_below_60.Rda"))
#load("2024-02-26_alb_hrs_incl_egfr_below_60.Rda")
#load("2024-02-26_egfr_hrs_incl_egfr_below_60.Rda")

# prep data frames with row for overall
overall <- all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$analysis == "adjusted",] %>% select(-model)
overall$contrast <- "overall"

alb_hrs <- rbind(overall, alb_hrs)
egfr_hrs <- rbind(overall, egfr_hrs)

alb_hrs$model <- paste0(alb_hrs$string, " [", alb_hrs$contrast, "]")
alb_hrs$model <- factor(alb_hrs$model, levels = unique(alb_hrs$model))

egfr_hrs$model <- paste0(egfr_hrs$string, " [", egfr_hrs$contrast, "]")
egfr_hrs$model <- factor(egfr_hrs$model, levels = unique(egfr_hrs$model))

# have to coerce HR and CI to class numeric as they sometimes default to character
class(alb_hrs$HR) <- class(alb_hrs$LB) <- class(alb_hrs$UB) <-
  class(egfr_hrs$HR) <- class(egfr_hrs$LB) <- class(egfr_hrs$UB) <- "numeric"

# create labels for forest plots
labels_plot3 <- alb_hrs

for (i in unique(alb_hrs$outcome)) {
    labels_temp <- labels %>% data.frame()
    labels_temp$outcome <- i
    labels_plot3 <- rbind(labels_temp, labels_plot3)
}

labels_plot3 <- labels_plot3 %>% mutate(model = ifelse(model == "Hazard Ratio (95% CI)", 
                                                       "Adjusted Hazard Ratio (95% CI)",
                                                       as.character(model)))

labels_plot3$model <- factor(labels_plot3$model, levels = rev(unique(labels_plot3$model)))

labels_plot4 <- egfr_hrs

for (m in unique(egfr_hrs$outcome)) {
    labels_temp <- labels %>% data.frame()
    labels_temp$outcome <- m
    labels_plot4 <- rbind(labels_temp, labels_plot4)
}

labels_plot4 <- labels_plot4 %>% mutate(model = ifelse(model == "Hazard Ratio (95% CI)", 
                                                       "Adjusted Hazard Ratio (95% CI)",
                                                       as.character(model)))

labels_plot4$model <- factor(labels_plot4$model, levels = rev(unique(labels_plot4$model)))

layout3 <- c(
  area(t = 0, l = 0, b = 6, r = 4), 
  area(t = 1, l = 5, b = 6, r = 9),
  area(t = 7, l = 0, b = 11, r = 4), 
  area(t = 7, l = 5, b = 11, r = 9)
)

# plot for albuminuria status
p_egfr40_alb <- 
  alb_hrs %>%
  filter(outcome == "ckd_egfr40") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
  coord_cartesian(ylim=c(1,length(unique(alb_hrs[alb_hrs$outcome == "ckd_egfr40",]$contrast)) + 1), 
                  xlim=c(0.33, 3)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(alb_hrs[alb_hrs$outcome == "ckd_egfr40",]$contrast)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.5, 
           y = length(unique(alb_hrs[alb_hrs$outcome == "ckd_egfr40",]$contrast)) + 1, 
           label = "Favours DPP4i/SU") +
  labs(x="", y="") +
  ggtitle("40% decline in eGFR or ESRD", subtitle = paste0(c("By albuminuria status (p = ", round(p_value_interaction_alb,2), " for interaction)"),collapse="")) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_alb <-
  labels_plot3 %>%
  filter(outcome == "ckd_egfr40") %>%
  ggplot(aes(y = (model))) + 
  geom_text(
    aes(x = 1, label = model),
    hjust = 0,
    fontface = ifelse(labels_plot3$
                        model == "Adjusted Hazard Ratio (95% CI)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))


# plot for egfr status
p_egfr40_egfr <- 
  egfr_hrs %>%
  filter(outcome == "ckd_egfr40") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.3, 2.0)) +
  coord_cartesian(ylim=c(1,length(unique(egfr_hrs[egfr_hrs$outcome == "ckd_egfr40",]$contrast)) + 1), 
                  xlim=c(0.33, 3)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(egfr_hrs[egfr_hrs$outcome == "ckd_egfr40",]$contrast)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.5, 
           y = length(unique(egfr_hrs[egfr_hrs$outcome == "ckd_egfr40",]$contrast)) + 1, 
           label = "Favours DPP4i/SU") +
  labs(x="Hazard Ratio", y="") +
  ggtitle(label = "", subtitle = paste0(c("By eGFR status (p = ", round(p_value_interaction_egfr,2), " for interaction)"),collapse="")) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_egfr <-
  labels_plot4 %>%
  filter(outcome == "ckd_egfr40") %>%
  ggplot(aes(y = (model))) + 
  geom_text(
    aes(x = 1, label = model),
    hjust = 0,
    fontface = ifelse(labels_plot4$
                        model == "Adjusted Hazard Ratio (95% CI)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

# final plot arrangement
p_left_alb + p_egfr40_alb + p_left_egfr + p_egfr40_egfr + plot_layout(design = layout3)

