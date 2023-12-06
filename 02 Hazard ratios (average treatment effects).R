## the question we are trying to answer in this script is
## whether the average treatment effects (HRs) for SGLT2 vs SU in CPRD are comparable to those from trials

## for reference, the trial meta-analysis HR reported in Lancet. 2022 Nov 19; 400(10365): 1788–1801 is:
# HR 0.62 (0.56 - 0.68)

## Here we will be calculating HRs for SGLT2 and DPP4 vs SU for:
## CKD: onset of CKD stage 3a-5
## CKD: fall in eGFR of <=40% from baseline or onset of CKD stage 5
## all-cause mortality

## originally, there were plans to add other outcomes, e.g. 
## prescription of medications for kidney disease (e.g. injections for renal anaemia, phosphate binders).
## however, due to the low number of subjects progressing to advanced CKD we are not looking into this at present.

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
load("2023-11-08_t2d_ckdpc_imputed_data.Rda")

n.imp <- 10
set.seed(123)
today <- as.character(Sys.Date(), format="%Y%m%d")

############################1 CALCULATE WEIGHTS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in them, then stack them again at the end
temp$IPTW <- temp$overlap <- NA

# take non-imputed dataset to append the imputed datasets to later on
temp2 <- temp[temp$.imp == 0,]

# propensity score formula
ps.formula <- formula(paste("studydrug ~ 

                               #sociodemographic characteristics:
                               dstartdate_age + malesex + imd2015_10 + 
                               ethnicity_qrisk2 + regstartdate + 
                               initiation_year +
                               
                               #laboratory and vital signs measurements:
                               preweight + prebmi + prehba1c + pretriglyceride + prehdl + preldl +
                               pretotalcholesterol + prealt + preegfr + uacr + presbp + predbp + 
                               ckdpc_egfr60_confirmed_score + ckdpc_40egfr_score +
                               
                               #medications:
                               ncurrtx + MFN + statin +
                               ACEi + ARB + BB + CCB + ThZD + loopD + MRA + 
                               steroids + immunosuppr +
                               
                               #Comorbidities
                               qrisk2_smoking_cat + dstartdate_dm_dur_all + 
                               predrug_hypertension + osteoporosis + predrug_dka +
                               predrug_falls + predrug_urinary_frequency + predrug_volume_depletion +
                               predrug_acutepancreatitis +
                               predrug_micturition_control + predrug_dementia + hosp_admission_prev_year"))


#calculate weights in each imputed dataset
for (i in 1:n.imp) {
  
  print(paste0("Calculating weights for imputed dataset number ", i))
  
  imp.x <- temp[temp$.imp == i,]
  
  w.overlap <- SumStat(ps.formula=ps.formula, 
                       data = imp.x, 
                       weight = c("IPW", "overlap"))
  
  ## Add weights to data frame
  weights <- w.overlap$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
  imp.x$IPTW <- weights$IPW
  imp.x$overlap <- weights$overlap
  
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
save(cohort, file=paste0(today, "_t2d_ckdpc_imputed_data_withweights.Rda"))
#load("2023-11-28_t2d_ckdpc_imputed_data_withweights.Rda")
############################2 CALCULATE HAZARD RATIOS################################################################

## 2 calculate hazard ratios (unadjusted, adjusted, weighted) and n events per study drug

#outcomes to be studied:
kf_outcomes <- c("ckd_345", "ckd_egfr40", "death", "ckd_345_pp", "ckd_egfr40_pp", "death_pp")

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

#create empty data frame to which we can append the hazard ratios once calculated
all_sglt2_hrs <- 
  all_dpp4_hrs <-
  all_sglt2vsdpp4_hrs <- 
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
                               ethnicity_qrisk2 + regstartdate + 
                               initiation_year +
                               
                               #laboratory and vital signs measurements:
                               preweight + prebmi + prehba1c + pretriglyceride + prehdl + preldl +
                               pretotalcholesterol + prealt + preegfr + uacr + presbp + predbp + 
                               ckdpc_egfr60_confirmed_score + ckdpc_40egfr_score +
                               
                               #medications:
                               ncurrtx + MFN + statin +
                               ACEi + ARB + BB + CCB + ThZD + loopD + MRA + 
                               steroids + immunosuppr +
                               
                               #medications:
                               ncurrtx + MFN + statin +
                               ACEi + ARB + BB + CCB + ThZD + loopD + MRA + 
                               steroids + immunosuppr +
                               
                               #Comorbidities
                               qrisk2_smoking_cat + dstartdate_dm_dur_all + 
                               predrug_hypertension + osteoporosis + predrug_dka +
                               predrug_falls + predrug_urinary_frequency + predrug_volume_depletion +
                               predrug_acutepancreatitis +
                               predrug_micturition_control + predrug_dementia + hosp_admission_prev_year"))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.SGLT2.unadj <- SE.SGLT2.unadj <-
    COEFS.DPP4.unadj <- SE.DPP4.unadj <-
    COEFS.SGLT2vsDPP4.unadj <- SE.SGLT2vsDPP4.unadj <-
  # for the adjusted survival models
    COEFS.SGLT2.adj <- SE.SGLT2.adj <-
    COEFS.DPP4.adj <- SE.DPP4.adj <-
    COEFS.SGLT2vsDPP4.adj <- SE.SGLT2vsDPP4.adj <-
  # for the overlap weighted survival models
    COEFS.SGLT2.ow <- SE.SGLT2.ow <-
    COEFS.DPP4.ow <- SE.DPP4.ow <-
    COEFS.SGLT2vsDPP4.ow <- SE.SGLT2vsDPP4.ow <-
  # for the inverse probability of treatment weighted survival models
    COEFS.SGLT2.iptw <- SE.SGLT2.iptw <-
    COEFS.DPP4.iptw <- SE.DPP4.iptw <-
    COEFS.SGLT2vsDPP4.iptw <- SE.SGLT2vsDPP4.iptw <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.DPP4.unadj[i] <- fit.unadj$coefficients[1]
    COEFS.SGLT2.unadj[i] <- fit.unadj$coefficients[2]
    SE.DPP4.unadj[i] <- sqrt(fit.unadj$var[1,1])
    SE.SGLT2.unadj[i] <- sqrt(fit.unadj$var[2,2])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted, cohort[cohort$.imp == i,])
    
    COEFS.DPP4.adj[i] <- fit.adj$coefficients[1]
    COEFS.SGLT2.adj[i] <- fit.adj$coefficients[2]
    SE.DPP4.adj[i] <- sqrt(fit.adj$var[1,1])
    SE.SGLT2.adj[i] <- sqrt(fit.adj$var[2,2])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = overlap)
    
    COEFS.DPP4.ow[i] <- fit.ow$coefficients[1]
    COEFS.SGLT2.ow[i] <- fit.ow$coefficients[2]
    SE.DPP4.ow[i] <- sqrt(fit.ow$var[1,1])
    SE.SGLT2.ow[i] <- sqrt(fit.ow$var[2,2])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = IPTW)
    
    COEFS.DPP4.iptw[i] <- fit.iptw$coefficients[1]
    COEFS.SGLT2.iptw[i] <- fit.iptw$coefficients[2]
    SE.DPP4.iptw[i] <- sqrt(fit.iptw$var[1,1])
    SE.SGLT2.iptw[i] <- sqrt(fit.iptw$var[2,2])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
    
    #calculate HRs for SGLT2 vs DPP4
    cohort$studydrug <- relevel(factor(cohort$studydrug), ref = "DPP4")
    
    #unadjusted analyses first
    fit.unadj <- coxph(f, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.SGLT2vsDPP4.unadj[i] <- fit.unadj$coefficients[2]
    SE.SGLT2vsDPP4.unadj[i] <- sqrt(fit.unadj$var[2,2])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted, cohort[cohort$.imp == i,])
    
    COEFS.SGLT2vsDPP4.adj[i] <- fit.adj$coefficients[2]
    SE.SGLT2vsDPP4.adj[i] <- sqrt(fit.adj$var[2,2])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = overlap)
    
    COEFS.SGLT2vsDPP4.ow[i] <- fit.ow$coefficients[2]
    SE.SGLT2vsDPP4.ow[i] <- sqrt(fit.ow$var[2,2])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = IPTW)
    
    COEFS.SGLT2vsDPP4.iptw[i] <- fit.iptw$coefficients[2]
    SE.SGLT2vsDPP4.iptw[i] <- sqrt(fit.iptw$var[2,2])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
    
    cohort$studydrug <- relevel(factor(cohort$studydrug), ref = "SU")
    
  }
  
  # pool hazard ratios
  unadjusted_sglt2 <- pool.rubin.HR(COEFS.SGLT2.unadj, SE.SGLT2.unadj, n.imp)
  unadjusted_dpp4 <- pool.rubin.HR(COEFS.DPP4.unadj, SE.DPP4.unadj, n.imp)
  unadjusted_sglt2vsdpp4 <- pool.rubin.HR(COEFS.SGLT2vsDPP4.unadj, SE.SGLT2vsDPP4.unadj, n.imp)
  
  adjusted_sglt2 <- pool.rubin.HR(COEFS.SGLT2.adj, SE.SGLT2.adj, n.imp)
  adjusted_dpp4 <- pool.rubin.HR(COEFS.DPP4.adj, SE.DPP4.adj, n.imp)  
  adjusted_sglt2vsdpp4 <- pool.rubin.HR(COEFS.SGLT2vsDPP4.adj, SE.SGLT2vsDPP4.adj, n.imp)
  
  overlapweighted_sglt2 <- pool.rubin.HR(COEFS.SGLT2.ow, SE.SGLT2.ow, n.imp)
  overlapweighted_dpp4 <- pool.rubin.HR(COEFS.DPP4.ow, SE.DPP4.ow, n.imp)
  overlapweighted_sglt2vsdpp4 <- pool.rubin.HR(COEFS.SGLT2vsDPP4.ow, SE.SGLT2vsDPP4.ow, n.imp)
  
  iptw_sglt2 <- pool.rubin.HR(COEFS.SGLT2.iptw, SE.SGLT2.iptw, n.imp)
  iptw_dpp4 <- pool.rubin.HR(COEFS.DPP4.iptw, SE.DPP4.iptw, n.imp)
  iptw_sglt2vsdpp4 <- pool.rubin.HR(COEFS.SGLT2vsDPP4.iptw, SE.SGLT2vsDPP4.iptw, n.imp)
  
  
  # save pooled HR and 95% confidence interval
  unadjusted_sglt2_string <- paste0(round(unadjusted_sglt2[1], 2), " (", round(unadjusted_sglt2[2], 2), ", ", round(unadjusted_sglt2[3], 2), ")")
  unadjusted_dpp4_string <- paste0(round(unadjusted_dpp4[1], 2), " (", round(unadjusted_dpp4[2], 2), ", ", round(unadjusted_dpp4[3], 2), ")")
  unadjusted_sglt2vsdpp4_string <- paste0(round(unadjusted_sglt2vsdpp4[1], 2), " (", round(unadjusted_sglt2vsdpp4[2], 2), ", ", round(unadjusted_sglt2vsdpp4[3], 2), ")")

  adjusted_sglt2_string <- paste0(round(adjusted_sglt2[1], 2), " (", round(adjusted_sglt2[2], 2), ", ", round(adjusted_sglt2[3], 2), ")")
  adjusted_dpp4_string <- paste0(round(adjusted_dpp4[1], 2), " (", round(adjusted_dpp4[2], 2), ", ", round(adjusted_dpp4[3], 2), ")")
  adjusted_sglt2vsdpp4_string <- paste0(round(adjusted_sglt2vsdpp4[1], 2), " (", round(adjusted_sglt2vsdpp4[2], 2), ", ", round(adjusted_sglt2vsdpp4[3], 2), ")")
  
  overlapweighted_sglt2_string <- paste0(round(overlapweighted_sglt2[1], 2), " (", round(overlapweighted_sglt2[2], 2), ", ", round(overlapweighted_sglt2[3], 2), ")")
  overlapweighted_dpp4_string <- paste0(round(overlapweighted_dpp4[1], 2), " (", round(overlapweighted_dpp4[2], 2), ", ", round(overlapweighted_dpp4[3], 2), ")")
  overlapweighted_sglt2vsdpp4_string <- paste0(round(overlapweighted_sglt2vsdpp4[1], 2), " (", round(overlapweighted_sglt2vsdpp4[2], 2), ", ", round(overlapweighted_sglt2vsdpp4[3], 2), ")")

  iptw_sglt2_string <- paste0(round(iptw_sglt2[1], 2), " (", round(iptw_sglt2[2], 2), ", ", round(iptw_sglt2[3], 2), ")")
  iptw_dpp4_string <- paste0(round(iptw_dpp4[1], 2), " (", round(iptw_dpp4[2], 2), ", ", round(iptw_dpp4[3], 2), ")")
  iptw_sglt2vsdpp4_string <- paste0(round(iptw_sglt2vsdpp4[1], 2), " (", round(iptw_sglt2vsdpp4[2], 2), ", ", round(iptw_sglt2vsdpp4[3], 2), ")")
  
  
  # combine in dataframe that we can tabulate
  outcome_sglt2_hr <- cbind(outcome=k, count[c(1,3)], followup[c(1,3)], events[c(1,3)], 
                            unadjusted_sglt2_string, adjusted_sglt2_string, overlapweighted_sglt2_string, iptw_sglt2_string)
  
  outcome_dpp4_hr <- cbind(outcome=k, count[1:2], followup[1:2], events[1:2], 
                            unadjusted_dpp4_string, adjusted_dpp4_string, overlapweighted_dpp4_string, iptw_dpp4_string)
  
  outcome_sglt2vsdpp4_hr <- cbind(outcome=k, count[2:3], followup[2:3], events[2:3], 
                            unadjusted_sglt2vsdpp4_string, adjusted_sglt2vsdpp4_string, overlapweighted_sglt2vsdpp4_string, iptw_sglt2vsdpp4_string)
  
  all_sglt2_hrs <- rbind(all_sglt2_hrs, outcome_sglt2_hr)
  all_dpp4_hrs <- rbind(all_dpp4_hrs, outcome_dpp4_hr)
  all_sglt2vsdpp4_hrs <- rbind(all_sglt2vsdpp4_hrs, outcome_sglt2vsdpp4_hr)
  
  outcome_hr <- rbind(
    cbind(outcome = k, contrast = "SGLT2 vs SU", analysis = "unadjusted", 
          HR = unadjusted_sglt2[1], LB = unadjusted_sglt2[2], UB = unadjusted_sglt2[3], string = unadjusted_sglt2_string),
    cbind(outcome = k, contrast = "SGLT2 vs SU", analysis = "adjusted", 
          HR = adjusted_sglt2[1], LB = adjusted_sglt2[2], UB = adjusted_sglt2[3], string = adjusted_sglt2_string),
    cbind(outcome = k, contrast = "SGLT2 vs SU", analysis = "overlap-weighted", 
          HR = overlapweighted_sglt2[1], LB = overlapweighted_sglt2[2], UB = overlapweighted_sglt2[3], string = overlapweighted_sglt2_string),
    cbind(outcome = k, contrast = "SGLT2 vs SU", analysis = "IPTW", 
          HR = iptw_sglt2[1], LB = iptw_sglt2[2], UB = iptw_sglt2[3], string = iptw_sglt2_string),
    cbind(outcome = k, contrast = "DPP4 vs SU", analysis = "unadjusted", 
          HR = unadjusted_dpp4[1], LB = unadjusted_dpp4[2], UB = unadjusted_dpp4[3], string = unadjusted_dpp4_string),
    cbind(outcome = k, contrast = "DPP4 vs SU", analysis = "adjusted", 
          HR = adjusted_dpp4[1], LB = adjusted_dpp4[2], UB = adjusted_dpp4[3], string = adjusted_dpp4_string),
    cbind(outcome = k, contrast = "DPP4 vs SU", analysis = "overlap-weighted", 
          HR = overlapweighted_dpp4[1], LB = overlapweighted_dpp4[2], UB = overlapweighted_dpp4[3], string = overlapweighted_dpp4_string),
    cbind(outcome = k, contrast = "DPP4 vs SU", analysis = "IPTW", 
          HR = iptw_dpp4[1], LB = iptw_dpp4[2], UB = iptw_dpp4[3], string = iptw_dpp4_string),
    cbind(outcome = k, contrast = "SGLT2 vs DPP4", analysis = "unadjusted", 
          HR = unadjusted_sglt2vsdpp4[1], LB = unadjusted_sglt2vsdpp4[2], UB = unadjusted_sglt2vsdpp4[3], string = unadjusted_sglt2vsdpp4_string),
    cbind(outcome = k, contrast = "SGLT2 vs DPP4", analysis = "adjusted", 
          HR = adjusted_sglt2vsdpp4[1], LB = adjusted_sglt2vsdpp4[2], UB = adjusted_sglt2vsdpp4[3], string = adjusted_sglt2vsdpp4_string),
    cbind(outcome = k, contrast = "SGLT2 vs DPP4", analysis = "overlap-weighted", 
          HR = overlapweighted_sglt2vsdpp4[1], LB = overlapweighted_sglt2vsdpp4[2], UB = overlapweighted_sglt2vsdpp4[3], string = overlapweighted_sglt2vsdpp4_string),
    cbind(outcome = k, contrast = "SGLT2 vs DPP4", analysis = "IPTW", 
          HR = iptw_sglt2vsdpp4[1], LB = iptw_sglt2vsdpp4[2], UB = iptw_sglt2vsdpp4[3], string = iptw_sglt2vsdpp4_string)
  )
  
  all_hrs <- rbind(all_hrs, outcome_hr)
  
}

# show table with events, follow up time, and hazard ratios
flextable(all_sglt2_hrs)
flextable(all_dpp4_hrs)
flextable(all_sglt2vsdpp4_hrs)


save(all_hrs, file=paste0(today, "_all_hrs.Rda"))
