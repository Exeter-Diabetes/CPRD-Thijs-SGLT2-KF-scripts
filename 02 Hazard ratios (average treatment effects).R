## In this script our aim is to see whether the average treatment for SGLT2i vs SU is similar
## in the literature and in this CPRD cohort
## we will compare hazard ratios for SGLT2i vs SU in CPRD to those from trials.

## for reference, the trial meta-analysis RR reported in Lancet (Nuffield Group 2022) is:
# HR 0.62 (0.56 - 0.68)

## Here we will be calculating HRs for SGLT2i and DPP4i vs SU for:
## primary outcome kidney disease progression: fall in eGFR of ≥50% from baseline or onset of CKD stage 5
## secondary outcomes: progression to macroalbuminuria, all-cause mortality, dka, amputation, and mycotic genital infections

## Contents:
# 0 setup
# 1 calculate weights (IPTW, overlap weights)
# 2 calculate hazard ratios
# 3 store hazard ratios
# 4 subgroup analyses (by albuminuria and eGFR status)
# 5 forest plots

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

n.imp <- 10
set.seed(123)
#today <- as.character(Sys.Date(), format="%Y%m%d")
today <- "2024-07-13"
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

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-07-13_t2d_ckdpc_imputed_data.Rda")

covariates <- c("dstartdate_age", "malesex", "imd2015_10", "ethnicity_4cat", "initiation_year", "prebmi", "prehba1c",
                "pretotalcholesterol", "preegfr", "uacr", "presbp", "ckdpc_50egfr_score", "ncurrtx", "statin", "INS", 
                "ACEi_or_ARB", "smoking_status", "dstartdate_dm_dur_all", "predrug_hypertension", "predrug_af", "hosp_admission_prev_year")

# we exclude initiation_year from the propensity score model
covariates_ps <- covariates[-5]

############################1 CALCULATE WEIGHTS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in them, then stack them again at the end
temp$IPTW <- temp$overlap <- 
temp$IPTW2 <- temp$overlap2 <- NA

# take non-imputed dataset to append the imputed datasets to later on
temp2 <- temp[temp$.imp == 0,]

# propensity score formula
ps.formula <- formula(paste("studydrug ~ ", paste(covariates_ps, collapse=" + ")))

ps.formula2 <- formula(paste("studydrug2 ~ ", paste(covariates_ps, collapse=" + ")))


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
  
  # truncate IPTW at 2nd and 98% percentile
  w.overlap$ps.weights$IPW <- ifelse(w.overlap$ps.weights$IPW < quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1], 
                                     quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                     ifelse(
                                       w.overlap$ps.weights$IPW > quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                       quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                       w.overlap$ps.weights$IPW
                                     ))
  
  
  w.overlap2$ps.weights$IPW <- ifelse(w.overlap2$ps.weights$IPW < quantile(w.overlap2$ps.weights$IPW, probs = c(0.02, 0.98))[1], 
                                      quantile(w.overlap2$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                      ifelse(
                                        w.overlap2$ps.weights$IPW > quantile(w.overlap2$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                        quantile(w.overlap2$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                        w.overlap2$ps.weights$IPW
                                      ))
  
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

#plot SMD in weighted / unweighted populations (love plot - shows balance by variable)
plot(w.overlap2)
#histogram of propensity scores by treatment group (I have amended the code from the package to fit with our colours)
plot.hist.ps.2drugs<-function(x, weighted.var=TRUE, threshold=0.1, metric="ASD", breaks=50,...){
#get object info
m<-length(names(x$ps.weights)[-c(1,2)])+1
zname<-names(x$ps.weights)[1]
Z<-unlist(x$ps.weights[zname])
ncate<-length(unique(Z))
metric<-toupper(metric)

#extract original treatment labels
dic0<-rep(NA,ncate)
for (i in 1:ncate){
  dic0[i]<-as.character(x$ps.weights[which(x$ps.weights$zindex==i)[1],1])
}
e<-x$propensity[,2]
z<-x$ps.weights$zindex
#enough room for 8 character-strings in legend
par(mar=c(5,4,4,10.1),xpd=TRUE)
he1<-hist(e[z==1],breaks=breaks,plot = FALSE)
ylim1<-max(he1$counts)+0.5
he2<-hist(e[z==2],breaks=breaks,plot = FALSE)
ylim2<-max(he2$counts+0.5)
ylims<-max(c(ylim1,ylim2))

hist(e[z==1],breaks=breaks,col="#0072B2",border="black",bty='L',
     main=NULL,ylab=NULL,xlim=c(max(min(e)-0.2,0),min(max(e)+0.2,1)),xlab="Estimated propensity score",
     freq=T,cex.lab = 1.5, cex.axis = 1.5 ,cex.main = 2, cex = 2,bty='L',ylim = c(0,ylims))
hist(e[z==2],breaks=breaks,add=TRUE,col=rgb(230/255, 159/255, 0, alpha = 0.9),freq=T)

legend("right",inset=c(-0.2, 0), title="",legend=dic0,col=c("#0072B2","#E69F00"),
       lty=1,lwd=1.5,bty='n',cex=1.5,seg.len = 0.65,xjust = 1)
}
plot.hist.ps.2drugs(w.overlap2)

summary(w.overlap, weighted.var = TRUE, metric = "ASD")

# save dataset with weights so this can be used in subsequent scripts
cohort <- temp
rm(temp)
cohort <- cohort %>% filter(!.imp == 0)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_imputed_data_withweights.Rda"))
#load("2024-07-13_t2d_ckdpc_imputed_data_withweights.Rda")
############################2 CALCULATE HAZARD RATIOS################################################################

## 2 calculate hazard ratios (unadjusted, adjusted, weighted) and n events per study drug

#outcomes to be studied:
outcomes_per_drugclass <- c("ckd_egfr50", "ckd_egfr50_pp")

kf_key_outcomes <- c("ckd_egfr40", "ckd_egfr50", "death", "macroalb", "dka", "amputation", "side_effect")

#create empty data frame to which we can append the hazard ratios once calculated
all_sglt2i_hrs <- 
  all_dpp4i_hrs <-
  all_SGLT2ivsDPP4i_hrs <- 
  all_hrs <- 
  data.frame()

for (k in outcomes_per_drugclass) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  # calculate number of subjects in each group
  count <- cohort %>%
    group_by(studydrug) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort %>%
    group_by(studydrug) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- cohort %>%
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
  
  f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug + ", paste(covariates, collapse=" + ")))
  
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
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "Unadjusted", 
          HR = unadjusted_SGLT2i[1], LB = unadjusted_SGLT2i[2], UB = unadjusted_SGLT2i[3], 
          string = unadjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "Adjusted", 
          HR = adjusted_SGLT2i[1], LB = adjusted_SGLT2i[2], UB = adjusted_SGLT2i[3], 
          string = adjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "Overlap-weighted", 
          HR = overlapweighted_SGLT2i[1], LB = overlapweighted_SGLT2i[2], UB = overlapweighted_SGLT2i[3], 
          string = overlapweighted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "IPTW", 
          HR = iptw_SGLT2i[1], LB = iptw_SGLT2i[2], UB = iptw_SGLT2i[3], 
          string = iptw_SGLT2i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "Unadjusted", 
          HR = unadjusted_DPP4i[1], LB = unadjusted_DPP4i[2], UB = unadjusted_DPP4i[3], 
          string = unadjusted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "Adjusted", 
          HR = adjusted_DPP4i[1], LB = adjusted_DPP4i[2], UB = adjusted_DPP4i[3], 
          string = adjusted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "Overlap-weighted", 
          HR = overlapweighted_DPP4i[1], LB = overlapweighted_DPP4i[2], UB = overlapweighted_DPP4i[3], 
          string = overlapweighted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "IPTW", 
          HR = iptw_DPP4i[1], LB = iptw_DPP4i[2], UB = iptw_DPP4i[3], 
          string = iptw_DPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "Unadjusted", 
          HR = unadjusted_SGLT2ivsDPP4i[1], LB = unadjusted_SGLT2ivsDPP4i[2], UB = unadjusted_SGLT2ivsDPP4i[3], 
          string = unadjusted_SGLT2ivsDPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "Adjusted", 
          HR = adjusted_SGLT2ivsDPP4i[1], LB = adjusted_SGLT2ivsDPP4i[2], UB = adjusted_SGLT2ivsDPP4i[3], 
          string = adjusted_SGLT2ivsDPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "Overlap-weighted", 
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

# store n/N for every drug class in all_hrs
all_hrs <- all_hrs %>% left_join (all_dpp4i_hrs %>% select(1:7)) %>% left_join(all_sglt2i_hrs %>% select(1,3,5,7))



## compare SGLT2i vs combined group of DPP4i/SU

#create empty data frame to which we can append the hazard ratios once calculated
all_SGLT2ivsDPP4iSU_hrs <- SGLT2ivsDPP4iSU_hrs <- data.frame()

## remove double overlapping entries for DPP4i and SU that overlap (take one only)
cohort <- cohort %>% group_by(.imp, patid) %>% filter(
  !duplicated(studydrug2)
) %>% ungroup()


for (k in kf_key_outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  if (k == "macroalb") {
    temp <- cohort[cohort$uacr < 30,]
  } else {
    temp <- cohort
  }
  
  # calculate number of subjects in each group
  count <- temp[temp$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- temp[temp$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- temp[temp$.imp > 0,] %>%
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
  
  f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", paste(covariates, collapse=" + ")))
  
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
    fit.unadj <- coxph(f2, temp[temp$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.SGLT2i.unadj[i] <- fit.unadj$coefficients[1]
    SE.SGLT2i.unadj[i] <- sqrt(fit.unadj$var[1,1])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted2, temp[temp$.imp == i,])
    
    COEFS.SGLT2i.adj[i] <- fit.adj$coefficients[1]
    SE.SGLT2i.adj[i] <- sqrt(fit.adj$var[1,1])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted2, temp[temp$.imp ==i,], weights = overlap2)
    
    COEFS.SGLT2i.ow[i] <- fit.ow$coefficients[1]
    SE.SGLT2i.ow[i] <- sqrt(fit.ow$var[1,1])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted2, temp[temp$.imp ==i,], weights = IPTW2)
    
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
  outcome_SGLT2ivsDPP4iSU_hr <- cbind(outcome=k, count[c(1:2)], followup[c(1:2)], events[c(1:2)],
                                      unadjusted_SGLT2i_string, adjusted_SGLT2i_string, overlapweighted_SGLT2i_string, iptw_SGLT2i_string)
  
  SGLT2ivsDPP4iSU_hrs <- rbind(SGLT2ivsDPP4iSU_hrs, outcome_SGLT2ivsDPP4iSU_hr)
  
  SGLT2ivsDPP4iSU_hr <- rbind(
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Unadjusted", 
          HR = unadjusted_SGLT2i[1], LB = unadjusted_SGLT2i[2], UB = unadjusted_SGLT2i[3], string = unadjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Adjusted", 
          HR = adjusted_SGLT2i[1], LB = adjusted_SGLT2i[2], UB = adjusted_SGLT2i[3], string = adjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Overlap-weighted", 
          HR = overlapweighted_SGLT2i[1], LB = overlapweighted_SGLT2i[2], UB = overlapweighted_SGLT2i[3], string = overlapweighted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "IPTW", 
          HR = iptw_SGLT2i[1], LB = iptw_SGLT2i[2], UB = iptw_SGLT2i[3], string = iptw_SGLT2i_string))
    
    all_SGLT2ivsDPP4iSU_hrs <- rbind(all_SGLT2ivsDPP4iSU_hrs, SGLT2ivsDPP4iSU_hr)
  
    rm(temp)
}
############################3 STORE AND DISPLAY HAZARD RATIOS################################################################

# save all_hrs table and SGLT2i vs DPP4i/su table
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")

save(all_hrs, file=paste0(today, "_all_hrs.Rda"))
save(all_SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_all_SGLT2ivsDPP4iSU_hrs.Rda"))
save(SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_SGLT2ivsDPP4iSU_hrs.Rda"))
# load("2024-07-13_all_hrs.Rda")
# load("2024-07-13_all_SGLT2ivsDPP4iSU_hrs.Rda")
# load("2024-07-13_SGLT2ivsDPP4iSU_hrs.Rda")

# show table with events, follow up time, and hazard ratios
flextable(all_sglt2i_hrs)
flextable(all_dpp4i_hrs)
flextable(all_SGLT2ivsDPP4i_hrs)
flextable(all_SGLT2ivsDPP4iSU_hrs)
flextable(SGLT2ivsDPP4iSU_hrs)

all_hrs$model <- paste0(all_hrs$string, " [", all_hrs$analysis, "]")
all_hrs$model <- factor(all_hrs$model, levels = unique(all_hrs$model))

all_SGLT2ivsDPP4iSU_hrs$model <- paste0(all_SGLT2ivsDPP4iSU_hrs$string, " [", all_SGLT2ivsDPP4iSU_hrs$analysis, "]")
all_SGLT2ivsDPP4iSU_hrs$model <- factor(all_SGLT2ivsDPP4iSU_hrs$model, levels = unique(all_SGLT2ivsDPP4iSU_hrs$model))

# have to coerce HR and CI to class numeric as they sometimes default to character
class(all_hrs$HR) <- class(all_hrs$LB) <- class(all_hrs$UB) <-
  class(all_SGLT2ivsDPP4iSU_hrs$HR) <- class(all_SGLT2ivsDPP4iSU_hrs$LB) <- class(all_SGLT2ivsDPP4iSU_hrs$UB) <- "numeric"

# remove unadjusted HR as we do not want to plot these
all_hrs <- all_hrs[!all_hrs$analysis == "Unadjusted",]
all_SGLT2ivsDPP4iSU_hrs <- all_SGLT2ivsDPP4iSU_hrs[!all_SGLT2ivsDPP4iSU_hrs$analysis == "Unadjusted",]

#add event count and total count to hr dataframe
all_SGLT2ivsDPP4iSU_hrs <- all_SGLT2ivsDPP4iSU_hrs %>% left_join(SGLT2ivsDPP4iSU_hrs %>% select(1:7))


all_hrs <- all_hrs %>%
  separate(`DPP4i_events`, into = c("DPP4i_events_number", "DPP4i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(`SU_events`, into = c("SU_events_number", "SU_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(SGLT2i_events, into = c("SGLT2i_events_number", "SGLT2i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `DPP4i_events_percentage` = str_replace(`DPP4i_events_percentage`, "\\)", ""),
    SU_events_percentage = str_replace(SU_events_percentage, "\\)", ""),
    SGLT2i_events_percentage = str_replace(SGLT2i_events_percentage, "\\)", ""),
    `DPP4i_nN` = paste0(`DPP4i_events_number`, "/", `DPP4i_count`),
    `SU_nN` = paste0(`SU_events_number`, "/", `SU_count`),
    SGLT2i_nN = paste0(SGLT2i_events_number, "/", SGLT2i_count)
  )


all_SGLT2ivsDPP4iSU_hrs <- all_SGLT2ivsDPP4iSU_hrs %>%
  separate(`DPP4i/SU_events`, into = c("DPP4i/SU_events_number", "DPP4i/SU_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(SGLT2i_events, into = c("SGLT2i_events_number", "SGLT2i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `DPP4i/SU_events_percentage` = str_replace(`DPP4i/SU_events_percentage`, "\\)", ""),
    SGLT2i_events_percentage = str_replace(SGLT2i_events_percentage, "\\)", ""),
    `DPP4i/SU_nN` = paste0(`DPP4i/SU_events_number`, "/", `DPP4i/SU_count`),
    SGLT2i_nN = paste0(SGLT2i_events_number, "/", SGLT2i_count)
  )

############################4 SUBGROUP ANALYSES################################################################

cohort <- cohort %>% mutate(
  risk_group = factor(risk_group, levels=c("eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                           "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                                           "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                           "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                           "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                           "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol")))


#create empty data frame to which we can append the hazard ratios once calculated
subgroup_SGLT2ivsDPP4iSU_hrs <- subgroup_hrs <- data.frame()


## analyses stratified by risk group
for (k in kf_key_outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  if (!k == "macroalb") {
    
    # calculate number of subjects in each group
    count <- cohort %>%
      mutate(risk_group = factor(risk_group, levels=c("eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                                      "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))) %>%
      group_by(studydrug2,risk_group) %>%
      summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_count",
                  values_from=count)
    
    # calculate median follow up time (years) per group
    followup <- cohort %>%
      mutate(risk_group = factor(risk_group, levels=c("eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                                      "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))) %>%
      group_by(studydrug2,risk_group) %>%
      summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_followup",
                  values_from=time)
    
    # summarise number of events per group
    events <- cohort %>%
      mutate(risk_group = factor(risk_group, levels=c("eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                                      "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))) %>%
      group_by(studydrug2,risk_group) %>%
      summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
                drug_count=round(n()/n.imp, 0)) %>%
      mutate(events_perc=round(event_count*100/drug_count, 1),
             events=paste0(event_count, " (", events_perc, "%)")) %>%
      select(studydrug2, risk_group, events) %>%
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_events",
                  values_from=events)
    
    
    # write formulas for adjusted and unadjusted analyses
    f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*risk_group"))
    
    f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*risk_group + ", paste(covariates, collapse=" + ")))
    
    # create empty vectors to store the hazard ratios from every imputed dataset
    # for the unadjusted survival models
    COEFS.presegfr.noalb.unadj <- SE.presegfr.noalb.unadj <-
      COEFS.presegfr.microalb.unadj <- SE.presegfr.microalb.unadj <-
      COEFS.presegfr.macroalb.unadj <- SE.presegfr.macroalb.unadj <-      
      COEFS.redegfr.noalb.unadj <- SE.redegfr.noalb.unadj <-
      COEFS.redegfr.microalb.unadj <- SE.redegfr.microalb.unadj <-
      COEFS.redegfr.macroalb.unadj <- SE.redegfr.macroalb.unadj <-
      # for the adjusted survival models
      COEFS.presegfr.noalb.adj <- SE.presegfr.noalb.adj <-
      COEFS.presegfr.microalb.adj <- SE.presegfr.microalb.adj <-
      COEFS.presegfr.macroalb.adj <- SE.presegfr.macroalb.adj <-      
      COEFS.redegfr.noalb.adj <- SE.redegfr.noalb.adj <-
      COEFS.redegfr.microalb.adj <- SE.redegfr.microalb.adj <-
      COEFS.redegfr.macroalb.adj <- SE.redegfr.macroalb.adj <-
      rep(NA,n.imp)
    
    for (i in 1:n.imp) {
      print(paste0("Analyses in imputed dataset number ", i))
      
      #unadjusted analyses first
      fit.unadj <- coxph(f2, cohort[cohort$.imp == i,])
      
      #store coefficients and standard errors from this model
      COEFS.presegfr.noalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"]
      SE.presegfr.noalb.unadj[i] <- sqrt(fit.unadj$var[1,1])
      
      COEFS.presegfr.microalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:risk_groupeGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol"]
      SE.presegfr.microalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var)-4,nrow(fit.unadj$var)-4]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)-4])
      
      COEFS.presegfr.macroalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:risk_groupeGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol"]
      SE.presegfr.macroalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var)-3,nrow(fit.unadj$var)-3]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)-3])
      
      COEFS.redegfr.noalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR <3mg/mmol"]
      SE.redegfr.noalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var)-2,nrow(fit.unadj$var)-2]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)-2])
      
      COEFS.redegfr.microalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR 3-30mg/mmol"]
      SE.redegfr.microalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var)-1,nrow(fit.unadj$var)-1]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)-1])
      
      COEFS.redegfr.macroalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"]
      SE.redegfr.macroalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var),nrow(fit.unadj$var)]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)])
      
      #adjusted analyses
      fit.adj <- coxph(f_adjusted2, cohort[cohort$.imp == i,])
      
      COEFS.presegfr.noalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"]
      SE.presegfr.noalb.adj[i] <- sqrt(fit.adj$var[1,1])
      
      COEFS.presegfr.microalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:risk_groupeGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol"]
      SE.presegfr.microalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var)-4,nrow(fit.adj$var)-4]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)-4])
      
      COEFS.presegfr.macroalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:risk_groupeGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol"]
      SE.presegfr.macroalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var)-3,nrow(fit.adj$var)-3]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)-3])
      
      COEFS.redegfr.noalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR <3mg/mmol"]
      SE.redegfr.noalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var)-2,nrow(fit.adj$var)-2]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)-2])
      
      COEFS.redegfr.microalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR 3-30mg/mmol"]
      SE.redegfr.microalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var)-1,nrow(fit.adj$var)-1]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)-1])
      
      COEFS.redegfr.macroalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"]
      SE.redegfr.macroalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var),nrow(fit.adj$var)]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)])
      
      if (k == "ckd_egfr50") {
        if (i == n.imp) {
          f_adjusted3 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + risk_group + ", paste(covariates, collapse=" + ")))
          fit.no_interaction <- coxph(f_adjusted3, cohort[cohort$.imp == i,])
          
          loglikelihood_test <- anova(fit.no_interaction, fit.adj, test = "Chisq")
          p_value_interaction <- loglikelihood_test$`Pr(>|Chi|)`[2]      
        }}
      
    }
    
    # pool hazard ratios
    unadjusted_presegfr_noalb <- pool.rubin.HR(COEFS.presegfr.noalb.unadj, SE.presegfr.noalb.unadj, n.imp)
    adjusted_presegfr_noalb <- pool.rubin.HR(COEFS.presegfr.noalb.adj, SE.presegfr.noalb.adj, n.imp)
    
    unadjusted_presegfr_microalb <- pool.rubin.HR(COEFS.presegfr.microalb.unadj, SE.presegfr.microalb.unadj, n.imp)
    adjusted_presegfr_microalb <- pool.rubin.HR(COEFS.presegfr.microalb.adj, SE.presegfr.microalb.adj, n.imp)
    
    unadjusted_presegfr_macroalb <- pool.rubin.HR(COEFS.presegfr.macroalb.unadj, SE.presegfr.macroalb.unadj, n.imp)
    adjusted_presegfr_macroalb <- pool.rubin.HR(COEFS.presegfr.macroalb.adj, SE.presegfr.macroalb.adj, n.imp)
    
    unadjusted_redegfr_noalb <- pool.rubin.HR(COEFS.redegfr.noalb.unadj, SE.redegfr.noalb.unadj, n.imp)
    adjusted_redegfr_noalb <- pool.rubin.HR(COEFS.redegfr.noalb.adj, SE.redegfr.noalb.adj, n.imp)
    
    unadjusted_redegfr_microalb <- pool.rubin.HR(COEFS.redegfr.microalb.unadj, SE.redegfr.microalb.unadj, n.imp)
    adjusted_redegfr_microalb <- pool.rubin.HR(COEFS.redegfr.microalb.adj, SE.redegfr.microalb.adj, n.imp)
    
    unadjusted_redegfr_macroalb <- pool.rubin.HR(COEFS.redegfr.macroalb.unadj, SE.redegfr.macroalb.unadj, n.imp)
    adjusted_redegfr_macroalb <- pool.rubin.HR(COEFS.redegfr.macroalb.adj, SE.redegfr.macroalb.adj, n.imp)
    
    # save pooled HR and 95% confidence interval
    unadjusted_presegfr_noalb_string <- paste0(round(unadjusted_presegfr_noalb[1], 2), " (", round(unadjusted_presegfr_noalb[2], 2), ", ", round(unadjusted_presegfr_noalb[3], 2), ")")
    adjusted_presegfr_noalb_string <- paste0(round(adjusted_presegfr_noalb[1], 2), " (", round(adjusted_presegfr_noalb[2], 2), ", ", round(adjusted_presegfr_noalb[3], 2), ")")
    
    unadjusted_presegfr_microalb_string <- paste0(round(unadjusted_presegfr_microalb[1], 2), " (", round(unadjusted_presegfr_microalb[2], 2), ", ", round(unadjusted_presegfr_microalb[3], 2), ")")
    adjusted_presegfr_microalb_string <- paste0(round(adjusted_presegfr_microalb[1], 2), " (", round(adjusted_presegfr_microalb[2], 2), ", ", round(adjusted_presegfr_microalb[3], 2), ")")
    
    unadjusted_presegfr_macroalb_string <- paste0(round(unadjusted_presegfr_macroalb[1], 2), " (", round(unadjusted_presegfr_macroalb[2], 2), ", ", round(unadjusted_presegfr_macroalb[3], 2), ")")
    adjusted_presegfr_macroalb_string <- paste0(round(adjusted_presegfr_macroalb[1], 2), " (", round(adjusted_presegfr_macroalb[2], 2), ", ", round(adjusted_presegfr_macroalb[3], 2), ")")
    
    unadjusted_redegfr_noalb_string <- paste0(round(unadjusted_redegfr_noalb[1], 2), " (", round(unadjusted_redegfr_noalb[2], 2), ", ", round(unadjusted_redegfr_noalb[3], 2), ")")
    adjusted_redegfr_noalb_string <- paste0(round(adjusted_redegfr_noalb[1], 2), " (", round(adjusted_redegfr_noalb[2], 2), ", ", round(adjusted_redegfr_noalb[3], 2), ")")
    
    unadjusted_redegfr_microalb_string <- paste0(round(unadjusted_redegfr_microalb[1], 2), " (", round(unadjusted_redegfr_microalb[2], 2), ", ", round(unadjusted_redegfr_microalb[3], 2), ")")
    adjusted_redegfr_microalb_string <- paste0(round(adjusted_redegfr_microalb[1], 2), " (", round(adjusted_redegfr_microalb[2], 2), ", ", round(adjusted_redegfr_microalb[3], 2), ")")
    
    unadjusted_redegfr_macroalb_string <- paste0(round(unadjusted_redegfr_macroalb[1], 2), " (", round(unadjusted_redegfr_macroalb[2], 2), ", ", round(unadjusted_redegfr_macroalb[3], 2), ")")
    adjusted_redegfr_macroalb_string <- paste0(round(adjusted_redegfr_macroalb[1], 2), " (", round(adjusted_redegfr_macroalb[2], 2), ", ", round(adjusted_redegfr_macroalb[3], 2), ")")
    
    # combine in dataframe that we can tabulate
    presegfr_noalb_hr <- cbind(outcome=k, count[1,c(2:3)], followup[1,c(2:3)], events[1,c(2:3)],
                               unadjusted=unadjusted_presegfr_noalb_string, adjusted=adjusted_presegfr_noalb_string
    )
    presegfr_microalb_hr <- cbind(outcome=k, count[2,c(2:3)], followup[2,c(2:3)], events[2,c(2:3)],
                                  unadjusted=unadjusted_presegfr_microalb_string, adjusted=adjusted_presegfr_microalb_string
    )
    presegfr_macroalb_hr <- cbind(outcome=k, count[3,c(2:3)], followup[3,c(2:3)], events[3,c(2:3)],
                                  unadjusted=unadjusted_presegfr_macroalb_string, adjusted=adjusted_presegfr_macroalb_string
    )
    redegfr_noalb_hr <- cbind(outcome=k, count[4,c(2:3)], followup[4,c(2:3)], events[4,c(2:3)],
                              unadjusted=unadjusted_redegfr_noalb_string, adjusted=adjusted_redegfr_noalb_string
    )
    redegfr_microalb_hr <- cbind(outcome=k, count[5,c(2:3)], followup[5,c(2:3)], events[5,c(2:3)],
                                 unadjusted=unadjusted_redegfr_microalb_string, adjusted=adjusted_redegfr_microalb_string
    )
    redegfr_macroalb_hr <- cbind(outcome=k, count[6,c(2:3)], followup[6,c(2:3)], events[6,c(2:3)],
                                 unadjusted=unadjusted_redegfr_macroalb_string, adjusted=adjusted_redegfr_macroalb_string
    )
    
    outcome_subgroup_SGLT2ivsDPP4iSU_hrs <- rbind(presegfr_noalb_hr, presegfr_microalb_hr, presegfr_macroalb_hr,
                                                  redegfr_noalb_hr, redegfr_microalb_hr, redegfr_macroalb_hr)
    
    temp <- rbind(
      cbind(outcome = k, contrast = "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", analysis = "Adjusted",
            HR = adjusted_presegfr_noalb[1], LB = adjusted_presegfr_noalb[2], UB = adjusted_presegfr_noalb[3], string = adjusted_presegfr_noalb_string),
      cbind(outcome = k, contrast = "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol", analysis = "Adjusted",
            HR = adjusted_presegfr_microalb[1], LB = adjusted_presegfr_microalb[2], UB = adjusted_presegfr_microalb[3], string = adjusted_presegfr_microalb_string),
      cbind(outcome = k, contrast = "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", analysis = "Adjusted",
            HR = adjusted_presegfr_macroalb[1], LB = adjusted_presegfr_macroalb[2], UB = adjusted_presegfr_macroalb[3], string = adjusted_presegfr_macroalb_string),
      cbind(outcome = k, contrast = "eGFR <60mL/min/1.73m2, uACR <3mg/mmol", analysis = "Adjusted",
            HR = adjusted_redegfr_noalb[1], LB = adjusted_redegfr_noalb[2], UB = adjusted_redegfr_noalb[3], string = adjusted_redegfr_noalb_string),
      cbind(outcome = k, contrast = "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol", analysis = "Adjusted",
            HR = adjusted_redegfr_microalb[1], LB = adjusted_redegfr_microalb[2], UB = adjusted_redegfr_microalb[3], string = adjusted_redegfr_microalb_string),
      cbind(outcome = k, contrast = "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol", analysis = "Adjusted",
            HR = adjusted_redegfr_macroalb[1], LB = adjusted_redegfr_macroalb[2], UB = adjusted_redegfr_macroalb[3], string = adjusted_redegfr_macroalb_string)
    )
    
    subgroup_hrs <- rbind(subgroup_hrs, temp)
    subgroup_SGLT2ivsDPP4iSU_hrs <- rbind(subgroup_SGLT2ivsDPP4iSU_hrs, outcome_subgroup_SGLT2ivsDPP4iSU_hrs)
    
  } 
  else 
  {
    
    temp2 <- cohort[cohort$uacr < 30,]
    
    # calculate number of subjects in each group
    count <- temp2 %>%
      mutate(risk_group = factor(risk_group, levels=c("eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                                      "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))) %>%
      group_by(studydrug2, risk_group) %>%
      summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_count",
                  values_from=count)
    
    # calculate median follow up time (years) per group
    followup <- temp2 %>%
      mutate(risk_group = factor(risk_group, levels=c("eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                                      "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))) %>%
      group_by(studydrug2,risk_group) %>%
      summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_followup",
                  values_from=time)
    
    # summarise number of events per group
    events <- temp2 %>%
      mutate(risk_group = factor(risk_group, levels=c("eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                                      "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                      "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))) %>%
      group_by(studydrug2,risk_group) %>%
      summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
                drug_count=round(n()/n.imp, 0)) %>%
      mutate(events_perc=round(event_count*100/drug_count, 1),
             events=paste0(event_count, " (", events_perc, "%)")) %>%
      select(studydrug2, risk_group, events) %>%
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_events",
                  values_from=events)
    
    
    # write formulas for adjusted and unadjusted analyses
    f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*risk_group"))
    
    f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*risk_group + ", paste(covariates, collapse=" + ")))
    
    # create empty vectors to store the hazard ratios from every imputed dataset
    # for the unadjusted survival models
    COEFS.presegfr.noalb.unadj <- SE.presegfr.noalb.unadj <-
      COEFS.presegfr.microalb.unadj <- SE.presegfr.microalb.unadj <-
      # COEFS.presegfr.macroalb.unadj <- SE.presegfr.macroalb.unadj <-      
      COEFS.redegfr.noalb.unadj <- SE.redegfr.noalb.unadj <-
      COEFS.redegfr.microalb.unadj <- SE.redegfr.microalb.unadj <-
      # COEFS.redegfr.macroalb.unadj <- SE.redegfr.macroalb.unadj <-
      # for the adjusted survival models
      COEFS.presegfr.noalb.adj <- SE.presegfr.noalb.adj <-
      COEFS.presegfr.microalb.adj <- SE.presegfr.microalb.adj <-
      # COEFS.presegfr.macroalb.adj <- SE.presegfr.macroalb.adj <-      
      COEFS.redegfr.noalb.adj <- SE.redegfr.noalb.adj <-
      COEFS.redegfr.microalb.adj <- SE.redegfr.microalb.adj <-
      # COEFS.redegfr.macroalb.adj <- SE.redegfr.macroalb.adj <-
      rep(NA,n.imp)
    
    for (i in 1:n.imp) {
      print(paste0("Analyses in imputed dataset number ", i))
      
      #unadjusted analyses first
      fit.unadj <- coxph(f2, temp2[temp2$.imp == i,])
      
      #store coefficients and standard errors from this model
      COEFS.presegfr.noalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"]
      SE.presegfr.noalb.unadj[i] <- sqrt(fit.unadj$var[1,1])
      
      COEFS.presegfr.microalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:risk_groupeGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol"]
      SE.presegfr.microalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var)-2,nrow(fit.unadj$var)-2]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)-2])
      
      COEFS.redegfr.noalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR <3mg/mmol"]
      SE.redegfr.noalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var)-1,nrow(fit.unadj$var)-1]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)-1])
      
      COEFS.redegfr.microalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR 3-30mg/mmol"]
      SE.redegfr.microalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var),nrow(fit.unadj$var)]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)])
      
      #adjusted analyses
      fit.adj <- coxph(f_adjusted2, temp2[temp2$.imp == i,])
      
      COEFS.presegfr.noalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"]
      SE.presegfr.noalb.adj[i] <- sqrt(fit.adj$var[1,1])
      
      COEFS.presegfr.microalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:risk_groupeGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol"]
      SE.presegfr.microalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var)-2,nrow(fit.adj$var)-2]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)-2])
      
      COEFS.redegfr.noalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR <3mg/mmol"]
      SE.redegfr.noalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var)-1,nrow(fit.adj$var)-1]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)-1])
      
      COEFS.redegfr.microalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:risk_groupeGFR <60mL/min/1.73m2, uACR 3-30mg/mmol"]
      SE.redegfr.microalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var),nrow(fit.adj$var)]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)])
      
    }
    
    # pool hazard ratios
    unadjusted_presegfr_noalb <- pool.rubin.HR(COEFS.presegfr.noalb.unadj, SE.presegfr.noalb.unadj, n.imp)
    adjusted_presegfr_noalb <- pool.rubin.HR(COEFS.presegfr.noalb.adj, SE.presegfr.noalb.adj, n.imp)
    
    unadjusted_presegfr_microalb <- pool.rubin.HR(COEFS.presegfr.microalb.unadj, SE.presegfr.microalb.unadj, n.imp)
    adjusted_presegfr_microalb <- pool.rubin.HR(COEFS.presegfr.microalb.adj, SE.presegfr.microalb.adj, n.imp)
    
    unadjusted_redegfr_noalb <- pool.rubin.HR(COEFS.redegfr.noalb.unadj, SE.redegfr.noalb.unadj, n.imp)
    adjusted_redegfr_noalb <- pool.rubin.HR(COEFS.redegfr.noalb.adj, SE.redegfr.noalb.adj, n.imp)
    
    unadjusted_redegfr_microalb <- pool.rubin.HR(COEFS.redegfr.microalb.unadj, SE.redegfr.microalb.unadj, n.imp)
    adjusted_redegfr_microalb <- pool.rubin.HR(COEFS.redegfr.microalb.adj, SE.redegfr.microalb.adj, n.imp)
    
    # save pooled HR and 95% confidence interval
    unadjusted_presegfr_noalb_string <- paste0(round(unadjusted_presegfr_noalb[1], 2), " (", round(unadjusted_presegfr_noalb[2], 2), ", ", round(unadjusted_presegfr_noalb[3], 2), ")")
    adjusted_presegfr_noalb_string <- paste0(round(adjusted_presegfr_noalb[1], 2), " (", round(adjusted_presegfr_noalb[2], 2), ", ", round(adjusted_presegfr_noalb[3], 2), ")")
    
    unadjusted_presegfr_microalb_string <- paste0(round(unadjusted_presegfr_microalb[1], 2), " (", round(unadjusted_presegfr_microalb[2], 2), ", ", round(unadjusted_presegfr_microalb[3], 2), ")")
    adjusted_presegfr_microalb_string <- paste0(round(adjusted_presegfr_microalb[1], 2), " (", round(adjusted_presegfr_microalb[2], 2), ", ", round(adjusted_presegfr_microalb[3], 2), ")")
    
    unadjusted_redegfr_noalb_string <- paste0(round(unadjusted_redegfr_noalb[1], 2), " (", round(unadjusted_redegfr_noalb[2], 2), ", ", round(unadjusted_redegfr_noalb[3], 2), ")")
    adjusted_redegfr_noalb_string <- paste0(round(adjusted_redegfr_noalb[1], 2), " (", round(adjusted_redegfr_noalb[2], 2), ", ", round(adjusted_redegfr_noalb[3], 2), ")")
    
    unadjusted_redegfr_microalb_string <- paste0(round(unadjusted_redegfr_microalb[1], 2), " (", round(unadjusted_redegfr_microalb[2], 2), ", ", round(unadjusted_redegfr_microalb[3], 2), ")")
    adjusted_redegfr_microalb_string <- paste0(round(adjusted_redegfr_microalb[1], 2), " (", round(adjusted_redegfr_microalb[2], 2), ", ", round(adjusted_redegfr_microalb[3], 2), ")")
    
    # combine in dataframe that we can tabulate
    presegfr_noalb_hr <- cbind(outcome=k, count[1,c(2:3)], followup[1,c(2:3)], events[1,c(2:3)],
                               unadjusted=unadjusted_presegfr_noalb_string, adjusted=adjusted_presegfr_noalb_string
    )
    presegfr_microalb_hr <- cbind(outcome=k, count[2,c(2:3)], followup[2,c(2:3)], events[2,c(2:3)],
                                  unadjusted=unadjusted_presegfr_microalb_string, adjusted=adjusted_presegfr_microalb_string
    )
    presegfr_macroalb_hr <- cbind(outcome=k, "DPP4i/SU_count"=NA, SGLT2i_count=NA, "DPP4i/SU_followup"=NA,
                                  SGLT2i_followup=NA, "DPP4i/SU_events"=NA, SGLT2i_events=NA,
                                  unadjusted=NA, adjusted=NA
    ) %>% data.frame(check.names=F)
    redegfr_noalb_hr <- cbind(outcome=k, count[3,c(2:3)], followup[3,c(2:3)], events[3,c(2:3)],
                              unadjusted=unadjusted_redegfr_noalb_string, adjusted=adjusted_redegfr_noalb_string
    )
    redegfr_microalb_hr <- cbind(outcome=k, count[4,c(2:3)], followup[4,c(2:3)], events[4,c(2:3)],
                                 unadjusted=unadjusted_redegfr_microalb_string, adjusted=adjusted_redegfr_microalb_string
    )
    redegfr_macroalb_hr <- cbind(outcome=k, "DPP4i/SU_count"=NA, SGLT2i_count=NA, "DPP4i/SU_followup"=NA,
                                 SGLT2i_followup=NA, "DPP4i/SU_events"=NA, SGLT2i_events=NA,
                                 unadjusted=NA, adjusted=NA
    ) %>% data.frame(check.names=F)
    
    outcome_subgroup_SGLT2ivsDPP4iSU_hrs <- rbind(presegfr_noalb_hr, presegfr_microalb_hr, presegfr_macroalb_hr,
                                                  redegfr_noalb_hr, redegfr_microalb_hr, redegfr_macroalb_hr)
    
    temp <- rbind(
      cbind(outcome = k, contrast = "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", analysis = "Adjusted",
            HR = adjusted_presegfr_noalb[1], LB = adjusted_presegfr_noalb[2], UB = adjusted_presegfr_noalb[3], string = adjusted_presegfr_noalb_string),
      cbind(outcome = k, contrast = "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol", analysis = "Adjusted",
            HR = adjusted_presegfr_microalb[1], LB = adjusted_presegfr_microalb[2], UB = adjusted_presegfr_microalb[3], string = adjusted_presegfr_microalb_string),
      cbind(outcome = k, contrast = "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", analysis = "Adjusted",
            HR = NA, LB = NA, UB = NA, string = NA),
      cbind(outcome = k, contrast = "eGFR <60mL/min/1.73m2, uACR <3mg/mmol", analysis = "Adjusted",
            HR = adjusted_redegfr_noalb[1], LB = adjusted_redegfr_noalb[2], UB = adjusted_redegfr_noalb[3], string = adjusted_redegfr_noalb_string),
      cbind(outcome = k, contrast = "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol", analysis = "Adjusted",
            HR = adjusted_redegfr_microalb[1], LB = adjusted_redegfr_microalb[2], UB = adjusted_redegfr_microalb[3], string = adjusted_redegfr_microalb_string),
      cbind(outcome = k, contrast = "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol", analysis = "Adjusted",
            HR = NA, LB = NA, UB = NA, string = NA)
    )
    
    subgroup_hrs <- rbind(subgroup_hrs, temp)
    subgroup_SGLT2ivsDPP4iSU_hrs <- rbind(subgroup_SGLT2ivsDPP4iSU_hrs, outcome_subgroup_SGLT2ivsDPP4iSU_hrs)
  }
}


## save HRs and display
flextable(subgroup_SGLT2ivsDPP4iSU_hrs)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
save(subgroup_hrs, file=paste0(today, "_subgroup_hrs.Rda"))
save(subgroup_SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_subgroup_SGLT2ivsDPP4iSU_hrs.Rda"))

############################5A FOREST PLOT FOR HRs PER DRUG CLASS (SUPPLEMENTAL FIGURE)################################################################

# create labels
labels <- data.frame(matrix("", nrow = 1, ncol = length(all_hrs)))
names(labels) <- names(all_hrs)
labels <- labels %>% mutate(contrast = "Overall (by analytic approach)", 
                  analysis = "Overall (by analytic approach)",
                  string = "Hazard Ratio (95% CI)",
                  SU_nN = "Events/subjects (SU)",
                  DPP4i_nN = "Events/subjects (DPP4i)",
                  SGLT2i_nN = "Events/subjects (SGLT2i)")

labels_plot <- all_hrs

for (k in unique(all_hrs$contrast)) {
  for (m in unique(all_hrs$outcome)) {
    labels_temp <- labels %>% data.frame()
    labels_temp$outcome <- m
    labels_temp$contrast <- k
    labels_plot <- rbind(labels_temp, labels_plot)
  }
}

labels_plot$analysis <- factor(labels_plot$analysis, levels = rev(unique(labels_plot$analysis)))



# plot
p_hr_1 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.3, 0.5, 0.75, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp",]$analysis)) + 1), 
                  xlim=c(0.3, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.25,
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp",]$analysis)) + 1, 
           label = "Favours SU") +
  labs(x="", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_1 <-
  labels_plot %>%
  mutate(analysis = ifelse(analysis == "Overall (by analytic approach)",
                           "SGLT2-inhibitors vs sulfonylureas",
                           as.character(analysis)),
         analysis = factor(analysis, levels = c(
           "IPTW",
           "Overlap-weighted",
           "Adjusted",           
           "SGLT2-inhibitors vs sulfonylureas"
         ))) %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = (analysis))) + 
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs SU",]$
                        analysis == "Overall (by analytic approach)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_1 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50" & labels_plot$contrast == "SGLT2i vs SU",]$
                        string == "Hazard Ratio (95% CI)", "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_counts_1 <- labels_plot %>% filter(outcome=="ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
            fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs SU",]$SGLT2i_nN == labels$SGLT2i_nN, "bold", "plain")) +
  geom_text(aes(x = 4, label = `SU_nN`), hjust = 1, fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs SU",]$SU_nN == labels$SU_nN, "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(-2, 5))

## same for SGLT2i vs DPP4

p_hr_2 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.3, 0.5, 0.75, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1), 
                  xlim=c(0.3, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.25,
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours DPP4i") +
  labs(x="", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 

p_left_2 <-
  labels_plot %>%
  mutate(analysis = ifelse(analysis == "Overall (by analytic approach)",
                           "SGLT2-inhibitors vs DPP4-inhibitors",
                           as.character(analysis)),
         analysis = factor(analysis, levels = c(
           "IPTW",
           "Overlap-weighted",
           "Adjusted",           
           "SGLT2-inhibitors vs DPP4-inhibitors"
         ))) %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = (analysis))) + 
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50" & labels_plot$contrast == "SGLT2i vs DPP4i",]$
                        analysis == "Overall (by analytic approach)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_2 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    colour = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50" & labels_plot$contrast == "SGLT2i vs DPP4i",]$
                      string == "Hazard Ratio (95% CI)", "white", "black")) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_counts_2 <- labels_plot %>% filter(outcome=="ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
            fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs DPP4i",]$SGLT2i_nN == labels$SGLT2i_nN, "bold", "plain")) +
  geom_text(aes(x = 4, label = `DPP4i_nN`), hjust = 1, fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs DPP4i",]$DPP4i_nN == labels$DPP4i_nN, "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(-2, 5))

## same for DPP4i vs SU

p_hr_3 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.3, 0.5, 0.75, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1), 
                  xlim=c(0.3, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours DPP4i") +
  annotate("text", x = 1.25,
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours SU") +
  labs(x="Hazard Ratio", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 

p_left_3 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  mutate(analysis = ifelse(analysis == "Overall (by analytic approach)",
                           "DPP4-inhibitors vs sulfonylureas",
                           as.character(analysis)),
         analysis = factor(analysis, levels = c(
           "IPTW",
           "Overlap-weighted",
           "Adjusted",           
           "DPP4-inhibitors vs sulfonylureas"
         ))) %>%
  ggplot(aes(y = (analysis))) + 
  # geom_text(aes(x = 0, label = analysis), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "DPP4i vs SU",]$
                        analysis == "Overall (by analytic approach)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_3 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    colour = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50" & labels_plot$contrast == "DPP4i vs SU",]$
                      string == "Hazard Ratio (95% CI)", "white", "black")) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_counts_3 <- labels_plot %>% filter(outcome=="ckd_egfr50_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(aes(x = 1, label = DPP4i_nN), hjust = 1, 
            fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "DPP4i vs SU",]$DPP4i_nN == labels$DPP4i_nN, "bold", "plain")) +
  geom_text(aes(x = 4, label = `SU_nN`), hjust = 1, fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "DPP4i vs SU",]$SU_nN == labels$SU_nN, "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(-2, 5))

layout <- c(
  area(t = 0, l = 0, b = 5, r = 4), 
  area(t = 0, l = 5, b = 5, r = 9), 
  area(t = 0, l = 10, b = 5, r = 14),
  area(t = 0, l = 15, b = 5, r = 19),
  area(t = 6, l = 0, b = 11, r = 4), 
  area(t = 6, l = 5, b = 11, r = 9),
  area(t = 6, l = 10, b = 11, r = 14),
  area(t = 6, l = 15, b = 11, r = 19),
  area(t = 12, l = 0, b = 17, r = 4), 
  area(t = 12, l = 5, b = 17, r = 9),
  area(t = 12, l = 10, b = 17, r = 14),
  area(t = 12, l = 15, b = 17, r = 19)
)

p_left_1 + p_counts_1 + p_hr_1 + p_right_1 + 
  p_left_2 + p_counts_2 + p_hr_2 + p_right_2 + 
  p_left_3 + p_counts_3 + p_hr_3 + p_right_3 + plot_layout(design = layout)

############################5B FOREST PLOT FOR HRs OVERALL (FIGURE 1)################################################################


# prep all_hrs dataframe for forest plot to show hazard ratios and add literature-reported HR
trial_hr <- data.frame(matrix("", nrow=1, ncol=length(all_SGLT2ivsDPP4iSU_hrs))) 
names(trial_hr) <- names(all_SGLT2ivsDPP4iSU_hrs)
trial_hr <- trial_hr %>% mutate(
  outcome = "ckd_egfr50", contrast = "SGLT2i vs SU", analysis = "Nuffield Group, 2022", 
                  HR = 0.62, LB = 0.56, UB = 0.68, string = paste0(HR, " (", LB, ", ", UB, ")"), 
                                                                               model = paste0(string, " [", analysis, "]"))

labels1 <- data.frame(matrix("", nrow = 1, ncol = length(trial_hr)))
names(labels1) <- names(trial_hr)
labels1 <- labels1 %>% mutate(contrast = "Meta-analysis of RCTs", 
                            analysis = "Meta-analysis of RCTs",
                            string = "Hazard Ratio (95% CI)")
labels_plot <- trial_hr

for (k in unique(trial_hr$contrast)) {
  for (m in unique(trial_hr$outcome)) {
    labels_temp <- labels1
    labels_temp$outcome <- m
    labels_temp$contrast <- k
    labels_plot <- rbind(labels_temp, labels_plot)
  }
}

# create labels for forest plot with studydrug2 (SGLT2i vs DPP4i/SU combined)
labels_plot2 <- all_SGLT2ivsDPP4iSU_hrs

labels2 <- data.frame(matrix("", nrow = 1, ncol = length(all_SGLT2ivsDPP4iSU_hrs)))
names(labels2) <- names(all_SGLT2ivsDPP4iSU_hrs)
labels2 <- labels2 %>% mutate(contrast = "CPRD (by analytic approach)", 
                              analysis = "CPRD (by analytic approach)",
                              string = "Hazard Ratio (95% CI)",
                              )



for (k in unique(all_SGLT2ivsDPP4iSU_hrs$contrast)) {
  for (m in unique(all_SGLT2ivsDPP4iSU_hrs$outcome)) {
    labels_temp <- labels2 
    labels_temp$outcome <- m
    labels_temp$contrast <- k
    labels_plot2 <- rbind(labels_temp, labels_plot2)
  }
}

labels_plot2$analysis <- factor(labels_plot2$analysis, levels = rev(unique(labels_plot2$analysis)))


class(trial_hr$HR) <- class(trial_hr$LB) <- class(trial_hr$UB) <-
  class(all_SGLT2ivsDPP4iSU_hrs$HR) <- class(all_SGLT2ivsDPP4iSU_hrs$LB) <- class(all_SGLT2ivsDPP4iSU_hrs$UB) <- "numeric"

# plot
p_hr_all <- 
  all_SGLT2ivsDPP4iSU_hrs %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.45, 0.6, 0.80, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$outcome == "ckd_egfr50",]$analysis)) + 1), 
                  xlim=c(0.45, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  labs(x="Hazard Ratio", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5)) 


p_left_all <-
  labels_plot2 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = (analysis))) + 
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot2[labels_plot2$outcome == "ckd_egfr50",]$
                        analysis == labels2$analysis, "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_all <-
  labels_plot2 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    fontface = "plain",
    colour = ifelse(labels_plot2[labels_plot2$outcome == "ckd_egfr50",]$string == labels2$string, "white", "black")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

# trial hr
p_hr_trial <- 
  trial_hr %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.45, 0.60, 0.80, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(trial_hr[trial_hr$outcome == "ckd_egfr50",]$analysis)) + 1), 
                  xlim=c(0.45, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3, colour = "#0072B2") +
  geom_linerange(aes(xmin=LB, xmax=UB), colour = "#0072B2") +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(trial_hr[trial_hr$outcome == "ckd_egfr50",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.25,
           y = length(unique(trial_hr[trial_hr$outcome == "ckd_egfr50",]$analysis)) + 1, 
           label = "Favours\nDPP4i/SU") +
  labs(x="", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


p_left_trial <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = rev(analysis))) + 
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50",]$
                        analysis == labels1$analysis, "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_trial <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50",]$
                        string == labels1$string, "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

# layout for plots below
layout <- c(
  area(t = 21, l = 0, b = 55, r = 4), 
  area(t = 21, l = 5, b = 55, r = 9),
  area(t = 21, l = 10, b = 55, r = 14),
  area(t = 0, l = 0, b = 23, r = 4), 
  area(t = 0, l = 5, b = 23, r = 9), 
  area(t = 0, l = 10, b = 23, r = 14)
  
)

# Final plot arrangement
p_left_all + p_hr_all + p_right_all + 
  p_left_trial + p_hr_trial + p_right_trial + plot_layout(design = layout)


############################5C FOREST PLOT FOR HRs BY SUBGROUP (SUPPLEMENTAL FIGURE)################################################################

# load("2024-07-13_subgroup_hrs.Rda")
# load("2024-07-13_subgroup_SGLT2ivsDPP4iSU_hrs.Rda")

subgroup_hrs <- subgroup_hrs %>% cbind(subgroup_SGLT2ivsDPP4iSU_hrs %>% select(-c(outcome, adjusted, unadjusted)))

subgroup_hrs <- subgroup_hrs %>%
  separate(`DPP4i/SU_events`, into = c("DPP4i/SU_events_number", "DPP4i/SU_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(SGLT2i_events, into = c("SGLT2i_events_number", "SGLT2i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `DPP4i/SU_events_percentage` = str_replace(`DPP4i/SU_events_percentage`, "\\)", ""),
    SGLT2i_events_percentage = str_replace(SGLT2i_events_percentage, "\\)", ""),
    `DPP4i/SU_nN` = paste0(`DPP4i/SU_events_number`, "/", `DPP4i/SU_count`),
    SGLT2i_nN = paste0(SGLT2i_events_number, "/", SGLT2i_count)
  )

# prep data frames with row for overall
overall <- all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$analysis == "Adjusted",]
overall$contrast <- "Entire cohort"

labels_plot5 <- overall

labels5 <- data.frame(matrix("", nrow = 1, ncol = length(overall)))
names(labels5) <- names(overall)
labels5 <- labels5 %>% mutate(contrast = "Overall", 
                 analysis = "Overall", 
                 string = "Hazard Ratio (95% CI)",
                 SGLT2i_nN = "Events/subjects (SGLT2i)",
                 `DPP4i/SU_nN` = "(DPP4i/SU)")

for (m in unique(overall$outcome)) {
  labels_temp <- labels5 
  labels_temp$outcome <- m
  labels_plot5 <- rbind(labels_temp, labels_plot5)
}

# have to coerce HR and CI to class numeric as they sometimes default to character
class(subgroup_hrs$HR) <- class(subgroup_hrs$LB) <- class(subgroup_hrs$UB) <- "numeric"

labels3 <- data.frame(matrix("", nrow = 1, ncol = length(subgroup_hrs)))
names(labels3) <- names(subgroup_hrs)
labels3 <- labels3 %>% mutate(contrast = paste0("By subgroup (p = ", round(p_value_interaction,2), " for trend)"), 
                  analysis = "By subgroup", 
                  string = "Hazard Ratio (95% CI)",
                  SGLT2i_nN = "Events/subjects",
                  `DPP4i/SU_nN` = "")

labels_plot3 <- subgroup_hrs

for (k in unique(subgroup_hrs$outcome)) {
  labels_temp <- labels3
  labels_temp$outcome <- k
  labels_plot3 <- rbind(labels_temp, labels_plot3)
}

labels_plot3$contrast <- factor(labels_plot3$contrast, levels = c(labels3$contrast,
                                                                  "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                                                  "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                                  "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                                                  "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                                                  "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                                                  "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))

# plot by risk group
p_counts_subgroup <- labels_plot3 %>% filter(outcome=="ckd_egfr50") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
            colour = ifelse(labels_plot3[labels_plot3$outcome == "ckd_egfr50",]$SGLT2i_nN == labels3$SGLT2i_nN, "white", "black")) +
  geom_text(aes(x = 3, label = `DPP4i/SU_nN`), hjust = 1, fontface = "plain") +
  theme_void() +
  coord_cartesian(xlim = c(-2, 5))

p_hr_subgroup <- 
  subgroup_hrs %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0)) +
  coord_cartesian(ylim=c(1,length(unique(subgroup_hrs[subgroup_hrs$outcome == "ckd_egfr50",]$contrast)) + 1), 
                  xlim=c(0.25, 2)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(subgroup_hrs[subgroup_hrs$outcome == "ckd_egfr50",]$contrast)) + 1, 
           label = "") +
  annotate("text", x = 1.25,
           y = length(unique(subgroup_hrs[subgroup_hrs$outcome == "ckd_egfr50",]$contrast)) + 1, 
           label = "") +
  labs(x="Hazard ratio", y="") + 
  geom_point(aes(x = ifelse(UB > 2, 2.2, NA),
                 y=1.01), 
             shape = 62, size = 5, color = "black") + # add symbol to indicate CI extends beyond axis limits
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_subgroup <-
  labels_plot3 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = rev(unique((contrast))))) + 
  geom_text(
    aes(x = 1, label = contrast),
    hjust = 0,
    fontface = ifelse(labels_plot3[labels_plot3$outcome == "ckd_egfr50",]$
                        contrast == labels3$contrast, "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_subgroup <-
  labels_plot3 %>%
  filter(outcome == "ckd_egfr50") %>%
  mutate(string = 
           ifelse(coalesce(string == lead(string), FALSE), paste0(string, " "), string)) %>% # if hazard ratios are the same, add a space to the end of one of them so they do not get taken as one
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    colour = ifelse(labels_plot3[labels_plot3$outcome == "ckd_egfr50",]$string == "Hazard Ratio (95% CI)", "white", "black")) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))


# plot for overall HR
p_counts_overall <- labels_plot5 %>% filter(outcome=="ckd_egfr50") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
            fontface = ifelse(labels_plot5[labels_plot5$outcome == "ckd_egfr50",]$SGLT2i_nN == labels5$SGLT2i_nN, "bold", "plain")) +
  geom_text(aes(x = 3, label = `DPP4i/SU_nN`), hjust = 1, 
            fontface = ifelse(labels_plot5[labels_plot5$outcome == "ckd_egfr50",]$`DPP4i/SU_nN` == labels5$`DPP4i/SU_nN`, "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(-2, 5))

p_hr_overall <- 
  overall %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.25, 0.5, 0.75, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast)) + 1), 
                  xlim=c(0.25, 2)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3, colour = "#D55E00") +
  geom_linerange(aes(xmin=LB, xmax=UB), colour = "#D55E00") +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.25,
           y = length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast)) + 1, 
           label = "Favours\nDPP4i/SU") +
  labs(x="Hazard Ratio", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_overall <-
  labels_plot5 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = (unique((contrast))))) + 
  geom_text(
    aes(x = 1, label = contrast),
    hjust = 0,
    fontface = ifelse(labels_plot5[labels_plot5$outcome == "ckd_egfr50",]$
                        contrast == labels5$contrast, "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_overall <-
  labels_plot5 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    fontface = ifelse(labels_plot5[labels_plot5$outcome == "ckd_egfr50",]$string == "Hazard Ratio (95% CI)", "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))


# layout for plots below
layout <- c(
  area(t = 11, l = 7, b = 55, r = 13),
  area(t = 11, l = 0, b = 55, r = 7), 
  area(t = 11, l = 12, b = 55, r = 18),
  area(t = 11, l = 19, b = 55, r = 24),
  area(t = 0, l = 7, b = 13, r = 13), 
  area(t = 0, l = 0, b = 13, r = 7),
  area(t = 0, l = 12, b = 13, r = 18),
  area(t = 0, l = 19, b = 13, r = 24)
  
)

# Final plot arrangement
p_counts_subgroup + p_left_subgroup + p_hr_subgroup + p_right_subgroup + 
  p_counts_overall + p_left_overall + p_hr_overall + p_right_overall + plot_layout(design = layout)


############################5D FOREST PLOT FOR HRs OF SECONDARY OUTCOMES (SUPPLEMENTAL FIGURE)################################################################

secondary <- all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$analysis == "Adjusted",]
secondary <- secondary %>% filter(!outcome %in% c("death")) # will not display this in current figure

labels_plot6 <- secondary

labels6 <- data.frame(matrix("", nrow = 1, ncol = length(secondary)))
names(labels6) <- names(secondary)
labels6 <- labels6 %>% mutate(analysis = "Overall", 
                              string = "Hazard Ratio (95% CI)",
                              SGLT2i_nN = "Events/subjects (SGLT2i)",
                              `DPP4i/SU_nN` = "(DPP4i/SU)")

for (m in unique(secondary$outcome)) {
  labels_temp <- labels6 
  labels_temp$outcome <- m
  labels_plot6 <- rbind(labels_temp, labels_plot6)
}

labels_plot6 <- labels_plot6 %>% mutate(
  contrast = ifelse(
    outcome == "ckd_egfr40", " 40% eGFR decline / ESKD", 
    ifelse(outcome == "ckd_egfr50", " 50% eGFR decline / ESKD",
           ifelse(outcome == "death", " All-cause mortality", ifelse(
             outcome == "macroalb", " Progression to significant albuminuria (≥30mg/mmol)", ifelse(
               outcome == "dka", " Diabetic keto-acidosis", ifelse(
                 outcome == "amputation", " Amputation", ifelse(
                   outcome == "side_effect", " Mycotic genital infection", NA 
                 )
               )
             )
           )
         )
       )
  ),
  contrast = ifelse(
    analysis == "Overall", "Outcome", as.character(contrast)
  )
)

# have to coerce HR and CI to class numeric as they sometimes default to character
class(secondary$HR) <- class(secondary$LB) <- class(secondary$UB) <- "numeric"

plot_expression <- ""

for (m in rev(unique(secondary$outcome))) {

    p_counts <- labels_plot6 %>% filter(outcome==m) %>%
    ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
    geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
              colour = ifelse(!m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$SGLT2i_nN == labels6$SGLT2i_nN,
                              "white", "black"),
              fontface = ifelse(m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$SGLT2i_nN == labels6$SGLT2i_nN,
                                "bold", "plain")) +
    geom_text(aes(x = 3, label = `DPP4i/SU_nN`), hjust = 1, 
              colour = ifelse(!m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$`DPP4i/SU_nN` == labels6$`DPP4i/SU_nN`,
                              "white", "black"),
              fontface = ifelse(m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$`DPP4i/SU_nN` == labels6$`DPP4i/SU_nN`,
                                "bold", "plain")) +
    theme_void() +
    coord_cartesian(xlim = c(-2, 5))
  
  p_hr <- 
    secondary %>%
    filter(outcome == m) %>%
    ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
    scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.5, 2.25, 3.25)) +
    coord_cartesian(ylim=c(1,length(unique(labels_plot6[labels_plot6$outcome == m,]$contrast)) + 1), 
                    xlim=c(0.5, 3.25)) +
    theme_classic() +
    geom_point(aes(x=HR), shape=15, size=3) +
    geom_linerange(aes(xmin=LB, xmax=UB)) +
    geom_vline(xintercept = 1, linetype="dashed") +
    annotate("text", x = .65, 
             y = length(unique(secondary[secondary$outcome == m,]$contrast)) + 2, 
             label = ifelse(m==unique(secondary$outcome)[1], "Favours\nSGLT2i", "")) +
    annotate("text", x = 1.5,
             y = length(unique(secondary[secondary$outcome == m,]$contrast)) + 2, 
             label = ifelse(m==unique(secondary$outcome)[1], "Favours\nDPP4i/SU", "")) +
    labs(x=ifelse(m==unique(secondary$outcome)[nlevels(as.factor(secondary$outcome))], "Hazard ratio", ""), y="") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y= element_blank(),
          axis.title.y= element_blank(),
          axis.line.x = if (!m==unique(secondary$outcome)[nlevels(as.factor(secondary$outcome))]) {element_blank()},
          axis.text.x = if (!m==unique(secondary$outcome)[nlevels(as.factor(secondary$outcome))]) {element_blank()},
          axis.ticks.x = if (!m==unique(secondary$outcome)[nlevels(as.factor(secondary$outcome))]) {element_blank()},
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) 
  
  p_left <-
    labels_plot6 %>%
    filter(outcome == m) %>%
    ggplot(aes(y = (unique((contrast))))) + 
    geom_text(
      aes(x = 1, label = contrast),
      hjust = 0,
      fontface = ifelse(labels_plot6[labels_plot6$outcome == m,]$
                          contrast == "Outcome", "bold", "plain"),
      colour = ifelse(m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$
                        contrast == "Outcome" | !labels_plot6[labels_plot6$outcome == m,]$
                        contrast == "Outcome", "black", "white")
    ) +
    theme_void() +
    coord_cartesian(xlim = c(0, 4))
  
  p_right <-
    labels_plot6 %>%
    filter(outcome == m) %>%
    ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
    geom_text(
      aes(x = 0, label = string),
      hjust = 0,
      fontface = ifelse(labels_plot6[labels_plot6$outcome == m,]$string == "Hazard Ratio (95% CI)", "bold", "plain"),
      colour = ifelse(labels_plot6[labels_plot6$outcome == m,]$string == "Hazard Ratio (95% CI)" & !m==unique(secondary$outcome)[1], 
                      "white", "black")) +
    theme_void() +
    coord_cartesian(xlim = c(0, 4))
  
  assign(paste0("p_counts_", m), p_counts)
  assign(paste0("p_hr_", m), p_hr)
  assign(paste0("p_left_", m), p_left)
  assign(paste0("p_right_", m), p_right)
  
  plot_expression <- paste0(plot_expression, "p_counts_", m, " + p_left_", m, " + p_hr_", m, " + p_right_", m, " + ")
}

n.plots <- nlevels(as.factor(secondary$outcome))

# layout for plots below

i <- 1
layout <- paste("area(t = ",(i-1)*10, ", l = 7, b = ",(i-1)*10+22,", r = 13), area(t = ",(i-1)*10, ", l = 0, b = ",(i-1)*10+22,", r = 7), area(t = ",(i-1)*10, ", l = 12, b = ",(i-1)*10+22,", r = 18), area(t = ",(i-1)*10, ", l = 19, b = ", (i-1)*10+22,", r = 24)")


for (i in 2:n.plots) {
  layout <- paste("area(t = ",(i-1)*10, ", l = 7, b = ",(i-1)*10+22,", r = 13), area(t = ",(i-1)*10, ", l = 0, b = ",(i-1)*10+22,", r = 7), area(t = ",(i-1)*10, ", l = 12, b = ",(i-1)*10+22,", r = 18), area(t = ",(i-1)*10, ", l = 19, b = ", (i-1)*10+22,", r = 24), ", layout)
}

layout <- paste0("c(", layout, ")")

layout <- eval(str2lang(layout))

# Final plot arrangement

plot_expression <- paste0(plot_expression, "plot_layout(design = layout)")

eval(str2lang(plot_expression))

