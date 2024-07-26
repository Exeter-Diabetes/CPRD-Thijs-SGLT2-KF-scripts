## In this script our aim is to predict baseline risk using the CKD-PC risk scores.
## we will be looking at the risk score for 50% decline in eGFR/ESKD (50egfr).
## which we will recalibrate in the CPRD cohort for endpoint 50% decline in eGFR/ESKD.

## Contents
# 0 setup
# 1 recalibration of risk score in eGFR ≥60
# 2 recalibration of risk score in eGFR <60
# 3 store dataset with recalibrated risk score

############################0 SETUP################################################################

# 0 Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(rms)
library(pROC)
library(riskRegression)

options(dplyr.dpp4isummarise.inform = FALSE)

# set seed
set.seed(123)

# number of imputations 
n.imp <- 10

# number of bootstraps for bootstrap validation
n.bootstrap <- 500

# number of quantiles (for risk scores later on)
n.quantiles <- 10

# today's date
#today <- as.character(Sys.Date(), format="%Y%m%d")
today <- "2024-07-13"
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


# load data
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-07-13_t2d_ckdpc_imputed_data_withweights.Rda")

# select imputed data only (ie. remove non-imputed data)
cohort <- cohort[cohort$.imp > 0,]

## remove double overlapping entries for DPP4i and SU that overlap (take one only)
cohort <- cohort %>% group_by(.imp, patid) %>% filter(
  !duplicated(studydrug2)
) %>% ungroup()

# check number of subjects
table(cohort$studydrug)
# SU  DPP4i SGLT2i 
# 389470 624700 559760   # 10 imputations therefore number of subjects per group appears 10 times larger

# set default colour-blind accessible colours for figures later on
cols <- c("SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00")
#in further analyses, the dpp4/su group will be combined, and we will use the dpp4 colour for this (strongest contrast)
cols <- c(cols, "DPP4i/SU" = "#0072B2")
cols <- cols[names(cols) %in% cohort$studydrug2]
cols <- cols[order(names(cols))]

############################1A RISK SCORE eGFR ≥60 ################################################################

# there are in effect 2 risk scores - one for people with an eGFR below 60mL/min/1.73m2, 
# and one for people with an eGFR above that.

# recalibrate risk score in those with preserved eGFR:

cohort_redegfr <- cohort %>% filter(preegfr < 60)
cohort_presegfr <- cohort %>% filter(preegfr >=60)



# make variable for risk deciles
cohort_presegfr$risk_decile <- ntile(cohort_presegfr$ckdpc_50egfr_score, n.quantiles)

## Get mean predicted probabilities by risk decile and studydrug
predicted <- cohort_presegfr %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score)/100)

# get mean predicted probabilities by risk decile (not by studydrug)
predicted_all <- cohort_presegfr %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score)/100)

## Find actual observed probabilities by risk score category and studydrug

EST.dpp4isu <- SE.dpp4isu <-
  EST.sglt2i <- SE.sglt2i <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_dpp4isu <- tibble() %>% mutate(
  observed_dpp4isu=NA,
  lower_ci_dpp4isu=NA,
  upper_ci_dpp4isu=NA,
  strata=NA
)

observed_sglt2i <- tibble() %>% mutate(
  observed_sglt2i=NA,
  lower_ci_sglt2i=NA,
  upper_ci_sglt2i=NA,
  strata=NA
)

observed_all <- tibble() %>% mutate(
  observed=NA,
  lower_ci=NA,
  upper_ci=NA,
  strata=NA
)

for (k in 1:n.quantiles) {
  for (i in 1:n.imp) {
    
    observed_dpp4isu_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                      data=cohort_presegfr[cohort_presegfr$.imp == i & 
                                                             cohort_presegfr$risk_decile == k &
                                                             cohort_presegfr$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                     data=cohort_presegfr[cohort_presegfr$.imp == i & 
                                                            cohort_presegfr$risk_decile == k &
                                                            cohort_presegfr$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                  data=cohort_presegfr[cohort_presegfr$.imp == i & 
                                                         cohort_presegfr$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd50$estimate
    SE.all[k,i] <- observed_all_ckd50$std.error
  }
  
  est.dpp4isu <- pool.rubin.KM(EST.dpp4isu[k,], SE.dpp4isu[k,], n.imp)
  observed_dpp4isu[k,] <- observed_dpp4isu[k,] %>% 
    mutate(
      observed_dpp4isu=est.dpp4isu[1],
      lower_ci_dpp4isu=est.dpp4isu[2],
      upper_ci_dpp4isu=est.dpp4isu[3],
      strata=k
    )
  
  est.sglt2i <- pool.rubin.KM(EST.sglt2i[k,], SE.sglt2i[k,], n.imp)
  observed_sglt2i[k,] <- observed_sglt2i[k,] %>% 
    mutate(
      observed_sglt2i=est.sglt2i[1],
      lower_ci_sglt2i=est.sglt2i[2],
      upper_ci_sglt2i=est.sglt2i[3],
      strata=k
    )
  
  est.all <- pool.rubin.KM(EST.all[k,], SE.all[k,], n.imp)
  observed_all[k,] <- observed_all[k,] %>% 
    mutate(
      observed=est.all[1],
      lower_ci=est.all[2],
      upper_ci=est.all[3],
      strata=k
    )
  
}


dpp4isu_events <- cohort_presegfr %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort_presegfr %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLTi=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug2=="DPP4i/SU")), observed_dpp4isu),
  cbind((predicted %>% filter(studydrug2=="SGLT2i")), observed_sglt2i)
) %>%
  mutate(observed=coalesce(observed_dpp4isu, observed_sglt2i),
         lower_ci=coalesce(lower_ci_dpp4isu, lower_ci_sglt2i),
         upper_ci=coalesce(upper_ci_dpp4isu, upper_ci_sglt2i))

events_table <- data.frame(t(dpp4isu_events %>% 
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd50_pred=NA, risk_decile=0)

## FINAL PLOT
p_presegfr_uncal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Raw CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = c(0,2.5,5,7.5,10,12.5))+
  scale_y_continuous(limits=c(-1,100)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Raw risk score, by risk decile") +
  coord_cartesian(xlim = c(0,13), ylim = c(0,12))


p_presegfr_uncal_bydeciles_dpp4isu

## C-stat
cohort_presegfr <- cohort_presegfr %>%
  mutate(ckdpc_50egfr_survival=(100-ckdpc_50egfr_score)/100)

raw_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival, data=cohort_presegfr, method="breslow")
cstat_presegfr_est <- summary(raw_mod)$concordance[1]
cstat_presegfr_est_ll <- summary(raw_mod)$concordance[1]-(1.96*summary(raw_mod)$concordance[2])
cstat_presegfr_est_ul <- summary(raw_mod)$concordance[1]+(1.96*summary(raw_mod)$concordance[2])
paste0("C statistic: ", round(cstat_presegfr_est, 4), ", 95% CI ", round(cstat_presegfr_est_ll, 4), "-", round(cstat_presegfr_est_ul,4))
# C statistic: 0.7262, 95% CI 0.7232-0.7291

## AUC
ROC_presegfr_raw <- roc(cohort_presegfr, ckd_egfr50_censvar, ckdpc_50egfr_score)
auc(ROC_presegfr_raw)
ci.auc(ROC_presegfr_raw)


## brier score for raw risk score
brier_presegfr_raw <- rep(NA, n.imp)
brier_presegfr_raw_se <- rep(NA, n.imp)


for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort_presegfr %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival)
  temp <- temp %>% # error if times = 3, therefore adding extra row with time beyond t=3
    rbind(    # adds rows below your dataset
      temp %>%
        slice(1) %>% # this selects the first patients in your dataset
        mutate(ckd_egfr50_censtime_yrs = 3.5)   # changes censored time to 3.5
    )
  raw_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival, 
                   data=temp, x=T)
  
  
  score_raw <- 
    Score(object = list(raw_mod), # need to pass cox model to Score() as a list in order for it to be processed
          formula = Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ 1, # null model
          data = temp,
          summary = "ibs", # statistic of interest is integrated brier score
          times = 3,  
          splitMethod = "bootcv",  # use bootstrapping for confidence intervals
          B = n.bootstrap,
          verbose = T) 
  
  brier_presegfr_raw[i] <- score_raw$Brier$score$Brier[2]
  brier_presegfr_raw_se[i] <- score_raw$Brier$score$se[2]
  
  rm(temp)
}

brier_presegfr_raw_se_pooled <- sqrt(mean(brier_presegfr_raw_se^2) + (1+1/n.imp)*var(brier_presegfr_raw))

#pool and print brier score
print(paste0("Brier score for raw risk score ", mean(brier_presegfr_raw), ", 95% CI ", mean(brier_presegfr_raw)-1.96*brier_presegfr_raw_se_pooled, "-", mean(brier_presegfr_raw)+1.96*brier_presegfr_raw_se_pooled))


# this risk score overestimates risk and needs to be recalibrated.

############################1B RECALIBRATION - BASELINE HAZARD UPDATE ################################################################

# initially we will assess calibration if we only update the baseline hazard
# this is done by fitting a cox proportional model with the linear predictor as the only variable as an offset.
# as this would only involve adjusting the overall event probability for this cohort_presegfr, no resampling (internal validation) is needed
recal_mod <- cph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ stats::offset(ckdpc_50egfr_lin_predictor), 
                 data = cohort_presegfr[!cohort_presegfr$studydrug=="SGLT2i",], x = TRUE, y = TRUE, surv = TRUE)

x <- summary(survfit(recal_mod),time=3)
bh_update_presegfr <- x$surv
bh_update_presegfr_se <- x$std.err
# print updated baseline hazard
print(paste0("Baseline hazard ", bh_update_presegfr, ", 95% CI ", bh_update_presegfr-1.96*bh_update_presegfr_se, "-", bh_update_presegfr+1.96*bh_update_presegfr_se))

cohort_presegfr <- cohort_presegfr %>% mutate(
  ckdpc_50egfr_survival_cal_bh=bh_update_presegfr^exp(ckdpc_50egfr_lin_predictor-mean(ckdpc_50egfr_lin_predictor)), # baseline hazard ^ e ^ centred linear predictor
  ckdpc_50egfr_score_cal_bh=(1-ckdpc_50egfr_survival_cal_bh)*100
)

## Show observed (estimated frequency binned per risk decile) vs predicted (risk-score predicted)
cohort_presegfr$risk_decile <- ntile(cohort_presegfr$ckdpc_50egfr_score_cal_bh, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- cohort_presegfr %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal_bh)/100)

predicted_all <- cohort_presegfr %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal_bh)/100)

### Find actual observed probabilities by risk score category and studydrug

EST.dpp4isu <- SE.dpp4isu <-
  EST.sglt2i <- SE.sglt2i <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_dpp4isu <- tibble() %>% mutate(
  observed_dpp4isu=NA,
  lower_ci_dpp4isu=NA,
  upper_ci_dpp4isu=NA,
  strata=NA
)

observed_sglt2i <- tibble() %>% mutate(
  observed_sglt2i=NA,
  lower_ci_sglt2i=NA,
  upper_ci_sglt2i=NA,
  strata=NA
)

observed_all <- tibble() %>% mutate(
  observed=NA,
  lower_ci=NA,
  upper_ci=NA,
  strata=NA
)

for (k in 1:n.quantiles) {
  for (i in 1:n.imp) {
    
    observed_dpp4isu_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                      data=cohort_presegfr[cohort_presegfr$.imp == i & 
                                                             cohort_presegfr$risk_decile == k &
                                                             cohort_presegfr$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                     data=cohort_presegfr[cohort_presegfr$.imp == i & 
                                                            cohort_presegfr$risk_decile == k &
                                                            cohort_presegfr$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                  data=cohort_presegfr[cohort_presegfr$.imp == i & 
                                                         cohort_presegfr$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd50$estimate
    SE.all[k,i] <- observed_all_ckd50$std.error
    
  }
  
  est.dpp4isu <- pool.rubin.KM(EST.dpp4isu[k,], SE.dpp4isu[k,], n.imp)
  observed_dpp4isu[k,] <- observed_dpp4isu[k,] %>% 
    mutate(
      observed_dpp4isu=est.dpp4isu[1],
      lower_ci_dpp4isu=est.dpp4isu[2],
      upper_ci_dpp4isu=est.dpp4isu[3],
      strata=k
    )
  
  est.sglt2i <- pool.rubin.KM(EST.sglt2i[k,], SE.sglt2i[k,], n.imp)
  observed_sglt2i[k,] <- observed_sglt2i[k,] %>% 
    mutate(
      observed_sglt2i=est.sglt2i[1],
      lower_ci_sglt2i=est.sglt2i[2],
      upper_ci_sglt2i=est.sglt2i[3],
      strata=k
    )
  
  est.all <- pool.rubin.KM(EST.all[k,], SE.all[k,], n.imp)
  observed_all[k,] <- observed_all[k,] %>% 
    mutate(
      observed=est.all[1],
      lower_ci=est.all[2],
      upper_ci=est.all[3],
      strata=k
    )
  
}


dpp4isu_events <- cohort_presegfr %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort_presegfr %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLTi=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug2=="DPP4i/SU")), observed_dpp4isu),
  cbind((predicted %>% filter(studydrug2=="SGLT2i")), observed_sglt2i)
) %>%
  mutate(observed=coalesce(observed_dpp4isu, observed_sglt2i),
         lower_ci=coalesce(lower_ci_dpp4isu,lower_ci_sglt2i),
         upper_ci=coalesce(upper_ci_dpp4isu, upper_ci_sglt2i))

events_table <- data.frame(t(dpp4isu_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd50_pred=NA, risk_decile=0)

## FINAL PLOT
p_presegfr_interim_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = c(0,2.5,5,7.5,10,12.5))+
  scale_y_continuous(limits=c(-1,100)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Recalibrated risk score (baseline hazard updated), by risk decile") +
  coord_cartesian(xlim = c(0,13), ylim = c(0,12))


p_presegfr_interim_bydeciles_dpp4isu

## C-stat

surv_mod_ckd50 <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival_cal_bh, data=cohort_presegfr, method="breslow")
cstat_presegfr_est <- summary(surv_mod_ckd50)$concordance[1]
cstat_presegfr_est_ll <- summary(surv_mod_ckd50)$concordance[1]-(1.96*summary(surv_mod_ckd50)$concordance[2])
cstat_presegfr_est_ul <- summary(surv_mod_ckd50)$concordance[1]+(1.96*summary(surv_mod_ckd50)$concordance[2])
paste0("C statistic: ", round(cstat_presegfr_est, 4), ", 95% CI ", round(cstat_presegfr_est_ll, 4), "-", round(cstat_presegfr_est_ul,4))
# C statistic: 0.7262, 95% CI 0.7232-0.7291

## AUC
ROC_presegfr_cal <- roc(cohort_presegfr, ckd_egfr50_censvar, ckdpc_50egfr_score_cal_bh)
auc(ROC_presegfr_cal)
ci.auc(ROC_presegfr_cal)

##brier score after baseline hazard updated
brier_presegfr_recal_bh <- rep(NA, n.imp)
brier_presegfr_recal_bh_se <- rep(NA, n.imp)
for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort_presegfr %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival_cal_bh)
  temp <- temp %>%
    rbind(    # adds rows below your dataset
      temp %>%
        slice(1) %>% # this selects the first patients in your dataset
        mutate(ckd_egfr50_censtime_yrs = 3.5)   # changes censored time to 3.5
    )
  recal_bh_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival_cal_bh, 
                        data=temp, x=T)
  
  score_recal_bh <- 
    Score(object = list(recal_bh_mod), # need to pass cox model to Score() as a list in order for it to be processed
          formula = Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ 1, # null model
          data = temp,
          summary = "ibs", # statistic of interest is integrated brier score
          times = 3,  
          splitMethod = "bootcv",  # use bootstrapping for confidence intervals
          B = n.bootstrap,
          verbose = T) 
  
  brier_presegfr_recal_bh[i] <- score_recal_bh$Brier$score$Brier[2]
  brier_presegfr_recal_bh_se[i] <- score_recal_bh$Brier$score$se[2]
  
  rm(temp)
  
}

brier_presegfr_recal_bh_se_pooled <- sqrt(mean(brier_presegfr_recal_bh_se^2) + (1+1/n.imp)*var(brier_presegfr_recal_bh))

#pool and print brier score
print(paste0("Brier score for risk score after baseline hazard updated ", mean(brier_presegfr_recal_bh), ", 95% CI ", mean(brier_presegfr_recal_bh)-1.96*brier_presegfr_recal_bh_se_pooled, "-", mean(brier_presegfr_recal_bh)+1.96*brier_presegfr_recal_bh_se_pooled))


############################1C RECALIBRATION - OVERALL SLOPE RECALIBRATION################################################################

#before recalibrating, the risk score systematically overestimated risk. However after recalibration, it underestimates risk in those at high risk.

bh_new <- rep(NA, n.imp)
se_bh_new <- rep(NA, n.imp)
cal_slope <- rep(NA, n.imp)
slope_optimism_presegfr <- rep(NA, n.imp)
var_slope <- rep(NA, n.imp)

for (i in 1:n.imp) {
  print(paste0("Calculations in imputation ", i))
  recal_mod2 <- cph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ (ckdpc_50egfr_lin_predictor), data = cohort_presegfr[cohort_presegfr$.imp == i,], x = TRUE, y = TRUE, surv = TRUE)
  
  # Obtain baseline survival estimates for model
  x <- summary(survfit(recal_mod2),time=3)
  bh_new[i] <- x$surv
  se_bh_new[i] <- x$std.err
  
  # bootstrap internal validation
  boot <- validate(recal_mod2)
  slope_optimism_presegfr[i] <- boot[3,5]
  
  # store calibration slope
  cal_slope[i] <- recal_mod2$coefficients * boot[3,5]
  var_slope[i] <- recal_mod2$var
}

bh_presegfr_recal <- mean(bh_new)
bh_presegfr_recal_se <- sqrt(mean(se_bh_new)^2 + (1+1/n.imp)*var(bh_new))
# print baseline hazard with 95% CI
print(paste0("Baseline hazard ", bh_presegfr_recal, ", 95% CI ", bh_presegfr_recal-1.96*bh_presegfr_recal_se, "-", bh_presegfr_recal+1.96*bh_presegfr_recal_se))

coef_presegfr_recal <- mean(cal_slope)
coef_presegfr_recal_se <- sqrt(mean(var_slope) + (1+1/n.imp)*var(cal_slope))
# print calibration slope with 95% CI
print(paste0("Calibration slope ", coef_presegfr_recal, ", 95% CI ", coef_presegfr_recal-1.96*coef_presegfr_recal_se, "-", coef_presegfr_recal+1.96*coef_presegfr_recal_se))

cohort_presegfr <- cohort_presegfr %>% mutate(
  ckdpc_50egfr_lin_predictor_cal=coef_presegfr_recal*ckdpc_50egfr_lin_predictor,
  ckdpc_50egfr_survival_cal=bh_presegfr_recal^exp(ckdpc_50egfr_lin_predictor_cal-mean(ckdpc_50egfr_lin_predictor_cal)), # baseline hazard ^ e ^ centred linear predictor
  ckdpc_50egfr_score_cal=(1-ckdpc_50egfr_survival_cal)*100
)


## Plot
cohort_presegfr$risk_decile <- ntile(cohort_presegfr$ckdpc_50egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- cohort_presegfr %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal)/100)

predicted_all <- cohort_presegfr %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal)/100)

### Find actual observed probabilities by risk score category and studydrug

EST.dpp4isu <- SE.dpp4isu <-
  EST.sglt2i <- SE.sglt2i <-
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_dpp4isu <- tibble() %>% mutate(
  observed_dpp4isu=NA,
  lower_ci_dpp4isu=NA,
  upper_ci_dpp4isu=NA,
  strata=NA
)

observed_sglt2i <- tibble() %>% mutate(
  observed_sglt2i=NA,
  lower_ci_sglt2i=NA,
  upper_ci_sglt2i=NA,
  strata=NA
)

observed_all <- tibble() %>% mutate(
  observed=NA,
  lower_ci=NA,
  upper_ci=NA,
  strata=NA
)

for (k in 1:n.quantiles) {
  for (i in 1:n.imp) {
    
    observed_dpp4isu_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile,
                                      data=cohort_presegfr[cohort_presegfr$.imp == i &
                                                             cohort_presegfr$risk_decile == k &
                                                             cohort_presegfr$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile,
                                     data=cohort_presegfr[cohort_presegfr$.imp == i &
                                                            cohort_presegfr$risk_decile == k &
                                                            cohort_presegfr$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile,
                                  data=cohort_presegfr[cohort_presegfr$.imp == i &
                                                         cohort_presegfr$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd50$estimate
    SE.all[k,i] <- observed_all_ckd50$std.error
    
  }
  
  est.dpp4isu <- pool.rubin.KM(EST.dpp4isu[k,], SE.dpp4isu[k,], n.imp)
  observed_dpp4isu[k,] <- observed_dpp4isu[k,] %>%
    mutate(
      observed_dpp4isu=est.dpp4isu[1],
      lower_ci_dpp4isu=est.dpp4isu[2],
      upper_ci_dpp4isu=est.dpp4isu[3],
      strata=k
    )
  
  est.sglt2i <- pool.rubin.KM(EST.sglt2i[k,], SE.sglt2i[k,], n.imp)
  observed_sglt2i[k,] <- observed_sglt2i[k,] %>%
    mutate(
      observed_sglt2i=est.sglt2i[1],
      lower_ci_sglt2i=est.sglt2i[2],
      upper_ci_sglt2i=est.sglt2i[3],
      strata=k
    )
  
  est.all <- pool.rubin.KM(EST.all[k,], SE.all[k,], n.imp)
  observed_all[k,] <- observed_all[k,] %>%
    mutate(
      observed=est.all[1],
      lower_ci=est.all[2],
      upper_ci=est.all[3],
      strata=k
    )
  
}


dpp4isu_events <- cohort_presegfr %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort_presegfr %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLTi=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug2=="DPP4i/SU")), observed_dpp4isu),
  cbind((predicted %>% filter(studydrug2=="SGLT2i")), observed_sglt2i)
) %>%
  mutate(observed=coalesce(observed_dpp4isu, observed_sglt2i),
         lower_ci=coalesce(lower_ci_dpp4isu,  lower_ci_sglt2i),
         upper_ci=coalesce(upper_ci_dpp4isu, upper_ci_sglt2i))

events_table <- data.frame(t(dpp4isu_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd50_pred=NA, risk_decile=0)

## FINAL PLOT
p_presegfr_cal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = c(0,2.5,5,7.5,10,12.5))+
  scale_y_continuous(limits=c(-1,100)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Recalibrated risk score (calibration slope applied), by risk decile") +
  coord_cartesian(xlim = c(0,13), ylim = c(0,12))


p_presegfr_cal_bydeciles_dpp4isu


## C-stat
recal_mod2 <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_lin_predictor_cal, data=cohort_presegfr, method="breslow")
cstat_presegfr_est <- summary(recal_mod2)$concordance[1]
cstat_presegfr_est_ll <- summary(recal_mod2)$concordance[1]-(1.96*summary(recal_mod2)$concordance[2])
cstat_presegfr_est_ul <- summary(recal_mod2)$concordance[1]+(1.96*summary(recal_mod2)$concordance[2])
paste0("C statistic: ", round(cstat_presegfr_est, 4), ", 95% CI ", round(cstat_presegfr_est_ll, 4), "-", round(cstat_presegfr_est_ul,4))
# C statistic: 0.7262, 95% CI 0.7232-0.7291

## AUC
ROC_presegfr_cal <- roc(cohort_presegfr, ckd_egfr50_censvar, ckdpc_50egfr_score_cal)
auc(ROC_presegfr_cal)
ci.auc(ROC_presegfr_cal)
# the C statistic and AUC have not changed - this is supposed to be the case as the ranking of cases stays the same!

##brier score after overall calibration slope applied
brier_presegfr_recal <- rep(NA, n.imp)
brier_presegfr_recal_se <- rep(NA, n.imp)

for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort_presegfr %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival_cal)
  temp <- temp %>%
    rbind(    # adds rows below your dataset
      temp %>%
        slice(1) %>% # this selects the first patients in your dataset
        mutate(ckd_egfr50_censtime_yrs = 3.5)   # changes censored time to 3.5
    )
  recal_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival_cal, 
                     data=temp, x=T)
  
  score_recal <- 
    Score(object = list(recal_mod), # need to pass cox model to Score() as a list in order for it to be processed
          formula = Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ 1, # null model
          data = temp,
          summary = "ibs", # statistic of interest is integrated brier score
          times = 3,  
          splitMethod = "bootcv",  # use bootstrapping for confidence intervals
          B = n.bootstrap,
          verbose = T) 
  
  brier_presegfr_recal[i] <- score_recal$Brier$score$Brier[2]
  brier_presegfr_recal_se[i] <- score_recal$Brier$score$se[2]
  
  rm(temp)
  
}
brier_presegfr_recal_se_pooled <- sqrt(mean(brier_presegfr_recal_se^2) + (1+1/n.imp)*var(brier_presegfr_recal))

#pool and print brier score
print(paste0("Brier score for risk score after calibration slope applied ", mean(brier_presegfr_recal), ", 95% CI ", mean(brier_presegfr_recal)-1.96*brier_presegfr_recal_se_pooled, "-", mean(brier_presegfr_recal)+1.96*brier_presegfr_recal_se_pooled))

############################2A RISK SCORE eGFR <60 ################################################################

# make variable for risk deciles
cohort_redegfr$risk_decile <- ntile(cohort_redegfr$ckdpc_50egfr_score, n.quantiles)

## Get mean predicted probabilities by risk decile and studydrug
predicted <- cohort_redegfr %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score)/100)

# get mean predicted probabilities by risk decile (not by studydrug)
predicted_all <- cohort_redegfr %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score)/100)

## Find actual observed probabilities by risk score category and studydrug

EST.dpp4isu <- SE.dpp4isu <-
  EST.sglt2i <- SE.sglt2i <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_dpp4isu <- tibble() %>% mutate(
  observed_dpp4isu=NA,
  lower_ci_dpp4isu=NA,
  upper_ci_dpp4isu=NA,
  strata=NA
)

observed_sglt2i <- tibble() %>% mutate(
  observed_sglt2i=NA,
  lower_ci_sglt2i=NA,
  upper_ci_sglt2i=NA,
  strata=NA
)

observed_all <- tibble() %>% mutate(
  observed=NA,
  lower_ci=NA,
  upper_ci=NA,
  strata=NA
)

for (k in 1:n.quantiles) {
  for (i in 1:n.imp) {
    
    observed_dpp4isu_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                      data=cohort_redegfr[cohort_redegfr$.imp == i & 
                                                            cohort_redegfr$risk_decile == k &
                                                            cohort_redegfr$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                     data=cohort_redegfr[cohort_redegfr$.imp == i & 
                                                           cohort_redegfr$risk_decile == k &
                                                           cohort_redegfr$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                  data=cohort_redegfr[cohort_redegfr$.imp == i & 
                                                        cohort_redegfr$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd50$estimate
    SE.all[k,i] <- observed_all_ckd50$std.error
  }
  
  est.dpp4isu <- pool.rubin.KM(EST.dpp4isu[k,], SE.dpp4isu[k,], n.imp)
  observed_dpp4isu[k,] <- observed_dpp4isu[k,] %>% 
    mutate(
      observed_dpp4isu=est.dpp4isu[1],
      lower_ci_dpp4isu=est.dpp4isu[2],
      upper_ci_dpp4isu=est.dpp4isu[3],
      strata=k
    )
  
  est.sglt2i <- pool.rubin.KM(EST.sglt2i[k,], SE.sglt2i[k,], n.imp)
  observed_sglt2i[k,] <- observed_sglt2i[k,] %>% 
    mutate(
      observed_sglt2i=est.sglt2i[1],
      lower_ci_sglt2i=est.sglt2i[2],
      upper_ci_sglt2i=est.sglt2i[3],
      strata=k
    )
  
  est.all <- pool.rubin.KM(EST.all[k,], SE.all[k,], n.imp)
  observed_all[k,] <- observed_all[k,] %>% 
    mutate(
      observed=est.all[1],
      lower_ci=est.all[2],
      upper_ci=est.all[3],
      strata=k
    )
  
}


dpp4isu_events <- cohort_redegfr %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort_redegfr %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLTi=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug2=="DPP4i/SU")), observed_dpp4isu),
  cbind((predicted %>% filter(studydrug2=="SGLT2i")), observed_sglt2i)
) %>%
  mutate(observed=coalesce(observed_dpp4isu, observed_sglt2i),
         lower_ci=coalesce(lower_ci_dpp4isu, lower_ci_sglt2i),
         upper_ci=coalesce(upper_ci_dpp4isu, upper_ci_sglt2i))

events_table <- data.frame(t(dpp4isu_events %>% 
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd50_pred=NA, risk_decile=0)

## FINAL PLOT
p_redegfr_uncal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Raw CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = c(0,2.5,5,7.5,10,12.5))+
  scale_y_continuous(limits=c(-1,100)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Raw risk score, by risk decile") +
  coord_cartesian(xlim = c(0,13), ylim = c(0,12))


p_redegfr_uncal_bydeciles_dpp4isu

## C-stat
cohort_redegfr <- cohort_redegfr %>%
  mutate(ckdpc_50egfr_survival=(100-ckdpc_50egfr_score)/100)

raw_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival, data=cohort_redegfr, method="breslow")
cstat_redegfr_est <- summary(raw_mod)$concordance[1]
cstat_redegfr_est_ll <- summary(raw_mod)$concordance[1]-(1.96*summary(raw_mod)$concordance[2])
cstat_redegfr_est_ul <- summary(raw_mod)$concordance[1]+(1.96*summary(raw_mod)$concordance[2])
paste0("C statistic: ", round(cstat_redegfr_est, 4), ", 95% CI ", round(cstat_redegfr_est_ll, 4), "-", round(cstat_redegfr_est_ul,4))
# C statistic: 0.7262, 95% CI 0.7232-0.7291

## AUC
ROC_redegfr <- roc(cohort_redegfr, ckd_egfr50_censvar, ckdpc_50egfr_score)
auc(ROC_redegfr)
ci.auc(ROC_redegfr)


## brier score for raw risk score
brier_redegfr_raw <- rep(NA, n.imp)
brier_redegfr_raw_se <- rep(NA, n.imp)


for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort_redegfr %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival)
  temp <- temp %>% # error if times = 3, therefore adding extra row with time beyond t=3
    rbind(    # adds rows below your dataset
      temp %>%
        slice(1) %>% # this selects the first patients in your dataset
        mutate(ckd_egfr50_censtime_yrs = 3.5)   # changes censored time to 3.5
    )
  raw_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival, 
                   data=temp, x=T)
  
  
  score_raw <- 
    Score(object = list(raw_mod), # need to pass cox model to Score() as a list in order for it to be processed
          formula = Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ 1, # null model
          data = temp,
          summary = "ibs", # statistic of interest is integrated brier score
          times = 3,  
          splitMethod = "bootcv",  # use bootstrapping for confidence intervals
          B = n.bootstrap,
          verbose = T) 
  
  brier_redegfr_raw[i] <- score_raw$Brier$score$Brier[2]
  brier_redegfr_raw_se[i] <- score_raw$Brier$score$se[2]
  
  rm(temp)
}

brier_redegfr_raw_se_pooled <- sqrt(mean(brier_redegfr_raw_se^2) + (1+1/n.imp)*var(brier_redegfr_raw))

#pool and print brier score
print(paste0("Brier score for raw risk score ", mean(brier_redegfr_raw), ", 95% CI ", mean(brier_redegfr_raw)-1.96*brier_redegfr_raw_se_pooled, "-", mean(brier_redegfr_raw)+1.96*brier_redegfr_raw_se_pooled))


# this risk score underestimates risk and needs to be recalibrated.

############################2B RECALIBRATION - BASELINE HAZARD UPDATE ################################################################

# initially we will assess calibration if we only update the baseline hazard
# this is done by fitting a cox proportional model with the linear predictor as the only variable as an offset.
# as this would only involve adjusting the overall event probability for this cohort_redegfr, no resampling (internal validation) is needed
recal_mod <- cph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ stats::offset(ckdpc_50egfr_lin_predictor), 
                 data = cohort_redegfr[!cohort_redegfr$studydrug=="SGLT2i",], x = TRUE, y = TRUE, surv = TRUE)

x <- summary(survfit(recal_mod),time=3)
bh_update_redegfr <- x$surv
bh_update_redegfr_se <- x$std.err
# print updated baseline hazard
print(paste0("Baseline hazard ", bh_update_redegfr, ", 95% CI ", bh_update_redegfr-1.96*bh_update_redegfr_se, "-", bh_update_redegfr+1.96*bh_update_redegfr_se))

cohort_redegfr <- cohort_redegfr %>% mutate(
  ckdpc_50egfr_survival_cal_bh=bh_update_redegfr^exp(ckdpc_50egfr_lin_predictor-mean(ckdpc_50egfr_lin_predictor)), # baseline hazard ^ e ^ centred linear predictor
  ckdpc_50egfr_score_cal_bh=(1-ckdpc_50egfr_survival_cal_bh)*100
)

## Show observed (estimated frequency binned per risk decile) vs predicted (risk-score predicted)
cohort_redegfr$risk_decile <- ntile(cohort_redegfr$ckdpc_50egfr_score_cal_bh, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- cohort_redegfr %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal_bh)/100)

predicted_all <- cohort_redegfr %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal_bh)/100)

### Find actual observed probabilities by risk score category and studydrug

EST.dpp4isu <- SE.dpp4isu <-
  EST.sglt2i <- SE.sglt2i <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_dpp4isu <- tibble() %>% mutate(
  observed_dpp4isu=NA,
  lower_ci_dpp4isu=NA,
  upper_ci_dpp4isu=NA,
  strata=NA
)

observed_sglt2i <- tibble() %>% mutate(
  observed_sglt2i=NA,
  lower_ci_sglt2i=NA,
  upper_ci_sglt2i=NA,
  strata=NA
)

observed_all <- tibble() %>% mutate(
  observed=NA,
  lower_ci=NA,
  upper_ci=NA,
  strata=NA
)

for (k in 1:n.quantiles) {
  for (i in 1:n.imp) {
    
    observed_dpp4isu_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                      data=cohort_redegfr[cohort_redegfr$.imp == i & 
                                                            cohort_redegfr$risk_decile == k &
                                                            cohort_redegfr$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                     data=cohort_redegfr[cohort_redegfr$.imp == i & 
                                                           cohort_redegfr$risk_decile == k &
                                                           cohort_redegfr$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                  data=cohort_redegfr[cohort_redegfr$.imp == i & 
                                                        cohort_redegfr$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd50$estimate
    SE.all[k,i] <- observed_all_ckd50$std.error
    
  }
  
  est.dpp4isu <- pool.rubin.KM(EST.dpp4isu[k,], SE.dpp4isu[k,], n.imp)
  observed_dpp4isu[k,] <- observed_dpp4isu[k,] %>% 
    mutate(
      observed_dpp4isu=est.dpp4isu[1],
      lower_ci_dpp4isu=est.dpp4isu[2],
      upper_ci_dpp4isu=est.dpp4isu[3],
      strata=k
    )
  
  est.sglt2i <- pool.rubin.KM(EST.sglt2i[k,], SE.sglt2i[k,], n.imp)
  observed_sglt2i[k,] <- observed_sglt2i[k,] %>% 
    mutate(
      observed_sglt2i=est.sglt2i[1],
      lower_ci_sglt2i=est.sglt2i[2],
      upper_ci_sglt2i=est.sglt2i[3],
      strata=k
    )
  
  est.all <- pool.rubin.KM(EST.all[k,], SE.all[k,], n.imp)
  observed_all[k,] <- observed_all[k,] %>% 
    mutate(
      observed=est.all[1],
      lower_ci=est.all[2],
      upper_ci=est.all[3],
      strata=k
    )
  
}


dpp4isu_events <- cohort_redegfr %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort_redegfr %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLTi=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug2=="DPP4i/SU")), observed_dpp4isu),
  cbind((predicted %>% filter(studydrug2=="SGLT2i")), observed_sglt2i)
) %>%
  mutate(observed=coalesce(observed_dpp4isu, observed_sglt2i),
         lower_ci=coalesce(lower_ci_dpp4isu,lower_ci_sglt2i),
         upper_ci=coalesce(upper_ci_dpp4isu, upper_ci_sglt2i))

events_table <- data.frame(t(dpp4isu_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd50_pred=NA, risk_decile=0)

## FINAL PLOT
p_redegfr_interim_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = c(0,2.5,5,7.5,10,12.5))+
  scale_y_continuous(limits=c(-1,100)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Recalibrated risk score (baseline hazard updated), by risk decile") +
  coord_cartesian(xlim = c(0,13), ylim = c(0,12))


p_redegfr_interim_bydeciles_dpp4isu

## C-stat

surv_mod_ckd50 <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival_cal_bh, data=cohort_redegfr, method="breslow")
cstat_redegfr_est <- summary(surv_mod_ckd50)$concordance[1]
cstat_redegfr_est_ll <- summary(surv_mod_ckd50)$concordance[1]-(1.96*summary(surv_mod_ckd50)$concordance[2])
cstat_redegfr_est_ul <- summary(surv_mod_ckd50)$concordance[1]+(1.96*summary(surv_mod_ckd50)$concordance[2])
paste0("C statistic: ", round(cstat_redegfr_est, 4), ", 95% CI ", round(cstat_redegfr_est_ll, 4), "-", round(cstat_redegfr_est_ul,4))
# C statistic: 0.7262, 95% CI 0.7232-0.7291

## AUC
ROC_redegfr_cal <- roc(cohort_redegfr, ckd_egfr50_censvar, ckdpc_50egfr_score_cal_bh)
auc(ROC_redegfr_cal)
ci.auc(ROC_redegfr_cal)

##brier score after baseline hazard updated
brier_redegfr_recal_bh <- rep(NA, n.imp)
brier_redegfr_recal_bh_se <- rep(NA, n.imp)
for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort_redegfr %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival_cal_bh)
  temp <- temp %>%
    rbind(    # adds rows below your dataset
      temp %>%
        slice(1) %>% # this selects the first patients in your dataset
        mutate(ckd_egfr50_censtime_yrs = 3.5)   # changes censored time to 3.5
    )
  recal_bh_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival_cal_bh, 
                        data=temp, x=T)
  
  score_recal_bh <- 
    Score(object = list(recal_bh_mod), # need to pass cox model to Score() as a list in order for it to be processed
          formula = Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ 1, # null model
          data = temp,
          summary = "ibs", # statistic of interest is integrated brier score
          times = 3,  
          splitMethod = "bootcv",  # use bootstrapping for confidence intervals
          B = n.bootstrap,
          verbose = T) 
  
  brier_redegfr_recal_bh[i] <- score_recal_bh$Brier$score$Brier[2]
  brier_redegfr_recal_bh_se[i] <- score_recal_bh$Brier$score$se[2]
  
  rm(temp)
  
}

brier_redegfr_recal_bh_se_pooled <- sqrt(mean(brier_redegfr_recal_bh_se^2) + (1+1/n.imp)*var(brier_redegfr_recal_bh))

#pool and print brier score
print(paste0("Brier score for risk score after baseline hazard updated ", mean(brier_redegfr_recal_bh), ", 95% CI ", mean(brier_redegfr_recal_bh)-1.96*brier_redegfr_recal_bh_se_pooled, "-", mean(brier_redegfr_recal_bh)+1.96*brier_redegfr_recal_bh_se_pooled))


############################2C RECALIBRATION - OVERALL SLOPE RECALIBRATION################################################################

#before recalibrating, the risk score systematically overestimated risk. However after recalibration, it underestimates risk in those at high risk.

bh_new <- rep(NA, n.imp)
se_bh_new <- rep(NA, n.imp)
cal_slope <- rep(NA, n.imp)
slope_optimism_redegfr <- rep(NA, n.imp)
var_slope <- rep(NA, n.imp)

for (i in 1:n.imp) {
  print(paste0("Calculations in imputation ", i))
  recal_mod2 <- cph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ (ckdpc_50egfr_lin_predictor), data = cohort_redegfr[cohort_redegfr$.imp == i,], x = TRUE, y = TRUE, surv = TRUE)
  
  # Obtain baseline survival estimates for model
  x <- summary(survfit(recal_mod2),time=3)
  bh_new[i] <- x$surv
  se_bh_new[i] <- x$std.err
  
  # bootstrap internal validation
  boot <- validate(recal_mod2)
  slope_optimism_redegfr[i] <- boot[3,5]
  
  # store calibration slope
  cal_slope[i] <- recal_mod2$coefficients * boot[3,5]
  var_slope[i] <- recal_mod2$var
}

bh_redegfr_recal <- mean(bh_new)
bh_redegfr_recal_se <- sqrt(mean(se_bh_new)^2 + (1+1/n.imp)*var(bh_new))
# print baseline hazard with 95% CI
print(paste0("Baseline hazard ", bh_redegfr_recal, ", 95% CI ", bh_redegfr_recal-1.96*bh_redegfr_recal_se, "-", bh_redegfr_recal+1.96*bh_redegfr_recal_se))

coef_redegfr_recal <- mean(cal_slope)
coef_redegfr_recal_se <- sqrt(mean(var_slope) + (1+1/n.imp)*var(cal_slope))
# print calibration slope with 95% CI
print(paste0("Calibration slope ", coef_redegfr_recal, ", 95% CI ", coef_redegfr_recal-1.96*coef_redegfr_recal_se, "-", coef_redegfr_recal+1.96*coef_redegfr_recal_se))

cohort_redegfr <- cohort_redegfr %>% mutate(
  ckdpc_50egfr_lin_predictor_cal=coef_redegfr_recal*ckdpc_50egfr_lin_predictor,
  ckdpc_50egfr_survival_cal=bh_redegfr_recal^exp(ckdpc_50egfr_lin_predictor_cal-mean(ckdpc_50egfr_lin_predictor_cal)), # baseline hazard ^ e ^ centred linear predictor
  ckdpc_50egfr_score_cal=(1-ckdpc_50egfr_survival_cal)*100
)


## Plot
cohort_redegfr$risk_decile <- ntile(cohort_redegfr$ckdpc_50egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- cohort_redegfr %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal)/100)

predicted_all <- cohort_redegfr %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal)/100)

### Find actual observed probabilities by risk score category and studydrug

EST.dpp4isu <- SE.dpp4isu <-
  EST.sglt2i <- SE.sglt2i <-
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_dpp4isu <- tibble() %>% mutate(
  observed_dpp4isu=NA,
  lower_ci_dpp4isu=NA,
  upper_ci_dpp4isu=NA,
  strata=NA
)

observed_sglt2i <- tibble() %>% mutate(
  observed_sglt2i=NA,
  lower_ci_sglt2i=NA,
  upper_ci_sglt2i=NA,
  strata=NA
)

observed_all <- tibble() %>% mutate(
  observed=NA,
  lower_ci=NA,
  upper_ci=NA,
  strata=NA
)

for (k in 1:n.quantiles) {
  for (i in 1:n.imp) {
    
    observed_dpp4isu_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile,
                                      data=cohort_redegfr[cohort_redegfr$.imp == i &
                                                            cohort_redegfr$risk_decile == k &
                                                            cohort_redegfr$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile,
                                     data=cohort_redegfr[cohort_redegfr$.imp == i &
                                                           cohort_redegfr$risk_decile == k &
                                                           cohort_redegfr$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile,
                                  data=cohort_redegfr[cohort_redegfr$.imp == i &
                                                        cohort_redegfr$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd50$estimate
    SE.all[k,i] <- observed_all_ckd50$std.error
    
  }
  
  est.dpp4isu <- pool.rubin.KM(EST.dpp4isu[k,], SE.dpp4isu[k,], n.imp)
  observed_dpp4isu[k,] <- observed_dpp4isu[k,] %>%
    mutate(
      observed_dpp4isu=est.dpp4isu[1],
      lower_ci_dpp4isu=est.dpp4isu[2],
      upper_ci_dpp4isu=est.dpp4isu[3],
      strata=k
    )
  
  est.sglt2i <- pool.rubin.KM(EST.sglt2i[k,], SE.sglt2i[k,], n.imp)
  observed_sglt2i[k,] <- observed_sglt2i[k,] %>%
    mutate(
      observed_sglt2i=est.sglt2i[1],
      lower_ci_sglt2i=est.sglt2i[2],
      upper_ci_sglt2i=est.sglt2i[3],
      strata=k
    )
  
  est.all <- pool.rubin.KM(EST.all[k,], SE.all[k,], n.imp)
  observed_all[k,] <- observed_all[k,] %>%
    mutate(
      observed=est.all[1],
      lower_ci=est.all[2],
      upper_ci=est.all[3],
      strata=k
    )
  
}


dpp4isu_events <- cohort_redegfr %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort_redegfr %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLTi=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug2=="DPP4i/SU")), observed_dpp4isu),
  cbind((predicted %>% filter(studydrug2=="SGLT2i")), observed_sglt2i)
) %>%
  mutate(observed=coalesce(observed_dpp4isu, observed_sglt2i),
         lower_ci=coalesce(lower_ci_dpp4isu,  lower_ci_sglt2i),
         upper_ci=coalesce(upper_ci_dpp4isu, upper_ci_sglt2i))

events_table <- data.frame(t(dpp4isu_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd50_pred=NA, risk_decile=0)

## FINAL PLOT
p_redegfr_cal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = c(0,2.5,5,7.5,10,12.5))+
  scale_y_continuous(limits=c(-1,100)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Recalibrated risk score (calibration slope applied), by risk decile") +
  coord_cartesian(xlim = c(0,13), ylim = c(0,12))


p_redegfr_cal_bydeciles_dpp4isu


## C-stat
recal_mod2 <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_lin_predictor_cal, data=cohort_redegfr, method="breslow")
cstat_redegfr_est <- summary(recal_mod2)$concordance[1]
cstat_redegfr_est_ll <- summary(recal_mod2)$concordance[1]-(1.96*summary(recal_mod2)$concordance[2])
cstat_redegfr_est_ul <- summary(recal_mod2)$concordance[1]+(1.96*summary(recal_mod2)$concordance[2])
paste0("C statistic: ", round(cstat_redegfr_est, 4), ", 95% CI ", round(cstat_redegfr_est_ll, 4), "-", round(cstat_redegfr_est_ul,4))
# C statistic: 0.7262, 95% CI 0.7232-0.7291

## AUC
ROC_redegfr_cal <- roc(cohort_redegfr, ckd_egfr50_censvar, ckdpc_50egfr_score_cal)
auc(ROC_redegfr_cal)
ci.auc(ROC_redegfr_cal)
# the C statistic and AUC have not changed - this is supposed to be the case as the ranking of cases stays the same!

##brier score after overall calibration slope applied
brier_redegfr_recal <- rep(NA, n.imp)
brier_redegfr_recal_se <- rep(NA, n.imp)

for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort_redegfr %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival_cal)
  temp <- temp %>%
    rbind(    # adds rows below your dataset
      temp %>%
        slice(1) %>% # this selects the first patients in your dataset
        mutate(ckd_egfr50_censtime_yrs = 3.5)   # changes censored time to 3.5
    )
  recal_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival_cal, 
                     data=temp, x=T)
  
  score_recal <- 
    Score(object = list(recal_mod), # need to pass cox model to Score() as a list in order for it to be processed
          formula = Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ 1, # null model
          data = temp,
          summary = "ibs", # statistic of interest is integrated brier score
          times = 3,  
          splitMethod = "bootcv",  # use bootstrapping for confidence intervals
          B = n.bootstrap,
          verbose = T) 
  
  brier_redegfr_recal[i] <- score_recal$Brier$score$Brier[2]
  brier_redegfr_recal_se[i] <- score_recal$Brier$score$se[2]
  
  rm(temp)
  
}
brier_redegfr_recal_se_pooled <- sqrt(mean(brier_redegfr_recal_se^2) + (1+1/n.imp)*var(brier_redegfr_recal))

#pool and print brier score
print(paste0("Brier score for risk score after calibration slope applied ", mean(brier_redegfr_recal), ", 95% CI ", mean(brier_redegfr_recal)-1.96*brier_redegfr_recal_se_pooled, "-", mean(brier_redegfr_recal)+1.96*brier_redegfr_recal_se_pooled))

#### model details summed up:

# eGFR ≥60:
auc(ROC_presegfr_cal)
ci.auc(ROC_presegfr_cal)
paste0("C statistic: ", round(cstat_presegfr_est, 4), ", 95% CI ", round(cstat_presegfr_est_ll, 4), "-", round(cstat_presegfr_est_ul,4))
print(paste0("Brier score for raw risk score ", mean(brier_presegfr_raw), ", 95% CI ", mean(brier_presegfr_raw)-1.96*brier_presegfr_raw_se_pooled, "-", mean(brier_presegfr_raw)+1.96*brier_presegfr_raw_se_pooled))
print(paste0("Brier score for risk score after baseline hazard updated ", mean(brier_presegfr_recal_bh), ", 95% CI ", mean(brier_presegfr_recal_bh)-1.96*brier_presegfr_recal_bh_se_pooled, "-", mean(brier_presegfr_recal_bh)+1.96*brier_presegfr_recal_bh_se_pooled))
print(paste0("Baseline hazard ", bh_update_presegfr, ", 95% CI ", bh_update_presegfr-1.96*bh_update_presegfr_se, "-", bh_update_presegfr+1.96*bh_update_presegfr_se))
print(paste0("Brier score for risk score after calibration slope applied ", mean(brier_presegfr_recal), ", 95% CI ", mean(brier_presegfr_recal)-1.96*brier_presegfr_recal_se_pooled, "-", mean(brier_presegfr_recal)+1.96*brier_presegfr_recal_se_pooled))
print(paste0("Baseline hazard ", bh_presegfr_recal, ", 95% CI ", bh_presegfr_recal-1.96*bh_presegfr_recal_se, "-", bh_presegfr_recal+1.96*bh_presegfr_recal_se))
print(paste0("Calibration slope ", coef_presegfr_recal, ", 95% CI ", coef_presegfr_recal-1.96*coef_presegfr_recal_se, "-", coef_presegfr_recal+1.96*coef_presegfr_recal_se))
print(paste0("Slope optimism ", mean(slope_optimism_presegfr)))

# eGFR <60:
auc(ROC_redegfr_cal)
ci.auc(ROC_redegfr_cal)
paste0("C statistic: ", round(cstat_redegfr_est, 4), ", 95% CI ", round(cstat_redegfr_est_ll, 4), "-", round(cstat_redegfr_est_ul,4))
print(paste0("Brier score for raw risk score ", mean(brier_redegfr_raw), ", 95% CI ", mean(brier_redegfr_raw)-1.96*brier_redegfr_raw_se_pooled, "-", mean(brier_redegfr_raw)+1.96*brier_redegfr_raw_se_pooled))
print(paste0("Brier score for risk score after baseline hazard updated ", mean(brier_redegfr_recal_bh), ", 95% CI ", mean(brier_redegfr_recal_bh)-1.96*brier_redegfr_recal_bh_se_pooled, "-", mean(brier_redegfr_recal_bh)+1.96*brier_redegfr_recal_bh_se_pooled))
print(paste0("Baseline hazard ", bh_update_redegfr, ", 95% CI ", bh_update_redegfr-1.96*bh_update_redegfr_se, "-", bh_update_redegfr+1.96*bh_update_redegfr_se))
print(paste0("Brier score for risk score after calibration slope applied ", mean(brier_redegfr_recal), ", 95% CI ", mean(brier_redegfr_recal)-1.96*brier_redegfr_recal_se_pooled, "-", mean(brier_redegfr_recal)+1.96*brier_redegfr_recal_se_pooled))
print(paste0("Baseline hazard ", bh_redegfr_recal, ", 95% CI ", bh_redegfr_recal-1.96*bh_redegfr_recal_se, "-", bh_redegfr_recal+1.96*bh_redegfr_recal_se))
print(paste0("Calibration slope ", coef_redegfr_recal, ", 95% CI ", coef_redegfr_recal-1.96*coef_redegfr_recal_se, "-", coef_redegfr_recal+1.96*coef_redegfr_recal_se))
print(paste0("Slope optimism ", mean(slope_optimism_redegfr)))

############################3 STORE RECALIBRATED SCORES################################################################
# save dataset with calibrated risk score so this can be used in the subsequent scripts
cohort <- bind_rows(cohort_presegfr, cohort_redegfr)
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated.Rda"))
