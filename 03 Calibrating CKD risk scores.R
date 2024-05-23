## In this script our aim is to predict baseline risk using the CKD-PC risk scores.
## we will be looking at 2 established risk scores for: 
##  - developing CKD stage 3 or more (egfr60 or ckd345).
##  - 40% decline in eGFR/ESKD (40egfr).
## which we will review regarding need for recalibration in the CPRD cohort.

## Contents
# 0 setup
# 1 review uncalibrated risk score (40% decline in eGFR)
# 2 recalibration of risk score (40% decline in eGFR)
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

options(dplyr.dpp4isummarise.inform = FALSE)

# set seed
set.seed(123)

# number of imputations 
n.imp <- 10

# number of bootstraps
n.bootstrap <- 500

# number of quantiles (for risk scores later on)
n.quantiles <- 10

# today's date
today <- as.character(Sys.Date(), format="%Y%m%d")

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
load("2024-04-30_t2d_ckdpc_imputed_data_withweights.Rda")

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

# select calibration cohort and non-calibration cohort 

# # assign random 20% (SU/DPP4) as recalibration cohort and remove from main cohort
# # we will do this in each imputation and combine 
# cal_cohort <- cohort %>% filter(!studydrug2 == "SGLT2i") %>% group_by(.imp) %>% slice_sample(prop=0.2)
# 
# # select those not in the calibration cohort
# cohort <- cohort %>%
#   anti_join(cal_cohort, by=c("patid", "dstartdate", "studydrug2", ".imp"))


# set default colour-blind accessible colours for figures later on
cols <- c("SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00")
#in further analyses, the dpp4/su group will be combined, and we will use the dpp4 colour for this (strongest contrast)
cols <- c(cols, "DPP4i/SU" = "#0072B2")
cols <- cols[names(cols) %in% cohort$studydrug2]
cols <- cols[order(names(cols))]

############################1 UNCALIBRATED SCORE - 40% DECLINE IN EGFR################################################################

# 40% decline in eGFR / ESKD

# make variable for risk deciles
cohort$ckd40_risk_decile <- ntile(cohort$ckdpc_40egfr_score, n.quantiles)

## Get mean predicted probabilities by risk decile and studydrug
predicted <- cohort %>%
  group_by(ckd40_risk_decile, studydrug2) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score)/100)

# get mean predicted probabilities by risk decile (not by studydrug)
predicted_all <- cohort %>%
  group_by(ckd40_risk_decile) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score)/100)

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
    
    observed_dpp4isu_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile, 
                                      data=cohort[cohort$.imp == i & 
                                                    cohort$ckd40_risk_decile == k &
                                                    cohort$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd40$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd40$std.error
    
    
    observed_sglt2i_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile, 
                                     data=cohort[cohort$.imp == i & 
                                                   cohort$ckd40_risk_decile == k &
                                                   cohort$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd40$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd40$std.error
    
    observed_all_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile, 
                                  data=cohort[cohort$.imp == i & 
                                                cohort$ckd40_risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd40$estimate
    SE.all[k,i] <- observed_all_ckd40$std.error
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


dpp4isu_events <- cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_decile) %>%
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
  filter(rowname!="ckd40_risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(ckd40_risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd40_pred=NA, ckd40_risk_decile=0)

## FINAL PLOT
p_ckd40_uncal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Uncalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,8))+
  scale_y_continuous(limits=c(-1,7)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Calibration plot of CKD-PC risk score for 40% decline in eGFR / ESKD", subtitle = "Uncalibrated risk score, binned by risk decile") +
  coord_cartesian(xlim = c(0,7.5), ylim = c(0,7.5))


p_ckd40_uncal_bydeciles_dpp4isu

## C-stat
raw_mod <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_lin_predictor, data=cohort, method="breslow")
cstat_est <- summary(raw_mod)$concordance[1]
cstat_est_ll <- summary(raw_mod)$concordance[1]-(1.96*summary(raw_mod)$concordance[2])
cstat_est_ul <- summary(raw_mod)$concordance[1]+(1.96*summary(raw_mod)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7693, 95% CI 0.7642-0.7744

# calibration slope:
raw_mod$coefficients

# brier_score <- brier(surv_mod_ckd40, times = 3)
# paste0("Brier score: ", round(brier_score$brier,4))
ROC <- roc(cohort, ckd_egfr40_censvar, ckdpc_40egfr_score)
auc(ROC)

# this risk score overestimates risk and needs to be recalibrated.

############################2A RECALIBRATION - BASELINE HAZARD UPDATE ################################################################

# no need for resampling - this would only be adjusting the overall event probability to this cohort
recal_mod <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ stats::offset(ckdpc_40egfr_lin_predictor), 
                 data = cohort[!cohort$studydrug=="SGLT2i",], x = TRUE, y = TRUE, surv = TRUE)

x <- summary(survfit(recal_mod),time=3)
bh_update <- x$surv
bh_update_se <- x$std.error
# print updated baseline hazard
print(paste0("Baseline hazard ", bh_update, ", 95% CI ", bh_update-1.96*bh_update_se, "-", bh_update+1.96*bh_update_se))

cohort <- cohort %>% mutate(
  ckdpc_40egfr_survival_cal=bh_update^exp(ckdpc_40egfr_lin_predictor-mean(ckdpc_40egfr_lin_predictor)), # baseline hazard ^ e ^ centred linear predictor
  ckdpc_40egfr_score_cal=(1-ckdpc_40egfr_survival_cal)*100
)

## Show observed (estimated frequency binned per risk decile) vs predicted (risk-score predicted)
cohort$ckd40_risk_decile <- ntile(cohort$ckdpc_40egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- cohort %>%
  group_by(ckd40_risk_decile, studydrug2) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

predicted_all <- cohort %>%
  group_by(ckd40_risk_decile) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

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
    
    observed_dpp4isu_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile, 
                                      data=cohort[cohort$.imp == i & 
                                                    cohort$ckd40_risk_decile == k &
                                                    cohort$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd40$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd40$std.error
    
    observed_sglt2i_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile, 
                                     data=cohort[cohort$.imp == i & 
                                                   cohort$ckd40_risk_decile == k &
                                                   cohort$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd40$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd40$std.error
    
    observed_all_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile, 
                                  data=cohort[cohort$.imp == i & 
                                                cohort$ckd40_risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd40$estimate
    SE.all[k,i] <- observed_all_ckd40$std.error
    
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


dpp4isu_events <- cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_decile) %>%
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
  filter(rowname!="ckd40_risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(ckd40_risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd40_pred=NA, ckd40_risk_decile=0)

## FINAL PLOT
p_ckd40_interim_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,8))+
  scale_y_continuous(limits=c(-1,7)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Calibration plot of CKD-PC risk score for 40% decline in eGFR / ESKD", subtitle = "Baseline hazard (intercept) recalibration of risk score, binned by risk decile") +
  coord_cartesian(xlim = c(0,7.5), ylim = c(0,7.5))


p_ckd40_interim_bydeciles_dpp4isu

## C-stat

surv_mod_ckd40 <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival_cal, data=cohort, method="breslow")
cstat_est <- summary(surv_mod_ckd40)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd40)$concordance[1]-(1.96*summary(surv_mod_ckd40)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd40)$concordance[1]+(1.96*summary(surv_mod_ckd40)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7693, 95% CI 0.7642-0.7744
# brier_score <- brier(surv_mod_ckd40, times = 3)
# paste0("Brier score: ", round(brier_score$brier,4))
ROC_cal <- roc(cohort, ckd_egfr40_censvar, ckdpc_40egfr_score_cal)
auc(ROC_cal)

############################2B RECALIBRATION - OVERALL SLOPE RECALIBRATION################################################################

#before recalibrating, the risk score systematically overestimated risk. However after recalibration, it underestimates risk in those at high risk.

bh_new <- rep(NA, n.imp)
se_bh_new <- rep(NA, n.imp)
cal_slope <- rep(NA, n.imp)
slope_optimism <- rep(NA, n.imp)
var_slope <- rep(NA, n.imp)
for (i in 1:n.imp) {
  print(paste0("Calculations in imputation ", i))
  recal_mod2 <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ (ckdpc_40egfr_lin_predictor), data = cohort[cohort$.imp == i,], x = TRUE, y = TRUE, surv = TRUE)
  
  # Obtain baseline survival estimates for model
  x <- summary(survfit(recal_mod2),time=3)
  bh_new[i] <- x$surv
  se_bh_new[i] <- x$std.err
  
  # bootstrap internal validation
  boot <- validate(recal_mod2)
  slope_optimism[i] <- boot[3,5]
  
  # store calibration slope
  cal_slope[i] <- recal_mod2$coefficients * boot[3,5]
  var_slope[i] <- recal_mod2$var
}

bh_recal <- mean(bh_update)
bh_recal_se <- sqrt(mean(se_bh_new)^2 + 1+1/n.imp*var(bh_recal))
# print baseline hazard with 95% CI
print(paste0("Baseline hazard ", bh_recal, ", 95% CI ", bh_recal-1.96*bh_recal_se, "-", bh_recal+1.96*bh_recal_se))

coef_recal <- mean(cal_slope)
coef_recal_se <- sqrt(mean(var_slope) + (1+1/n.imp)*var(cal_slope))
# print calibration slope with 95% CI
print(paste0("Calibration slope ", coef_recal, ", 95% CI ", coef_recal-1.96*coef_recal_se, "-", coef_recal+1.96*coef_recal_se))

cohort <- cohort %>% mutate(
  ckdpc_40egfr_lin_predictor_cal=coef_recal*ckdpc_40egfr_lin_predictor,
  ckdpc_40egfr_survival_cal=bh_recal^exp(ckdpc_40egfr_lin_predictor_cal-mean(ckdpc_40egfr_lin_predictor_cal)), # baseline hazard ^ e ^ centred linear predictor
  ckdpc_40egfr_score_cal=(1-ckdpc_40egfr_survival_cal)*100
)


## Plot
cohort$ckd40_risk_decile <- ntile(cohort$ckdpc_40egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- cohort %>%
  group_by(ckd40_risk_decile, studydrug2) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

predicted_all <- cohort %>%
  group_by(ckd40_risk_decile) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

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
    
    observed_dpp4isu_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile,
                                      data=cohort[cohort$.imp == i &
                                                    cohort$ckd40_risk_decile == k &
                                                    cohort$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd40$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd40$std.error
    
    observed_sglt2i_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile,
                                     data=cohort[cohort$.imp == i &
                                                   cohort$ckd40_risk_decile == k &
                                                   cohort$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd40$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd40$std.error
    
    observed_all_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile,
                                  data=cohort[cohort$.imp == i &
                                                cohort$ckd40_risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd40$estimate
    SE.all[k,i] <- observed_all_ckd40$std.error
    
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


dpp4isu_events <- cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_decile) %>%
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
  filter(rowname!="ckd40_risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(ckd40_risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd40_pred=NA, ckd40_risk_decile=0)

## FINAL PLOT
p_ckd40_cal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,6))+
  scale_y_continuous(limits=c(-1,7)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Calibration plot of CKD-PC risk score for 40% decline in eGFR / ESKD", subtitle = "Recalibrated risk score, binned by risk decile") +
  coord_cartesian(xlim = c(0,6), ylim = c(0,6))


p_ckd40_cal_bydeciles_dpp4isu


## C-stat
recal_mod2 <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_lin_predictor_cal, data=cohort, method="breslow")
cstat_est <- summary(recal_mod2)$concordance[1]
cstat_est_ll <- summary(recal_mod2)$concordance[1]-(1.96*summary(recal_mod2)$concordance[2])
cstat_est_ul <- summary(recal_mod2)$concordance[1]+(1.96*summary(recal_mod2)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7693, 95% CI 0.7642-0.7744
# brier_score <- brier(surv_mod_ckd40, times = 3)
#  paste0("Brier score: ", round(brier_score$brier,4))
ROC_cal <- roc(cohort, ckd_egfr40_censvar, ckdpc_40egfr_score_cal)
auc(ROC_cal)


############################3 STORE RECALIBRATED SCORES################################################################
# save dataset with calibrated risk score so this can be used in the subsequent scripts
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated.Rda"))
