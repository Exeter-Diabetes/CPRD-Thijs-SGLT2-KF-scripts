## In this script our aim is to predict baseline risk using the CKD-PC risk scores.
## we will be looking at the risk score for 50% decline in eGFR/ESKD (50egfr).
## which we will recalibrate in the CPRD cohort for endpoint 50% decline in eGFR/ESKD.

## Contents
# 0 setup
# 1 uncalibrated risk score
# 2 baseline hazard update
# 3 calibration slope update
# 4 store dataset with recalibrated risk score
# 5 spline plot of HR by baseline risk

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
today <- "2024-10-29"
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
load("2024-10-29_t2d_ckdpc_imputed_data_withweights.Rda")

# select imputed data only (ie. remove non-imputed data)
cohort <- cohort[cohort$.imp > 0,]

## remove double overlapping entries for DPP4i and SU that overlap (take one only)
cohort <- cohort %>% group_by(.imp, patid) %>% filter(
  !duplicated(studydrug2)
) %>% ungroup()

# check number of subjects
table(cohort$studydrug)
# SU  DPP4i SGLT2i 
# 363910 563110 550020   # 10 imputations therefore number of subjects per group appears 10 times larger

# set default colour-blind accessible colours for figures later on
cols <- c("SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00")
#in further analyses, the dpp4/su group will be combined, and we will use the dpp4 colour for this (strongest contrast)
cols <- c(cols, "DPP4i/SU" = "#0072B2")
cols <- cols[names(cols) %in% cohort$studydrug2]
cols <- cols[order(names(cols))]

# covariates for multivariable adjustment
covariates <- c("dstartdate_age", "malesex", "imd2015_10", "ethnicity_4cat", "initiation_year", "prebmi", "prehba1c",
                "pretotalcholesterol", "preegfr", "uacr", "presbp", "ckdpc_50egfr_score_cal", "ncurrtx", "statin", "INS", 
                "ACEi_or_ARB", "smoking_status", "dstartdate_dm_dur_all", "predrug_hypertension", "predrug_af", "hosp_admission_prev_year")


############################1 UNCALIBRATED RISK SCORE################################################################


# make variable for risk deciles
cohort$risk_decile <- ntile(cohort$ckdpc_50egfr_score, n.quantiles)

## Get mean predicted probabilities by risk decile and studydrug
predicted <- cohort %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score)/100)

# get mean predicted probabilities by risk decile (not by studydrug)
predicted_all <- cohort %>%
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
                                      data=cohort[cohort$.imp == i & 
                                                             cohort$risk_decile == k &
                                                             cohort$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                     data=cohort[cohort$.imp == i & 
                                                            cohort$risk_decile == k &
                                                            cohort$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                  data=cohort[cohort$.imp == i & 
                                                         cohort$risk_decile == k,]) %>%
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


dpp4isu_events <- cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort %>%
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
p_uncal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Uncalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = seq(0, 5, 1))+
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
  ggtitle("Uncalibrated risk score, by risk decile") +
  coord_cartesian(xlim = c(0,5), ylim = c(0,5))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_uncalibrated_risk_score_calibration.tiff"), width=6, height=5.5, units = "in", res=800) 
p_uncal_bydeciles_dpp4isu
dev.off()

## C-stat
cohort <- cohort %>%
  mutate(ckdpc_50egfr_survival=(100-ckdpc_50egfr_score)/100)

raw_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival, data=cohort[!cohort$studydrug2 == "SGLT2i",], method="breslow")
cstat_est <- summary(raw_mod)$concordance[1]
cstat_est_ll <- summary(raw_mod)$concordance[1]-(1.96*summary(raw_mod)$concordance[2])
cstat_est_ul <- summary(raw_mod)$concordance[1]+(1.96*summary(raw_mod)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7262, 95% CI 0.7232-0.7291

## AUC
ROC_raw <- roc(cohort, ckd_egfr50_censvar, ckdpc_50egfr_score)
auc(ROC_raw)
ci.auc(ROC_raw)


## brier score for raw risk score
brier_raw <- rep(NA, n.imp)
brier_raw_se <- rep(NA, n.imp)


for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival)
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
  
  brier_raw[i] <- score_raw$Brier$score$Brier[2]
  brier_raw_se[i] <- score_raw$Brier$score$se[2]
  
  rm(temp)
}

brier_raw_se_pooled <- sqrt(mean(brier_raw_se^2) + (1+1/n.imp)*var(brier_raw))

#pool and print brier score
print(paste0("Brier score for raw risk score ", mean(brier_raw), ", 95% CI ", mean(brier_raw)-1.96*brier_raw_se_pooled, "-", mean(brier_raw)+1.96*brier_raw_se_pooled))


# this risk score overestimates risk and needs to be recalibrated.

############################2 RECALIBRATION - BASELINE HAZARD UPDATE ################################################################

# initially we will assess calibration if we only update the baseline hazard
# this is done by fitting a cox proportional model with the linear predictor as the only variable as an offset.
# as this would only involve adjusting the overall event probability for this cohort, no resampling (internal validation) is needed
recal_mod <- cph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ stats::offset(ckdpc_50egfr_lin_predictor), 
                 data = cohort, x = TRUE, y = TRUE, surv = TRUE)

x <- summary(survfit(recal_mod),time=3)
bh_update_presegfr <- x$surv
bh_update_se <- x$std.err
# print updated baseline hazard
print(paste0("Baseline hazard ", bh_update_presegfr, ", 95% CI ", bh_update_presegfr-1.96*bh_update_se, "-", bh_update_presegfr+1.96*bh_update_se))

cohort <- cohort %>% mutate(
  ckdpc_50egfr_survival_cal_bh=bh_update_presegfr^exp(ckdpc_50egfr_lin_predictor-mean(ckdpc_50egfr_lin_predictor)), # baseline hazard ^ e ^ centred linear predictor
  ckdpc_50egfr_score_cal_bh=(1-ckdpc_50egfr_survival_cal_bh)*100
)

## Show observed (estimated frequency binned per risk decile) vs predicted (risk-score predicted)
cohort$risk_decile <- ntile(cohort$ckdpc_50egfr_score_cal_bh, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- cohort %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal_bh)/100)

predicted_all <- cohort %>%
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
                                      data=cohort[cohort$.imp == i & 
                                                             cohort$risk_decile == k &
                                                             cohort$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                     data=cohort[cohort$.imp == i & 
                                                            cohort$risk_decile == k &
                                                            cohort$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile, 
                                  data=cohort[cohort$.imp == i & 
                                                         cohort$risk_decile == k,]) %>%
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


dpp4isu_events <- cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort %>%
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
p_interim_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = seq(0, 5, 1))+
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
  coord_cartesian(xlim = c(0,5), ylim = c(0,5))


p_interim_bydeciles_dpp4isu

##brier score after baseline hazard updated
brier_recal_bh <- rep(NA, n.imp)
brier_recal_bh_se <- rep(NA, n.imp)
for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival_cal_bh)
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
  
  brier_recal_bh[i] <- score_recal_bh$Brier$score$Brier[2]
  brier_recal_bh_se[i] <- score_recal_bh$Brier$score$se[2]
  
  rm(temp)
  
}

brier_recal_bh_se_pooled <- sqrt(mean(brier_recal_bh_se^2) + (1+1/n.imp)*var(brier_recal_bh))

#pool and print brier score
print(paste0("Brier score for risk score after baseline hazard updated ", mean(brier_recal_bh), ", 95% CI ", mean(brier_recal_bh)-1.96*brier_recal_bh_se_pooled, "-", mean(brier_recal_bh)+1.96*brier_recal_bh_se_pooled))


############################3 RECALIBRATION - OVERALL SLOPE RECALIBRATION################################################################

#before recalibrating, the risk score systematically overestimated risk. However after recalibration, it underestimates risk in those at high risk.

bh_new <- rep(NA, n.imp)
se_bh_new <- rep(NA, n.imp)
cal_slope <- rep(NA, n.imp)
slope_optimism_presegfr <- rep(NA, n.imp)
var_slope <- rep(NA, n.imp)

for (i in 1:n.imp) {
  print(paste0("Calculations in imputation ", i))
  recal_mod2 <- cph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ (ckdpc_50egfr_lin_predictor), data = cohort[cohort$.imp == i,], x = TRUE, y = TRUE, surv = TRUE)
  
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

bh_recal <- mean(bh_new)
bh_recal_se <- sqrt(mean(se_bh_new)^2 + (1+1/n.imp)*var(bh_new))
# print baseline hazard with 95% CI
print(paste0("Baseline hazard ", bh_recal, ", 95% CI ", bh_recal-1.96*bh_recal_se, "-", bh_recal+1.96*bh_recal_se))

coef_recal <- mean(cal_slope)
coef_recal_se <- sqrt(mean(var_slope) + (1+1/n.imp)*var(cal_slope))
# print calibration slope with 95% CI
print(paste0("Calibration slope ", coef_recal, ", 95% CI ", coef_recal-1.96*coef_recal_se, "-", coef_recal+1.96*coef_recal_se))

cohort <- cohort %>% mutate(
  ckdpc_50egfr_lin_predictor_cal=coef_recal*ckdpc_50egfr_lin_predictor,
  ckdpc_50egfr_survival_cal=bh_recal^exp(ckdpc_50egfr_lin_predictor_cal-mean(ckdpc_50egfr_lin_predictor_cal)), # baseline hazard ^ e ^ centred linear predictor
  ckdpc_50egfr_score_cal=(1-ckdpc_50egfr_survival_cal)*100
)


## Plot
cohort$risk_decile <- ntile(cohort$ckdpc_50egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- cohort %>%
  group_by(risk_decile, studydrug2) %>%
  summarise(mean_ckd50_pred=mean(ckdpc_50egfr_score_cal)/100)

predicted_all <- cohort %>%
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
                                      data=cohort[cohort$.imp == i &
                                                             cohort$risk_decile == k &
                                                             cohort$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd50$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd50$std.error
    
    observed_sglt2i_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile,
                                     data=cohort[cohort$.imp == i &
                                                            cohort$risk_decile == k &
                                                            cohort$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd50$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd50$std.error
    
    observed_all_ckd50 <- survfit(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ risk_decile,
                                  data=cohort[cohort$.imp == i &
                                                         cohort$risk_decile == k,]) %>%
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


dpp4isu_events <- cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr50_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- cohort %>%
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
p_cal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd50_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100))+
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
  ggtitle("Calibrated risk score, by risk decile") +
  coord_cartesian(xlim = c(0,5), ylim = c(0,5))


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_calibrated_risk_score_calibration.tiff"), width=6, height=5.5, units = "in", res=800) 
p_cal_bydeciles_dpp4isu
dev.off()

##brier score after overall calibration slope applied
brier_recal <- rep(NA, n.imp)
brier_recal_se <- rep(NA, n.imp)

for (i in 1:n.imp) {
  print(paste("Imputation ", i))
  
  temp <- cohort %>% filter(.imp == i) %>% select(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar, ckdpc_50egfr_survival_cal)
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
  
  brier_recal[i] <- score_recal$Brier$score$Brier[2]
  brier_recal_se[i] <- score_recal$Brier$score$se[2]
  
  rm(temp)
  
}
brier_recal_se_pooled <- sqrt(mean(brier_recal_se^2) + (1+1/n.imp)*var(brier_recal))

#pool and print brier score
print(paste0("Brier score for risk score after calibration slope applied ", mean(brier_recal), ", 95% CI ", mean(brier_recal)-1.96*brier_recal_se_pooled, "-", mean(brier_recal)+1.96*brier_recal_se_pooled))





############################4 STORE RECALIBRATED SCORES################################################################

#### model details summed up:
auc(ROC_raw)
ci.auc(ROC_raw)
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
print(paste0("Brier score for raw risk score ", mean(brier_raw), ", 95% CI ", mean(brier_raw)-1.96*brier_raw_se_pooled, "-", mean(brier_raw)+1.96*brier_raw_se_pooled))
print(paste0("Brier score for risk score after baseline hazard updated ", mean(brier_recal_bh), ", 95% CI ", mean(brier_recal_bh)-1.96*brier_recal_bh_se_pooled, "-", mean(brier_recal_bh)+1.96*brier_recal_bh_se_pooled))
print(paste0("Baseline hazard ", bh_update_presegfr, ", 95% CI ", bh_update_presegfr-1.96*bh_update_se, "-", bh_update_presegfr+1.96*bh_update_se))
print(paste0("Brier score for risk score after calibration slope applied ", mean(brier_recal), ", 95% CI ", mean(brier_recal)-1.96*brier_recal_se_pooled, "-", mean(brier_recal)+1.96*brier_recal_se_pooled))
print(paste0("Baseline hazard ", bh_recal, ", 95% CI ", bh_recal-1.96*bh_recal_se, "-", bh_recal+1.96*bh_recal_se))
print(paste0("Calibration slope ", coef_recal, ", 95% CI ", coef_recal-1.96*coef_recal_se, "-", coef_recal+1.96*coef_recal_se))
print(paste0("Slope optimism ", mean(slope_optimism_presegfr)))

# save dataset with calibrated risk score so this can be used in the subsequent scripts
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated.Rda"))

############################5 HR BY RISK SCORE################################################################
## check whether there is evidence of treatment heterogeneity by baseline risk (figure 1B)

# fit model using interaction term of treatment with risk score, 
# modelled with restricted cubic splines [rcs()] with 3-5 knots

ddist <- cohort %>% datadist()
options(datadist='ddist')


# Define the range of knots to test
k_range <- 3:5

# Initialize empty vectors to store results
aic_values <- numeric(length(k_range))
bic_values <- numeric(length(k_range))

# Loop over each value of k
for (i in seq_along(k_range)) {
  k <- k_range[i]
  
  # Fit the model with k knots
  model <- cph(
    as.formula(paste0(
      "Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ studydrug2*rcs(ckdpc_50egfr_score_cal,", k, ") + ",
      paste(setdiff(covariates, "ckdpc_50egfr_score_cal"), collapse=" + ") # need to remove risk score from covariate list as already specified in interaction term
    )),
    data = cohort %>% filter(.imp == n.imp), x = TRUE, y = TRUE
  )
  
  # Store the AIC and BIC values
  aic_values[i] <- AIC(model)
  bic_values[i] <- BIC(model)
}

# Find the optimal k based on minimum AIC or BIC
optimal_k_aic <- k_range[which.min(aic_values)]
optimal_k_bic <- k_range[which.min(bic_values)]

# Fit the final model using the optimal number of knots based on AIC
final_model <- cph(
  as.formula(paste0(
    "Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ studydrug2*rcs(ckdpc_50egfr_score_cal,", optimal_k_bic, ") + ",
    paste(setdiff(covariates, "ckdpc_50egfr_score_cal"), collapse=" + ")
  )),
  data = cohort %>% filter(.imp == n.imp), x = TRUE, y = TRUE
)

# Print optimal k values
cat("Optimal number of knots based on AIC:", optimal_k_aic, "\n")
cat("Optimal number of knots based on BIC:", optimal_k_bic, "\n")

anova(final_model)
anova(final_model)[2,3] # p value for non-linear interaction term

# create data frame with range of scores by study drug
contrast_spline <- contrast(final_model, 
                            list(studydrug2 = "SGLT2i", ckdpc_50egfr_score_cal = seq(quantile(cohort$ckdpc_50egfr_score_cal, .01, na.rm=TRUE), quantile(cohort$ckdpc_50egfr_score_cal, .99, na.rm=TRUE), by=0.05)), 
                            list(studydrug2 = "DPP4i/SU", ckdpc_50egfr_score_cal = seq(quantile(cohort$ckdpc_50egfr_score_cal, .01, na.rm=TRUE), quantile(cohort$ckdpc_50egfr_score_cal, .99, na.rm=TRUE), by=0.05))
)

contrast_spline_df <- as.data.frame(contrast_spline[c('ckdpc_50egfr_score_cal','Contrast','Lower','Upper')])
# plot
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_HR_by_ckd_egfr50_risk.tiff"), width=10, height=4, units = "in", res=800) 
ggplot(data=contrast_spline_df, aes(x=ckdpc_50egfr_score_cal, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=ckdpc_50egfr_score_cal, y=exp(Contrast)), size=1) +
  xlab(expression(paste("Predicted 3-year risk of kidney disease progression"))) +
  ylab("Hazard ratio") +
  coord_trans(y = "log10") +
  scale_x_continuous(breaks = seq(0,20,.5)) +
  scale_y_continuous(breaks = c(seq(0, 0.8, 0.1), seq(0.8, 1.6, 0.2))) +
  geom_ribbon(data=contrast_spline_df, aes(x=ckdpc_50egfr_score_cal, ymin=exp(Lower), ymax=exp(Upper)), alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  geom_hline(aes(yintercept = 0.57, linetype = "hr", size="hr"), color="#D55E00")  +
  geom_hline(aes(yintercept = 0.67, linetype = "hr_95", size="hr_95"), color="#D55E00")  +
  geom_hline(aes(yintercept = 0.49, linetype = "hr_95", size="hr_95"), color="#D55E00")  +
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="bottom",
        legend.title = element_text(size=14, face = "italic"),
        legend.text = element_text(face="italic")) +
  scale_linetype_manual(values = c(hr = "twodash", hr_95 = "twodash"), labels = c(hr = "Overall hazard ratio", hr_95 = "95% CI"), name="") +
  scale_size_manual(values = c(hr = 1, hr_95 = 0.5), labels = c(hr = "Overall hazard ratio", hr_95 = "95% CI"), name="")
dev.off()

options(datadist = NULL)
## there is no significant treatment heterogeneity by baseline risk score