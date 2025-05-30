## In this script our aim is to predict baseline risk using the CKD-PC risk score.
## we will be looking at the risk score for 50% decline in eGFR/ESKD (50egfr).
## which we will recalibrate in the CPRD cohort for endpoint 50% decline in eGFR/ESKD.

## Contents
# 0 setup
# 1 calibration of risk score `as is`
# 2 baseline hazard update (redundant)
# 3 calibration slope
# 4 store dataset with recalibrated risk score as required

############################0 SETUP################################################################

# 0 Setup
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/")
source("00 Setup.R")

# load data
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
load(paste0(today, "_t2d_ckdpc_imputed_data_withweights.Rda"))

cols <- cols[names(cols) %in% cohort$studydrug2]
cols <- cols[order(names(cols))]

# select imputed data only (ie. remove non-imputed data)
cohort <- cohort[cohort$.imp > 0,]

## remove double overlapping entries for DPP4i and SU that overlap (take one only)
cohort <- cohort %>% group_by(.imp, patid) %>% filter(
  !duplicated(studydrug2)
) %>% ungroup()

# check number of subjects
table(cohort$studydrug2)

############################1 CALIBRATION OF RISK SCORE `AS IS`################################################################



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
  xlab("Predicted 3-year risk of kidney disease progression (%)") + ylab("Observed risk (%)")+
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
  #ggtitle("Predicted risk of kidney disease progression") +
  coord_cartesian(xlim = c(0,5), ylim = c(0,5))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_uncalibrated_risk_score_calibration.tiff"), width=6, height=5.5, units = "in", res=800) 
print(p_uncal_bydeciles_dpp4isu)
dev.off()

## C-stat
cohort <- cohort %>%
  mutate(ckdpc_50egfr_survival=(100-ckdpc_50egfr_score)/100)

raw_mod <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ ckdpc_50egfr_survival, data=cohort[!cohort$studydrug2 == "SGLT2i",], method="breslow")
cstat_est <- summary(raw_mod)$concordance[1]
cstat_est_ll <- summary(raw_mod)$concordance[1]-(1.96*summary(raw_mod)$concordance[2])
cstat_est_ul <- summary(raw_mod)$concordance[1]+(1.96*summary(raw_mod)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))

## compare this to c-stat for albuminuria as a continuous variable:
mod_alb <- coxph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ uacr, data=cohort[!cohort$studydrug2 == "SGLT2i",], method="breslow")
cstat_est_alb <- summary(mod_alb)$concordance[1]
cstat_est_alb_ll <- summary(mod_alb)$concordance[1]-(1.96*summary(mod_alb)$concordance[2])
cstat_est_alb_ul <- summary(mod_alb)$concordance[1]+(1.96*summary(mod_alb)$concordance[2])
paste0("C statistic for albuminuria (continuous variable): ", round(cstat_est_alb, 4), ", 95% CI ", round(cstat_est_alb_ll, 4), "-", round(cstat_est_alb_ul,4))


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
        slice(1) %>% # this selects the first patients in dataset
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


# as can be seen from the figure above, the uncalibrated risk score (`as is`) predicts the outcome well.
# further down, we will calculate a formal calibration slope. 
# the next bit of code would be the first step if needing to recalibrate the risk score, ie. updating the baseline hazard.
# as this is not required in this case, this code is redundant.

############################2 CALIBRATION SLOPE################################################################

# in order to calculate an overall calibration slope, we are going to fit a cox model with the linear predictor as the only variable.
# this would be the same step we would take if needing to recalibrate the model beyond just updating the baseline hazard, as this would only require 
# recalculation of 2 variables (baseline hazard + overall calibration slope)
# this is why bootstrap internal validation is being used, to correct for optimism (this is a more efficient alternative to k-fold cross-validation or split-sample validation)
# in this case we will not need to use the recalibrated risk score, but rather will use the calibration slope estimated here
# to formally evidence calibration accuracy of the risk score `as is`.

# create empty vectors to store values of each calculation in multiple imputations:
bh_new <- rep(NA, n.imp)
se_bh_new <- rep(NA, n.imp)
cal_slope <- rep(NA, n.imp)
slope_optimism_presegfr <- rep(NA, n.imp)
var_slope <- rep(NA, n.imp)

options(datadist=NULL)

for (i in 1:n.imp) {
  print(paste0("Calculations in imputation ", i))
  # fit model with linear predictor as only variable
  recal_mod2 <- cph(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ (ckdpc_50egfr_lin_predictor), 
                    data = cohort %>% filter(.imp == i & studydrug2 == "DPP4i/SU"), x = TRUE, y = TRUE, surv = TRUE)
  
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

# pool baseline hazard results
bh_recal <- mean(bh_new)
bh_recal_se <- sqrt(mean(se_bh_new)^2 + (1+1/n.imp)*var(bh_new))
# print baseline hazard with 95% CI
print(paste0("Baseline hazard ", bh_recal, ", 95% CI ", bh_recal-1.96*bh_recal_se, "-", bh_recal+1.96*bh_recal_se))

# pool calibration slope results
coef_recal <- mean(cal_slope)
coef_recal_se <- sqrt(mean(var_slope) + (1+1/n.imp)*var(cal_slope))
# print calibration slope with 95% CI
print(paste0("Calibration slope ", coef_recal, ", 95% CI ", coef_recal-1.96*coef_recal_se, "-", coef_recal+1.96*coef_recal_se))


############################3 MODEL DETAILS SUMMED UP################################################################

#### model details summed up:
auc(ROC_raw)
ci.auc(ROC_raw)
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
print(paste0("Brier score for risk score `as is`", mean(brier_raw), ", 95% CI ", mean(brier_raw)-1.96*brier_raw_se_pooled, "-", mean(brier_raw)+1.96*brier_raw_se_pooled))
print(paste0("Baseline hazard ", bh_recal, ", 95% CI ", bh_recal-1.96*bh_recal_se, "-", bh_recal+1.96*bh_recal_se))
print(paste0("Calibration slope ", coef_recal, ", 95% CI ", coef_recal-1.96*coef_recal_se, "-", coef_recal+1.96*coef_recal_se))
print(paste0("Slope optimism ", mean(slope_optimism_presegfr)))

options(datadist=NULL)