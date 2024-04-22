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

options(dplyr.dpp4isummarise.inform = FALSE)

# set seed
set.seed(123)

# number of imputations 
n.imp <- 10

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
load("2024-03-06_t2d_ckdpc_imputed_data_withweights_incl_egfr_below_60.Rda")

# select imputed data only (ie. remove non-imputed data)
cohort <- cohort[cohort$.imp > 0,]

# set reference level for variable studydrug
cohort$studydrug <- relevel(as.factor(cohort$studydrug), ref = "SU")

## remove double overlapping entries for DPP4i and SU that overlap (take one only)
cohort <- cohort %>% group_by(.imp, patid) %>% filter(
  !duplicated(studydrug2)
) %>% ungroup()

# check number of subjects
table(cohort$studydrug)
# SU  DPP4i SGLT2i 
# 389470 624700 559760   # 10 imputations therefore number of subjects per group appears 10 times larger

# select calibration cohort and non-calibration cohort 

# assign random 20% (SU/DPP4) as recalibration cohort and remove from main cohort
# we will do this in each imputation and combine 
cal_cohort <- cohort %>% filter(!studydrug2 == "SGLT2i") %>% group_by(.imp) %>% slice_sample(prop=0.2)

# select those not in the calibration cohort
noncal_cohort <- cohort %>%
  anti_join(cal_cohort, by=c("patid", "dstartdate", "studydrug2", ".imp"))

# check number of subjects
table(noncal_cohort$studydrug)
# SU  DPP4i SGLT2i 
# 308945 475755 499270 


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

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=ckd40_risk_decile, group=studydrug2, color=studydrug2)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-10,10)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated risk score vs 40% eGFR decline/ESKD (3 year)") +
  coord_cartesian(ylim = c(0,7.5))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p_ckd40_uncal_bydeciles_bydrug <- p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))
p_ckd40_uncal_bydeciles_bydrug

# plot with obs vs predicted with loess smoother

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  #geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed*100, group=studydrug2, color=studydrug2), position=dodge, shape=17, size=2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab("CKD-PC risk score (%)") + ylab("Estimated risk (%)")+
  scale_x_continuous(limits=c(0,6))+
  scale_y_continuous(limits=c(-5,10)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Predicted vs observed incidence of 40% eGFR decline/ESKD (3 year)")+
  coord_cartesian(ylim = c(0,3.5))

p_ckd40_uncal_bydrug <- p1 + geom_smooth(aes(y=observed*100, x=mean_ckd40_pred*100), colour = "black", data = obs_v_pred)
p_ckd40_uncal_bydrug

## first plot with all observed combined (not grouped by studydrug) vs predicted

obs_v_pred2 <- rbind(
  (predicted_all %>% mutate(
    lower_ci = NA,
    upper_ci = NA,
    risk_type = "predicted"
  ) %>%
    select(strata=ckd40_risk_decile, estimate=mean_ckd40_pred, lower_ci, upper_ci, risk_type)
  ),
  (observed_all %>% mutate(
    risk_type = "observed"
  ) %>% 
    relocate(strata, .before = observed) %>%
    relocate(risk_type, .after = upper_ci) %>%
    select(strata, estimate=observed, lower_ci, upper_ci, risk_type)
  )
)

p_ckd40_uncal_bydeciles <- ggplot(data=obs_v_pred2, aes(x=strata, y=estimate*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.1,size=.75) +
  geom_point(aes(fill = risk_type, size = risk_type), shape = 21, colour = "black", stroke = 1.5) +
  scale_fill_manual("",
                    breaks = c("predicted", "observed"),
                    labels = c("Uncalibrated risk score", "Incidence of 40% eGFR decline/ESKD (Kaplan-Meier estimate)"),
                    values = c("white", "cadetblue3")) +
  scale_size_manual(breaks = c("predicted", "observed"),
                    values = c(4,3),
                    guide = "none") +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +  
  theme(legend.position = c(0.05,0.95),
        legend.justification = c(0,1),
        legend.key.size = unit(1, "cm"),
        legend.text=element_text(size=rel(1))) +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated risk score vs incidence of 40% eGFR decline/ESKD (3-year)") +
  coord_cartesian(ylim = c(0,7.5))

p_ckd40_uncal_bydeciles

## FINAL PLOT
p_ckd40_uncal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Uncalibrated CKD-PC risk score (%)") + ylab("Estimated risk (%)")+
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
cohort <- cohort %>%
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score)/100)

surv_mod_ckd40 <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival, data=cohort, method="breslow")
cstat_est <- summary(surv_mod_ckd40)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd40)$concordance[1]-(1.96*summary(surv_mod_ckd40)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd40)$concordance[1]+(1.96*summary(surv_mod_ckd40)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7693, 95% CI 0.7643-0.7744
# brier_score <- brier(surv_mod_ckd40, times = 3)
# paste0("Brier score: ", round(brier_score$brier,4))

# similarly, this risk score also overestimates risk.
# both will need some recalibration.

############################2A RECALIBRATION - INTERCEPT ONLY ################################################################
# 40% decline in eGFR / CKD 5
# risk equation uses a logistic model. There are several ways to update these models.
# the most straightforward one is to adjust the intercept ("calibration in the large")
#see: Journal of Clinical Epidemiology 61 (2008) 76e86. Janssen KJM, Moons KGM, Kalkman CJ, Grobbee DE, Vergouwe Y. Updating methods improved the performance of a clinical prediction model in new patients


#original intercept:
EHRBiomarkr::ckdpc40EgfrRiskConstants$intercept %>% as.numeric()

#calculate correction factor for intercept = ln(observed outcome frequency/1-observed outcome frequency / predicted outcome frequency/1-predicted outcome frequency)
correction_factor <- rep(NA, n.imp)
for (i in 1:n.imp) {
  correction_factor[i] <- ((mean(cal_cohort[cal_cohort$.imp == i,]$ckd_egfr40_censvar) / 
                              (1-mean(cal_cohort[cal_cohort$.imp == i,]$ckd_egfr40_censvar))) /
                             (mean(cal_cohort[cal_cohort$.imp == i,]$ckdpc_40egfr_score*0.01) / 
                                (1-mean(cal_cohort[cal_cohort$.imp == i,]$ckdpc_40egfr_score*0.01))) ) %>% 
    log()
}

#new intercept = original intercept + correction factor
new_intercept <- EHRBiomarkr::ckdpc40EgfrRiskConstants$intercept %>% as.numeric() + mean(correction_factor)

## Recalculate scores in rest of cohort
noncal_cohort <- noncal_cohort %>%
  mutate(centred_40egfr_lin_predictor=ckdpc_40egfr_lin_predictor-mean(ckdpc_40egfr_lin_predictor)) %>%
  ungroup() %>%
  mutate(ckdpc_40egfr_score_cal=100 * (exp(centred_40egfr_lin_predictor + new_intercept)/(1 + exp(centred_40egfr_lin_predictor + new_intercept))))


## Show observed (estimated frequency binned per risk decile) vs predicted (risk-score predicted)
noncal_cohort$ckd40_risk_decile <- ntile(noncal_cohort$ckdpc_40egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- noncal_cohort %>%
  group_by(ckd40_risk_decile, studydrug2) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

predicted_all <- noncal_cohort %>%
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
                                 data=noncal_cohort[noncal_cohort$.imp == i & 
                                                      noncal_cohort$ckd40_risk_decile == k &
                                                      noncal_cohort$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd40$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd40$std.error
    
    observed_sglt2i_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile, 
                                    data=noncal_cohort[noncal_cohort$.imp == i & 
                                                         noncal_cohort$ckd40_risk_decile == k &
                                                         noncal_cohort$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd40$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd40$std.error
    
    observed_all_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile, 
                                  data=noncal_cohort[noncal_cohort$.imp == i & 
                                                       noncal_cohort$ckd40_risk_decile == k,]) %>%
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


dpp4isu_events <- noncal_cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- noncal_cohort %>%
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

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=ckd40_risk_decile, group=studydrug2, color=studydrug2)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Recalibrated risk score vs 40% eGFR decline/ESKD (3 year)\n(Calibration-in-the-large)") +
  coord_cartesian(ylim = c(0,7.5))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p_ckd40_interim_bydeciles_bydrug <- p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))
p_ckd40_interim_bydeciles_bydrug

# obs vs predicted (not decile on x-axis but predicted risk) with loess smoother

p1  <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  #geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed*100, group=studydrug2, color=studydrug2), position=dodge, shape=17, size=2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab("CKD-PC risk score (%)") + ylab("Estimated risk (%)")+
  scale_x_continuous(limits=c(0,3.5))+
  scale_y_continuous(limits=c(-1,7)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Predicted vs observed incidence of 40% eGFR decline/ESKD (3 year)") +
  coord_cartesian(ylim = c(0,3.5))

p_ckd40_interim_bydrug <- p1 + geom_smooth(aes(y=observed*100, x=mean_ckd40_pred*100), colour = "black", data = obs_v_pred)
p_ckd40_interim_bydrug

## plot with all observed combined vs predicted

obs_v_pred_all <- rbind(
  (predicted_all %>% mutate(
    lower_ci = NA,
    upper_ci = NA,
    risk_type = "predicted"
  ) %>%
    select(strata=ckd40_risk_decile, estimate=mean_ckd40_pred, lower_ci, upper_ci, risk_type)
  ),
  (observed_all %>% mutate(
    risk_type = "observed"
  ) %>% 
    relocate(strata, .before = observed) %>%
    relocate(risk_type, .after = upper_ci) %>%
    select(strata, estimate=observed, lower_ci, upper_ci, risk_type)
  )
)

p_ckd40_interim_bydeciles <- ggplot(data=obs_v_pred_all, aes(x=strata, y=estimate*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.1,size=.75) +
  geom_point(aes(fill = risk_type, size = risk_type), shape = 21, colour = "black", stroke = 1.5) +
  scale_fill_manual("",
                    breaks = c("predicted", "observed"),
                    labels = c("Recalibrated risk score", "Incidence of 40% eGFR decline/ESKD (Kaplan-Meier estimate)"),
                    values = c("white", "cadetblue3")) +
  scale_size_manual(breaks = c("predicted", "observed"),
                    values = c(4,3),
                    guide = "none") +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +  
  theme(legend.position = c(0.05,0.95),
        legend.justification = c(0,1),
        legend.key.size = unit(1, "cm"),
        legend.text=element_text(size=rel(1))) +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("CKD-PC risk score vs estimated risk of 40% decline in eGFR / ESKD", subtitle = "Intercept recalibration of CKD-PC risk score, binned by risk decile") +  
  coord_cartesian(ylim = c(0,7.5))

p_ckd40_interim_bydeciles

## FINAL PLOT
p_ckd40_interim_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Estimated risk (%)")+
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
  ggtitle("Calibration plot of CKD-PC risk score for 40% decline in eGFR / ESKD", subtitle = "Intercept recalibration of risk score, binned by risk decile") +
  coord_cartesian(xlim = c(0,7.5), ylim = c(0,7.5))


p_ckd40_interim_bydeciles_dpp4isu

## C-stat
noncal_cohort <- noncal_cohort %>%
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score_cal)/100)

surv_mod_ckd40 <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival, data=noncal_cohort, method="breslow")
cstat_est <- summary(surv_mod_ckd40)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd40)$concordance[1]-(1.96*summary(surv_mod_ckd40)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd40)$concordance[1]+(1.96*summary(surv_mod_ckd40)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7515, 95% CI 0.7454-0.7577
# brier_score <- brier(surv_mod_ckd40, times = 3)
# paste0("Brier score: ", round(brier_score$brier,4))

############################2B RECALIBRATION - LOGISTIC RECALIBRATION################################################################

#before recalibrating, the risk score systematically overestimated risk. However after recalibration, it underestimates risk in those at high risk.
#an alternative method is logistic calibration:
#a logistic regression model is fitted with the linear predictor as the only covariate in the updating set.
#The calibration slope is used to recalibrate the original regression coefficients (multiplied with the calibration slope).
#And we will add the new calibration intercept

cal_cohort <- cal_cohort %>%
  mutate(temp_lin_predictor=ckdpc_40egfr_lin_predictor + EHRBiomarkr::ckdpc40EgfrRiskConstants$intercept %>% as.numeric(),
         centred_lin_predictor=temp_lin_predictor-mean(temp_lin_predictor))

intercept_cf <- coef_cf <- rep(NA,n.imp)

for (i in 1:n.imp) {
  recal_mod <- glm(ckd_egfr40_censvar ~ centred_lin_predictor,
                   data = cal_cohort[cal_cohort$.imp == i,], family = "binomial"
  )
  
  new_intercept[i] <- recal_mod$coef[1]
  coef_cf[i] <- recal_mod$coef[2]
  
}

new_intercept <- mean(new_intercept)

## Recalculate scores in rest of cohort
noncal_cohort <- noncal_cohort %>%
  mutate(
    ckdpc_40egfr_score_cal=
      100 * 
      (exp(mean(coef_cf)*ckdpc_40egfr_lin_predictor + new_intercept))/(1 + exp(mean(coef_cf)*ckdpc_40egfr_lin_predictor + new_intercept))
  )

## Plot
noncal_cohort$ckd40_risk_decile <- ntile(noncal_cohort$ckdpc_40egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- noncal_cohort %>%
  group_by(ckd40_risk_decile, studydrug2) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

predicted_all <- noncal_cohort %>%
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
                                      data=noncal_cohort[noncal_cohort$.imp == i &
                                                           noncal_cohort$ckd40_risk_decile == k &
                                                           noncal_cohort$studydrug2=="DPP4i/SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4isu[k,i] <- observed_dpp4isu_ckd40$estimate
    SE.dpp4isu[k,i] <- observed_dpp4isu_ckd40$std.error
    
    observed_sglt2i_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile,
                                     data=noncal_cohort[noncal_cohort$.imp == i &
                                                          noncal_cohort$ckd40_risk_decile == k &
                                                          noncal_cohort$studydrug2=="SGLT2i",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2i[k,i] <- observed_sglt2i_ckd40$estimate
    SE.sglt2i[k,i] <- observed_sglt2i_ckd40$std.error
    
    observed_all_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckd40_risk_decile,
                                  data=noncal_cohort[noncal_cohort$.imp == i &
                                                       noncal_cohort$ckd40_risk_decile == k,]) %>%
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


dpp4isu_events <- noncal_cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_decile) %>%
  summarise(DPP4iSU=round(n()/n.imp, 0))

sglt2_events <- noncal_cohort %>%
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

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=ckd40_risk_decile, group=studydrug2, color=studydrug2)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Recalibrated risk score vs 40% eGFR decline/ESKD (3 year)\n(Logistic recalibration)") +
  coord_cartesian(ylim = c(0,7.5))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p_ckd40_cal_bydeciles_bydrug <- p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))
p_ckd40_cal_bydeciles_bydrug

## obs v predicted plot (not deciles but predicted risk on x-axis)
p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  #geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed*100, group=studydrug2, color=studydrug2), position=dodge, shape=17, size=2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab("CKD-PC risk score (%)") + ylab("Estimated risk (%)")+
  scale_x_continuous(limits=c(0,3.5))+
  scale_y_continuous(limits=c(-1,7)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Predicted vs observed incidence of 40% eGFR decline/ESKD (3 year)\n(Logistic recalibration)") +
  coord_cartesian(ylim = c(0,3.5))

p_ckd40_cal_bydrug <- p1 + geom_smooth(aes(y=observed*100, x=mean_ckd40_pred*100), colour = "black", data = obs_v_pred)
p_ckd40_cal_bydrug

## plot with all observed combined (like first plot but not grouped by drug) vs predicted


obs_v_pred2 <- rbind(
  (predicted_all %>% mutate(
    lower_ci = NA,
    upper_ci = NA,
    risk_type = "predicted"
  ) %>%
    select(strata=ckd40_risk_decile, estimate=mean_ckd40_pred, lower_ci, upper_ci, risk_type)
  ),
  (observed_all %>% mutate(
    risk_type = "observed"
  ) %>%
    relocate(strata, .before = observed) %>%
    relocate(risk_type, .after = upper_ci) %>%
    select(strata, estimate=observed, lower_ci, upper_ci, risk_type)
  )
)

p_ckd40_cal_bydeciles <- ggplot(data=obs_v_pred2, aes(x=strata, y=estimate*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.1,size=.75) +
  geom_point(aes(fill = risk_type, size = risk_type), shape = 21, colour = "black", stroke = 1.5) +
  scale_fill_manual("",
                    breaks = c("predicted", "observed"),
                    labels = c("Recalibrated risk score", "Incidence of 40% eGFR decline/ESKD (Kaplan-Meier estimate)"),
                    values = c("white", "cadetblue3")) +
  scale_size_manual(breaks = c("predicted", "observed"),
                    values = c(4,3),
                    guide = "none") +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.position = c(0.05,0.95),
        legend.justification = c(0,1),
        legend.key.size = unit(1, "cm"),
        legend.text=element_text(size=rel(1))) +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Recalibrated risk score vs incidence of 40% eGFR decline/ESKD (3-year)\n(Logistic recalibration)") +
  coord_cartesian(ylim = c(0,7.5))

p_ckd40_cal_bydeciles

## FINAL PLOT
p_ckd40_cal_bydeciles_dpp4isu <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_dpp4isu*100,ymin=lower_ci_dpp4isu*100, color=studydrug2),width=0.1,size=1) +
  geom_point(aes(y = observed_dpp4isu*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Recalibrated CKD-PC risk score (%)") + ylab("Estimated risk (%)")+
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
  ggtitle("Calibration plot of CKD-PC risk score for 40% decline in eGFR / ESKD", subtitle = "Logistic recalibration of risk score, binned by risk decile") +
  coord_cartesian(xlim = c(0,6), ylim = c(0,6))


p_ckd40_cal_bydeciles_dpp4isu

# with SGLT2i included (but not to show as these are overpredicted)
# p_ckd40_cal_bydeciles <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
#   geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100, color=studydrug2),width=0.1,size=1) +
#   geom_point(aes(y = observed*100, group=studydrug2, color=studydrug2), shape=18, size=3) +
#   geom_abline(intercept = 0, slope = 1, lty = 2) +
#   theme_bw() +
#   xlab("CKD-PC risk score (%)") + ylab("Estimated risk (%)")+
#   scale_x_continuous(limits=c(0,6))+
#   scale_y_continuous(limits=c(-1,7)) +
#   scale_colour_manual(values = cols) +
#   theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
#         axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
#         plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
#   theme(axis.text=element_text(size=rel(1.5)),
#         axis.title=element_text(size=rel(1.5)),
#         plot.title=element_text(hjust = 0.5),
#         plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
#         legend.position = c(0.2,0.8),
#         legend.title = element_text(colour = "white")) +
#   ggtitle("CKD-PC risk score vs estimated risk of 40% decline in eGFR / ESKD", subtitle = "Recalibrated score using logistic recalibration, binned by risk decile") +
#   coord_cartesian(xlim = c(0,6), ylim = c(0,6))
# p_ckd40_cal_bydeciles

## C-stat
noncal_cohort <- noncal_cohort %>%
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score_cal)/100)

surv_mod_ckd40 <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival, data=noncal_cohort, method="breslow")
cstat_est <- summary(surv_mod_ckd40)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd40)$concordance[1]-(1.96*summary(surv_mod_ckd40)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd40)$concordance[1]+(1.96*summary(surv_mod_ckd40)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7656, 95% CI 0.76-0.7713
# brier_score <- brier(surv_mod_ckd40, times = 3)
#  paste0("Brier score: ", round(brier_score$brier,4))

############################3 STORE RECALIBRATED SCORES################################################################
# save dataset with calibrated risk score so this can be used in the subsequent scripts
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(noncal_cohort, file=paste0(today, "_t2d_ckdpc_recalibrated_incl_egfr_below_60.Rda"))
