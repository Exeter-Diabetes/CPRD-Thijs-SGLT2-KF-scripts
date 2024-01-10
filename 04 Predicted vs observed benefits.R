
## In this script our aim is to compare risk-score predicted SGLT2i benefit with (pseudo-)observed in our data.
# we will estimate the (pseudo-)observed or counterfactual absolute risk reduction (ARR or "benefit") with SGLT2i,
# and compare this to the risk-score predicted SGLT2i benefit (risk-score predicted baseline risk * overal relative risk [HR]).
# we will also make graphs as display items.

## Contents
# 0 setup
# 1 ARR, NNT, and calibration plots of observed vs predicted ARR with SGLT2i
# 2 survival curves
# 3 (in progress) adjusted curves


############################0 SETUP################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(rms)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())

set.seed(123)
n.imp <- 10
n.quantiles <- 10

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2023-12-06_t2d_ckdpc_recalibrated.Rda")

############################1 COMPARE PREDICTED/OBSERVED ARR AND NNT################################################################

# 1 compare ARR from SGLT2i vs SU as predicted by risk scores (multiplied by relative risk) versus (pseudo-)observed in data
# hazard ratio from trials meta-analysis: Lancet. 2022 Nov 19; 400(10365): 1788–1801

trial_hr_kf_sglt2 <- 0.62

# Add HR from trials
noncal_cohort <- noncal_cohort %>%
  mutate(ckdpc_egfr60_cal=(100-ckdpc_egfr60_confirmed_score_cal)/100,
         ckdpc_egfr60_cal_sglt2=ckdpc_egfr60_cal^trial_hr_kf_sglt2,
         ckdpc_egfr60_5yrisk_cal_sglt2=100-(ckdpc_egfr60_cal_sglt2*100))
         
noncal_cohort <- noncal_cohort %>%
  mutate(ckdpc_40egfr_cal=(100-ckdpc_40egfr_score_cal)/100,
         ckdpc_40egfr_cal_sglt2=ckdpc_40egfr_cal^trial_hr_kf_sglt2,
         ckdpc_40egfr_3yrisk_cal_sglt2=100-(ckdpc_40egfr_cal_sglt2*100))


# Distribution of predicted benefits
noncal_cohort <- noncal_cohort %>%
  mutate(ckdpc_egfr60_sglt2_benefit=ckdpc_egfr60_cal_sglt2 - ckdpc_egfr60_cal)

noncal_cohort <- noncal_cohort %>%
         mutate(ckdpc_40egfr_sglt2_benefit=ckdpc_40egfr_cal_sglt2 - ckdpc_40egfr_cal)

describe(noncal_cohort$ckdpc_egfr60_sglt2_benefit)
#mean=1.9%

describe(noncal_cohort$ckdpc_40egfr_sglt2_benefit)
#mean=0.3%



# Compare predicted and actual benefits by decile (predicted = mean of decile, as previous)

## Unadjusted

## eGFR <60

noncal_cohort$sglt2_ckd60_benefit_decile <- ntile(noncal_cohort$ckdpc_egfr60_sglt2_benefit, n.quantiles)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_ckd60_benefit_decile=as.factor(sglt2_ckd60_benefit_decile))

pred <- noncal_cohort %>%
  group_by(sglt2_ckd60_benefit_decile) %>%
  summarise(mean_ckdpc_egfr60_sglt2_benefit=mean(ckdpc_egfr60_sglt2_benefit, na.rm=T))


noncal_cohort <- noncal_cohort %>% mutate(studydrug=as.vector(studydrug))
# ddist <- datadist(noncal_cohort)


# to get the pooled benefits, we need to pool the imputations for every single decile per group
# first we will create a list of length n.imp with the benefit per decile per group for each imputation in a separate matrix
survival_est <- list()

for (i in 1:n.imp) {
model <- cph(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ studydrug*sglt2_ckd60_benefit_decile, 
             data=noncal_cohort[noncal_cohort$.imp == i,], x=TRUE, y=TRUE, surv=TRUE)


survival_est[[i]] <- survest(model, 
                        newdata=expand.grid(studydrug=c("SU", "DPP4", "SGLT2"), 
                                            sglt2_ckd60_benefit_decile=c(1:n.quantiles)), times=5)
}

# now we will create empty vectors for the pooled benefits per decile

surv <- rep(NA, n.quantiles*nlevels(as.factor(noncal_cohort$studydrug)))
se <- rep(NA, n.quantiles*nlevels(as.factor(noncal_cohort$studydrug)))

# for every decile per group we will extract the n.imp imputed values and pool these, then store these in the vectors created above

for (k in 1:(n.quantiles*nlevels(as.factor(noncal_cohort$studydrug)))) {

  #create empty vectors to store the imputed values for every individual decile 
  surv_temp=paste0("surv_", k)
  se_temp=paste0("surv_se_", k)
  
  assign(surv_temp, rep(NA, n.imp))
  assign(se_temp, rep(NA, n.imp))

  #store the imputed values in this vector
  for (i in 1:n.imp) {
    
    #surv_temp[i] <- survival_est[[i]]$surv[k]
    eval(str2lang(paste0(surv_temp, "[", i, "]", " <- survival_est", "[[", i, "]]$surv[", k, "]")))
    
    #se_temp[i] <- survival_est[[i]]$std.error[k]
    eval(str2lang(paste0(se_temp, "[", i, "]", " <- survival_est", "[[", i, "]]$std.err[", k, "]")))

    }
  # pool the imputed values and store the pooled value in the main vector
  surv[k] <- mean(eval(parse(text = surv_temp)))
  w <- mean(eval(parse(text = se_temp))^2)
  b <- var(eval(parse(text = surv_temp)))
  t.var <- w + (1+1/n.imp)*b
  se[k] <- sqrt(t.var)
}

obs <- data.frame(surv=surv, studydrug=rep(c("SU","DPP4","SGLT2"),n.quantiles), 
                  sglt2_ckd60_benefit_decile=rep(1:n.quantiles, rep_len(nlevels(as.factor(noncal_cohort$studydrug)), n.quantiles)))

obs <- obs %>%
  pivot_wider(id_cols=sglt2_ckd60_benefit_decile, values_from=surv, names_from=studydrug) %>%
  mutate(surv_diff=SGLT2-SU) %>% 
  select(sglt2_ckd60_benefit_decile, surv_diff)


se <- data.frame(se=se, studydrug=rep(c("SU","DPP4","SGLT2"),n.quantiles), 
                 sglt2_ckd60_benefit_decile=rep(1:n.quantiles, rep_len(nlevels(as.factor(noncal_cohort$studydrug)), n.quantiles)))

se <- se %>%
  pivot_wider(id_cols=sglt2_ckd60_benefit_decile, values_from=se, names_from=studydrug) %>%
  mutate(se=sqrt((SGLT2^2)+(SU^2))) %>%
  select(sglt2_ckd60_benefit_decile, se)

obs <- obs %>%
  inner_join(se, by="sglt2_ckd60_benefit_decile") %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

su_events <- noncal_cohort %>%
  filter(studydrug=="SU" & ckd_345_censvar==1) %>%
  group_by(sglt2_ckd60_benefit_decile) %>%
  summarise(SU=round(n()/n.imp,0))

su_events <- su_events[!is.na(su_events$sglt2_ckd60_benefit_decile),]

sglt2_events <- noncal_cohort %>%
  filter(studydrug=="SGLT2" & ckd_345_censvar==1) %>%
  group_by(sglt2_ckd60_benefit_decile) %>%
  summarise(SGLT2=round(n()/n.imp,0))

sglt2_events <- sglt2_events[!is.na(sglt2_events$sglt2_ckd60_benefit_decile),]


events_table <- data.frame(t(su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="sglt2_ckd60_benefit_decile")


obs <- obs %>% mutate(sglt2_ckd60_benefit_decile=as.factor(sglt2_ckd60_benefit_decile))
obs_v_pred <- pred %>% inner_join(obs, by="sglt2_ckd60_benefit_decile")


empty_tick <- obs_v_pred %>%
  filter(sglt2_ckd60_benefit_decile==1) %>%
  mutate(mean_predicted_benefit=NA, surv_diff=NA, lower_ci=NA, upper_ci=NA, sglt2_ckd60_benefit_decile=as.factor(0))

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckdpc_egfr60_sglt2_benefit*100)) + 
  geom_point(aes(y = surv_diff*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1) +
  theme_bw() +
  geom_abline(slope=1, intercept=0, linetype=2) +
  xlab("Predicted ARR with SGLT2i") + ylab("Observed ARR with SGLT2i") +
  ggtitle("Absolute risk reduction (ARR) of new CKD binned by risk score decile\n(SGLT2-inhibitors compared with sulfonylureas)") +
    scale_x_continuous(breaks=c(seq(0,10,by=1)), limits=c(0,10)) +
  scale_y_continuous(breaks=c(seq(-1,20,by=1)), limits=c(-1,20)) +
  theme(axis.text   = element_text(size=rel(1.5)),
        axis.title  = element_text(size=rel(1.5)),
        plot.title  = element_text(size=rel(1.5), face="bold")) +
  coord_cartesian(ylim=c(0,15))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

#NNT
obs_v_pred <- obs_v_pred %>%
  mutate(nnt_predicted=1/(mean_ckdpc_egfr60_sglt2_benefit),
         nnt_observed=1/(surv_diff))
obs_v_pred


## 40% eGFR decline/ESKD

noncal_cohort$sglt2_egfr40_benefit_decile <- ntile(noncal_cohort$ckdpc_40egfr_sglt2_benefit, n.quantiles)
noncal_cohort <- noncal_cohort %>% mutate(sglt2_egfr40_benefit_decile=as.factor(sglt2_egfr40_benefit_decile))

pred <- noncal_cohort %>%
  group_by(sglt2_egfr40_benefit_decile) %>%
  summarise(mean_ckdpc_40egfr_predicted_benefit=mean(ckdpc_40egfr_sglt2_benefit, na.rm=T))


noncal_cohort <- noncal_cohort %>% mutate(studydrug=as.vector(studydrug))
# ddist <- datadist(noncal_cohort)

# to get the pooled benefits, we need to pool the imputations for every single decile per group
# first we will create a list of length n.imp with the benefit per decile per group for each imputation in a separate matrix
survival_est <- list()

for (i in 1:n.imp) {
  model <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug*sglt2_egfr40_benefit_decile, 
               data=noncal_cohort[noncal_cohort$.imp == i,], x=TRUE, y=TRUE, surv=TRUE)
  
  
  survival_est[[i]] <- survest(model, 
                               newdata=expand.grid(studydrug=c("SU", "DPP4", "SGLT2"), 
                                                   sglt2_egfr40_benefit_decile=c(1:n.quantiles)), times=5)
}

# now we will create empty vectors for the pooled benefits per decile

surv <- rep(NA, n.quantiles*nlevels(as.factor(noncal_cohort$studydrug)))
se <- rep(NA, n.quantiles*nlevels(as.factor(noncal_cohort$studydrug)))

# for every decile per group we will extract the n.imp imputed values and pool these, then store these in the vectors created above

for (k in 1:(n.quantiles*nlevels(as.factor(noncal_cohort$studydrug)))) {
  
  #create empty vectors to store the imputed values for every individual decile 
  surv_temp=paste0("surv_", k)
  se_temp=paste0("surv_se_", k)
  
  assign(surv_temp, rep(NA, n.imp))
  assign(se_temp, rep(NA, n.imp))
  
  #store the imputed values in this vector
  for (i in 1:n.imp) {
    
    #surv_temp[i] <- survival_est[[i]]$surv[k]
    eval(str2lang(paste0(surv_temp, "[", i, "]", " <- survival_est", "[[", i, "]]$surv[", k, "]")))
    
    #se_temp[i] <- survival_est[[i]]$std.error[k]
    eval(str2lang(paste0(se_temp, "[", i, "]", " <- survival_est", "[[", i, "]]$std.err[", k, "]")))
    
  }
  # pool the imputed values and store the pooled value in the main vector
  surv[k] <- mean(eval(parse(text = surv_temp)))
  w <- mean(eval(parse(text = se_temp))^2)
  b <- var(eval(parse(text = surv_temp)))
  t.var <- w + (1+1/n.imp)*b
  se[k] <- sqrt(t.var)
}

obs <- data.frame(surv=surv, studydrug=rep(c("SU","DPP4","SGLT2"),n.quantiles), 
                  sglt2_egfr40_benefit_decile=rep(1:n.quantiles, rep_len(nlevels(as.factor(noncal_cohort$studydrug)), n.quantiles)))

obs <- obs %>%
  pivot_wider(id_cols=sglt2_egfr40_benefit_decile, values_from=surv, names_from=studydrug) %>%
  mutate(surv_diff=SGLT2-SU) %>%
  select(sglt2_egfr40_benefit_decile, surv_diff)

se <- data.frame(se=se, studydrug=rep(c("SU","DPP4","SGLT2"),n.quantiles), 
                 sglt2_egfr40_benefit_decile=rep(1:n.quantiles, rep_len(nlevels(as.factor(noncal_cohort$studydrug)), n.quantiles)))

se <- se %>%
  pivot_wider(id_cols=sglt2_egfr40_benefit_decile, values_from=se, names_from=studydrug) %>%
  mutate(se=sqrt((SGLT2^2)+(SU^2))) %>%
  select(sglt2_egfr40_benefit_decile, se)

obs <- obs %>%
  inner_join(se, by="sglt2_egfr40_benefit_decile") %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

su_events <- noncal_cohort %>%
  filter(studydrug=="SU" & ckd_egfr40_censvar==1) %>%
  group_by(sglt2_egfr40_benefit_decile) %>%
  summarise(SU=round(n()/n.imp,0))

sglt2_events <- noncal_cohort %>%
  filter(studydrug=="SGLT2" & ckd_egfr40_censvar==1) %>%
  group_by(sglt2_egfr40_benefit_decile) %>%
  summarise(SGLT2=round(n()/n.imp,0))

events_table <- data.frame(t(su_events %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="sglt2_egfr40_benefit_decile")

obs <- obs %>% mutate(sglt2_egfr40_benefit_decile=as.factor(sglt2_egfr40_benefit_decile))
obs_v_pred <- pred %>% inner_join(obs, by="sglt2_egfr40_benefit_decile")


empty_tick <- obs_v_pred %>%
  filter(sglt2_egfr40_benefit_decile==1) %>%
  mutate(mean_predicted_benefit=NA, surv_diff=NA, lower_ci=NA, upper_ci=NA, sglt2_egfr40_benefit_decile=as.factor(0))

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckdpc_40egfr_predicted_benefit*100)) + 
  geom_point(aes(y = surv_diff*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),alpha=1) +
  theme_bw() +
  geom_abline(slope=1, intercept=0, linetype=2) +
  xlab("Predicted ARR with SGLT2i") + ylab("Observed ARR with SGLT2i") +
  ggtitle("Absolute risk reduction (ARR) of 40% decline in eGFR/ESKD\nbinned by risk score decile (SGLT2-inhibitors compared with sulfonylureas)") +
    scale_x_continuous(breaks=c(seq(0,6,by=.5)), limits=c(0,2.5)) +
  scale_y_continuous(breaks=c(seq(0,8,by=.5)), limits=c(-2,5)) +
  theme(axis.text   = element_text(size=rel(1.5)),
        axis.title  = element_text(size=rel(1.5)),
        plot.title  = element_text(size=rel(1.5), face="bold")) +
  coord_cartesian(ylim=c(-.5,2.1), xlim=c(0,1.1))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

#NNT
obs_v_pred <- obs_v_pred %>%
  mutate(nnt_predicted=1/(mean_ckdpc_40egfr_predicted_benefit),
         nnt_observed=1/(surv_diff))
obs_v_pred

############################2 SURVIVAL CURVES################################################################


## 2 survival curves

## survival curves per study drug with panel for different risk deciles.

#split each risk score into 3 groups (visually decided based on graphs above)
noncal_cohort <- noncal_cohort %>% mutate(
  ckd60_risk_group = ifelse(ckd60_risk_decile %in% c(9,10), "Risk decile 9-10", ifelse(ckd60_risk_decile %in% c(6,7,8), "Risk decile 6-8", "Risk decile 1-5")),
  ckd40_risk_group = ifelse(ckd40_risk_decile == 10, "Risk decile 10", ifelse(ckd60_risk_decile %in% c(4,5,6,7,8,9), "Risk decile 4-9", "Risk decile 1-3")),
  ckd40_risk_group = factor(ckd40_risk_group, levels = c("Risk decile 1-3", "Risk decile 4-9", "Risk decile 10"))
)


surv_ckd345 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ studydrug,
                     data = noncal_cohort,
                     conf.type = "log", conf.int = 0.95)

ggsurvplot(
  fit = surv_ckd345,
  fun = "event",
  data = noncal_cohort[noncal_cohort$.imp == 1,],
  facet.by = "ckd60_risk_group",
  palette = c("darkgoldenrod1", "darkmagenta", "darkolivegreen2"),
  color = "studydrug",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors", "SGLT2 inhibitors", "Sulfonylureas"),
  panel.labs = list(ckd_risk_group = c("Risk decile 1-5", "Risk decile 6-8", "Risk decile 9-10")),
  short.panel.labs = T,
  legend = c(0.12, 0.88),
  break.y.by = 0.05,
  censor.size = 1,
  surv.scale = "percent",
  xlab = "Years",
  ylab = "Incidence of new CKD (eGFR <60 mL/min)"
    )

surv_egfr40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug,
                     data = noncal_cohort,
                     conf.type = "log", conf.int = 0.95)

ggsurvplot(
  fit = surv_egfr40,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1,],
  facet.by = "ckd40_risk_group",
  palette = c("darkgoldenrod1", "darkmagenta", "darkolivegreen2"),
  color = "studydrug",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors", "SGLT2 inhibitors", "Sulfonylureas"),
  panel.labs = list(ckd_risk_group = c("Risk decile 1-3", "Risk decile 4-9", "Risk decile 10")),
  short.panel.labs = T,
  legend = c(0.12, 0.88),
  break.y.by = 0.01,
  censor.size = 1,
  surv.scale = "percent",
  xlab = "Years",
  ylab = "Incidence of 40% decline in eGFR or ESKD"
)

# alternative is to use plot.survfit below but this doesn't allow for shaded confidence intervals unfortunately

par(mfrow=c(1,3))
for (i in 1:3) {
  surv_temp <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ studydrug,
                       data = noncal_cohort[noncal_cohort$ckd60_risk_group == i,],
                   conf.type = "log", conf.int = 0.95)
  plot(surv_temp, xlim = c(0,5), xscale = 1, ylim = c(0, 0.25),
       xlab = ifelse(i == 2, "Years", ""),
       ylab = ifelse(i == 1, "CKD incidence", ""),
       col = c(1:3), lty = c(2:1), lwd = 2,
       cumhaz = T)
  legend(x=0, y=ifelse(i == 1, 0.05, 100), c("DPP4i", "SGLT2i", "SU"),
         cex = 1, col = c(1:3), lty = 2:1, bty = 'n', lwd = 2)
  title(paste0("Risk decile ", ifelse(
    i == 1, "1-5", ifelse(
      i == 2, "6-8", "9-10"
    )
  )))
}

par(mfrow=c(1,3))
for (i in 1:3) {
  surv_temp <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug,
                       data = noncal_cohort[noncal_cohort$ckd40_risk_group == i,],
                       conf.type = "log", conf.int = 0.95)
  plot(surv_temp, xlim = c(0,3), xscale = 1, ylim = c(0, 0.05),
       xlab = ifelse(i == 2, "Years", ""),
       ylab = ifelse(i == 1, "40% decline in eGFR/ESKD", ""),
       col = c(1:3), lty = c(2:1), lwd = 2,
       cumhaz = T)
  legend(x=0, y=ifelse(i == 1, 0.05, 100), c("DPP4i", "SGLT2i", "SU"),
           cex = 1, col = c(1:3), lty = 2:1, bty = 'n', lwd = 2)
  title(paste0("Risk decile ", ifelse(
    i == 1, "1-3", ifelse(
      i == 2, "4-9", "10"
    )
  )))
}

############################3 ADJUSTED CURVES################################################################
