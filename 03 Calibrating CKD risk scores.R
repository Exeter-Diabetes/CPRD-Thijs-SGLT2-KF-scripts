## the question we try to answer here is whether we can predict baseline risk with the CKDPC risk scores
# we will be looking at 2 risk scores for: 
#  - developing CKD stage 3 or more (egfr60 or ckd345) 
#  - 40% decline in eGFR/ESKD (40egfr)

## 1 load data

## 2 uncalibrated risk scores

## 3 re-calibration


############################0 SETUP################################################################

# 0 Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(rms)

options(dplyr.summarise.inform = FALSE)

# set seed
set.seed(123)

# number of imputations 
n.imp <- 10

# number of quantiles (for risk scores later on)
n.quantiles <- 10
############################1 LOAD DATA################################################################

# 1 load data

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2023-11-08_t2d_ckdpc_imputed_data.Rda")

cohort <- temp[temp$.imp > 0,]

cohort$studydrug <- relevel(as.factor(cohort$studydrug), ref = "SU")

table(cohort[cohort$.imp>0,]$studydrug)
#    SU   DPP4  SGLT2 
#463810  459650  209640   # 10 imputations therefore number of subjects per group appears 10 times larger


ckd60_cohort <- cohort %>%
  filter(!is.na(ckdpc_egfr60_confirmed_score)) %>%
  select(-c(contains("40"), contains("risk_total")))

ckd40_cohort <- cohort %>%
  filter(!is.na(ckdpc_40egfr_score)) %>%
  select(-c(contains("60"), contains("ckd_345")))


############################2A UNCALIBRATED SCORES - NEW CKD################################################################

# 2 How well do uncalibrated risk scores predict incidence?


# CKD <60
# note that there is a "total" and "confirmed score" - the confirmed score means that in the original study the lower eGFR was confirmed by a second measurement (which is what we do here) therefore we need to stick with the confirmed score.

ckd60_cohort$risk_decile <- ntile(ckd60_cohort$ckdpc_egfr60_confirmed_score, n.quantiles)

## Get mean predicted probabilities by studydrug
predicted <- ckd60_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd60_pred=mean(ckdpc_egfr60_confirmed_score)/100)

predicted_all <- noncal_ckd60_cohort %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd60_pred=mean(ckdpc_egfr60_confirmed_score_cal)/100)

## Find actual observed probabilities by risk score category and studydrug

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

EST.su <- SE.su <-
  EST.dpp4 <- SE.dpp4 <-
  EST.sglt2 <- SE.sglt2 <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_su <- tibble() %>% mutate(
  observed_su=NA,
  lower_ci_su=NA,
  upper_ci_su=NA,
  strata=NA
)

observed_dpp4 <- tibble() %>% mutate(
  observed_dpp4=NA,
  lower_ci_dpp4=NA,
  upper_ci_dpp4=NA,
  strata=NA
)

observed_sglt2 <- tibble() %>% mutate(
  observed_sglt2=NA,
  lower_ci_sglt2=NA,
  upper_ci_sglt2=NA,
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
  
  observed_su_ckd60 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, 
                               data=ckd60_cohort[ckd60_cohort$.imp == i & 
                                                   ckd60_cohort$risk_decile == k &
                                                   ckd60_cohort$studydrug=="SU",]) %>%
    tidy() %>%
    # group_by(strata) %>%
    filter(time==max(time))
  
  EST.su[k,i] <- observed_su_ckd60$estimate
  SE.su[k,i] <- observed_su_ckd60$std.error
  
  observed_dpp4_ckd60 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, 
                               data=ckd60_cohort[ckd60_cohort$.imp == i & 
                                                   ckd60_cohort$risk_decile == k &
                                                   ckd60_cohort$studydrug=="DPP4",]) %>%
    tidy() %>%
    # group_by(strata) %>%
    filter(time==max(time))
  
  EST.dpp4[k,i] <- observed_dpp4_ckd60$estimate
  SE.dpp4[k,i] <- observed_dpp4_ckd60$std.error
  
  observed_sglt2_ckd60 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, 
                                  data=ckd60_cohort[ckd60_cohort$.imp == i & 
                                                      ckd60_cohort$risk_decile == k &
                                                      ckd60_cohort$studydrug=="SGLT2",]) %>%
    tidy() %>%
    # group_by(strata) %>%
    filter(time==max(time))
  
  EST.sglt2[k,i] <- observed_sglt2_ckd60$estimate
  SE.sglt2[k,i] <- observed_sglt2_ckd60$std.error
  
  observed_all_ckd60 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, 
                                  data=ckd60_cohort[ckd60_cohort$.imp == i & 
                                                      ckd60_cohort$risk_decile == k,]) %>%
    tidy() %>%
    # group_by(strata) %>%
    filter(time==max(time))
  
  EST.all[k,i] <- observed_all_ckd60$estimate
  SE.all[k,i] <- observed_all_ckd60$std.error
}
  
  est.su <- pool.rubin.KM(EST.su[k,], SE.su[k,], n.imp)
  observed_su[k,] <- observed_su[k,] %>% 
    mutate(
      observed_su=est.su[1],
      lower_ci_su=est.su[2],
      upper_ci_su=est.su[3],
      strata=k
    )
  
  est.dpp4 <- pool.rubin.KM(EST.dpp4[k,], SE.dpp4[k,], n.imp)
  observed_dpp4[k,] <- observed_dpp4[k,] %>% 
    mutate(
      observed_dpp4=est.dpp4[1],
      lower_ci_dpp4=est.dpp4[2],
      upper_ci_dpp4=est.dpp4[3],
      strata=k
    )
  
  est.sglt2 <- pool.rubin.KM(EST.sglt2[k,], SE.sglt2[k,], n.imp)
  observed_sglt2[k,] <- observed_sglt2[k,] %>% 
    mutate(
      observed_sglt2=est.sglt2[1],
      lower_ci_sglt2=est.sglt2[2],
      upper_ci_sglt2=est.sglt2[3],
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


su_events <- ckd60_cohort %>%
  filter(studydrug=="SU" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SU=round(n()/n.imp, 0))

dpp4_events <- ckd60_cohort %>%
  filter(studydrug=="DPP4" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4=round(n()/n.imp, 0))

sglt2_events <- ckd60_cohort %>%
  filter(studydrug=="SGLT2" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=round(n()/n.imp, 0))

obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="SU")), observed_su),
  cbind((predicted %>% filter(studydrug=="DPP4")), observed_dpp4),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_su, observed_dpp4, observed_sglt2),
         lower_ci=coalesce(lower_ci_su, lower_ci_dpp4, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_su, upper_ci_dpp4, upper_ci_sglt2))

events_table <- data.frame(t(su_events %>%
  inner_join(sglt2_events) %>% 
    inner_join(dpp4_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd60_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd60_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,60,by=10)), limits=c(-10,75)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated risk score vs CKD incidence (5 year)") +
  coord_cartesian(ylim = c(0,60))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))



## obs vs pred calibration plot:
p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd60_pred*100)) +
  #geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed*100, group=studydrug, color=studydrug), position=dodge, shape=17, size=2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab("Predicted risk (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,33))+
  scale_y_continuous(limits=c(-10,40)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Predicted versus observed CKD incidence (5 year)") +
  coord_cartesian(ylim = c(0,33))

p1 + geom_smooth(aes(y=observed*100, x=mean_ckd60_pred*100), colour = "black", data = obs_v_pred)



## plot with all observed combined vs predicted with loess smoother

obs_v_pred <- cbind(predicted_all, observed_all)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.1,size=.75, position=dodge) +
  geom_point(aes(y = observed*100), shape = 21, position=dodge, size = 4, colour = "black", fill = "cadetblue3", stroke = 1.5) +
  geom_point(aes(y = mean_ckd60_pred*100), position=dodge, shape=21, size=3, colour = "black", fill = "white", stroke = 1.5) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,60,by=10)), limits=c(-2,75)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated risk score vs CKD incidence (5 year)") +
  coord_cartesian(ylim = c(0,60))

p1

## C-stat
ckd60_cohort <- ckd60_cohort %>%
  mutate(ckdpc_egfr60_confirmed_survival=(100-ckdpc_egfr60_confirmed_score)/100)

surv_mod_ckd60 <- coxph(Surv(ckd_345_censtime_yrs, ckd_345_censvar)~ckdpc_egfr60_confirmed_survival, data=ckd60_cohort, method="breslow")

cstat_est <- summary(surv_mod_ckd60)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd60)$concordance[1]-(1.96*summary(surv_mod_ckd60)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd60)$concordance[1]+(1.96*summary(surv_mod_ckd60)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.8927, 95% CI 0.8909-0.8944

# the above clearly shows the risk score overestimates risk


############################2B UNCALIBRATED SCORES - CKD PROGRESSION################################################################

# 40% decline in eGFR / CKD stage 5

ckd40_cohort$risk_decile <- ntile(ckd40_cohort$ckdpc_40egfr_score, n.quantiles)

## Get mean predicted probabilities by studydrug
predicted <- ckd40_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score)/100)

predicted_all <- ckd40_cohort %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score)/100)

## Find actual observed probabilities by risk score category and studydrug


EST.su <- SE.su <-
  EST.dpp4 <- SE.dpp4 <-
  EST.sglt2 <- SE.sglt2 <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_su <- tibble() %>% mutate(
  observed_su=NA,
  lower_ci_su=NA,
  upper_ci_su=NA,
  strata=NA
)

observed_dpp4 <- tibble() %>% mutate(
  observed_dpp4=NA,
  lower_ci_dpp4=NA,
  upper_ci_dpp4=NA,
  strata=NA
)

observed_sglt2 <- tibble() %>% mutate(
  observed_sglt2=NA,
  lower_ci_sglt2=NA,
  upper_ci_sglt2=NA,
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
    
    observed_su_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                 data=ckd40_cohort[ckd40_cohort$.imp == i & 
                                                     ckd40_cohort$risk_decile == k &
                                                     ckd40_cohort$studydrug=="SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.su[k,i] <- observed_su_ckd40$estimate
    SE.su[k,i] <- observed_su_ckd40$std.error
    
    observed_dpp4_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                   data=ckd40_cohort[ckd40_cohort$.imp == i & 
                                                       ckd40_cohort$risk_decile == k &
                                                       ckd40_cohort$studydrug=="DPP4",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4[k,i] <- observed_dpp4_ckd40$estimate
    SE.dpp4[k,i] <- observed_dpp4_ckd40$std.error
    
    observed_sglt2_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                    data=ckd40_cohort[ckd40_cohort$.imp == i & 
                                                        ckd40_cohort$risk_decile == k &
                                                        ckd40_cohort$studydrug=="SGLT2",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2[k,i] <- observed_sglt2_ckd40$estimate
    SE.sglt2[k,i] <- observed_sglt2_ckd40$std.error
    
    observed_all_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                    data=ckd40_cohort[ckd40_cohort$.imp == i & 
                                                        ckd40_cohort$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd40$estimate
    SE.all[k,i] <- observed_all_ckd40$std.error
  }
  
  est.su <- pool.rubin.KM(EST.su[k,], SE.su[k,], n.imp)
  observed_su[k,] <- observed_su[k,] %>% 
    mutate(
      observed_su=est.su[1],
      lower_ci_su=est.su[2],
      upper_ci_su=est.su[3],
      strata=k
    )
  
  est.dpp4 <- pool.rubin.KM(EST.dpp4[k,], SE.dpp4[k,], n.imp)
  observed_dpp4[k,] <- observed_dpp4[k,] %>% 
    mutate(
      observed_dpp4=est.dpp4[1],
      lower_ci_dpp4=est.dpp4[2],
      upper_ci_dpp4=est.dpp4[3],
      strata=k
    )
  
  est.sglt2 <- pool.rubin.KM(EST.sglt2[k,], SE.sglt2[k,], n.imp)
  observed_sglt2[k,] <- observed_sglt2[k,] %>% 
    mutate(
      observed_sglt2=est.sglt2[1],
      lower_ci_sglt2=est.sglt2[2],
      upper_ci_sglt2=est.sglt2[3],
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


su_events <- ckd40_cohort %>%
  filter(studydrug=="SU" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SU=round(n()/n.imp, 0))

dpp4_events <- ckd40_cohort %>%
  filter(studydrug=="DPP4" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4=round(n()/n.imp, 0))

sglt2_events <- ckd40_cohort %>%
  filter(studydrug=="SGLT2" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="SU")), observed_su),
  cbind((predicted %>% filter(studydrug=="DPP4")), observed_dpp4),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_su, observed_dpp4, observed_sglt2),
         lower_ci=coalesce(lower_ci_su, lower_ci_dpp4, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_su, upper_ci_dpp4, upper_ci_sglt2))

events_table <- data.frame(t(su_events %>% inner_join(dpp4_events) %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd40_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-10,10)) +
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

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

# plot with obs vs predicted

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  #geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed*100, group=studydrug, color=studydrug), position=dodge, shape=17, size=2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab("Predicted risk (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,3.5))+
  scale_y_continuous(limits=c(-5,10)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Predicted vs observed incidence of 40% eGFR decline/ESKD (3 year)")+
  coord_cartesian(ylim = c(0,3.5))

p1 + geom_smooth(aes(y=observed*100, x=mean_ckd40_pred*100), colour = "black", data = obs_v_pred)

## plot with all observed combined vs predicted

obs_v_pred <- cbind(predicted_all, observed_all)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.1,size=.75, position=dodge) +
  geom_point(aes(y = observed*100), shape = 21, position=dodge, size = 4, colour = "black", fill = "cadetblue3", stroke = 1.5) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=21, size=3, colour = "black", fill = "white", stroke = 1.5) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Uncalibrated risk score vs 40% eGFR decline/ESKD (3 year)") +
  coord_cartesian(ylim = c(0,7.5))

p1

## C-stat
ckd40_cohort <- ckd40_cohort %>%
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score)/100)

surv_mod_ckd40 <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival, data=ckd40_cohort, method="breslow")
cstat_est <- summary(surv_mod_ckd40)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd40)$concordance[1]-(1.96*summary(surv_mod_ckd40)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd40)$concordance[1]+(1.96*summary(surv_mod_ckd40)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7198, 95% CI 0.7122-0.7274


# similarly, this risk score also overestimates risk.
# both will need some recalibration.

############################3A RECALIBRATION - NEW CKD################################################################

# 3 Recalibration

# CKD < 60

## Assign random 20% as recalibration cohort and remove from main cohort

cal_ckd60_cohort <- data.frame()

for (i in 1:n.imp) {
  temp <- ckd60_cohort[ckd60_cohort$.imp == i & !ckd60_cohort$studydrug == "SGLT2",] %>% slice_sample(prop=0.2)
  cal_ckd60_cohort <- rbind(cal_ckd60_cohort, temp)
  rm(temp)
}

#184690

noncal_ckd60_cohort <- ckd60_cohort %>%
  anti_join(cal_ckd60_cohort, by=c("patid", "dstartdate", "studydrug", ".imp"))
table(noncal_ckd60_cohort$studydrug)
#    SU   DPP4  SGLT2 
#371095 367675 209640 

## Original survival constant: 
#baseline hazard:
EHRBiomarkr::ckdpcEgfr60RiskConstants$new_surv_confirmed %>% as.numeric()
#0.01221876
#Baseline hazard (BH) = exp(-5^survival constant)
#therefore original survival constant was (log(-log(BH)))/log(5)
as.numeric(EHRBiomarkr::ckdpcEgfr60RiskConstants$new_surv_confirmed) %>% log() %>% abs() %>% log(base=5) 
#0.9212477

#Now we will recalculate a survival constant in the calibration cohort that we will use to recalibrate in the non-calibration cohort.
new_ckd60_surv <- rep(NA, n.imp)
for (i in 1:n.imp) {
recal_mod <- coxph(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ offset(ckdpc_egfr60_confirmed_lin_predictor), 
                   data=cal_ckd60_cohort[cal_ckd60_cohort$.imp == i,])
new_ckd60_surv[i] <- summary(survfit(recal_mod),time=5)$surv
sprintf("%.15f", new_ckd60_surv[i])
}
new_ckd60_surv2 <- mean(new_ckd60_surv)
# 0.970392
## new surv constant=exp(-5^new_ckd60_surv2)

## Recalculate scores in rest of cohort
noncal_ckd60_cohort <- noncal_ckd60_cohort %>%
  mutate(centred_egfr60_lin_predictor=ckdpc_egfr60_confirmed_lin_predictor-mean(ckdpc_egfr60_confirmed_lin_predictor)) %>%
  ungroup() %>%
  mutate(ckdpc_egfr60_confirmed_score_cal=(1-(new_ckd60_surv2^exp(centred_egfr60_lin_predictor)))*100)

cohort <- cohort %>%
  mutate(centred_egfr60_lin_predictor=ckdpc_egfr60_confirmed_lin_predictor-mean(ckdpc_egfr60_confirmed_lin_predictor, na.rm = T)) %>%
  ungroup() %>%
  mutate(ckdpc_egfr60_confirmed_score_cal=(1-(new_ckd60_surv^exp(centred_egfr60_lin_predictor)))*100)

## Plot
noncal_ckd60_cohort$risk_decile <- ntile(noncal_ckd60_cohort$ckdpc_egfr60_confirmed_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- noncal_ckd60_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd60_pred=mean(ckdpc_egfr60_confirmed_score_cal)/100)

predicted_all <- noncal_ckd60_cohort %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd60_pred=mean(ckdpc_egfr60_confirmed_score_cal)/100)

### Find actual observed probabilities by risk score category and studydrug

EST.su <- SE.su <-
  EST.dpp4 <- SE.dpp4 <-
  EST.sglt2 <- SE.sglt2 <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_su <- tibble() %>% mutate(
  observed_su=NA,
  lower_ci_su=NA,
  upper_ci_su=NA,
  strata=NA
)

observed_dpp4 <- tibble() %>% mutate(
  observed_dpp4=NA,
  lower_ci_dpp4=NA,
  upper_ci_dpp4=NA,
  strata=NA
)

observed_sglt2 <- tibble() %>% mutate(
  observed_sglt2=NA,
  lower_ci_sglt2=NA,
  upper_ci_sglt2=NA,
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
    
    observed_su_ckd60 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, 
                                 data=noncal_ckd60_cohort[noncal_ckd60_cohort$.imp == i & 
                                                     noncal_ckd60_cohort$risk_decile == k &
                                                     noncal_ckd60_cohort$studydrug=="SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.su[k,i] <- observed_su_ckd60$estimate
    SE.su[k,i] <- observed_su_ckd60$std.error
    
    observed_dpp4_ckd60 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, 
                                   data=noncal_ckd60_cohort[noncal_ckd60_cohort$.imp == i & 
                                                       noncal_ckd60_cohort$risk_decile == k &
                                                       noncal_ckd60_cohort$studydrug=="DPP4",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4[k,i] <- observed_dpp4_ckd60$estimate
    SE.dpp4[k,i] <- observed_dpp4_ckd60$std.error
    
    observed_sglt2_ckd60 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, 
                                    data=noncal_ckd60_cohort[noncal_ckd60_cohort$.imp == i & 
                                                               noncal_ckd60_cohort$risk_decile == k &
                                                               noncal_ckd60_cohort$studydrug=="SGLT2",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2[k,i] <- observed_sglt2_ckd60$estimate
    SE.sglt2[k,i] <- observed_sglt2_ckd60$std.error
    
    observed_all_ckd60 <- survfit(Surv(ckd_345_censtime_yrs, ckd_345_censvar) ~ risk_decile, 
                                    data=noncal_ckd60_cohort[noncal_ckd60_cohort$.imp == i & 
                                                               noncal_ckd60_cohort$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd60$estimate
    SE.all[k,i] <- observed_all_ckd60$std.error
  }
  
  est.su <- pool.rubin.KM(EST.su[k,], SE.su[k,], n.imp)
  observed_su[k,] <- observed_su[k,] %>% 
    mutate(
      observed_su=est.su[1],
      lower_ci_su=est.su[2],
      upper_ci_su=est.su[3],
      strata=k
    )
  
  est.dpp4 <- pool.rubin.KM(EST.dpp4[k,], SE.dpp4[k,], n.imp)
  observed_dpp4[k,] <- observed_dpp4[k,] %>% 
    mutate(
      observed_dpp4=est.dpp4[1],
      lower_ci_dpp4=est.dpp4[2],
      upper_ci_dpp4=est.dpp4[3],
      strata=k
    )
  
  est.sglt2 <- pool.rubin.KM(EST.sglt2[k,], SE.sglt2[k,], n.imp)
  observed_sglt2[k,] <- observed_sglt2[k,] %>% 
    mutate(
      observed_sglt2=est.sglt2[1],
      lower_ci_sglt2=est.sglt2[2],
      upper_ci_sglt2=est.sglt2[3],
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


su_events <- noncal_ckd60_cohort %>%
  filter(studydrug=="SU" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SU=round(n()/n.imp, 0))

dpp4_events <- noncal_ckd60_cohort %>%
  filter(studydrug=="DPP4" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4=round(n()/n.imp, 0))

sglt2_events <- noncal_ckd60_cohort %>%
  filter(studydrug=="SGLT2" & ckd_345_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=round(n()/n.imp, 0))

obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="SU")), observed_su),
  cbind((predicted %>% filter(studydrug=="DPP4")), observed_dpp4),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_su, observed_dpp4, observed_sglt2),
         lower_ci=coalesce(lower_ci_su, lower_ci_dpp4, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_su, upper_ci_dpp4, upper_ci_sglt2))

events_table <- data.frame(t(su_events %>%
                               inner_join(sglt2_events) %>% 
                               inner_join(dpp4_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd60_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd60_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,60,by=10)), limits=c(-2,75)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Recalibrated risk score vs CKD incidence (5 year)") +
  coord_cartesian(ylim = c(0,60))

p2 <- gridExtra::tableGrob(events_table, rows = NULL, cols = NULL)
p2$widths <- unit(rep(1, ncol(p2)), "null")
p2$heights <- unit(rep(1, nrow(p2)), "null")

p3 <- ggplot() +
  annotation_custom(p2)

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))


## obs vs pred calibration plot:
p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd60_pred*100)) +
  #geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed*100, group=studydrug, color=studydrug), position=dodge, shape=17, size=2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab("Predicted risk (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,33))+
  scale_y_continuous(limits=c(0,50)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Predicted versus observed CKD incidence (5 year)") +
  coord_cartesian(ylim = c(0,33))

p1 + geom_smooth(aes(y=observed*100, x=mean_ckd60_pred*100), colour = "black", data = obs_v_pred)

## plot with all observed combined vs predicted

obs_v_pred <- cbind(predicted_all, observed_all)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.1,size=.75, position=dodge) +
  geom_point(aes(y = observed*100), shape = 21, position=dodge, size = 4, colour = "black", fill = "cadetblue3", stroke = 1.5) +
  geom_point(aes(y = mean_ckd60_pred*100), position=dodge, shape=21, size=3, colour = "black", fill = "white", stroke = 1.5) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,60,by=10)), limits=c(-2,75)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Recalibrated risk score vs CKD incidence (5 year)") +
  coord_cartesian(ylim = c(0,60))

p1


## C-stat
noncal_ckd60_cohort <- noncal_ckd60_cohort %>%
  mutate(ckdpc_egfr60_confirmed_survival=(100-ckdpc_egfr60_confirmed_score_cal)/100)

surv_mod_ckd60 <- coxph(Surv(ckd_345_censtime_yrs, ckd_345_censvar)~ckdpc_egfr60_confirmed_survival, data=noncal_ckd60_cohort, method="breslow")

cstat_est <- summary(surv_mod_ckd60)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd60)$concordance[1]-(1.96*summary(surv_mod_ckd60)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd60)$concordance[1]+(1.96*summary(surv_mod_ckd60)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.892, 95% CI 0.89-0.8939



############################3B RECALIBRATION - CKD PROGRESSION################################################################

# 40% decline in eGFR / CKD 5
# risk equation uses a logistic model. There are several ways to update these models.
# the most straightforward one is to adjust the intercept ("calibration in the large")
#see: Journal of Clinical Epidemiology 61 (2008) 76e86. Janssen KJM, Moons KGM, Kalkman CJ, Grobbee DE, Vergouwe Y. Updating methods improved the performance of a clinical prediction model in new patients


## Assign random 20% as calibration cohort and remove from main cohort

cal_ckd40_cohort <- data.frame()

for (i in 1:n.imp) {
  temp <- ckd40_cohort[ckd40_cohort$.imp == i & !ckd40_cohort$studydrug == "SGLT2",] %>% slice_sample(prop=0.2)
  cal_ckd40_cohort <- rbind(cal_ckd40_cohort, temp)
  rm(temp)
}

noncal_ckd40_cohort <- ckd40_cohort %>%
  anti_join(cal_ckd40_cohort, by=c("patid", "dstartdate", "studydrug", ".imp"))
table(noncal_ckd40_cohort$studydrug)
#    SU   DPP4  SGLT2 
#371030 367740 209640                                                                                                                                                                         371030 367740 209640

#original intercept:
EHRBiomarkr::ckdpc40EgfrRiskConstants$intercept %>% as.numeric()

#calculate correction factor for intercept = ln(observed outcome frequency/1-observed outcome frequency / predicted outcome frequency/1-predicted outcome frequency)
correction_factor <- rep(NA, n.imp)
for (i in 1:n.imp) {
correction_factor[i] <- ((mean(cal_ckd40_cohort[cal_ckd40_cohort$.imp == i,]$ckd_egfr40_censvar) / 
                         (1-mean(cal_ckd40_cohort[cal_ckd40_cohort$.imp == i,]$ckd_egfr40_censvar))) /
                        (mean(cal_ckd40_cohort[cal_ckd40_cohort$.imp == i,]$ckdpc_40egfr_score*0.01) / 
                           (1-mean(cal_ckd40_cohort[cal_ckd40_cohort$.imp == i,]$ckdpc_40egfr_score*0.01))) ) %>% 
  log()
}

#new intercept = original intercept + correction factor
new_intercept <- EHRBiomarkr::ckdpc40EgfrRiskConstants$intercept %>% as.numeric() + mean(correction_factor)

## Recalculate scores in rest of cohort
noncal_ckd40_cohort <- noncal_ckd40_cohort %>%
  mutate(centred_40egfr_lin_predictor=ckdpc_40egfr_lin_predictor-mean(ckdpc_40egfr_lin_predictor)) %>%
  ungroup() %>%
  mutate(ckdpc_40egfr_score_cal=100 * (exp(centred_40egfr_lin_predictor + new_intercept)/(1 + exp(centred_40egfr_lin_predictor + new_intercept))))


## Show observed (estimated frequency binned per risk decile) vs predicted (risk-score predicted)
noncal_ckd40_cohort$risk_decile <- ntile(noncal_ckd40_cohort$ckdpc_40egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- noncal_ckd40_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

predicted_all <- noncal_ckd40_cohort %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

### Find actual observed probabilities by risk score category and studydrug

EST.su <- SE.su <-
  EST.dpp4 <- SE.dpp4 <-
  EST.sglt2 <- SE.sglt2 <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_su <- tibble() %>% mutate(
  observed_su=NA,
  lower_ci_su=NA,
  upper_ci_su=NA,
  strata=NA
)

observed_dpp4 <- tibble() %>% mutate(
  observed_dpp4=NA,
  lower_ci_dpp4=NA,
  upper_ci_dpp4=NA,
  strata=NA
)

observed_sglt2 <- tibble() %>% mutate(
  observed_sglt2=NA,
  lower_ci_sglt2=NA,
  upper_ci_sglt2=NA,
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
    
    observed_su_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                 data=noncal_ckd40_cohort[noncal_ckd40_cohort$.imp == i & 
                                                            noncal_ckd40_cohort$risk_decile == k &
                                                            noncal_ckd40_cohort$studydrug=="SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.su[k,i] <- observed_su_ckd40$estimate
    SE.su[k,i] <- observed_su_ckd40$std.error
    
    observed_dpp4_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                   data=noncal_ckd40_cohort[noncal_ckd40_cohort$.imp == i & 
                                                              noncal_ckd40_cohort$risk_decile == k &
                                                              noncal_ckd40_cohort$studydrug=="DPP4",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4[k,i] <- observed_dpp4_ckd40$estimate
    SE.dpp4[k,i] <- observed_dpp4_ckd40$std.error
    
    observed_sglt2_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                    data=noncal_ckd40_cohort[noncal_ckd40_cohort$.imp == i & 
                                                               noncal_ckd40_cohort$risk_decile == k &
                                                               noncal_ckd40_cohort$studydrug=="SGLT2",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2[k,i] <- observed_sglt2_ckd40$estimate
    SE.sglt2[k,i] <- observed_sglt2_ckd40$std.error
    
    observed_all_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                    data=noncal_ckd40_cohort[noncal_ckd40_cohort$.imp == i & 
                                                               noncal_ckd40_cohort$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd40$estimate
    SE.all[k,i] <- observed_all_ckd40$std.error
    
  }
  
  est.su <- pool.rubin.KM(EST.su[k,], SE.su[k,], n.imp)
  observed_su[k,] <- observed_su[k,] %>% 
    mutate(
      observed_su=est.su[1],
      lower_ci_su=est.su[2],
      upper_ci_su=est.su[3],
      strata=k
    )
  
  est.dpp4 <- pool.rubin.KM(EST.dpp4[k,], SE.dpp4[k,], n.imp)
  observed_dpp4[k,] <- observed_dpp4[k,] %>% 
    mutate(
      observed_dpp4=est.dpp4[1],
      lower_ci_dpp4=est.dpp4[2],
      upper_ci_dpp4=est.dpp4[3],
      strata=k
    )
  
  est.sglt2 <- pool.rubin.KM(EST.sglt2[k,], SE.sglt2[k,], n.imp)
  observed_sglt2[k,] <- observed_sglt2[k,] %>% 
    mutate(
      observed_sglt2=est.sglt2[1],
      lower_ci_sglt2=est.sglt2[2],
      upper_ci_sglt2=est.sglt2[3],
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


su_events <- noncal_ckd40_cohort %>%
  filter(studydrug=="SU" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SU=round(n()/n.imp, 0))

dpp4_events <- noncal_ckd40_cohort %>%
  filter(studydrug=="DPP4" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4=round(n()/n.imp, 0))

sglt2_events <- noncal_ckd40_cohort %>%
  filter(studydrug=="SGLT2" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="SU")), observed_su),
  cbind((predicted %>% filter(studydrug=="DPP4")), observed_dpp4),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_su, observed_dpp4, observed_sglt2),
         lower_ci=coalesce(lower_ci_su, lower_ci_dpp4, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_su, upper_ci_dpp4, upper_ci_sglt2))

events_table <- data.frame(t(su_events %>% inner_join(dpp4_events) %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd40_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
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

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))


## plot with all observed combined vs predicted

obs_v_pred <- cbind(predicted_all, observed_all)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.1,size=.75, position=dodge) +
  geom_point(aes(y = observed*100), shape = 21, position=dodge, size = 4, colour = "black", fill = "cadetblue3", stroke = 1.5) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=21, size=3, colour = "black", fill = "white", stroke = 1.5) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Recalibrated risk score vs 40% eGFR decline/ESKD (3 year)\n(Calibration-in-the-large)") +
  coord_cartesian(ylim = c(0,7.5))

p1

## C-stat
noncal_ckd40_cohort <- noncal_ckd40_cohort %>%
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score_cal)/100)

surv_mod_ckd40 <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival, data=noncal_ckd40_cohort, method="breslow")
cstat_est <- summary(surv_mod_ckd40)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd40)$concordance[1]-(1.96*summary(surv_mod_ckd40)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd40)$concordance[1]+(1.96*summary(surv_mod_ckd40)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7192, 95% CI 0.7107-0.7278

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  #geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed*100, group=studydrug, color=studydrug), position=dodge, shape=17, size=2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab("Predicted risk (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,3.5))+
  scale_y_continuous(limits=c(-1,7)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Predicted vs observed incidence of 40% eGFR decline/ESKD (3 year)") +
  coord_cartesian(ylim = c(0,3.5))

p1 + geom_smooth(aes(y=observed*100, x=mean_ckd40_pred*100), colour = "black", data = obs_v_pred)

#before recalibrating, the risk score systematically overestimated risk. However after recalibration, it underestimates risk in those at high risk.
#an alternative method is logistic calibration:
#a logistic regression model is fitted with the linear predictor as the only covariate in the updating set.
#The calibration slope is used to recalibrate the original regression coefficients the regression coefficients of the original model (multiplied with the calibration slope).
#The intercept of the original prediction model is adjusted by adding the calibration intercept

intercept_cf <- coef_cf <- rep(NA,n.imp)

for (i in 1:n.imp) {
recal_mod <- glm(ckd_egfr40_censvar ~ (ckdpc_40egfr_lin_predictor), 
                 data = cal_ckd40_cohort[cal_ckd40_cohort$.imp == i,], family = "binomial")

intercept_cf[i] <- recal_mod$coef[1]
coef_cf[i] <- recal_mod$coef[2]

}

# new intercept:
new_intercept <- EHRBiomarkr::ckdpc40EgfrRiskConstants$intercept %>% as.numeric() + mean(intercept_cf)


## Recalculate scores in rest of cohort
noncal_ckd40_cohort <- noncal_ckd40_cohort %>%
  mutate(centred_40egfr_lin_predictor=ckdpc_40egfr_lin_predictor-mean(ckdpc_40egfr_lin_predictor)) %>%
  mutate(recal_40egfr_lin_predictor=centred_40egfr_lin_predictor*mean(coef_cf)) %>%
  ungroup() %>%
  mutate(ckdpc_40egfr_score_cal= 
           100 *      100 *
           (exp(recal_40egfr_lin_predictor + new_intercept)/(1 + exp(recal_40egfr_lin_predictor + new_intercept))))

## Plot
noncal_ckd40_cohort$risk_decile <- ntile(noncal_ckd40_cohort$ckdpc_40egfr_score_cal, n.quantiles)

### Get mean predicted probabilities by studydrug
predicted <- noncal_ckd40_cohort %>%
  group_by(risk_decile, studydrug) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

predicted_all <- noncal_ckd40_cohort %>%
  group_by(risk_decile) %>%
  summarise(mean_ckd40_pred=mean(ckdpc_40egfr_score_cal)/100)

### Find actual observed probabilities by risk score category and studydrug

EST.su <- SE.su <-
  EST.dpp4 <- SE.dpp4 <-
  EST.sglt2 <- SE.sglt2 <- 
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_su <- tibble() %>% mutate(
  observed_su=NA,
  lower_ci_su=NA,
  upper_ci_su=NA,
  strata=NA
)

observed_dpp4 <- tibble() %>% mutate(
  observed_dpp4=NA,
  lower_ci_dpp4=NA,
  upper_ci_dpp4=NA,
  strata=NA
)

observed_sglt2 <- tibble() %>% mutate(
  observed_sglt2=NA,
  lower_ci_sglt2=NA,
  upper_ci_sglt2=NA,
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
    
    observed_su_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                 data=noncal_ckd40_cohort[noncal_ckd40_cohort$.imp == i & 
                                                            noncal_ckd40_cohort$risk_decile == k &
                                                            noncal_ckd40_cohort$studydrug=="SU",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.su[k,i] <- observed_su_ckd40$estimate
    SE.su[k,i] <- observed_su_ckd40$std.error
    
    observed_dpp4_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                   data=noncal_ckd40_cohort[noncal_ckd40_cohort$.imp == i & 
                                                              noncal_ckd40_cohort$risk_decile == k &
                                                              noncal_ckd40_cohort$studydrug=="DPP4",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.dpp4[k,i] <- observed_dpp4_ckd40$estimate
    SE.dpp4[k,i] <- observed_dpp4_ckd40$std.error
    
    observed_sglt2_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                    data=noncal_ckd40_cohort[noncal_ckd40_cohort$.imp == i & 
                                                               noncal_ckd40_cohort$risk_decile == k &
                                                               noncal_ckd40_cohort$studydrug=="SGLT2",]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.sglt2[k,i] <- observed_sglt2_ckd40$estimate
    SE.sglt2[k,i] <- observed_sglt2_ckd40$std.error
    
    observed_all_ckd40 <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                  data=noncal_ckd40_cohort[noncal_ckd40_cohort$.imp == i & 
                                                             noncal_ckd40_cohort$risk_decile == k,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.all[k,i] <- observed_all_ckd40$estimate
    SE.all[k,i] <- observed_all_ckd40$std.error
    
  }
  
  est.su <- pool.rubin.KM(EST.su[k,], SE.su[k,], n.imp)
  observed_su[k,] <- observed_su[k,] %>% 
    mutate(
      observed_su=est.su[1],
      lower_ci_su=est.su[2],
      upper_ci_su=est.su[3],
      strata=k
    )
  
  est.dpp4 <- pool.rubin.KM(EST.dpp4[k,], SE.dpp4[k,], n.imp)
  observed_dpp4[k,] <- observed_dpp4[k,] %>% 
    mutate(
      observed_dpp4=est.dpp4[1],
      lower_ci_dpp4=est.dpp4[2],
      upper_ci_dpp4=est.dpp4[3],
      strata=k
    )
  
  est.sglt2 <- pool.rubin.KM(EST.sglt2[k,], SE.sglt2[k,], n.imp)
  observed_sglt2[k,] <- observed_sglt2[k,] %>% 
    mutate(
      observed_sglt2=est.sglt2[1],
      lower_ci_sglt2=est.sglt2[2],
      upper_ci_sglt2=est.sglt2[3],
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


su_events <- noncal_ckd40_cohort %>%
  filter(studydrug=="SU" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SU=round(n()/n.imp, 0))

dpp4_events <- noncal_ckd40_cohort %>%
  filter(studydrug=="DPP4" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(DPP4=round(n()/n.imp, 0))

sglt2_events <- noncal_ckd40_cohort %>%
  filter(studydrug=="SGLT2" & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(SGLT2=round(n()/n.imp, 0))


obs_v_pred <- rbind(
  cbind((predicted %>% filter(studydrug=="SU")), observed_su),
  cbind((predicted %>% filter(studydrug=="DPP4")), observed_dpp4),
  cbind((predicted %>% filter(studydrug=="SGLT2")), observed_sglt2)
) %>%
  mutate(observed=coalesce(observed_su, observed_dpp4, observed_sglt2),
         lower_ci=coalesce(lower_ci_su, lower_ci_dpp4, lower_ci_sglt2),
         upper_ci=coalesce(upper_ci_su, upper_ci_dpp4, upper_ci_sglt2))

events_table <- data.frame(t(su_events %>% inner_join(dpp4_events) %>%
                               inner_join(sglt2_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_ckd40_pred=NA, risk_decile=0)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile, group=studydrug, color=studydrug)) +
  geom_point(aes(y = observed*100), position=dodge) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=17, size=2) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
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

p1 + p3 + plot_layout(ncol = 1, heights=c(5,1))

## obs v predicted plot
p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_ckd40_pred*100)) +
  #geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.25,size=1, position=dodge) +
  geom_point(aes(y = observed*100, group=studydrug, color=studydrug), position=dodge, shape=17, size=2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab("Predicted risk (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,3.5))+
  scale_y_continuous(limits=c(-1,7)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Predicted vs observed incidence of 40% eGFR decline/ESKD (3 year)\n(Logistic recalibration)") +
  coord_cartesian(ylim = c(0,3.5))

p1 + geom_smooth(aes(y=observed*100, x=mean_ckd40_pred*100), colour = "black", data = obs_v_pred)

## plot with all observed combined vs predicted

obs_v_pred <- cbind(predicted_all, observed_all)

p1 <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=risk_decile)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100),width=0.1,size=.75, position=dodge) +
  geom_point(aes(y = observed*100), shape = 21, position=dodge, size = 4, colour = "black", fill = "cadetblue3", stroke = 1.5) +
  geom_point(aes(y = mean_ckd40_pred*100), position=dodge, shape=21, size=3, colour = "black", fill = "white", stroke = 1.5) +
  theme_bw() +
  xlab("Risk score decile") + ylab("Risk (%)")+
  scale_x_continuous(breaks=c(seq(0,10,by=1)))+
  scale_y_continuous(breaks=c(seq(0,7,by=1)), limits=c(-1,10)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5))) +
  ggtitle("Recalibrated risk score vs 40% eGFR decline/ESKD (3 year)\n(Logistic recalibration)") +
  coord_cartesian(ylim = c(0,7.5))

p1

## C-stat
noncal_ckd40_cohort <- noncal_ckd40_cohort %>%
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score_cal)/100)

surv_mod_ckd40 <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival, data=noncal_ckd40_cohort, method="breslow")
cstat_est <- summary(surv_mod_ckd40)$concordance[1]
cstat_est_ll <- summary(surv_mod_ckd40)$concordance[1]-(1.96*summary(surv_mod_ckd40)$concordance[2])
cstat_est_ul <- summary(surv_mod_ckd40)$concordance[1]+(1.96*summary(surv_mod_ckd40)$concordance[2])
paste0("C statistic: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))
# C statistic: 0.7192, 95% CI 0.7107-0.7278

