############################0 SETUP################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(rms)
library(tableone)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())

set.seed(123)
n.imp <- 10
n.quantiles <- 10
#set default colours (colour-blind accessible from the Okabe-Ito palette) for different drug classes
#this is the colour allocation that will be used across the department
cols <- c("SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00")
#in further analyses, the dpp4/su group will be combined, and we will use the dpp4 colour for this (strongest contrast)
cols <- c(cols, "DPP4i/SU" = "#0072B2")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-03-06_t2d_ckdpc_recalibrated_incl_egfr_below_60.Rda")

noncal_cohort$studydrug2 <- as.factor(noncal_cohort$studydrug2)

# calculate predicted sglt2 benefit (absolute risk reduction = ARR):
# ARR = S0(t)^HR - S0(t)
trial_hr_kf_sglt2i <- 0.62

noncal_cohort <- noncal_cohort %>% 
  mutate(ckdpc_40egfr_cal=(100-ckdpc_40egfr_score_cal)/100,
         ckdpc_40egfr_cal_sglt2i=ckdpc_40egfr_cal^trial_hr_kf_sglt2i,
         ckdpc_40egfr_sglt2i_benefit=ckdpc_40egfr_cal_sglt2i - ckdpc_40egfr_cal)

# calculate predicted NNT = 1/ARR
noncal_cohort  <- noncal_cohort %>%
  mutate(nnt_predicted = 1/(ckdpc_40egfr_sglt2i_benefit))

# make separate datasets for those with macroalbuminuria and those with eGFR below 60
noncal_cohort_macroalb <- noncal_cohort %>% filter(macroalbuminuria == T)
noncal_cohort <- noncal_cohort %>% filter(macroalbuminuria == F)

noncal_cohort_lowegfr <- noncal_cohort %>% filter(egfr_below_60 == T)
noncal_cohort <- noncal_cohort %>% filter(egfr_below_60 == F)

# noncal_cohort_lowegfr <- noncal_cohort %>% mutate(egfr_below_60=ifelse(preegfr >=60, F, T)) %>% filter(egfr_below_60 == T)
# noncal_cohort <- noncal_cohort %>% mutate(egfr_below_60=ifelse(preegfr >=60, F, T)) %>%  filter(egfr_below_60 == F)
# 
# noncal_cohort_lowegfr <- noncal_cohort %>% mutate(egfr_below_60=ifelse(preckdstage %in% c("stage_1", "stage_2"), F, T)) %>% filter(egfr_below_60 == T)
# noncal_cohort <- noncal_cohort %>% mutate(egfr_below_60=ifelse(preckdstage %in% c("stage_1", "stage_2"), F, T)) %>%  filter(egfr_below_60 == F)


today <- as.character(Sys.Date(), format="%Y%m%d")

# set parameters for all graphs so they don't have to be specified in individual graph calls
cols_fig <- cols[names(cols) %in% noncal_cohort$studydrug2]
cols_fig <- cols_fig[order(names(cols_fig))]
font_subtitle <- 14
font_legend <- 10
font_x <- font_y <- 12
font_risktable <- 3.5
font_risktable_title <- 12
risktable_height <- 0.15
font_risktable_text <- 3.5
limit_y <- c(0,0.5)
zoom_y <- c(0,0.05)
zoom_y_lowegfr <- c(0,0.1)
zoom_y_macroalb <- c(0,0.2)
zoom_x <- c(0,3)
break_y <- 0.02
legend_in <- c(0.40, 0.75)
legend_out <- c(100,-100)

# choose threshold for acceptable NNT for treatment for display items later on
# based on the literature, a risk threshold of 5-year NNT = 30 would be acceptable to on average 71% of patients
# this would work out as a 3-year NNT of around 50.
# this would be correspond to a risk threshold of around 5.35431%
# however this would be too few persons in the high-risk albumin group
# nnt_threshold <- 50
# risk_threshold_1 <- 2.5

# therefore we are trying it with a different risk threshold
# a risk of 2.654% would correspond to a 3-year NNT of around 100.
# nnt_threshold <- 100
# risk_threshold_1 <- 2.654

# a NNT of 75 would correspond to a risk of 3.548455%
# nnt_threshold <- 75
# risk_threshold_1 <- 3.548455

# We will now define risk group by top 5% in their category

risk_threshold <- .90

risk_threshold_1 <- noncal_cohort[noncal_cohort$albuminuria == T,]$ckdpc_40egfr_score_cal %>% quantile(risk_threshold) %>% as.numeric()
risk_threshold_2 <- noncal_cohort[noncal_cohort$albuminuria == F,]$ckdpc_40egfr_score_cal %>% quantile(risk_threshold) %>% as.numeric()

high_risk_cat <- paste(c("CKD-PC risk score ≥", risk_threshold*100, "th percentile"), collapse = "")
low_risk_cat <- paste(c("CKD-PC risk score 0-", risk_threshold*100, "th percentile"), collapse = "")

# high_risk_cat <- paste(c("risk ≥ ", risk_threshold_1, "%"), collapse = "")
# low_risk_cat <- paste(c("risk < ", risk_threshold_1, "%"), collapse = "")
# high_risk_cat <- paste(c("NNT ≤ ", nnt_threshold, ""), collapse = "")
# low_risk_cat <- paste(c("NNT > ", nnt_threshold, ""), collapse = "")

# alb_and_high_risk <- paste(c("Microalbuminuria and high predicted risk"), collapse = "")
# alb_and_low_risk <- paste(c("Microalbuminuria and low predicted risk"), collapse = "")
# nonalb_low_risk <- paste(c("No albuminuria and low predicted risk"), collapse = "")

# make risk groups based on risk cutoff
noncal_cohort <- noncal_cohort %>% mutate(
  ckd40_risk_group = ifelse(albuminuria == T, ifelse(
    ckdpc_40egfr_score_cal < (risk_threshold_1), low_risk_cat, high_risk_cat),
    ifelse(ckdpc_40egfr_score_cal < (risk_threshold_2), low_risk_cat, high_risk_cat))
  ,
  ckd40_risk_group = factor(ckd40_risk_group, levels = c(low_risk_cat, high_risk_cat))
)

# make risk groups in low egfr group
risk_threshold_3 <- noncal_cohort_lowegfr[noncal_cohort_lowegfr$albuminuria == F,]$ckdpc_40egfr_score_cal %>% quantile(risk_threshold) %>% as.numeric()

noncal_cohort_lowegfr <- noncal_cohort_lowegfr %>% mutate(
  ckd40_risk_group = ifelse(albuminuria == T, high_risk_cat,
                            ifelse(ckdpc_40egfr_score_cal < (risk_threshold_3), low_risk_cat, high_risk_cat))
  ,
  ckd40_risk_group = factor(ckd40_risk_group, levels = c(low_risk_cat, high_risk_cat))
)
# noncal_cohort <- noncal_cohort %>% mutate(risk_group=ifelse(albuminuria==T, ifelse(
#   ckd40_risk_group==low_risk_cat,
#   alb_and_low_risk, alb_and_high_risk),
#   ifelse(
#     ckd40_risk_group==low_risk_cat,
#     nonalb_low_risk, "X"
#   )))
# print(paste0(c("Only ", nrow(noncal_cohort[noncal_cohort$risk_group == "X",]), " subjects had no albuminuria and predicted NNT ≤ ", nnt_threshold, ". These subjects are therefore not shown (removed from current dataset)"), collapse = ""))

# 
# noncal_cohort <- noncal_cohort %>% filter(!risk_group == "X") %>% mutate(
#   risk_group = factor(risk_group, levels = c(alb_and_low_risk,
#                                              alb_and_high_risk,
#                                              nonalb_low_risk))
# )

############################1 INCIDENCE CURVES IN PRESERVED EGFR################################################################

# add number per group and events
pred40 <- noncal_cohort %>%
  group_by(albuminuria, ckd40_risk_group 
  ) %>%
  summarise(mean_predicted_benefit = mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            mean_predicted_background_risk = mean(1-ckdpc_40egfr_cal, na.rm=T),
            event_count=round(sum(ckd_egfr40_censvar)/n.imp, 0),
            drug_count=round(n()/n.imp, 0),
            time=round(median(ckd_egfr40_censtime_yrs), 2)) %>%
  mutate(events_perc=round(event_count*100/drug_count, 1),
         events=paste0(event_count, " (", events_perc, "%)"))  %>%
  mutate(nnt_predicted=1/(mean_predicted_benefit)) %>%
  select(-c(event_count, events_perc))

# we also want to add the observed background risk 
# to do that we need to run cox models and take these estimates (as counterfactuals)

n.groups <- nlevels(as.factor(noncal_cohort$albuminuria))*nlevels(noncal_cohort$ckd40_risk_group)
# run cox model in each imputed dataset, extract estimates, and pool results
survival_est <- list()

for (i in 1:n.imp) {
  model <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2*ckd40_risk_group*albuminuria, 
               data=noncal_cohort[noncal_cohort$.imp == i,], x=TRUE, y=TRUE, surv=TRUE)
  
  
  survival_est[[i]] <- survest(model, 
                               newdata=expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
                                                   ckd40_risk_group=levels(noncal_cohort$ckd40_risk_group),
                                                   albuminuria=c(T,F)),
                               times=5)
}

# now we will create empty vectors for the pooled benefits per category

surv <- rep(NA, n.groups*nlevels(as.factor(noncal_cohort$studydrug2)))
se <- rep(NA, n.groups*nlevels(as.factor(noncal_cohort$studydrug2)))

# for every decile per group we will extract the n.imp imputed values and pool these, then store these in the vectors created above

for (k in 1:(n.groups*nlevels(as.factor(noncal_cohort$studydrug2)))) {
  
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

# create data frame
obs <- cbind(surv, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
                               ckd40_risk_group=levels(noncal_cohort$ckd40_risk_group),
                               albuminuria=c(T,F)))

obs <- obs %>%
  pivot_wider(id_cols=c(ckd40_risk_group,albuminuria), values_from=surv, names_from=studydrug2) %>%
  mutate(surv_diff=SGLT2i-`DPP4i/SU`,
         surv=1-`DPP4i/SU`) %>%
  select(albuminuria,ckd40_risk_group, surv, surv_diff)

se <- cbind(se, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
                            ckd40_risk_group=levels(noncal_cohort$ckd40_risk_group),
                            albuminuria=c(T,F)))

se <- se %>%
  pivot_wider(id_cols=c(ckd40_risk_group,albuminuria), values_from=se, names_from=studydrug2) %>%
  mutate(se=sqrt((SGLT2i^2)+(`DPP4i/SU`^2))) %>%
  select(albuminuria, ckd40_risk_group, se)

obs <- obs %>%
  inner_join(se, by=c("ckd40_risk_group", "albuminuria")) %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

dpp4isu_events <- noncal_cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_group, albuminuria) %>%
  summarise(`DPP4i/SU`=round(n()/n.imp,0))

sglt2i_events <- noncal_cohort %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
  group_by(ckd40_risk_group, albuminuria) %>%
  summarise(SGLT2i=round(n()/n.imp,0))

events_table <- data.frame(t(dpp4isu_events %>%
                               inner_join(sglt2i_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="ckd40_risk_group")

obs40 <- obs %>% mutate(ckd40_risk_group=as.factor(ckd40_risk_group))
obs_v_pred40 <- pred40 %>% inner_join(obs40, by=c("ckd40_risk_group", "albuminuria")) %>%
  mutate(nnt_observed=1/surv_diff) %>% select(-se)


#### unadjusted curves:


surv_alb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                            data = noncal_cohort[noncal_cohort$.imp == 1 &
                                                   noncal_cohort$albuminuria == T & 
                                                   noncal_cohort$ckd40_risk_group == low_risk_cat,],
                            conf.type = "log", conf.int = 0.95)

surv_alb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                             data = noncal_cohort[noncal_cohort$.imp == 1 &
                                                    noncal_cohort$albuminuria == T &
                                                    noncal_cohort$ckd40_risk_group == high_risk_cat,],
                             conf.type = "log", conf.int = 0.95)

surv_noalb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                              data = noncal_cohort[noncal_cohort$.imp == 1 &
                                                     noncal_cohort$albuminuria == F &
                                                     noncal_cohort$ckd40_risk_group == low_risk_cat,],
                              conf.type = "log", conf.int = 0.95)

surv_noalb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                               data = noncal_cohort[noncal_cohort$.imp == 1 &
                                                      noncal_cohort$albuminuria == F &
                                                      noncal_cohort$ckd40_risk_group == high_risk_cat,],
                               conf.type = "log", conf.int = 0.95)


list_plots_unadjusted <- list()

list_plots_unadjusted[[1]] <- ggsurvplot(
  fit = surv_alb_lowrisk,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1 &
                         noncal_cohort$albuminuria == T &
                         noncal_cohort$ckd40_risk_group == low_risk_cat,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  
  legend = legend_in,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "",
  ylab = "40% decline in eGFR / ESKD",
  title = paste0(c("Microalbuminuria and ", low_risk_cat), collapse="")
)

list_plots_unadjusted[[3]] <- ggsurvplot(
  fit = surv_alb_highrisk,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1 &
                         noncal_cohort$albuminuria == T &
                         noncal_cohort$ckd40_risk_group == high_risk_cat,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "",
  ylab = "",
  title = paste0(c("Microalbuminuria and ", high_risk_cat), collapse="")
) 

list_plots_unadjusted[[2]] <- ggsurvplot(
  fit = surv_noalb_lowrisk,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1 &
                         noncal_cohort$albuminuria == F &
                         noncal_cohort$ckd40_risk_group == low_risk_cat,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "40% decline in eGFR / ESKD",
  title = paste0(c("No albuminuria and ", low_risk_cat), collapse="")
) 


list_plots_unadjusted[[4]] <- ggsurvplot(
  fit = surv_noalb_highrisk,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1 &
                         noncal_cohort$albuminuria == F &
                         noncal_cohort$ckd40_risk_group == high_risk_cat,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "40% decline in eGFR / ESKD",
  title = paste0(c("No albuminuria and ", high_risk_cat), collapse="")
) 


list_plots_unadjusted[[1]]$plot <- list_plots_unadjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_unadjusted[[1]]$table <- list_plots_unadjusted[[1]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[1]]$cumevents <- list_plots_unadjusted[[1]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))

list_plots_unadjusted[[2]]$plot <- list_plots_unadjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_unadjusted[[2]]$table <- list_plots_unadjusted[[2]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[2]]$cumevents <- list_plots_unadjusted[[2]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))

list_plots_unadjusted[[3]]$plot <- list_plots_unadjusted[[3]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_unadjusted[[3]]$table <- list_plots_unadjusted[[3]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[3]]$cumevents <- list_plots_unadjusted[[3]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))

list_plots_unadjusted[[4]]$plot <- list_plots_unadjusted[[4]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_unadjusted[[4]]$table <- list_plots_unadjusted[[4]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[4]]$cumevents <- list_plots_unadjusted[[4]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))

# arrange_ggsurvplots(list_plots_unadjusted, print = T, nrow = 1, ncol = 3)

# weighted curves:
surv_alb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                            data = noncal_cohort[noncal_cohort$.imp == 1 &
                                                   noncal_cohort$albuminuria == T & 
                                                   noncal_cohort$ckd40_risk_group == low_risk_cat,],
                            weights = overlap2,
                            conf.type = "log", conf.int = 0.95)

surv_alb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                             data = noncal_cohort[noncal_cohort$.imp == 1 &
                                                    noncal_cohort$albuminuria == T &
                                                    noncal_cohort$ckd40_risk_group == high_risk_cat,],
                             weights = overlap2,
                             conf.type = "log", conf.int = 0.95)

surv_noalb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                              data = noncal_cohort[noncal_cohort$.imp == 1 &
                                                     noncal_cohort$albuminuria == F &
                                                     noncal_cohort$ckd40_risk_group == low_risk_cat,],
                              weights = overlap2,
                              conf.type = "log", conf.int = 0.95)

surv_noalb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                               data = noncal_cohort[noncal_cohort$.imp == 1 &
                                                      noncal_cohort$albuminuria == F &
                                                      noncal_cohort$ckd40_risk_group == high_risk_cat,],
                               weights = overlap2,
                               conf.type = "log", conf.int = 0.95)

list_plots_adjusted <- list()

list_plots_adjusted[[1]] <- ggsurvplot(
  fit = surv_alb_lowrisk,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1 &
                         noncal_cohort$albuminuria == T &
                         noncal_cohort$ckd40_risk_group == low_risk_cat,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  
  legend = legend_in,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "",
  ylab = "40% decline in eGFR / ESKD",
  title = paste0(c("eGFR ≥60mL/min and microalbuminuria"), collapse=""),
  subtitle = paste0(c(
    low_risk_cat, 
    #" (mean ", round(obs_v_pred40[obs_v_pred40$albuminuria == T &
    #                                        obs_v_pred40$ckd40_risk_group == low_risk_cat
    #                                               ,]$mean_predicted_background_risk*100, 1), 
    "; 3-year NNT ",
    round(1/(summary(surv_alb_lowrisk, times = 3)$surv[2]-summary(surv_alb_lowrisk, times = 3)$surv[1]), 0)
  ),collapse="")
)

list_plots_adjusted[[3]] <- ggsurvplot(
  fit = surv_alb_highrisk,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1 &
                         noncal_cohort$albuminuria == T &
                         noncal_cohort$ckd40_risk_group == high_risk_cat,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "",
  ylab = "",
  title = ""
  # paste0(c("Microalbuminuria"), collapse="")
  ,
  subtitle = paste0(c(
    high_risk_cat, 
    #" (mean ", round(obs_v_pred40[obs_v_pred40$albuminuria == T &
    #                                               obs_v_pred40$ckd40_risk_group == high_risk_cat
    #                                             ,]$mean_predicted_background_risk*100, 1), 
   "; 3-year NNT ",
    round(1/(summary(surv_alb_highrisk, times = 3)$surv[2]-summary(surv_alb_highrisk, times = 3)$surv[1]), 0)
  ),collapse="")
) 

list_plots_adjusted[[2]] <- ggsurvplot(
  fit = surv_noalb_lowrisk,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1 &
                         noncal_cohort$albuminuria == F &
                         noncal_cohort$ckd40_risk_group == low_risk_cat,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "40% decline in eGFR / ESKD",
  title = paste0(c("eGFR ≥60mL/min and no albuminuria"), collapse=""),
  subtitle = paste0(c(
    low_risk_cat, 
    #" (mean ", round(obs_v_pred40[obs_v_pred40$albuminuria == F &
    #                                            obs_v_pred40$ckd40_risk_group == low_risk_cat
    #                                            ,]$mean_predicted_background_risk*100, 1), "%);
    "; 3-year NNT ",
    round(1/(summary(surv_noalb_lowrisk, times = 3)$surv[2]-summary(surv_noalb_lowrisk, times = 3)$surv[1]), 0)
  ),collapse="")
) 

list_plots_adjusted[[4]] <- ggsurvplot(
  fit = surv_noalb_highrisk,
  fun = "cumhaz",
  data = noncal_cohort[noncal_cohort$.imp == 1 &
                         noncal_cohort$albuminuria == F &
                         noncal_cohort$ckd40_risk_group == high_risk_cat,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "",
  title = ""
  # paste0(c("No albuminuria"), collapse="")
  ,
  subtitle = paste0(c(
    high_risk_cat, 
    #" (mean ", round(obs_v_pred40[obs_v_pred40$albuminuria == F &
    #                                               obs_v_pred40$ckd40_risk_group == high_risk_cat
    #                                             ,]$mean_predicted_background_risk*100, 1), 
    "; 3-year NNT ",
    round(1/(summary(surv_noalb_highrisk, times = 3)$surv[2]-summary(surv_noalb_highrisk, times = 3)$surv[1]), 0)
    
  ),collapse="")
) 


list_plots_adjusted[[1]]$plot <- list_plots_adjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_adjusted[[1]]$table <- list_plots_unadjusted[[1]]$table 
list_plots_adjusted[[1]]$cumevents <- list_plots_unadjusted[[1]]$cumevents 

list_plots_adjusted[[2]]$plot <- list_plots_adjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_adjusted[[2]]$table <- list_plots_unadjusted[[2]]$table
list_plots_adjusted[[2]]$cumevents <- list_plots_unadjusted[[2]]$cumevents

list_plots_adjusted[[3]]$plot <- list_plots_adjusted[[3]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_adjusted[[3]]$table <- list_plots_unadjusted[[3]]$table
list_plots_adjusted[[3]]$cumevents <- list_plots_unadjusted[[3]]$cumevents

list_plots_adjusted[[4]]$plot <- list_plots_adjusted[[4]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_adjusted[[4]]$table <- list_plots_unadjusted[[4]]$table
list_plots_adjusted[[4]]$cumevents <- list_plots_unadjusted[[4]]$cumevents


arrange_ggsurvplots(list_plots_adjusted, print = T, nrow = 2, ncol = 2)

############################2 INCIDENCE CURVE IN MACROALBUMINURIA################################################################


# # do same as above for macroalbuminuria group (2 subgroups)
# 
# noncal_cohort_macroalb <- noncal_cohort_macroalb %>% mutate(egfr_below_60=ifelse(preckdstage %in% c("stage_1", "stage_2")& preegfr >=60, F, T))
# 
# 
# pred40_macroalb <- noncal_cohort_macroalb %>%
#   group_by(egfr_below_60) %>%
#   summarise(mean_predicted_benefit = mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
#             mean_predicted_background_risk = mean(1-ckdpc_40egfr_cal, na.rm=T),
#             event_count=round(sum(ckd_egfr40_censvar)/n.imp, 0),
#             drug_count=round(n()/n.imp, 0),
#             time=round(median(ckd_egfr40_censtime_yrs), 2)) %>%
#   mutate(events_perc=round(event_count*100/drug_count, 1),
#          events=paste0(event_count, " (", events_perc, "%)"))  %>%
#   mutate(nnt_predicted=1/(mean_predicted_benefit)) %>%
#   select(-c(event_count, events_perc))
# 
# # we also want to add the observed background risk
# 
# # run cox model in each imputed dataset, extract estimates, and pool results
# survival_est <- list()
# 
# for (i in 1:n.imp) {
#   model <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2*egfr_below_60, 
#                data=noncal_cohort_macroalb[noncal_cohort_macroalb$.imp == i,], x=TRUE, y=TRUE, surv=TRUE)
#   
#   
#   survival_est[[i]] <- survest(model, 
#                                newdata=expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"),
#                                                    egfr_below_60=c(T,F)),
#                                times=5)
# }
# 
# # now we will create empty vectors for the pooled benefits per category
# 
# surv <- rep(NA, nlevels(as.factor(noncal_cohort_macroalb$egfr_below_60))* nlevels(as.factor(noncal_cohort$studydrug2)))
# se <- rep(NA, nlevels(as.factor(noncal_cohort_macroalb$egfr_below_60))* nlevels(as.factor(noncal_cohort$studydrug2)))
# 
# # for every decile per group we will extract the n.imp imputed values and pool these, then store these in the vectors created above
# 
# for (k in 1:(nlevels(as.factor(noncal_cohort_macroalb$egfr_below_60))*nlevels(as.factor(noncal_cohort$studydrug2)))) {
#   
#   #create empty vectors to store the imputed values for every individual decile 
#   surv_temp=paste0("surv_", k)
#   se_temp=paste0("surv_se_", k)
#   
#   assign(surv_temp, rep(NA, n.imp))
#   assign(se_temp, rep(NA, n.imp))
#   
#   #store the imputed values in this vector
#   for (i in 1:n.imp) {
#     
#     #surv_temp[i] <- survival_est[[i]]$surv[k]
#     eval(str2lang(paste0(surv_temp, "[", i, "]", " <- survival_est", "[[", i, "]]$surv[", k, "]")))
#     
#     #se_temp[i] <- survival_est[[i]]$std.error[k]
#     eval(str2lang(paste0(se_temp, "[", i, "]", " <- survival_est", "[[", i, "]]$std.err[", k, "]")))
#     
#   }
#   # pool the imputed values and store the pooled value in the main vector
#   surv[k] <- mean(eval(parse(text = surv_temp)))
#   w <- mean(eval(parse(text = se_temp))^2)
#   b <- var(eval(parse(text = surv_temp)))
#   t.var <- w + (1+1/n.imp)*b
#   se[k] <- sqrt(t.var)
# }
# 
# # create data frame
# obs <- cbind(surv, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"),
#                                egfr_below_60=c(T,F)))
# 
# obs <- obs %>%
#   pivot_wider(id_cols=egfr_below_60, values_from=surv, names_from=studydrug2) %>%
#   mutate(surv_diff=SGLT2i-`DPP4i/SU`,
#          surv=1-`DPP4i/SU`) %>%
#   select(egfr_below_60, surv, surv_diff)
# 
# se <- cbind(se, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"),
#                             egfr_below_60=c(T,F)))
# 
# se <- se %>%
#   pivot_wider(id_cols=egfr_below_60, values_from=se, names_from=studydrug2) %>%
#   mutate(se=sqrt((SGLT2i^2)+(`DPP4i/SU`^2))) %>%
#   select(egfr_below_60, se)
# 
# obs <- obs %>%
#   inner_join(se, by = "egfr_below_60") %>%
#   mutate(lower_ci=surv_diff-(1.96*se),
#          upper_ci=surv_diff+(1.96*se))
# 
# dpp4isu_events <- noncal_cohort_macroalb %>%
#   filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
#   group_by(egfr_below_60) %>%
#   summarise(`DPP4i/SU`=round(n()/n.imp,0))
# 
# sglt2i_events <- noncal_cohort_macroalb %>%
#   filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
#   group_by(egfr_below_60) %>%
#   summarise(SGLT2i=round(n()/n.imp,0))
# 
# events_table <- data.frame(t(dpp4isu_events %>%
#                                cbind(sglt2i_events))) %>%
#   rownames_to_column() %>%
#   filter(rowname!="egfr_below_60")
# 
# obs40_macroalb <- obs 
# obs_v_pred40_macroalb <- pred40_macroalb %>% inner_join(obs40_macroalb, by="egfr_below_60") %>%
#   mutate(nnt_observed=1/surv_diff) %>% select(-se)
# 
# ## unadjusted curve for macroalbuminuria
# 
# surv_macroalb_presegfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
#                                   data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp ==1 &
#                                                                   noncal_cohort_macroalb$egfr_below_60 == F,],
#                                   conf.type = "log", conf.int = 0.95)
# 
# surv_macroalb_lowegfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
#                                  data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp ==1 &
#                                                                  noncal_cohort_macroalb$egfr_below_60 == T,],
#                                  conf.type = "log", conf.int = 0.95)
# 
# 
# list_plots_macroalb_unadjusted <- list()
# 
# list_plots_macroalb_unadjusted[[1]] <- ggsurvplot(
#   fit = surv_macroalb_presegfr,
#   fun = "cumhaz",
#   data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp == 1 &
#                                   noncal_cohort_macroalb$egfr_below_60 == F,],
#   palette = unname(cols_fig),
#   color = "studydrug2",
#   conf.int = T,
#   legend.title = "",
#   legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
#   font.legend = font_legend,
#   font.title = "bold",
#   font.subtitle = font_subtitle,
#   font.x = font_x,
#   font.y = font_y,
#   legend = legend_out,
#   break.y.by = break_y,
#   ylim = limit_y,
#   censor.size = 1,
#   surv.scale = "percent",
#   xlab = "Years",
#   ylab = "",
#   risk.table = T,
#   cumevents = T,
#   tables.height = risktable_height,
#   fontsize = font_risktable,
#   tables.y.text = F,
#   tables.y.text.col = T,
#   tables.theme = theme_cleantable(),
#   title = paste0(c("Macroalbuminuria and eGFR ≥ 60 mL/min; 3-year NNT ", 
#                    round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == F,]$nnt_predicted)), collapse=""),
#   subtitle = paste0(c("Mean background risk ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == F,]$surv*100, 1), 
#                       "% (predicted ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == F,]$mean_predicted_background_risk*100, 1), 
#                       "%)"),collapse="")
# )
# 
# list_plots_macroalb_unadjusted[[2]] <- ggsurvplot(
#   fit = surv_macroalb_lowegfr,
#   fun = "cumhaz",
#   data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp == 1 &
#                                   noncal_cohort_macroalb$egfr_below_60 == T,],
#   palette = unname(cols_fig),
#   color = "studydrug2",
#   conf.int = T,
#   legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
#   font.legend = font_legend,
#   font.title = "bold",
#   font.subtitle = font_subtitle,
#   font.x = font_x,
#   font.y = font_y,
#   legend = legend_in,
#   break.y.by = break_y,
#   ylim = limit_y,
#   censor.size = 1,
#   surv.scale = "percent",
#   xlab = "Years",
#   ylab = "40% decline in eGFR / ESKD",
#   risk.table = T,
#   cumevents = T,
#   tables.height = risktable_height,
#   fontsize = font_risktable,
#   tables.y.text = F,
#   tables.y.text.col = T,
#   tables.theme = theme_cleantable(),
#   title = paste0(c("Macroalbuminuria and eGFR <60mL/min; 3-year NNT ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == T,]$nnt_predicted)), collapse=""),
#   subtitle = paste0(c("Mean background risk ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == T,]$surv*100, 1), 
#                       "% (predicted ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == T,]$mean_predicted_background_risk*100, 1), 
#                       "%)"),collapse="")
# )
# 
# list_plots_macroalb_unadjusted[[1]]$plot <- list_plots_macroalb_unadjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_macroalb)
# list_plots_macroalb_unadjusted[[1]]$table <- list_plots_macroalb_unadjusted[[1]]$table + theme(plot.title = element_text(size = font_risktable_title))
# list_plots_macroalb_unadjusted[[1]]$cumevents <- list_plots_macroalb_unadjusted[[1]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))
# 
# list_plots_macroalb_unadjusted[[2]]$plot <- list_plots_macroalb_unadjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_macroalb)
# list_plots_macroalb_unadjusted[[2]]$table <- list_plots_macroalb_unadjusted[[2]]$table + theme(plot.title = element_text(size = font_risktable_title))
# list_plots_macroalb_unadjusted[[2]]$cumevents <- list_plots_macroalb_unadjusted[[2]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))
# 
# # arrange_ggsurvplots(list_plots_macroalb_unadjusted, print = T, nrow = 1, ncol = 2)
# 
# 
# ## adjusted curves for macroalbuminuria
# 
# surv_macroalb_lowegfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
#                                  data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp ==1 &
#                                                                  noncal_cohort_macroalb$egfr_below_60 == T,],
#                                  weights = overlap2,
#                                  conf.type = "log", conf.int = 0.95)
# 
# surv_macroalb_presegfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
#                                   data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp ==1 &
#                                                                   noncal_cohort_macroalb$egfr_below_60 == F,],
#                                   weights = overlap2,
#                                   conf.type = "log", conf.int = 0.95)
# 
# list_plots_macroalb_adjusted <- list()
# 
# list_plots_macroalb_adjusted[[1]] <- ggsurvplot(
#   fit = surv_macroalb_presegfr,
#   fun = "cumhaz",
#   data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp == 1 &
#                                   noncal_cohort_macroalb$egfr_below_60 == F,],
#   palette = unname(cols_fig),
#   color = "studydrug2",
#   conf.int = T,  
#   legend.title = "",
#   legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
#   font.legend = font_legend,
#   font.title = "bold",
#   font.subtitle = font_subtitle,
#   font.x = font_x,
#   font.y = font_y,
#   legend = legend_in,
#   break.y.by = break_y,
#   ylim = limit_y,
#   censor.size = 1,
#   surv.scale = "percent",
#   xlab = "Years",
#   ylab = "40% decline in eGFR / ESKD",
#   risk.table = T,
#   cumevents = T,
#   tables.height = risktable_height,
#   fontsize = font_risktable,
#   tables.y.text = F,
#   tables.y.text.col = T,
#   tables.theme = theme_cleantable(),
#   title = paste0(c("Macroalbuminuria and eGFR ≥ 60 mL/min"), collapse=""),
#   subtitle = paste0(c("NNT ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == F,]$nnt_predicted), ", mean background risk ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == F,]$surv*100, 1), 
#                       "% (predicted ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == F,]$mean_predicted_background_risk*100, 1), 
#                       "%)"),collapse="")
# )
# 
# list_plots_macroalb_adjusted[[2]] <- ggsurvplot(
#   fit = surv_macroalb_lowegfr,
#   fun = "cumhaz",
#   data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp == 1 &
#                                   noncal_cohort_macroalb$egfr_below_60 == T,],
#   palette = unname(cols_fig),
#   color = "studydrug2",
#   conf.int = T,
#   legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
#   font.legend = font_legend,
#   font.title = "bold",
#   font.subtitle = font_subtitle,
#   font.x = font_x,
#   font.y = font_y,
#   legend = legend_out,
#   break.y.by = break_y,
#   ylim = limit_y,
#   censor.size = 1,
#   surv.scale = "percent",
#   xlab = "Years",
#   ylab = "",
#   risk.table = T,
#   cumevents = T,
#   tables.height = risktable_height,
#   fontsize = font_risktable,
#   tables.y.text = F,
#   tables.y.text.col = T,
#   tables.theme = theme_cleantable(),
#   title = paste0(c("Macroalbuminuria and eGFR <60mL/min"), collapse=""),
#   subtitle = paste0(c("NNT ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == T,]$nnt_predicted), ", mean background risk ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == T,]$surv*100, 1), 
#                       "% (predicted ", round(obs_v_pred40_macroalb[obs_v_pred40_macroalb$egfr_below_60 == T,]$mean_predicted_background_risk*100, 1), 
#                       "%)"),collapse="")
# )
# 
# list_plots_macroalb_adjusted[[1]]$plot <- list_plots_macroalb_adjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_macroalb)
# list_plots_macroalb_adjusted[[1]]$table <- list_plots_macroalb_unadjusted[[1]]$table
# list_plots_macroalb_adjusted[[1]]$cumevents <- list_plots_macroalb_unadjusted[[1]]$cumevents
# 
# list_plots_macroalb_adjusted[[2]]$plot <- list_plots_macroalb_adjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_macroalb)
# list_plots_macroalb_adjusted[[2]]$table <- list_plots_macroalb_unadjusted[[2]]$table
# list_plots_macroalb_adjusted[[2]]$cumevents <- list_plots_macroalb_unadjusted[[2]]$cumevents
# 
# arrange_ggsurvplots(list_plots_macroalb_adjusted, print = T, nrow = 1, ncol = 2)


# do same as above for macroalbuminuria as 1 group

pred40_macroalb2 <- noncal_cohort_macroalb %>%
  summarise(mean_predicted_benefit = mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            mean_predicted_background_risk = mean(1-ckdpc_40egfr_cal, na.rm=T),
            event_count=round(sum(ckd_egfr40_censvar)/n.imp, 0),
            drug_count=round(n()/n.imp, 0),
            time=round(median(ckd_egfr40_censtime_yrs), 2)) %>%
  mutate(events_perc=round(event_count*100/drug_count, 1),
         events=paste0(event_count, " (", events_perc, "%)"))  %>%
  mutate(nnt_predicted=1/(mean_predicted_benefit)) %>%
  select(-c(event_count, events_perc))

# we also want to add the observed background risk

# run cox model in each imputed dataset, extract estimates, and pool results
survival_est <- list()

for (i in 1:n.imp) {
  model <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2, 
               data=noncal_cohort_macroalb[noncal_cohort_macroalb$.imp == i,], x=TRUE, y=TRUE, surv=TRUE)
  
  
  survival_est[[i]] <- survest(model, 
                               newdata=expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i")),
                               times=5)
}

# now we will create empty vectors for the pooled benefits per category

surv <- rep(NA, nlevels(as.factor(noncal_cohort$studydrug2)))
se <- rep(NA, nlevels(as.factor(noncal_cohort$studydrug2)))

# for every decile per group we will extract the n.imp imputed values and pool these, then store these in the vectors created above

for (k in 1:(nlevels(as.factor(noncal_cohort$studydrug2)))) {
  
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

# create data frame
obs <- cbind(surv, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i")))

obs <- obs %>%
  pivot_wider(values_from=surv, names_from=studydrug2) %>%
  mutate(surv_diff=SGLT2i-`DPP4i/SU`,
         surv=1-`DPP4i/SU`) %>%
  select(surv, surv_diff)

se <- cbind(se, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i")))

se <- se %>%
  pivot_wider(values_from=se, names_from=studydrug2) %>%
  mutate(se=sqrt((SGLT2i^2)+(`DPP4i/SU`^2))) %>%
  select(se)

obs <- obs %>%
  cbind(se) %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

dpp4isu_events <- noncal_cohort_macroalb %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  summarise(`DPP4i/SU`=round(n()/n.imp,0))

sglt2i_events <- noncal_cohort_macroalb %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
  summarise(SGLT2i=round(n()/n.imp,0))

events_table <- data.frame(t(dpp4isu_events %>%
                               cbind(sglt2i_events))) 

obs40_macroalb2 <- obs 
obs_v_pred40_macroalb2 <- pred40_macroalb2 %>% cbind(obs40_macroalb2) %>%
  mutate(nnt_observed=1/surv_diff) %>% select(-se)

## curve for macroalbuminuria - single group

surv_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                         data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp ==1,],
                         conf.type = "log", conf.int = 0.95)
plot_macroalb_unadjusted <- ggsurvplot(
  fit = surv_macroalb,
  fun = "cumhaz",
  data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp == 1,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  legend = legend_in,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  xlab = "Years",
  ylab = "40% decline in eGFR / ESKD",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  title = paste0(c("Macroalbuminuria; 3-year NNT ", 
                   round(obs_v_pred40_macroalb2$nnt_predicted)), collapse=""),
  subtitle = paste0(c("Mean background risk ", round(obs_v_pred40_macroalb2$surv*100, 1), 
                      "% (predicted ", round(obs_v_pred40_macroalb2$mean_predicted_background_risk*100, 1), 
                      "%)"),collapse="")
)
plot_macroalb_unadjusted$plot <- plot_macroalb_unadjusted$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
plot_macroalb_unadjusted$table <- plot_macroalb_unadjusted$table + theme(plot.title = element_text(size = font_risktable_title))
plot_macroalb_unadjusted$cumevents <- plot_macroalb_unadjusted$cumevents + theme(plot.title = element_text(size = font_risktable_title))


surv_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                         data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp ==1,],
                         weights = overlap2,
                         conf.type = "log", conf.int = 0.95)
plot_macroalb_adjusted <- ggsurvplot(
  fit = surv_macroalb,
  fun = "cumhaz",
  data = noncal_cohort_macroalb[noncal_cohort_macroalb$.imp == 1,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  xlab = "Years",
  ylab = "",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  title = paste0(c("Macroalbuminuria and any eGFR"), collapse=""),
  subtitle = paste0(c(
    "3-year NNT ", 
    round(1/(summary(surv_macroalb, times=3)$surv[2]-summary(surv_macroalb, times=3)$surv[1])
          , 0)
    # "Mean background risk ", round(obs_v_pred40_macroalb2$surv*100, 1), 
    #                   "% (predicted ", round(obs_v_pred40_macroalb2$mean_predicted_background_risk*100, 1), 
    #                   "%)"
  ),collapse="")
)
plot_macroalb_adjusted$plot <- plot_macroalb_adjusted$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
plot_macroalb_adjusted$table <- plot_macroalb_unadjusted$table + theme(plot.title = element_text(size = font_risktable_title))
plot_macroalb_adjusted$cumevents <- plot_macroalb_unadjusted$cumevents + theme(plot.title = element_text(size = font_risktable_title))

# print(plot_macroalb_adjusted)

############################3 INCIDENCE CURVE IN REDUCED EGFR################################################################


## same for low egfr
pred40_lowegfr <- noncal_cohort_lowegfr %>%
  group_by(albuminuria
  ) %>%
  summarise(mean_predicted_benefit = mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            mean_predicted_background_risk = mean(1-ckdpc_40egfr_cal, na.rm=T),
            event_count=round(sum(ckd_egfr40_censvar)/n.imp, 0),
            drug_count=round(n()/n.imp, 0),
            time=round(median(ckd_egfr40_censtime_yrs), 2)) %>%
  mutate(events_perc=round(event_count*100/drug_count, 1),
         events=paste0(event_count, " (", events_perc, "%)"))  %>%
  mutate(nnt_predicted=1/(mean_predicted_benefit)) %>%
  select(-c(event_count, events_perc))

#
survival_est <- list()

for (i in 1:n.imp) {
  model <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2*albuminuria, 
               data=noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == i,], x=TRUE, y=TRUE, surv=TRUE)
  
  
  survival_est[[i]] <- survest(model, 
                               newdata=expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
                                                   albuminuria=c(T,F)),
                               times=5)
}

# now we will create empty vectors for the pooled benefits per category

surv <- rep(NA, nlevels(as.factor(noncal_cohort_lowegfr$albuminuria))*nlevels(as.factor(noncal_cohort$studydrug2)))
se <- rep(NA, nlevels(as.factor(noncal_cohort_lowegfr$albuminuria))*nlevels(as.factor(noncal_cohort$studydrug2)))

# for every decile per group we will extract the n.imp imputed values and pool these, then store these in the vectors created above

for (k in 1:(nlevels(as.factor(noncal_cohort_lowegfr$albuminuria))*nlevels(as.factor(noncal_cohort_lowegfr$studydrug2)))) {
  
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

# create data frame
obs <- cbind(surv, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
                               albuminuria=c(T,F)))

obs <- obs %>%
  pivot_wider(id_cols=c(albuminuria), values_from=surv, names_from=studydrug2) %>%
  mutate(surv_diff=SGLT2i-`DPP4i/SU`,
         surv=1-`DPP4i/SU`) %>%
  select(albuminuria, surv, surv_diff)

se <- cbind(se, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
                            albuminuria=c(T,F)))

se <- se %>%
  pivot_wider(id_cols=c(albuminuria), values_from=se, names_from=studydrug2) %>%
  mutate(se=sqrt((SGLT2i^2)+(`DPP4i/SU`^2))) %>%
  select(albuminuria, se)

obs <- obs %>%
  inner_join(se, by=c("albuminuria")) %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

dpp4isu_events <- noncal_cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  group_by(albuminuria) %>%
  summarise(`DPP4i/SU`=round(n()/n.imp,0))

sglt2i_events <- noncal_cohort %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
  group_by(albuminuria) %>%
  summarise(SGLT2i=round(n()/n.imp,0))

events_table <- data.frame(t(dpp4isu_events %>%
                               inner_join(sglt2i_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="albuminuria")

obs40_lowegfr <- obs 
obs_v_pred40_lowegfr <- pred40_lowegfr %>% inner_join(obs40_lowegfr, by=c("albuminuria")) %>%
  mutate(nnt_observed=1/surv_diff) %>% select(-se)


## unadjusted curves for low egfr

surv_lowegfr_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
                                                               noncal_cohort_lowegfr$albuminuria == F,],
                                conf.type = "log", conf.int = 0.95)

surv_lowegfr_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                 data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
                                                                noncal_cohort_lowegfr$albuminuria == T,],
                                 conf.type = "log", conf.int = 0.95)

list_plots_lowegfr_unadjusted <- list()

list_plots_lowegfr_unadjusted[[1]] <- ggsurvplot(
  fit = surv_lowegfr_lowrisk,
  fun = "cumhaz",
  data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
                                 noncal_cohort_lowegfr$albuminuria == F,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  legend = legend_in,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "40% decline in eGFR / ESKD",
  title = paste0(c("eGFR <60mL/min and no albuminuria"), collapse=""),
  subtitle = paste0(c(
    "Mean background risk ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == F
                                                        ,]$surv*100, 1), 
    "% (predicted ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == F 
                                                ,]$mean_predicted_background_risk*100, 1), "%)"),collapse=""))

list_plots_lowegfr_unadjusted[[1]]$plot <- list_plots_lowegfr_unadjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_lowegfr_unadjusted[[1]]$table <- list_plots_lowegfr_unadjusted[[1]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_lowegfr_unadjusted[[1]]$cumevents <- list_plots_lowegfr_unadjusted[[1]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))



list_plots_lowegfr_unadjusted[[2]] <- ggsurvplot(
  fit = surv_lowegfr_highrisk,
  fun = "cumhaz",
  data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 &
                                 noncal_cohort_lowegfr$albuminuria == T,],
  #  facet.by = "ckd40_risk_group",
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",  
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "",
  title = paste0(c("eGFR <60mL/min and microalbuminuria"), collapse=""),
  subtitle = paste0(c( 
    "Mean background risk ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == T 
                                                        ,]$surv*100, 1), 
    "% (predicted ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == T 
                                                ,]$mean_predicted_background_risk*100, 1), "%)"),collapse=""))

list_plots_lowegfr_unadjusted[[2]]$plot <- list_plots_lowegfr_unadjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_lowegfr_unadjusted[[2]]$table <- list_plots_lowegfr_unadjusted[[2]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_lowegfr_unadjusted[[2]]$cumevents <- list_plots_lowegfr_unadjusted[[2]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


# arrange_ggsurvplots(list_plots_lowegfr_unadjusted, print = T, nrow = 1, ncol = 2)
## adjusted curves for egfr < 60


surv_lowegfr_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
                                                               noncal_cohort_lowegfr$albuminuria == F,],
                                weights = overlap2,
                                conf.type = "log", conf.int = 0.95)

surv_lowegfr_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                 data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
                                                                noncal_cohort_lowegfr$albuminuria == T,],
                                 weights = overlap2,
                                 conf.type = "log", conf.int = 0.95)

list_plots_lowegfr_adjusted <- list()

list_plots_lowegfr_adjusted[[1]] <- ggsurvplot(
  fit = surv_lowegfr_lowrisk,
  fun = "cumhaz",
  data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
                                 noncal_cohort_lowegfr$albuminuria == F,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  # short.panel.labs = T,
  legend = legend_in,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "40% decline in eGFR / ESKD",
  title = paste0(c("eGFR <60mL/min and no albuminuria"), collapse=""),
  subtitle = paste0(c(
    "3-year NNT ",
    round(1/(summary(surv_lowegfr_lowrisk, times=3)$surv[2]-summary(surv_lowegfr_lowrisk, times = 3)$surv[1]), 0)
    # "Mean background risk ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == F
    #                                                     ,]$surv*100, 1), 
    # "% (predicted ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == F 
    #                                             ,]$mean_predicted_background_risk*100, 1), "%)"
  ),collapse=""))

list_plots_lowegfr_adjusted[[2]] <- ggsurvplot(
  fit = surv_lowegfr_highrisk,
  fun = "cumhaz",
  data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 &
                                 noncal_cohort_lowegfr$albuminuria == T,],
  #  facet.by = "ckd40_risk_group",
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
  font.legend = font_legend,
  font.title = "bold",
  font.subtitle = font_subtitle,
  font.x = font_x,
  font.y = font_y,
  # short.panel.labs = T,
  legend = legend_out,
  break.y.by = break_y,
  ylim = limit_y,
  censor.size = 1,
  surv.scale = "percent",  
  risk.table = T,
  cumevents = T,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "",
  title = paste0(c("eGFR <60mL/min and microalbuminuria"), collapse="")
  ,
  subtitle = paste0(c(
    "3-year NNT ",
    round(1/(summary(surv_lowegfr_highrisk, times=3)$surv[2]-summary(surv_lowegfr_highrisk, times = 3)$surv[1])
          , 0)
    # "Mean background risk ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == T 
    #                                                     ,]$surv*100, 1), 
    # "% (predicted ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == T 
    #                                             ,]$mean_predicted_background_risk*100, 1), "%)"
  ),collapse=""))


list_plots_lowegfr_adjusted[[1]]$plot <- list_plots_lowegfr_adjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_lowegfr_adjusted[[1]]$table <- list_plots_lowegfr_unadjusted[[1]]$table
list_plots_lowegfr_adjusted[[1]]$cumevents <- list_plots_lowegfr_unadjusted[[1]]$cumevents

list_plots_lowegfr_adjusted[[2]]$plot <- list_plots_lowegfr_adjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_lowegfr_adjusted[[2]]$table <- list_plots_lowegfr_unadjusted[[2]]$table
list_plots_lowegfr_adjusted[[2]]$cumevents <- list_plots_lowegfr_unadjusted[[2]]$cumevents


# arrange_ggsurvplots(list_plots_lowegfr_adjusted, print = T, nrow = 1, ncol = 2)

# show all 3
list_plots_lowegfr_adjusted[[3]] <- plot_macroalb_adjusted
arrange_ggsurvplots(list_plots_lowegfr_adjusted, print = T, nrow = 1, ncol = 3)


############################4 INCIDENCE CURVE IN REDUCED EGFR WITHOUT ALBUMINURIA SPLIT BY RISK################################################################
# 
# ## low egfr group without albuminuria split according to risk:
# pred40_lowegfr2 <- noncal_cohort_lowegfr[noncal_cohort_lowegfr$albuminuria == F,] %>%
#   group_by(ckd40_risk_group
#   ) %>%
#   summarise(mean_predicted_benefit = mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
#             mean_predicted_background_risk = mean(1-ckdpc_40egfr_cal, na.rm=T),
#             event_count=round(sum(ckd_egfr40_censvar)/n.imp, 0),
#             drug_count=round(n()/n.imp, 0),
#             time=round(median(ckd_egfr40_censtime_yrs), 2)) %>%
#   mutate(events_perc=round(event_count*100/drug_count, 1),
#          events=paste0(event_count, " (", events_perc, "%)"))  %>%
#   mutate(nnt_predicted=1/(mean_predicted_benefit)) %>%
#   select(-c(event_count, events_perc))
# 
# #
# survival_est <- list()
# 
# for (i in 1:n.imp) {
#   model <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2*ckd40_risk_group, 
#                data=noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == i,], x=TRUE, y=TRUE, surv=TRUE)
#   
#   
#   survival_est[[i]] <- survest(model, 
#                                newdata=expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
#                                                    ckd40_risk_group=levels(noncal_cohort_lowegfr$ckd40_risk_group)),
#                                times=5)
# }
# 
# # now we will create empty vectors for the pooled benefits per category
# 
# surv <- rep(NA, nlevels(as.factor(noncal_cohort_lowegfr$ckd40_risk_group))*nlevels(as.factor(noncal_cohort$studydrug2)))
# se <- rep(NA, nlevels(as.factor(noncal_cohort_lowegfr$ckd40_risk_group))*nlevels(as.factor(noncal_cohort$studydrug2)))
# 
# # for every decile per group we will extract the n.imp imputed values and pool these, then store these in the vectors created above
# 
# for (k in 1:(nlevels(as.factor(noncal_cohort_lowegfr$ckd40_risk_group))*nlevels(as.factor(noncal_cohort_lowegfr$studydrug2)))) {
#   
#   #create empty vectors to store the imputed values for every individual decile 
#   surv_temp=paste0("surv_", k)
#   se_temp=paste0("surv_se_", k)
#   
#   assign(surv_temp, rep(NA, n.imp))
#   assign(se_temp, rep(NA, n.imp))
#   
#   #store the imputed values in this vector
#   for (i in 1:n.imp) {
#     
#     #surv_temp[i] <- survival_est[[i]]$surv[k]
#     eval(str2lang(paste0(surv_temp, "[", i, "]", " <- survival_est", "[[", i, "]]$surv[", k, "]")))
#     
#     #se_temp[i] <- survival_est[[i]]$std.error[k]
#     eval(str2lang(paste0(se_temp, "[", i, "]", " <- survival_est", "[[", i, "]]$std.err[", k, "]")))
#     
#   }
#   # pool the imputed values and store the pooled value in the main vector
#   surv[k] <- mean(eval(parse(text = surv_temp)))
#   w <- mean(eval(parse(text = se_temp))^2)
#   b <- var(eval(parse(text = surv_temp)))
#   t.var <- w + (1+1/n.imp)*b
#   se[k] <- sqrt(t.var)
# }
# 
# # create data frame
# obs <- cbind(surv, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
#                                ckd40_risk_group=levels(noncal_cohort_lowegfr$ckd40_risk_group)))
# 
# obs <- obs %>%
#   pivot_wider(id_cols=c(ckd40_risk_group), values_from=surv, names_from=studydrug2) %>%
#   mutate(surv_diff=SGLT2i-`DPP4i/SU`,
#          surv=1-`DPP4i/SU`) %>%
#   select(ckd40_risk_group, surv, surv_diff)
# 
# se <- cbind(se, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
#                             ckd40_risk_group=levels(noncal_cohort_lowegfr$ckd40_risk_group)))
# 
# se <- se %>%
#   pivot_wider(id_cols=c(ckd40_risk_group), values_from=se, names_from=studydrug2) %>%
#   mutate(se=sqrt((SGLT2i^2)+(`DPP4i/SU`^2))) %>%
#   select(ckd40_risk_group, se)
# 
# obs <- obs %>%
#   inner_join(se, by=c("ckd40_risk_group")) %>%
#   mutate(lower_ci=surv_diff-(1.96*se),
#          upper_ci=surv_diff+(1.96*se))
# 
# dpp4isu_events <- noncal_cohort %>%
#   filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
#   group_by(ckd40_risk_group) %>%
#   summarise(`DPP4i/SU`=round(n()/n.imp,0))
# 
# sglt2i_events <- noncal_cohort %>%
#   filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
#   group_by(ckd40_risk_group) %>%
#   summarise(SGLT2i=round(n()/n.imp,0))
# 
# events_table <- data.frame(t(dpp4isu_events %>%
#                                inner_join(sglt2i_events))) %>%
#   rownames_to_column() %>%
#   filter(rowname!="ckd40_risk_group")
# 
# obs40_lowegfr <- obs 
# obs_v_pred40_lowegfr2 <- pred40_lowegfr2 %>% inner_join(obs40_lowegfr, by=c("ckd40_risk_group")) %>%
#   mutate(nnt_observed=1/surv_diff) %>% select(-se)
# 
# 
# ## unadjusted curves for low egfr
# 
# surv_lowegfr_noalb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
#                                       data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
#                                                                      noncal_cohort_lowegfr$albuminuria == F &
#                                                                      noncal_cohort_lowegfr$ckd40_risk_group == low_risk_cat,],
#                                       conf.type = "log", conf.int = 0.95)
# 
# surv_lowegfr_noalb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
#                                        data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
#                                                                       noncal_cohort_lowegfr$albuminuria == F &
#                                                                       noncal_cohort_lowegfr$ckd40_risk_group == high_risk_cat,],
#                                        conf.type = "log", conf.int = 0.95)
# 
# list_plots_lowegfr_noalb_unadjusted <- list()
# 
# list_plots_lowegfr_noalb_unadjusted[[1]] <- ggsurvplot(
#   fit = surv_lowegfr_noalb_lowrisk,
#   fun = "cumhaz",
#   data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
#                                  noncal_cohort_lowegfr$albuminuria == F &
#                                  noncal_cohort_lowegfr$ckd40_risk_group == low_risk_cat,],
#   palette = unname(cols_fig),
#   color = "studydrug2",
#   conf.int = T,
#   legend.title = "",
#   legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
#   font.legend = font_legend,
#   font.title = "bold",
#   font.subtitle = font_subtitle,
#   font.x = font_x,
#   font.y = font_y,
#   legend = legend_in,
#   break.y.by = break_y,
#   ylim = limit_y,
#   censor.size = 1,
#   surv.scale = "percent",
#   risk.table = T,
#   cumevents = T,
#   tables.height = risktable_height,
#   fontsize = font_risktable,
#   tables.y.text = F,
#   tables.y.text.col = T,
#   tables.theme = theme_cleantable(),
#   xlab = "Years",
#   ylab = "40% decline in eGFR / ESKD",
#   title = paste0(c("eGFR <60mL/min and no albuminuria"), collapse=""),
#   # subtitle = paste0(c(
#   #   "Mean background risk ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == F
#   #                                                       ,]$surv*100, 1), 
#   #   "% (predicted ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == F 
#   #                                               ,]$mean_predicted_background_risk*100, 1), "%)"),collapse="")
# )
# 
# list_plots_lowegfr_noalb_unadjusted[[1]]$plot <- list_plots_lowegfr_noalb_unadjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
# list_plots_lowegfr_noalb_unadjusted[[1]]$table <- list_plots_lowegfr_noalb_unadjusted[[1]]$table + theme(plot.title = element_text(size = font_risktable_title))
# list_plots_lowegfr_noalb_unadjusted[[1]]$cumevents <- list_plots_lowegfr_noalb_unadjusted[[1]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))
# 
# 
# 
# list_plots_lowegfr_noalb_unadjusted[[2]] <- ggsurvplot(
#   fit = surv_lowegfr_noalb_highrisk,
#   fun = "cumhaz",
#   data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
#                                  noncal_cohort_lowegfr$albuminuria == F &
#                                  noncal_cohort_lowegfr$ckd40_risk_group == high_risk_cat,],
#   #  facet.by = "ckd40_risk_group",
#   palette = unname(cols_fig),
#   color = "studydrug2",
#   conf.int = T,
#   legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
#   font.legend = font_legend,
#   font.title = "bold",
#   font.subtitle = font_subtitle,
#   font.x = font_x,
#   font.y = font_y,
#   legend = legend_out,
#   break.y.by = break_y,
#   ylim = limit_y,
#   censor.size = 1,
#   surv.scale = "percent",  
#   risk.table = T,
#   cumevents = T,
#   tables.height = risktable_height,
#   fontsize = font_risktable,
#   tables.y.text = F,
#   tables.y.text.col = T,
#   tables.theme = theme_cleantable(),
#   xlab = "Years",
#   ylab = "",
#   title = paste0(c(""), collapse=""),
#   # subtitle = paste0(c( 
#   #   "Mean background risk ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == T 
#   #                                                       ,]$surv*100, 1), 
#   #   "% (predicted ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == T 
#   #                                               ,]$mean_predicted_background_risk*100, 1), "%)"),collapse="")
# )
# 
# list_plots_lowegfr_noalb_unadjusted[[2]]$plot <- list_plots_lowegfr_noalb_unadjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
# list_plots_lowegfr_noalb_unadjusted[[2]]$table <- list_plots_lowegfr_noalb_unadjusted[[2]]$table + theme(plot.title = element_text(size = font_risktable_title))
# list_plots_lowegfr_noalb_unadjusted[[2]]$cumevents <- list_plots_lowegfr_noalb_unadjusted[[2]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))
# 
# 
# # arrange_ggsurvplots(list_plots_lowegfr_noalb_unadjusted, print = T, nrow = 1, ncol = 2)
# ## adjusted curves for egfr < 60
# 
# 
# surv_lowegfr_noalb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
#                                       data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
#                                                                      noncal_cohort_lowegfr$albuminuria == F &
#                                                                      noncal_cohort_lowegfr$ckd40_risk_group == low_risk_cat,],
#                                       weights = overlap2,
#                                       conf.type = "log", conf.int = 0.95)
# 
# surv_lowegfr_noalb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
#                                        data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
#                                                                       noncal_cohort_lowegfr$albuminuria == F &
#                                                                       noncal_cohort_lowegfr$ckd40_risk_group == high_risk_cat,],
#                                        weights = overlap2,
#                                        conf.type = "log", conf.int = 0.95)
# 
# list_plots_lowegfr_noalb_adjusted <- list()
# 
# list_plots_lowegfr_noalb_adjusted[[1]] <- ggsurvplot(
#   fit = surv_lowegfr_noalb_lowrisk,
#   fun = "cumhaz",
#   data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
#                                  noncal_cohort_lowegfr$albuminuria == F &
#                                  noncal_cohort_lowegfr$ckd40_risk_group == low_risk_cat,],
#   palette = unname(cols_fig),
#   color = "studydrug2",
#   conf.int = T,
#   legend.title = "",
#   legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
#   font.legend = font_legend,
#   font.title = "bold",
#   font.subtitle = font_subtitle,
#   font.x = font_x,
#   font.y = font_y,
#   # short.panel.labs = T,
#   legend = legend_in,
#   break.y.by = break_y,
#   ylim = limit_y,
#   censor.size = 1,
#   surv.scale = "percent",
#   risk.table = T,
#   cumevents = T,
#   tables.height = risktable_height,
#   fontsize = font_risktable,
#   tables.y.text = F,
#   tables.y.text.col = T,
#   tables.theme = theme_cleantable(),
#   xlab = "Years",
#   ylab = "40% decline in eGFR / ESKD",
#   title = paste0(c("eGFR <60mL/min and no albuminuria"), collapse=""),
#   subtitle = paste0(c(
#     low_risk_cat,
#     "; 3-year NNT ",
#     round(1/(summary(surv_lowegfr_noalb_lowrisk, times=3)$surv[2]-summary(surv_lowegfr_noalb_lowrisk, times = 3)$surv[1]), 0)
#     # "Mean background risk ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == F
#     #                                                     ,]$surv*100, 1), 
#     # "% (predicted ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == F 
#     #                                             ,]$mean_predicted_background_risk*100, 1), "%)"
#   ),collapse=""))
# 
# list_plots_lowegfr_noalb_adjusted[[2]] <- ggsurvplot(
#   fit = surv_lowegfr_noalb_highrisk,
#   fun = "cumhaz",
#   data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$.imp == 1 & 
#                                  noncal_cohort_lowegfr$albuminuria == F &
#                                  noncal_cohort_lowegfr$ckd40_risk_group == high_risk_cat,],
#   #  facet.by = "ckd40_risk_group",
#   palette = unname(cols_fig),
#   color = "studydrug2",
#   conf.int = T,
#   legend.labs = c("DPP4 inhibitors/sulfonylureas", "SGLT2 inhibitors"),
#   font.legend = font_legend,
#   font.title = "bold",
#   font.subtitle = font_subtitle,
#   font.x = font_x,
#   font.y = font_y,
#   # short.panel.labs = T,
#   legend = legend_out,
#   break.y.by = break_y,
#   ylim = limit_y,
#   censor.size = 1,
#   surv.scale = "percent",  
#   risk.table = T,
#   cumevents = T,
#   tables.height = risktable_height,
#   fontsize = font_risktable,
#   tables.y.text = F,
#   tables.y.text.col = T,
#   tables.theme = theme_cleantable(),
#   xlab = "Years",
#   ylab = "",
#   title = paste0(c("eGFR <60mL/min and no albuminuria"), collapse="")
#   ,
#   subtitle = paste0(c(
#     high_risk_cat,
#     "; 3-year NNT ",
#     round(1/(summary(surv_lowegfr_noalb_highrisk, times=3)$surv[2]-summary(surv_lowegfr_noalb_highrisk, times = 3)$surv[1])
#           , 0)
#     # "Mean background risk ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == T 
#     #                                                     ,]$surv*100, 1), 
#     # "% (predicted ", round(obs_v_pred40_lowegfr[obs_v_pred40_lowegfr$albuminuria == T 
#     #                                             ,]$mean_predicted_background_risk*100, 1), "%)"
#   ),collapse=""))
# 
# 
# list_plots_lowegfr_noalb_adjusted[[1]]$plot <- list_plots_lowegfr_noalb_adjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
# list_plots_lowegfr_noalb_adjusted[[1]]$table <- list_plots_lowegfr_noalb_unadjusted[[1]]$table
# list_plots_lowegfr_noalb_adjusted[[1]]$cumevents <- list_plots_lowegfr_noalb_unadjusted[[1]]$cumevents
# 
# list_plots_lowegfr_noalb_adjusted[[2]]$plot <- list_plots_lowegfr_noalb_adjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
# list_plots_lowegfr_noalb_adjusted[[2]]$table <- list_plots_lowegfr_noalb_unadjusted[[2]]$table
# list_plots_lowegfr_noalb_adjusted[[2]]$cumevents <- list_plots_lowegfr_noalb_unadjusted[[2]]$cumevents
# 
# 
# arrange_ggsurvplots(list_plots_lowegfr_noalb_adjusted, print = T, nrow = 1, ncol = 2)
############################5 SUBGROUP TABLES################################################################


vars <- c("dstartdate_age", "malesex", "ethnicity_5cat", "imd2015_10",             # sociodemographic variables
          "prebmi", "preegfr", "uacr", "albuminuria",                              # vital signs and laboratory measurements
          "preldl", "prehba1c", "presbp", "predbp",   
          "dstartdate_dm_dur_all", "qrisk2_smoking_cat", "predrug_hypertension",   # comorbidities
          "predrug_af", "predrug_dka", "osteoporosis", 
          "predrug_acutepancreatitis", "predrug_falls", 
          "predrug_urinary_frequency", "predrug_volume_depletion", 
          "predrug_micturition_control", "predrug_dementia", "hosp_admission_prev_year",
          "initiation_year",
          "ncurrtx", "MFN", "TZD", "INS", "ACEi_or_ARB",                            # medications
          "cv_high_risk", "qrisk2_above_10_pct"                                     # CV risk
)

#categorical variables
factors <- c("malesex", "ethnicity_qrisk2", "qrisk2_smoking_cat", "predrug_hypertension", 
             "predrug_af", "predrug_dka", "osteoporosis", "predrug_acutepancreatitis", 
             "predrug_falls", "predrug_urinary_frequency", "predrug_volume_depletion", 
             "predrug_micturition_control", "predrug_dementia", "hosp_admission_prev_year",
             "initiation_year", "albuminuria", 
             "ncurrtx", "MFN", "TZD", "INS", "ACEi_or_ARB",
             "cv_high_risk", "qrisk2_above_10_pct")

nonnormal <- c("imd2015_10", "uacr", "dstartdate_dm_dur_all")

noncal_cohort <- noncal_cohort %>% mutate(risk_group = ifelse(
  albuminuria == T, ifelse(
    ckd40_risk_group == high_risk_cat, "alb_highrisk", "alb_lowrisk"),
  ifelse(ckd40_risk_group == high_risk_cat, "noalb_highrisk", "noalb_lowrisk")
)
)

table <- CreateTableOne(vars = vars, strata = "risk_group", data = noncal_cohort, 
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv
# write.csv2(tabforprint, file = paste0(today, "_ckd40_highrisk_table.csv"))

table1 <- obs_v_pred40 %>% select(albuminuria, ckd40_risk_group, 
                                  risk_predicted=mean_predicted_background_risk, 
                                  risk_observed=surv, 
                                  arr_predicted=mean_predicted_benefit, 
                                  arr_observed=surv_diff, nnt_predicted, nnt_observed) %>%
  mutate(
    risk_predicted = round(risk_predicted*100,1),
    risk_observed = round(risk_observed*100,1),
    arr_predicted =round(arr_predicted*100,1),
    arr_observed = round(arr_observed*100,1),
    nnt_predicted = round(nnt_predicted),
    nnt_observed = round(nnt_observed)
  )
# write.csv2(table1, file = paste0(today, "_preservedegfr_arr_nnt_table.csv"))


table <- CreateTableOne(vars = vars, data = noncal_cohort_macroalb, 
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)
# write.csv2(tabforprint, file = paste0(today, "_macroalb_table.csv"))

table <- CreateTableOne(vars = vars, strata = "albuminuria", data = noncal_cohort_lowegfr, 
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)
# write.csv2(tabforprint, file = paste0(today, "_lowegfr_table.csv"))

# table <- CreateTableOne(vars = vars, strata = "ckd40_risk_group", data = noncal_cohort_lowegfr[noncal_cohort_lowegfr$albuminuria == F,], 
#                         factorVars = factors, test = F)
# 
# tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)
# # write.csv2(tabforprint, file = paste0(today, "_lowegfr_table.csv"))

