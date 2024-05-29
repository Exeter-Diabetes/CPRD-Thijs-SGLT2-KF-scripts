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

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-04-30_t2d_ckdpc_recalibrated.Rda")

cohort$studydrug2 <- as.factor(cohort$studydrug2)

# calculate predicted sglt2 benefit (absolute risk reduction = ARR):
# ARR = S0(t)^HR - S0(t)
trial_hr_kf_sglt2i <- 0.62

cohort <- cohort %>% 
  mutate(ckdpc_40egfr_cal=(100-ckdpc_40egfr_score_cal)/100,
         ckdpc_40egfr_cal_sglt2i=ckdpc_40egfr_cal^trial_hr_kf_sglt2i,
         ckdpc_40egfr_sglt2i_benefit=ckdpc_40egfr_cal_sglt2i - ckdpc_40egfr_cal)

# calculate predicted NNT = 1/ARR
cohort  <- cohort %>%
  mutate(nnt_predicted = 1/(ckdpc_40egfr_sglt2i_benefit))

# decile plots indicate only top 10% in people with preserved eGFR and low-level albuminuria have a significant risk

risk_threshold <- .90

risk_threshold_1 <- cohort[cohort$risk_group == "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",]$ckdpc_40egfr_score_cal %>% quantile(risk_threshold) %>% as.numeric()

high_risk_cat <- paste(c("CKD-PC risk score ≥", risk_threshold*100, "th percentile"), collapse = "")
low_risk_cat <- paste(c("CKD-PC risk score 0-", risk_threshold*100, "th percentile"), collapse = "")

# make risk groups based on risk cutoff
cohort <- cohort %>% mutate(
  risk_group = ifelse(risk_group == "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol", 
                      ifelse(ckdpc_40egfr_score_cal < risk_threshold_1,
                             paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat),
                             paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat)),
                      as.character(risk_group)),
  risk_group = factor(risk_group, levels=c("eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
                                           paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat),
                                           paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat),
                                           "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                           "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                           "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                                           "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))
)

# covariates for multivariable adjustment
covariates <- "dstartdate_age + malesex + imd2015_10 + ethnicity_5cat + initiation_year + prebmi + prehba1c + pretotalcholesterol + preegfr + uacr + presbp + ckdpc_40egfr_score + ncurrtx + MFN + statin + INS + ACEi_or_ARB + smoking_status + dstartdate_dm_dur_all + predrug_hypertension + predrug_af + hosp_admission_prev_year"

today <- as.character(Sys.Date(), format="%Y%m%d")
save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated_with_riskgroup.Rda"))

# set parameters for all graphs so they don't have to be specified in individual graph calls
cols_fig <- cols[names(cols) %in% cohort$studydrug2]
cols_fig <- cols_fig[order(names(cols_fig))]
font_subtitle <- 14
font_legend <- 12
font_x <- font_y <- 12
font_risktable <- 3.5
font_risktable_title <- 12
risktable_height <- 0.15
font_risktable_text <- 3.5
limit_y <- c(0,0.5)
zoom_y <- c(0,0.05)
zoom_y_lowegfr <- c(0,0.2)
zoom_y_macroalb <- c(0,0.2)
zoom_x <- c(0,3)
break_y <- 0.02
legend_in <- c(0.40, 0.75)
legend_out <- c(100,-100)

############################1 CALIBRATION PLOTS OF PREDICTED VS OBSERVED BENEFITS################################################################

# save dataset for later (performing analysis in entire dataset exceeds memory limit)
cohort <- cohort %>%
  select(patid, .imp, risk_group, studydrug2, ckd_egfr40_censtime_yrs, 
         ckd_egfr40_censvar, macroalb_censtime_yrs, macroalb_censvar, all_of(unlist(strsplit(covariates, " \\+ "))))
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_all_imputations_calibrated.Rda"))

rm(list = setdiff(ls(), c("n.imp", "covariates")))

k <- "ckd_egfr40"

for (i in 1:n.imp) {
  print(paste0("Survival estimates for imputation ", i, " (SGLT2i)"))
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
  load("2024-04-30_all_imputations_calibrated.Rda")
  cohort <- cohort[cohort$.imp == i,]
  gc()
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")  
  f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", covariates))
  
  model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE)
  
  obs_SGLT2 <- cohort %>%
    mutate(studydrug_original=studydrug2,
           studydrug2="SGLT2i",
           rowno=row_number())
  
  print(paste0(c("Extracting survival estimates")))
  
  observed_sglt2 <- survfit(model, newdata=as.data.frame(obs_SGLT2)) %>%
    tidy() %>%
    filter(time==3) %>%
    pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
    select(group, estimate_sglt2=estimate, se_sglt2=std.error) %>%
    mutate(group=as.numeric(group)) %>%
    inner_join(obs_SGLT2, by=c("group"="rowno"))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
  
  today <- as.character(Sys.Date(), format="%Y%m%d")
  save(observed_sglt2, file=paste0(today, "_adjusted_surv_SGLT2i_imp.", i, ".Rda"))
  
  rm(list = setdiff(ls(), c("n.imp", "covariates", "k")))
}


for (i in 1:n.imp) {
  print(paste0("Survival estimates for imputation ", i, " (DPP4i/SU)"))
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
  load("2024-04-30_all_imputations_calibrated.Rda")
  cohort <- cohort[cohort$.imp == i,]
  gc()
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")  
  f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", covariates))
  
  model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE)
  
  obs_DPP4SU <- cohort %>%
    mutate(studydrug_original=studydrug2,
           studydrug2="DPP4i/SU",
           rowno=row_number())
  
  print(paste0(c("Extracting survival estimates")))
  
  observed_dpp4su <- survfit(model, newdata=as.data.frame(obs_DPP4SU)) %>%
    tidy() %>%
    filter(time==3) %>%
    pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
    select(group, estimate_dpp4su=estimate, se_dpp4su=std.error) %>%
    mutate(group=as.numeric(group)) %>%
    inner_join(obs_DPP4SU, by=c("group"="rowno"))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
  today <- as.character(Sys.Date(), format="%Y%m%d")
  save(observed_dpp4su, file=paste0(today, "_adjusted_surv_DPP4iSU_imp.", i, ".Rda"))
  
  rm(list = setdiff(ls(), c("n.imp", "covariates", "k")))
}

temp_sglt2 <- temp_dpp4su <- data.frame()

for (i in 1:n.imp) {
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
  load(paste0("2024-04-30_adjusted_surv_SGLT2i_imp.", i, ".Rda"))
  load(paste0("2024-04-30_adjusted_surv_DPP4iSU_imp.", i, ".Rda"))
  temp_sglt2 <- temp_sglt2 %>% rbind(observed_sglt2)
  temp_dpp4su <- temp_dpp4su %>% rbind(observed_dpp4su)
  rm(observed_sglt2)
  rm(observed_dpp4su)
}

benefits <- temp_sglt2 %>%
  select(group, .imp, estimate_sglt2, se_sglt2) %>%
  inner_join(temp_dpp4su, by=c("group", ".imp")) %>%
  select(-studydrug2) %>%
  mutate(survdiff=estimate_sglt2-estimate_dpp4su,
         studydrug2=studydrug_original) %>%
  select(-studydrug_original)

### Predicted same as previous
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-04-30_t2d_ckdpc_recalibrated_with_riskgroup.Rda")

cohort <- cohort %>% inner_join(benefits %>% 
                                  select(.imp, estimate_sglt2, se_sglt2, estimate_dpp4su, se_dpp4su, patid, 
                                         studydrug2, survdiff), 
                                by=c(".imp", "patid", "studydrug2"))

cohort$studydrug2 <- as.factor(cohort$studydrug2)

cohort$benefit_decile <- ntile(cohort$ckdpc_40egfr_sglt2i_benefit, n.quantiles)

obs_v_pred_for_plot <- cohort %>%
  group_by(benefit_decile) %>%
  summarise(median_predicted_benefit=median(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            mean_predicted_benefit=mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            mean_benefit=mean(survdiff),
            median_benefit=median(survdiff),
            lq_benefit=quantile(survdiff, prob=c(.25)),
            uq_benefit=quantile(survdiff, prob=c(.75)),
            se_benefit=mean(1/2*(se_sglt2 + se_dpp4su)),
            upper_ci=mean_benefit + 1.96*se_benefit,
            lower_ci=mean_benefit - 1.96*se_benefit)



empty_tick <- obs_v_pred_for_plot %>% 
  filter(benefit_decile==1) %>%
  mutate(benefit_decile=0, estimate=NA, lower_ci=NA, upper_ci=NA, predicted=NA)

## SGLT2i benefit predicted vs observed - in all patients
p_benefit_bydeciles <- ggplot(data=bind_rows(empty_tick,obs_v_pred_for_plot), aes(x=mean_predicted_benefit*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100, color= "#0072B2"),width=0.1,size=1) +
  geom_point(aes(y = mean_benefit*100, color="#0072B2"), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Risk-score predicted SGLT2-inhibitor benefit (%)") + ylab("Adjusted observed benefit* (%)")+
  scale_colour_manual(values = "#0072B2") +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Predicted versus observed SGLT2-inhibitor benefit", subtitle = "Binned by predicted benefit decile") +
  coord_cartesian(xlim = c(0,4), ylim = c(0,4))

p_benefit_bydeciles

####
pred <- cohort %>%
  group_by(risk_group) %>%
  summarise(median_predicted_benefit=median(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            mean_predicted_benefit=mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T))

observed <- benefits %>%
  group_by(risk_group, .imp) %>%
  summarise(mean_benefit=mean(survdiff),
            median_benefit=median(survdiff),
            lq_benefit=quantile(survdiff, prob=c(.25)),
            uq_benefit=quantile(survdiff, prob=c(.75))) %>%
  group_by(risk_group) %>%
  summarise(mean_benefit=mean(mean_benefit),
            median_benefit=mean(median_benefit),
            lq_benefit=mean(lq_benefit),
            uq_benefit=mean(uq_benefit))

obs_v_pred <- pred %>% inner_join(observed, by="risk_group") %>% mutate(
  nnt_predicted_mean = 1/mean_predicted_benefit,
  nnt_predicted_median = 1/median_predicted_benefit,
  nnt_obs_adj_mean = 1/mean_benefit,
  nnt_obs_adj_lq = 1/uq_benefit,
  nnt_obs_adj_uq = 1/lq_benefit,
  nnt_obs_adj_median = 1/median_benefit
)

today <- as.character(Sys.Date(), format="%Y%m%d")
save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated_with_adjsurv.Rda"))

############################2 INCIDENCE CURVES################################################################

## plot for eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, low risk ====

#fit unadjusted curve so we can get risk table with number at risk and number of events:

surv_presegfr_microalb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                          data = cohort %>% filter(risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat) &
                                                                            .imp == n.imp),
                                          conf.type = "log", conf.int = 0.95)

list_plots_unadjusted <- list()

list_plots_unadjusted[[1]] <- ggsurvplot(
  fit = surv_presegfr_microalb_lowrisk,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat) &
                                    .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  ylab = "40% decline in eGFR / ESKD"
)


list_plots_unadjusted[[1]]$plot <- list_plots_unadjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_unadjusted[[1]]$table <- list_plots_unadjusted[[1]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[1]]$cumevents <- list_plots_unadjusted[[1]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


surv_presegfr_microalb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                          data = cohort %>% filter(risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat) &
                                                                            .imp == n.imp),
                                          weights = overlap2,
                                          conf.type = "log", conf.int = 0.95)


list_plots_adjusted <- list()

list_plots_adjusted[[1]] <- ggsurvplot(
  fit = surv_presegfr_microalb_lowrisk,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat) &
                                    .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  ylab = "40% decline in eGFR / ESKD",
  title = "(A) eGFR ≥60mL/min, uACR 3-30mg/mmol",
  subtitle = low_risk_cat
)


list_plots_adjusted[[1]]$plot <- list_plots_adjusted[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_adjusted[[1]]$table <- list_plots_unadjusted[[1]]$table 
list_plots_adjusted[[1]]$cumevents <- list_plots_unadjusted[[1]]$cumevents 


## plot for eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, high risk =====
surv_presegfr_microalb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                           data = cohort %>% filter(risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat) &
                                                                             .imp == n.imp),
                                           conf.type = "log", conf.int = 0.95)

list_plots_unadjusted[[3]] <- ggsurvplot(
  fit = surv_presegfr_microalb_highrisk,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat) &
                                    .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  ylab = "40% decline in eGFR / ESKD"
)


list_plots_unadjusted[[3]]$plot <- list_plots_unadjusted[[3]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_unadjusted[[3]]$table <- list_plots_unadjusted[[3]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[3]]$cumevents <- list_plots_unadjusted[[3]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


surv_presegfr_microalb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                           data = cohort %>% filter(risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat) &
                                                                             .imp == n.imp),
                                           weights = overlap2,
                                           conf.type = "log", conf.int = 0.95)


list_plots_adjusted[[3]] <- ggsurvplot(
  fit = surv_presegfr_microalb_highrisk,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat) &
                                    .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  title = "(B) eGFR ≥60mL/min, uACR 3-30mg/mmol",
  subtitle = high_risk_cat
)


list_plots_adjusted[[3]]$plot <- list_plots_adjusted[[3]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_adjusted[[3]]$table <- list_plots_unadjusted[[3]]$table 
list_plots_adjusted[[3]]$cumevents <- list_plots_unadjusted[[3]]$cumevents 


## plot for eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol ====

surv_presegfr_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                  data = cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol" & .imp == n.imp),
                                  conf.type = "log", conf.int = 0.95)

list_plots_unadjusted[[5]] <- ggsurvplot(
  fit = surv_presegfr_macroalb,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol" & .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  ylab = ""
) 

surv_presegfr_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                  data = cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol" & .imp == n.imp),
                                  weights = overlap2,
                                  conf.type = "log", conf.int = 0.95)

list_plots_adjusted[[5]] <- ggsurvplot(
  fit = surv_presegfr_macroalb,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol" & .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  title ="(C) eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
  subtitle = ""
) 

list_plots_adjusted[[5]]$plot <- list_plots_adjusted[[5]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_adjusted[[5]]$table <- list_plots_unadjusted[[5]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_adjusted[[5]]$cumevents <- list_plots_unadjusted[[5]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


## plot for eGFR 20-60mL/min/1.73m2 and uACR <3mg/mmol ====

surv_redegfr_noalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                              data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR <3mg/mmol" & .imp == n.imp),
                              conf.type = "log", conf.int = 0.95)

list_plots_unadjusted[[2]] <- ggsurvplot(
  fit = surv_redegfr_noalb,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR <3mg/mmol" & .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  ylab = "40% decline in eGFR / ESKD"
)

list_plots_unadjusted[[2]]$plot <- list_plots_unadjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_unadjusted[[2]]$table <- list_plots_unadjusted[[2]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[2]]$cumevents <- list_plots_unadjusted[[2]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


surv_redegfr_noalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                              data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR <3mg/mmol" & .imp == n.imp),
                              weights = overlap2,
                              conf.type = "log", conf.int = 0.95)

list_plots_adjusted[[2]] <- ggsurvplot(
  fit = surv_redegfr_noalb,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR <3mg/mmol" & .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  title = paste0(c("(D) eGFR 20-60mL/min/1.73m2, uACR <3mg/mmol"), collapse=""),
  subtitle = " "
)

list_plots_adjusted[[2]]$plot <- list_plots_adjusted[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_adjusted[[2]]$table <- list_plots_unadjusted[[2]]$table
list_plots_adjusted[[2]]$cumevents <- list_plots_unadjusted[[2]]$cumevents


## plot for eGFR 20-60mL/min/1.73m2 and uACR 3-30mg/mmol ====

surv_redegfr_microalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                 data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol" & .imp == n.imp),
                                 conf.type = "log", conf.int = 0.95)


list_plots_unadjusted[[4]] <- ggsurvplot(
  fit = surv_redegfr_microalb,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol" & .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  title = paste0(c("eGFR 20-60mL/min, uACR 3-30mg/mmol"), collapse=""),
  subtitle = " " 
)

list_plots_unadjusted[[4]]$plot <- list_plots_unadjusted[[4]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_unadjusted[[4]]$table <- list_plots_unadjusted[[4]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[4]]$cumevents <- list_plots_unadjusted[[4]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


surv_redegfr_microalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                 data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol" & .imp == n.imp),
                                 weights = overlap2,
                                 conf.type = "log", conf.int = 0.95)


list_plots_adjusted[[4]] <- ggsurvplot(
  fit = surv_redegfr_microalb,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol" & .imp == n.imp),
  
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  title = paste0(c("(E) eGFR 20-60mL/min/1.73m2, uACR 3-30mg/mmol"), collapse=""),
  subtitle = " "
)


list_plots_adjusted[[4]]$plot <- list_plots_adjusted[[4]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_adjusted[[4]]$table <- list_plots_unadjusted[[4]]$table
list_plots_adjusted[[4]]$cumevents <- list_plots_unadjusted[[4]]$cumevents


## plot for eGFR 20-60mL/min/1.73m2, uACR ≥30mg/mmol ====


surv_redegfr_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                 data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol" & .imp == n.imp),
                                 conf.type = "log", conf.int = 0.95)

list_plots_unadjusted[[6]] <- ggsurvplot(
  fit = surv_redegfr_macroalb,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol" & .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  tables.theme = theme_cleantable()
)

list_plots_unadjusted[[6]]$plot <- list_plots_unadjusted[[6]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
list_plots_unadjusted[[6]]$table <- list_plots_unadjusted[[6]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_unadjusted[[6]]$cumevents <- list_plots_unadjusted[[6]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


surv_redegfr_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                 data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol" & .imp == n.imp),
                                 weights = overlap2,
                                 conf.type = "log", conf.int = 0.95)

list_plots_adjusted[[6]] <- ggsurvplot(
  fit = surv_redegfr_macroalb,
  fun = "cumhaz",
  data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol" & .imp == n.imp),
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("DPP4-inhibitors/sulfonylureas", "SGLT2-inhibitors"),
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
  title = "(F) eGFR 20-60mL/min/1.73m2, uACR ≥30mg/mmol",
  subtitle = " "
)

list_plots_adjusted[[6]]$plot <- list_plots_adjusted[[6]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y_lowegfr)
list_plots_adjusted[[6]]$table <- list_plots_unadjusted[[6]]$table + theme(plot.title = element_text(size = font_risktable_title))
list_plots_adjusted[[6]]$cumevents <- list_plots_unadjusted[[6]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


## plot 6-panel plot ====
arrange_ggsurvplots(list_plots_adjusted, print = T, nrow = 2, ncol = 3)



############################3 TABLES WITH NNTs################################################################

## calculate ARR and NNT for eGFR 60mL/min
# add number per group and events
pred40 <- cohort %>%  group_by(risk_group
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

n.groups <- nlevels(as.factor(cohort$risk_group))

# run cox model in each imputed dataset, extract estimates, and pool results
survival_est <- list()

for (i in 1:n.imp) {
  
  ddist <- datadist(cohort[cohort$.imp == i,])
  options(datadist=ddist)
  
  model <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2*risk_group, 
               data=cohort[cohort$.imp == i,], x=TRUE, y=TRUE, surv=TRUE)
  
  
  survival_est[[i]] <- survest(model, 
                               newdata=expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
                                                   risk_group=levels(as.factor(cohort$risk_group))),
                               times=5)
}

# now we will create empty vectors for the pooled benefits per category

surv <- rep(NA, n.groups*nlevels(as.factor(cohort$studydrug2)))
se <- rep(NA, n.groups*nlevels(as.factor(cohort$studydrug2)))

# for every group we will extract the n.imp imputed values and pool these, then store these in the vectors created above

for (k in 1:(n.groups*nlevels(as.factor(cohort$studydrug2)))) {
  
  #create empty vectors to store the imputed values for every individual risk group
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
                               risk_group=levels(as.factor(cohort$risk_group))))

obs <- obs %>%
  pivot_wider(id_cols=c(risk_group), values_from=surv, names_from=studydrug2) %>%
  mutate(surv_diff=SGLT2i-`DPP4i/SU`,
         surv=1-`DPP4i/SU`) %>%
  select(risk_group, surv, surv_diff)

se <- cbind(se, expand.grid(studydrug2=c("DPP4i/SU", "SGLT2i"), 
                            risk_group=levels(as.factor(cohort$risk_group))))

se <- se %>%
  pivot_wider(id_cols=c(risk_group), values_from=se, names_from=studydrug2) %>%
  mutate(se=sqrt((SGLT2i^2)+(`DPP4i/SU`^2))) %>%
  select(risk_group, se)

obs <- obs %>%
  inner_join(se, by=c("risk_group")) %>%
  mutate(lower_ci=surv_diff-(1.96*se),
         upper_ci=surv_diff+(1.96*se))

dpp4isu_events <- cohort %>%
  filter(studydrug2=="DPP4i/SU" & ckd_egfr40_censvar==1) %>%
  group_by(risk_group) %>%
  summarise(`DPP4i/SU`=round(n()/n.imp,0))

sglt2i_events <- cohort %>%
  filter(studydrug2=="SGLT2i" & ckd_egfr40_censvar==1) %>%
  group_by(risk_group) %>%
  summarise(SGLT2i=round(n()/n.imp,0))

events_table <- data.frame(t(dpp4isu_events %>%
                               inner_join(sglt2i_events))) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_group")

obs40 <- obs %>% mutate(risk_group=as.factor(risk_group))
obs_v_pred40 <- pred40 %>% inner_join(obs40, by=c("risk_group")) %>%
  mutate(nnt_observed=1/surv_diff) #%>% select(-se)

nnt_observed_weighted <-   matrix(data = NA, nrow = nlevels(cohort$risk_group), ncol = n.imp)

for (i in 1:n.imp) {
  surv_presegfr_noalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                        data = cohort[cohort$.imp == i &
                                               cohort$risk_group == "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol",],
                        weights = overlap2,
                        conf.type = "log", conf.int = 0.95)
  
  surv_presegfr_microalb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                              data = cohort[cohort$.imp == i &
                                                     cohort$risk_group ==
                                                     paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat),],
                              weights = overlap2,
                              conf.type = "log", conf.int = 0.95)
  
  surv_presegfr_microalb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                               data = cohort[cohort$.imp == i &
                                                      cohort$risk_group ==
                                                      paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat),],
                               weights = overlap2,
                               conf.type = "log", conf.int = 0.95)
  
  surv_presegfr_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                    data = cohort[cohort$.imp == i &
                                                           cohort$risk_group== "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",],
                                    weights = overlap2,
                                    conf.type = "log", conf.int = 0.95)
  
  surv_redegfr_noalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                data = cohort[cohort$.imp == i & 
                                                               cohort$risk_group == "eGFR <60mL/min/1.73m2, uACR <3mg/mmol",],
                                weights = overlap2,
                                conf.type = "log", conf.int = 0.95)
  
  surv_redegfr_microalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                              data = cohort[cohort$.imp == i & 
                                                             cohort$risk_group == "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",],
                              weights = overlap2,
                              conf.type = "log", conf.int = 0.95)
  
  surv_redegfr_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                           data = cohort[cohort$.imp == i &
                                                  cohort$risk_group == "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol",],
                           weights = overlap2,
                           conf.type = "log", conf.int = 0.95)
  
  nnt_observed_weighted[1,i] <- summary(surv_presegfr_noalb, times = 3)$surv[2]-summary(surv_presegfr_noalb, times = 3)$surv[1]
  nnt_observed_weighted[2,i] <- summary(surv_presegfr_microalb_lowrisk, times = 3)$surv[2]-summary(surv_presegfr_microalb_lowrisk, times = 3)$surv[1]
  nnt_observed_weighted[3,i] <- summary(surv_presegfr_microalb_highrisk, times = 3)$surv[2]-summary(surv_presegfr_microalb_highrisk, times = 3)$surv[1]
  nnt_observed_weighted[4,i] <- summary(surv_presegfr_macroalb, times = 3)$surv[2]-summary(surv_presegfr_macroalb, times = 3)$surv[1]
  nnt_observed_weighted[5,i] <- summary(surv_redegfr_noalb, times = 3)$surv[2]-summary(surv_redegfr_noalb, times = 3)$surv[1]
  nnt_observed_weighted[6,i] <- summary(surv_redegfr_microalb, times = 3)$surv[2]-summary(surv_redegfr_microalb, times = 3)$surv[1]
  nnt_observed_weighted[7,i] <- summary(surv_redegfr_macroalb, times = 3)$surv[2]-summary(surv_redegfr_macroalb, times = 3)$surv[1]
}

nnt_observed_weighted <- nnt_observed_weighted %>% rowMeans()

nnt_observed_weighted <- (1/nnt_observed_weighted) %>% round(0)

nnt_table <- obs_v_pred40 %>% select(risk_group, 
                                           risk_predicted=mean_predicted_background_risk, 
                                           risk_observed=surv, 
                                           arr_predicted=mean_predicted_benefit, 
                                           arr_observed=surv_diff, 
                                           arr_lowerci=lower_ci,
                                           arr_upperci=upper_ci,
                                           nnt_predicted, 
                                           nnt_observed) %>%
  mutate(
    risk_predicted = round(risk_predicted*100,3),
    risk_observed = round(risk_observed*100,3),
    arr_predicted =arr_predicted,
    arr_observed = arr_observed,
    nnt_predicted = round(nnt_predicted),
    nnt_observed = round(nnt_observed),
    nnt_observed_lowerci = round(1/arr_upperci),
    nnt_observed_upperci = round(1/arr_lowerci)
  ) %>% cbind(nnt_observed_weighted)

# setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
write.csv2(nnt_table, file = paste0(today, "_nnt_table.csv"))


############################4 SUBGROUP TABLES################################################################


vars <- c("dstartdate_age", "malesex", "ethnicity_5cat", "imd2015_10",             # sociodemographic variables
          "prebmi", "preegfr", "uacr", "albuminuria",                              # vital signs and laboratory measurements
          "preldl", "prehba1c", "presbp", "predbp",   
          "dstartdate_dm_dur_all", "smoking_status", "predrug_hypertension",       # comorbidities
          "predrug_af", "predrug_dka", "osteoporosis", 
          "predrug_acutepancreatitis", "predrug_falls", 
          "predrug_urinary_frequency", "predrug_volume_depletion", 
          "predrug_micturition_control", "predrug_af", "hosp_admission_prev_year",
          "initiation_year",
          "ncurrtx", "MFN", "TZD", "INS", "ACEi_or_ARB",                            # medications
          "cv_high_risk", "qrisk2_above_10_pct"                                     # CV risk
)

#categorical variables
factors <- c("malesex", "ethnicity_qrisk2", "imd2015_10", "smoking_status", "predrug_hypertension", 
             "predrug_af", "predrug_dka", "osteoporosis", "predrug_acutepancreatitis", 
             "predrug_falls", "predrug_urinary_frequency", "predrug_volume_depletion", 
             "predrug_micturition_control", "predrug_af", "hosp_admission_prev_year",
             "initiation_year", "albuminuria", 
             "ncurrtx", "MFN", "TZD", "INS", "ACEi_or_ARB",
             "cv_high_risk", "qrisk2_above_10_pct")

nonnormal <- c("uacr", "dstartdate_dm_dur_all")

table <- CreateTableOne(vars = vars, strata = "risk_group", data = cohort, 
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv
# write.csv2(tabforprint, file = paste0(today, "_ckd40_highrisk_table.csv"))
