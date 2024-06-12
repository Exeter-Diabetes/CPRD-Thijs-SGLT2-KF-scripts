############################0 SETUP################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(rms)
library(broom)
library(adjustedCurves)
library(tableone)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())

# set seed
set.seed(123)
n.imp <- 10
n.quantiles <- 10

#today <- as.character(Sys.Date(), format="%Y%m%d")
today <- "2024-06-06"

#set default colours (colour-blind accessible from the Okabe-Ito palette) for different drug classes
#this is the colour allocation that will be used across the department
cols <- c("SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00")
#in further analyses, the dpp4/su group will be combined, and we will use the dpp4 colour for this (strongest contrast)
cols <- c(cols, "DPP4i/SU" = "#0072B2")

# covariates for multivariable adjustment
covariates <- "dstartdate_age + malesex + imd2015_10 + ethnicity_4cat + prebmi + prehba1c + preegfr + uacr + presbp + ckdpc_40egfr_score + INS + statin + ACEi_or_ARB + smoking_status + dstartdate_dm_dur_all + predrug_hypertension + hosp_admission_prev_year"
# these are slightly pruned compared with the covariates in previous codes - in order for models to converge

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

############################1 CALCULATE ADJUSTED OBSERVED BENEFITS################################################################


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-06-06_t2d_ckdpc_recalibrated.Rda")

cohort$studydrug2 <- as.factor(cohort$studydrug2)

# calculate predicted sglt2 benefit (absolute risk reduction = ARR):
# ARR = S0(t)^HR - S0(t)
trial_hr_kf_sglt2i <- 0.62

cohort <- cohort %>% 
  mutate(ckdpc_40egfr_cal=(100-ckdpc_40egfr_score_cal)/100,
         ckdpc_40egfr_cal_sglt2i=ckdpc_40egfr_cal^trial_hr_kf_sglt2i,
         ckdpc_40egfr_sglt2i_benefit=ckdpc_40egfr_cal_sglt2i - ckdpc_40egfr_cal)

cohort$benefit_decile <- ntile(cohort$ckdpc_40egfr_sglt2i_benefit, n.quantiles)

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

cohort <- cohort %>% mutate(
  ethnicity_4cat = ifelse(!ethnicity_4cat %in% c("White", "South Asian", "Black"), "Other", as.character(ethnicity_4cat)),
  initiation_year = ifelse(initiation_year %in% c("2019", "2020"), "2019-2020", as.character(initiation_year)),
  ncurrtx = ifelse(ncurrtx %in% c("3.", "4+"), "3+", as.character(ncurrtx))
)


save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated_with_riskgroup.Rda"))

# centre predictors and save dataset for later (performing analysis in entire dataset exceeds memory limit)

centre_and_reference <- function(df, covariates) {
  df %>%
    mutate(across(all_of(covariates), ~ if(is.numeric(.)) . - mean(., na.rm = TRUE) else . ),  # Center numeric variables
           across(all_of(covariates), ~ if(is.logical(.)) . == names(sort(table(.), decreasing = TRUE))[1] else . ) ,  # Set most frequent level as reference for logical
           across(all_of(covariates), ~ if(is.factor(.)) relevel(., ref = names(sort(table(.), decreasing = TRUE))[1]) else . ))  # Set most frequent level as reference for factor
}

# define outcomes to be analysed
outcomes <- c("ckd_egfr40", "death", "macroalb", "dka", "amputation", "side_effect")
# create regex pattern of censoring variables to select
outcome_variables <- paste0("(", paste(outcomes, collapse = "|"), ")(?!.*(5y|pp)).*(_censtime_yrs|_censvar)$")

cohort <- cohort %>%
  mutate(across(contains("predrug_"), as.logical),
         hosp_admission_prev_year=as.logical(hosp_admission_prev_year),
         INS=as.logical(INS),
         MFN=as.logical(MFN),
         malesex=as.factor(malesex),
         initiation_year=as.factor(initiation_year)) %>%
  select(patid, .imp, risk_group, studydrug2, benefit_decile,
         matches(outcome_variables, perl = TRUE), 
         all_of(unlist(strsplit(covariates, " \\+ ")))) %>% 
  centre_and_reference(unlist(strsplit(covariates, " \\+ ")))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_recalibrated_data_centred_predictors.Rda"))

rm(list = setdiff(ls(), c("n.imp", "covariates", "today", "outcomes")))


for (k in outcomes) {
  
  for (i in 1:n.imp) {
    print(paste0("Survival estimates for imputation ", i, " (SGLT2i)"))
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
    load("2024-06-06_recalibrated_data_centred_predictors.Rda")
    
    if (k == "macroalb") {
      cohort <- cohort %>% filter(.imp == i & !risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol") & macroalb_censtime_yrs >= 0)
    } else {
      cohort <- cohort %>% filter(.imp == i)
    }
    
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")  
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*benefit_decile + ", covariates)) ## add interaction ##
    
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE)
    
    obs_SGLT2 <- cohort %>%
      mutate(studydrug_original=studydrug2,
             studydrug2="SGLT2i",
             rowno=row_number())
    
    print(paste0("Extracting survival estimates (outcome ", k, ")"))
    
    observed_sglt2 <- survfit(model, newdata=as.data.frame(obs_SGLT2)) %>%
      tidy() %>%
      filter(time==3) %>%
      pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
      select(group, estimate_sglt2=estimate, se_sglt2=std.error) %>%
      mutate(group=as.numeric(group)) %>%
      inner_join(obs_SGLT2, by=c("group"="rowno"))
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
    
    save(observed_sglt2, file=paste0(today, "_adjusted_surv_",k,"_SGLT2i_imp.", i, ".Rda"))
    
    rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "outcomes")))
  }
  
  
  for (i in 1:n.imp) {
    print(paste0("Survival estimates for imputation ", i, " (DPP4i/SU)"))
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
    load("2024-06-06_recalibrated_data_centred_predictors.Rda")
    
    if (k == "macroalb") {
      cohort <- cohort %>% filter(.imp == i & !risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol") & macroalb_censtime_yrs >= 0)
    } else {
      cohort <- cohort %>% filter(.imp == i)
    }
    
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")  
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*benefit_decile + ", covariates))
    
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE)
    
    obs_DPP4SU <- cohort %>%
      mutate(studydrug_original=studydrug2,
             studydrug2="DPP4i/SU",
             rowno=row_number())
    
    print(paste0("Extracting survival estimates (outcome ", k, ")"))
    
    observed_dpp4su <- survfit(model, newdata=as.data.frame(obs_DPP4SU)) %>%
      tidy() %>%
      filter(time==3) %>%
      pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
      select(group, estimate_dpp4su=estimate, se_dpp4su=std.error) %>%
      mutate(group=as.numeric(group)) %>%
      inner_join(obs_DPP4SU, by=c("group"="rowno"))
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
    save(observed_dpp4su, file=paste0(today, "_adjusted_surv_",k,"_DPP4iSU_imp.", i, ".Rda"))
    
    rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "outcomes")))
  }
  
  temp_sglt2 <- temp_dpp4su <- data.frame()
  
  for (i in 1:n.imp) {
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
    load(paste0("2024-06-06_adjusted_surv_",k,"_SGLT2i_imp.", i, ".Rda"))
    load(paste0("2024-06-06_adjusted_surv_",k,"_DPP4iSU_imp.", i, ".Rda"))
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
           se_survdiff=sqrt(se_sglt2^2 + se_dpp4su^2),
           studydrug2=studydrug_original) %>%
    select(-studydrug_original) 
  
  #rename variables
  benefits <- benefits %>% rename_with(
    ~ paste0(.x, "_", k),
    contains("survdiff")
  )
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
  save(benefits, file=paste0(today, "_adjusted_surv_",k,".Rda"))
  rm(benefits)
}

### add adjusted observed survival probabilities to main dataset
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-06-06_t2d_ckdpc_recalibrated_with_riskgroup.Rda")

for (k in outcomes) {
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
  load(paste0("2024-06-05_adjusted_surv_", k, ".Rda"))
  cohort <- cohort %>% left_join(benefits %>%
                                   select(.imp, patid, studydrug2, contains("survdiff")), 
                                 by=c(".imp", "patid", "studydrug2"))
  rm(benefits)
}

cohort$studydrug2 <- as.factor(cohort$studydrug2)

# save dataset with adjusted survival probabilities

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated_with_adjsurv.Rda"))


############################2 CALIBRATION PLOTS OF PREDICTED VS OBSERVED BENEFITS################################################################

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-06-06_t2d_ckdpc_recalibrated_with_adjsurv.Rda")

obs_v_pred_for_plot <- cohort %>%
  group_by(benefit_decile) %>%
  summarise(median_predicted_benefit=median(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            mean_predicted_benefit=mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            mean_benefit=mean(survdiff_ckd_egfr40),
            se_benefit=mean(se_survdiff_ckd_egfr40),
            median_benefit=median(survdiff_ckd_egfr40),
            lq_benefit=quantile(survdiff_ckd_egfr40, prob=c(.25)),
            uq_benefit=quantile(survdiff_ckd_egfr40, prob=c(.75)),
            upper_ci=mean_benefit + 1.96*se_benefit,
            lower_ci=mean_benefit - 1.96*se_benefit)


empty_tick <- data.frame(matrix(NA, nrow = 1, ncol = length(obs_v_pred_for_plot)))
names(empty_tick) <- names(obs_v_pred_for_plot)
empty_tick <- empty_tick %>%
  mutate(benefit_decile=0)

## SGLT2i benefit predicted vs observed - mean
p_benefit_bydeciles_mean <- ggplot(data=bind_rows(empty_tick,obs_v_pred_for_plot), aes(x=mean_predicted_benefit*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100, color= "#0072B2"),width=0.1,size=1) +
  geom_point(aes(y = mean_benefit*100, color="#0072B2"), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Risk-score predicted SGLT2-inhibitor benefit (%)") + ylab("Observed benefit* (%)")+
  scale_colour_manual(values = "#0072B2") +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Mean predicted versus observed SGLT2-inhibitor benefit", subtitle = "By predicted benefit decile") +
  coord_cartesian(xlim = c(0,3), ylim = c(0,3))

p_benefit_bydeciles_mean

## SGLT2i benefit predicted vs observed - median
p_benefit_bydeciles_median <- ggplot(data=bind_rows(empty_tick,obs_v_pred_for_plot), aes(x=median_predicted_benefit*100)) +
  geom_errorbar(aes(ymax=uq_benefit*100,ymin=lq_benefit*100, color= "#0072B2"),width=0.1,size=1) +
  geom_point(aes(y = median_benefit*100, color="#0072B2"), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Risk-score predicted SGLT2-inhibitor benefit (%)") + ylab("Observed benefit* (%)")+
  scale_colour_manual(values = "#0072B2") +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Median predicted versus observed SGLT2-inhibitor benefit", subtitle = "By predicted benefit decile") +
  coord_cartesian(xlim = c(0,3), ylim = c(0,3))

p_benefit_bydeciles_median


####
pred <- cohort %>%
  group_by(risk_group, .imp) %>%
  summarise(median_predicted_benefit=median(ckdpc_40egfr_sglt2i_benefit, na.rm=T),
            lq_predicted_benefit=quantile(ckdpc_40egfr_sglt2i_benefit, prob=c(.25)),
            uq_predicted_benefit=quantile(ckdpc_40egfr_sglt2i_benefit, prob=c(.75)),
            mean_predicted_benefit=mean(ckdpc_40egfr_sglt2i_benefit, na.rm=T)) %>%
  group_by(risk_group) %>%
  summarise(median_predicted_benefit=mean(median_predicted_benefit),
            lq_predicted_benefit=mean(lq_predicted_benefit),
            uq_predicted_benefit=mean(uq_predicted_benefit),
            mean_predicted_benefit=mean(mean_predicted_benefit))

observed <- cohort %>%
  group_by(risk_group, .imp) %>%
  summarise(mean_benefit=mean(survdiff_ckd_egfr40),
            se_benefit=mean(se_survdiff_ckd_egfr40),
            median_benefit=median(survdiff_ckd_egfr40),
            lq_benefit=quantile(survdiff_ckd_egfr40, prob=c(.25)),
            uq_benefit=quantile(survdiff_ckd_egfr40, prob=c(.75))) %>%
  group_by(risk_group) %>%
  summarise(mean_benefit=mean(mean_benefit),
            se_benefit=mean(se_benefit),
            median_benefit=mean(median_benefit),
            lq_benefit=mean(lq_benefit),
            uq_benefit=mean(uq_benefit),
            lc_benefit=mean_benefit-1.96*se_benefit,
            uc_benefit=mean_benefit+1.96*se_benefit)

obs_v_pred <- pred %>% inner_join(observed, by="risk_group") %>% mutate(
  nnt_predicted_mean = 1/mean_predicted_benefit,
  nnt_predicted_median = 1/median_predicted_benefit,
  nnt_predicted_lq = 1/uq_predicted_benefit,
  nnt_predicted_uq = 1/lq_predicted_benefit,
  nnt_obs_adj_mean = 1/mean_benefit,
  nnt_obs_adj_lc = 1/uc_benefit,
  nnt_obs_adj_uc = 1/lc_benefit,
  nnt_obs_adj_uc = ifelse(nnt_obs_adj_uc < 0, Inf, nnt_obs_adj_uc),
  nnt_obs_adj_median = 1/median_benefit,
  nnt_obs_adj_lq = 1/uq_benefit,
  nnt_obs_adj_uq = 1/lq_benefit
)

# observed_macroalb <- cohort %>%
#   filter(!risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol")) %>%
#   group_by(risk_group, .imp) %>%
#   summarise(mean_benefit=mean(survdiff_macroalb, na.rm = T),
#             se_benefit=mean(se_survdiff_macroalb, na.rm = T),
#             median_benefit=median(survdiff_macroalb, na.rm = T),
#             lq_benefit=quantile(survdiff_macroalb, prob=c(.25), na.rm = T),
#             uq_benefit=quantile(survdiff_macroalb, prob=c(.75)), na.rm = T) %>%
#   group_by(risk_group) %>%
#   summarise(mean_benefit=mean(mean_benefit, na.rm = T),
#             se_benefit=mean(se_benefit, na.rm = T),
#             median_benefit=mean(median_benefit, na.rm = T),
#             lq_benefit=mean(lq_benefit, na.rm = T),
#             uq_benefit=mean(uq_benefit, na.rm = T),
#             lc_benefit=mean_benefit-1.96*se_benefit,
#             uc_benefit=mean_benefit+1.96*se_benefit) %>% 
#   mutate(nnt_obs_adj_mean = 1/mean_benefit,
#          nnt_obs_adj_lc = 1/uc_benefit,
#          nnt_obs_adj_uc = 1/lc_benefit,
#          nnt_obs_adj_uc = ifelse(nnt_obs_adj_uc < 0, Inf, nnt_obs_adj_uc),
#          nnt_obs_adj_median = 1/median_benefit,
#          nnt_obs_adj_lq = 1/uq_benefit,
#          nnt_obs_adj_uq = 1/lq_benefit
#          )


result_list <- list()

for (k in outcomes) {
  if (k == "macroalb") {
    temp <- cohort %>%
      filter(!risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol")) %>%
      group_by(risk_group, .imp) %>%
      summarise(across(starts_with(paste0("survdiff_", k)), 
                       list(mean_benefit = ~mean(.),
                            se_benefit = ~mean(get(paste0("se_survdiff_", k))),
                            median_benefit = ~median(.),
                            lq_benefit = ~quantile(., prob = .25),
                            uq_benefit = ~quantile(., prob = .75)),
                       .names = "{.fn}_{.col}")) %>%
      group_by(risk_group) %>%
      summarise(across(contains("mean_benefit"), ~mean(.), .names = "{sub('mean_benefit', '', .col)}"),
                across(contains("se_benefit"), ~mean(.), .names = "{sub('se_benefit','', .col)}"),
                across(contains("median_benefit"), ~mean(.), .names = "{sub('median_benefit','', .col)}"),
                across(contains("lq_benefit"), ~mean(.), .names = "{sub('lq_benefit', '',.col)}"),
                across(contains("uq_benefit"), ~mean(.), .names = "{sub('uq_benefit','', .col)}")) %>%
      mutate(across(contains("mean_benefit"),
                    list(lc_benefit = ~. - 1.96 * get(sub("mean_benefit", "se_benefit", cur_column())),
                         uc_benefit = ~. + 1.96 * get(sub("mean_benefit", "se_benefit", cur_column()))),
                    .names = "{.col}_{.fn}"))
  } else {
    temp <- cohort %>%
      group_by(risk_group, .imp) %>%
      summarise(across(starts_with(paste0("survdiff_", k)), 
                       list(mean_benefit = ~mean(.),
                            se_benefit = ~mean(get(paste0("se_survdiff_", k))),
                            median_benefit = ~median(.),
                            lq_benefit = ~quantile(., prob = .25),
                            uq_benefit = ~quantile(., prob = .75)),
                       .names = "{.fn}_{.col}")) %>%
      group_by(risk_group) %>%
      summarise(across(contains("_mean_benefit"), ~mean(.), .names = "{sub('_mean_benefit', '', .col)}"),
                across(contains("_se_benefit"), ~mean(.), .names = "{sub('_se_benefit', 'se_benefit', .col)}"),
                across(contains("_median_benefit"), ~mean(.), .names = "{sub('_median_benefit', 'median_benefit', .col)}"),
                across(contains("_lq_benefit"), ~mean(.), .names = "{sub('_lq_benefit', 'lq_benefit', .col)}"),
                across(contains("_uq_benefit"), ~mean(.), .names = "{sub('_uq_benefit', 'uq_benefit', .col)}")) %>%
      mutate(across(contains("_benefit"),
                    list(lc_benefit = ~. - 1.96 * get(sub("_benefit", "se_benefit", cur_column())),
                         uc_benefit = ~. + 1.96 * get(sub("_benefit", "se_benefit", cur_column()))),
                    .names = "{.col}_{.fn}"))
  }
  result_list[[k]] <- temp
}

observed <- reduce(result_list, full_join, by = "risk_group")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
write.csv2(obs_v_pred, file = paste0(today, "_nnt_table.csv"))

############################3 INCIDENCE CURVES################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-05-31_t2d_ckdpc_recalibrated_with_adjsurv.Rda")

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
limit_y <- c(0,0.38)
zoom_y <- c(0,0.25)
zoom_x <- c(0,3)
break_y <- 0.02
legend_in <- c(0.35, 0.55)
legend_out <- c(100,-100)


#fit unadjusted curve so we can get risk table with number at risk and number of events:
k <- "ckd_egfr40"
censvar_var=paste0(k, "_censvar")
censtime_var=paste0(k, "_censtime_yrs") 
f <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2"))
f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", covariates))

fit_presegfr_noalb <- coxph(f_adjusted,
                            data = cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol" &
                                                       .imp == n.imp) %>% as.data.frame(),
                            x=T)

plot0 <- adjustedsurv(data=cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol" &
                                               .imp == n.imp) %>% as.data.frame(),
                      variable="studydrug2",
                      ev_time="ckd_egfr40_censtime_yrs",
                      event="ckd_egfr40_censvar",
                      method="direct",
                      outcome_model=fit_presegfr_noalb,
                      conf_int=TRUE)

panel_0 <- plot(plot0, 
                conf_int = T, 
                cif = T, 
                title = "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol",
                subtitle = "",
                xlab = "Years",
                ylab = "Cumulative incidence of 40% decline in eGFR / ESKD",
                custom_colors = cols_fig,
                ylim = limit_y,
                legend.title = "",
                legend.position = legend_in,
                additional_layers = c(coord_cartesian(xlim = zoom_x, ylim = zoom_y)),
                risk_table = T,
                risk_table_stratify = T,
                risk_table_type = "n_at_risk",
                risk_table_xlab = NULL,
                risk_table_ylab = NULL,
                risk_table_size = 3.5,
                risk_table_title_size = 12,
                risk_table_warn = F
) 


fit_presegfr_microalb_lowrisk <- coxph(f_adjusted,
                                       data = cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score 0-90th percentile" &
                                                                  .imp == n.imp) %>% as.data.frame(),
                                       x=T)

plot1 <- adjustedsurv(data=cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score 0-90th percentile" &
                                               .imp == n.imp) %>% as.data.frame(),
                      variable="studydrug2",
                      ev_time="ckd_egfr40_censtime_yrs",
                      event="ckd_egfr40_censvar",
                      method="direct",
                      outcome_model=fit_presegfr_microalb_lowrisk,
                      conf_int=TRUE)

panel_1 <- plot(plot1, 
                conf_int = T, 
                cif = T, 
                title = "(A) eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                subtitle = "CKD-PC risk score 0-90th percentile",
                xlab = "",
                ylab = "Cumulative incidence of 40% decline in eGFR / ESKD",
                custom_colors = cols_fig,
                ylim = limit_y,
                legend.title = "",
                legend.position = legend_in,
                additional_layers = c(coord_cartesian(xlim = zoom_x, ylim = zoom_y)),
                risk_table = T,
                risk_table_stratify = T,
                risk_table_type = "n_at_risk",
                risk_table_xlab = NULL,
                risk_table_ylab = NULL,
                risk_table_size = 3.5,
                risk_table_title_size = 12,
                risk_table_warn = F
) 

fit_presegfr_microalb_highrisk <- coxph(f_adjusted,
                                        data = cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score ≥90th percentile" &
                                                                   .imp == n.imp),
                                        x=T)

plot2 <- adjustedsurv(data=cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score ≥90th percentile" &
                                               .imp == n.imp),
                      variable="studydrug2",
                      ev_time="ckd_egfr40_censtime_yrs",
                      event="ckd_egfr40_censvar",
                      method="direct",
                      outcome_model=fit_presegfr_microalb_highrisk,
                      conf_int=TRUE)

panel_2 <- plot(plot2, 
                conf_int = T, 
                cif = T, 
                xlab = "",
                title = "(B) eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",
                subtitle = "CKD-PC risk score ≥90th percentile",
                ylab = "",
                custom_colors = cols_fig,
                ylim = limit_y,
                legend.title = "",
                legend.position = legend_out,
                additional_layers = c(coord_cartesian(xlim = zoom_x, ylim = zoom_y)),
                risk_table = T,
                risk_table_stratify = T,
                risk_table_type = "n_at_risk",
                risk_table_xlab = NULL,
                risk_table_ylab = NULL,
                risk_table_size = 3.5,
                risk_table_title_size = 12,
                risk_table_warn = F
) 

fit_presegfr_macroalb <- coxph(f_adjusted,
                               data = cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol" &
                                                          .imp == n.imp),
                               x=T)

plot3 <- adjustedsurv(data=cohort %>% filter(risk_group == "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol" &
                                               .imp == n.imp),
                      variable="studydrug2",
                      ev_time="ckd_egfr40_censtime_yrs",
                      event="ckd_egfr40_censvar",
                      method="direct",
                      outcome_model=fit_presegfr_macroalb,
                      conf_int=TRUE)

panel_3 <- plot(plot3, 
                conf_int = T, 
                cif = T, 
                xlab = "",
                title = "(C) eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                subtitle = "",
                ylab = "",
                custom_colors = cols_fig,
                ylim = limit_y,
                legend.title = "",
                legend.position = legend_out,
                additional_layers = c(coord_cartesian(xlim = zoom_x, ylim = zoom_y)),
                risk_table = T,
                risk_table_stratify = T,
                risk_table_type = "n_at_risk",
                risk_table_xlab = NULL,
                risk_table_ylab = NULL,
                risk_table_size = 3.5,
                risk_table_title_size = 12,
                risk_table_warn = F
) 

fit_redegfr_noalb <- coxph(f_adjusted,
                           data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR <3mg/mmol" &
                                                      .imp == n.imp),
                           x=T)

plot4 <- adjustedsurv(data=cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR <3mg/mmol" &
                                               .imp == n.imp),
                      variable="studydrug2",
                      ev_time="ckd_egfr40_censtime_yrs",
                      event="ckd_egfr40_censvar",
                      method="direct",
                      outcome_model=fit_redegfr_noalb,
                      conf_int=TRUE)

panel_4 <- plot(plot4, 
                conf_int = T, 
                cif = T, 
                xlab = "Years",
                title = "(D) eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                subtitle = "",
                ylab = "Cumulative incidence of 40% decline in eGFR / ESKD",
                custom_colors = cols_fig,
                ylim = limit_y,
                legend.title = "",
                legend.position = legend_out,
                additional_layers = c(coord_cartesian(xlim = zoom_x, ylim = zoom_y)),
                risk_table = T,
                risk_table_stratify = T,
                risk_table_type = "n_at_risk",
                risk_table_xlab = NULL,
                risk_table_ylab = NULL,
                risk_table_size = 3.5,
                risk_table_title_size = 12,
                risk_table_warn = F
) 

fit_redegfr_microalb <- coxph(f_adjusted,
                              data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol" &
                                                         .imp == n.imp),
                              x=T)

plot5 <- adjustedsurv(data=cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol" &
                                               .imp == n.imp),
                      variable="studydrug2",
                      ev_time="ckd_egfr40_censtime_yrs",
                      event="ckd_egfr40_censvar",
                      method="direct",
                      outcome_model=fit_presegfr_macroalb,
                      conf_int=TRUE)

panel_5 <- plot(plot5, 
                conf_int = T, 
                cif = T, 
                xlab = "Years",
                title = "(E) eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol",
                subtitle = "",
                ylab = "",
                custom_colors = cols_fig,
                ylim = limit_y,
                legend.title = "",
                legend.position = legend_out,
                additional_layers = c(coord_cartesian(xlim = zoom_x, ylim = zoom_y)),
                risk_table = T,
                risk_table_stratify = T,
                risk_table_type = "n_at_risk",
                risk_table_xlab = NULL,
                risk_table_ylab = NULL,
                risk_table_size = 3.5,
                risk_table_title_size = 12,
                risk_table_warn = F
) 

fit_redegfr_macroalb <- coxph(f_adjusted,
                              data = cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol" &
                                                         .imp == n.imp),
                              x=T)

plot6 <- adjustedsurv(data=cohort %>% filter(risk_group == "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol" &
                                               .imp == n.imp),
                      variable="studydrug2",
                      ev_time="ckd_egfr40_censtime_yrs",
                      event="ckd_egfr40_censvar",
                      method="direct",
                      outcome_model=fit_redegfr_macroalb,
                      conf_int=TRUE)

panel_6 <- plot(plot6, 
                conf_int = T, 
                cif = T, 
                xlab = "Years",
                title = "(F) eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol",
                subtitle = "",
                ylab = "",
                custom_colors = cols_fig,
                ylim = limit_y,
                legend.title = "",
                legend.position = legend_out,
                additional_layers = c(coord_cartesian(xlim = zoom_x, ylim = zoom_y)),
                risk_table = T,
                risk_table_stratify = T,
                risk_table_type = "n_at_risk",
                risk_table_xlab = NULL,
                risk_table_ylab = NULL,
                risk_table_size = 3.5,
                risk_table_title_size = 12,
                risk_table_warn = F
) 


plot <- plot_grid(panel_1, panel_2, panel_3, panel_4, panel_5, panel_6, ncol = 3, nrow = 2)

plot

############################4 SUBGROUP TABLES################################################################


vars <- c("dstartdate_age", "malesex", "ethnicity_4cat", "imd2015_10",             # sociodemographic variables
          "prebmi", "preegfr", "uacr", "pretotalcholesterol",
          "prehba1c", "presbp", "predbp",   
          "dstartdate_dm_dur_all", "smoking_status", "predrug_hypertension",       # comorbidities
          "predrug_af", "predrug_dka", "hosp_admission_prev_year",
          "initiation_year",
          "ncurrtx", "MFN", "INS", "ACEi_or_ARB",                            # medications
          "cv_high_risk", "qrisk2_above_10_pct"                                     # CV risk
)

#categorical variables
factors <- c("malesex", "ethnicity_4cat", "imd2015_10", "smoking_status", "predrug_hypertension", 
             "predrug_af", "predrug_dka", "hosp_admission_prev_year",
             "initiation_year", 
             "ncurrtx", "MFN", "INS", "ACEi_or_ARB",
             "cv_high_risk", "qrisk2_above_10_pct")

nonnormal <- c("uacr", "dstartdate_dm_dur_all")

table <- CreateTableOne(vars = vars, strata = "risk_group", data = cohort, 
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv
# write.csv2(tabforprint, file = paste0(today, "_ckd40_highrisk_table.csv"))
