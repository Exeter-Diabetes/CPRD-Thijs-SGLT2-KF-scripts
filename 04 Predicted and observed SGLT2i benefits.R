
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
today <- "2024-07-13"

# covariates for multivariable adjustment
covariates <- c("dstartdate_age", "malesex", "imd2015_10", "ethnicity_4cat", "initiation_year", "prebmi", "prehba1c",
                "pretotalcholesterol", "preegfr", "uacr", "presbp", "ckdpc_40egfr_score", "ncurrtx", "statin", "INS", 
                "ACEi_or_ARB", "smoking_status", "dstartdate_dm_dur_all", "predurg_hypertension", "predrug_af", "hosp_admission_prev_year")

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
load("2024-07-13_t2d_ckdpc_recalibrated.Rda")

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

cohort <- cohort %>% mutate(
  ethnicity_4cat = ifelse(!ethnicity_4cat %in% c("White", "South Asian", "Black"), "Other", as.character(ethnicity_4cat)),
  initiation_year = ifelse(initiation_year %in% c("2019", "2020"), "2019-2020", as.character(initiation_year)),
  ncurrtx = ifelse(ncurrtx %in% c("3.", "4+"), "3+", as.character(ncurrtx))
)

cohort <- cohort %>% filter(!.imp > n.imp)

save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated_with_riskgroup.Rda"))

# centre predictors and save dataset for later (performing analysis in entire dataset exceeds memory limit)

centre_and_reference <- function(df, covariates) {
  df %>%
    mutate(across(all_of(covariates), ~ if(is.numeric(.)) . - mean(., na.rm = TRUE) else . ),  # Center numeric variables
           across(all_of(covariates), ~ if(is.logical(.)) . == names(sort(table(.), decreasing = TRUE))[1] else . ) ,  # Set most frequent level as reference for logical
           across(all_of(covariates), ~ if(is.factor(.)) relevel(., ref = names(sort(table(.), decreasing = TRUE))[1]) else . ))  # Set most frequent level as reference for factor
}

# define outcomes to be analysed
outcomes <- c("ckd_egfr40", "macroalb", "dka", "side_effect"#, "death", "amputation")
)

# create regex pattern of censoring variables to select
outcome_variables <- paste0("(", paste(outcomes, collapse = "|"), ")(?!.*(5y|pp)).*(_censtime_yrs|_censvar)$")

cohort <- cohort %>%
  mutate(across(contains("predrug_"), as.logical),
         hosp_admission_prev_year=as.logical(hosp_admission_prev_year),
         INS=as.logical(INS),
         MFN=as.logical(MFN),
         malesex=as.factor(malesex),
         initiation_year=as.factor(initiation_year),
         benefit_decile=as.factor(benefit_decile)) %>%
  select(patid, .imp, risk_group, studydrug2, benefit_decile,
         matches(outcome_variables, perl = TRUE), 
         all_of(covariates)) %>% 
  centre_and_reference(covariates)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_recalibrated_data_centred_predictors.Rda"))

rm(list = setdiff(ls(), c("n.imp", "covariates", "today", "outcomes")))


for (k in outcomes) {
  
  for (i in 1:n.imp) {
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
    load("2024-07-13_recalibrated_data_centred_predictors.Rda")
    
    if (k == "macroalb") {
      cohort <- cohort %>% filter(.imp == i & !risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol") & macroalb_censtime_yrs >= 0)
    } else {
      cohort <- cohort %>% filter(.imp == i)
    }
    
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")  
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*benefit_decile + ", paste(covariates, collapse=" + "))) ## add interaction ##
    
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE)
    
    obs_SGLT2 <- cohort %>%
      mutate(studydrug_original=studydrug2,
             studydrug2="SGLT2i",
             rowno=row_number())
    print(paste0("Survival estimates for SGLT2i in imputation ", i, "  (outcome ", k, ")"))
    
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
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
    load("2024-07-13_recalibrated_data_centred_predictors.Rda")
    
    if (k == "macroalb") {
      cohort <- cohort %>% filter(.imp == i & !risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol") & macroalb_censtime_yrs >= 0)
    } else {
      cohort <- cohort %>% filter(.imp == i)
    }
    
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")  
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*benefit_decile + ", paste(covariates, collapse=" + ")))
    
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE)
    
    obs_DPP4SU <- cohort %>%
      mutate(studydrug_original=studydrug2,
             studydrug2="DPP4i/SU",
             rowno=row_number())
    
    print(paste0("Survival estimates for DPP4i/SU in imputation ", i, "  (outcome ", k, ")"))
    
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
    load(paste0("2024-07-13_adjusted_surv_",k,"_SGLT2i_imp.", i, ".Rda"))
    load(paste0("2024-07-13_adjusted_surv_",k,"_DPP4iSU_imp.", i, ".Rda"))
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

rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "outcomes")))

### add adjusted observed survival probabilities to main dataset
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-07-13_t2d_ckdpc_recalibrated_with_riskgroup.Rda")

for (k in outcomes) {
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
  load(paste0("2024-07-13_adjusted_surv_", k, ".Rda"))
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
load("2024-07-13_t2d_ckdpc_recalibrated_with_adjsurv.Rda")

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
  xlab("Predicted SGLT2-inhibitor benefit (%)") + ylab("Observed benefit (%)")+
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
  coord_cartesian(xlim = c(0,4), ylim = c(0,4))

p_benefit_bydeciles_mean

## SGLT2i benefit predicted vs observed - median
p_benefit_bydeciles_median <- ggplot(data=bind_rows(empty_tick,obs_v_pred_for_plot), aes(x=median_predicted_benefit*100)) +
  geom_errorbar(aes(ymax=uq_benefit*100,ymin=lq_benefit*100, color= "#0072B2"),width=0.1,size=1) +
  geom_point(aes(y = median_benefit*100, color="#0072B2"), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Predicted SGLT2-inhibitor benefit (%)") + ylab("Observed benefit (%)")+
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
  coord_cartesian(xlim = c(0,4), ylim = c(0,4))

p_benefit_bydeciles_median


#### NNT TABLES

## make risk groups based on risk cutoff
risk_threshold <- .80

risk_threshold_1 <- cohort[cohort$risk_group == "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol",]$ckdpc_40egfr_score_cal %>% quantile(risk_threshold) %>% as.numeric()

high_risk_cat <- paste(c("CKD-PC risk score ≥", risk_threshold*100, "th percentile"), collapse = "")
low_risk_cat <- paste(c("CKD-PC risk score 0-", risk_threshold*100, "th percentile"), collapse = "")




pred <- cohort %>% mutate(
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
  ) %>%
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

pred_nnt <- pred %>%
  mutate(across(contains("_benefit"), ~ 1 / ., .names = "{str_replace(.col, '_benefit', '_nnt')}")) %>%
  select(risk_group, contains("_nnt"))

observed <- cohort %>% mutate(
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
) %>%
  select(risk_group, .imp, contains("survdiff")) %>%
  group_by(risk_group, .imp) %>%
  summarise(
    across(starts_with("survdiff_"), list(
      mean_benefit = ~ mean(.x, na.rm = TRUE),
      median_benefit = ~ median(.x, na.rm = TRUE),
      lq_benefit = ~ quantile(.x, probs = 0.25, na.rm = TRUE),
      uq_benefit = ~ quantile(.x, probs = 0.75, na.rm = TRUE)
    ), .names = "{.col}_{.fn}"),
    across(starts_with("se_survdiff_"), list(
      se_benefit = ~ mean(.x, na.rm = TRUE)
    ), .names = "{.col}_{.fn}")
  ) %>%
  ungroup()

# Rename the columns to match the desired format
observed <- observed %>%
  rename_with(~ gsub("survdiff_", "", .), starts_with("survdiff_")) %>%
  rename_with(~ gsub("se_survdiff_", "", .), starts_with("se_survdiff_"))


observed <- observed %>%
  group_by(risk_group) %>%
  summarise(
    across(contains("mean_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(contains("se_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(contains("median_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(contains("lq_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(contains("uq_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}")
  ) %>%
  mutate(
    across(ends_with("_mean_benefit"), ~ . - 1.96 * get(str_replace(cur_column(), "mean_", "se_")), .names = "{.col}_lc"),
    across(ends_with("_mean_benefit"), ~ . + 1.96 * get(str_replace(cur_column(), "mean_", "se_")), .names = "{.col}_uc")
  ) %>%
  ungroup()

observed <- observed %>%
  rename_with(~ str_replace(., "_mean_benefit_lc", "_lc_benefit"), ends_with("_mean_benefit_lc")) %>%
  rename_with(~ str_replace(., "_mean_benefit_uc", "_uc_benefit"), ends_with("_mean_benefit_uc"))

nnt <- observed %>%
  select(-contains("se_")) %>%
  mutate(across(contains("_benefit"), ~ 1 / ., .names = "{str_replace(.col, '_benefit', '_nnt')}"))

nnt <- nnt %>%
  select(risk_group, unlist(lapply(outcomes, function(outcome) {
    grep(paste0("^", outcome), names(nnt), value = TRUE)
  }))) %>%
  select(risk_group, contains("_nnt")) %>%
  inner_join(pred_nnt, by="risk_group")

## by risk group but with categories merged

pred2 <- cohort %>% mutate(
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
) %>%
  mutate(risk_group=
           ifelse(risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", 
                                    "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"),
                  "uACR ≥30mg/mmol", ifelse(
                    risk_group %in% c("eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol"),
                    "eGFR <60mL/min/1.73m2", as.character(risk_group)
                  )),
         risk_group=factor(risk_group, levels=c(
           "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
           paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat),
           paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat),
           "eGFR <60mL/min/1.73m2",
           "uACR ≥30mg/mmol"
         ))) %>%
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

pred_nnt2 <- pred2 %>%
  mutate(across(contains("_benefit"), ~ 1 / ., .names = "{str_replace(.col, '_benefit', '_nnt')}")) %>%
  select(risk_group, contains("_nnt"))

observed2 <- cohort %>% mutate(
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
) %>%
  mutate(risk_group=
           ifelse(risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", 
                                    "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"),
                  "uACR ≥30mg/mmol", ifelse(
                    risk_group %in% c("eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol"),
                    "eGFR <60mL/min/1.73m2", as.character(risk_group)
                  )),
         risk_group=factor(risk_group, levels=c(
           "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
           paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat),
           paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat),
           "eGFR <60mL/min/1.73m2",
           "uACR ≥30mg/mmol"
         ))) %>%
  select(risk_group, .imp, contains("survdiff")) %>%
  group_by(risk_group, .imp) %>%
  summarise(
    across(starts_with("survdiff_"), list(
      mean_benefit = ~ mean(.x, na.rm = TRUE),
      median_benefit = ~ median(.x, na.rm = TRUE),
      lq_benefit = ~ quantile(.x, probs = 0.25, na.rm = TRUE),
      uq_benefit = ~ quantile(.x, probs = 0.75, na.rm = TRUE)
    ), .names = "{.col}_{.fn}"),
    across(starts_with("se_survdiff_"), list(
      se_benefit = ~ mean(.x, na.rm = TRUE)
    ), .names = "{.col}_{.fn}")
  ) %>%
  ungroup()

observed2 <- observed2 %>%
  rename_with(~ gsub("survdiff_", "", .), starts_with("survdiff_")) %>%
  rename_with(~ gsub("se_survdiff_", "", .), starts_with("se_survdiff_"))


observed2 <- observed2 %>%
  group_by(risk_group) %>%
  summarise(
    across(contains("mean_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(contains("se_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(contains("median_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(contains("lq_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(contains("uq_benefit"), ~ mean(.x, na.rm = TRUE), .names = "{.col}")
  ) %>%
  mutate(
    across(ends_with("_mean_benefit"), ~ . - 1.96 * get(str_replace(cur_column(), "mean_", "se_")), .names = "{.col}_lc"),
    across(ends_with("_mean_benefit"), ~ . + 1.96 * get(str_replace(cur_column(), "mean_", "se_")), .names = "{.col}_uc")
  ) %>%
  ungroup()

observed2 <- observed2 %>%
  rename_with(~ str_replace(., "_mean_benefit_lc", "_lc_benefit"), ends_with("_mean_benefit_lc")) %>%
  rename_with(~ str_replace(., "_mean_benefit_uc", "_uc_benefit"), ends_with("_mean_benefit_uc"))

nnt2 <- observed2 %>%
  select(-contains("se_")) %>%
  mutate(across(contains("_benefit"), ~ 1 / ., .names = "{str_replace(.col, '_benefit', '_nnt')}"))

nnt2 <- nnt2 %>%
  select(risk_group, unlist(lapply(outcomes, function(outcome) {
    grep(paste0("^", outcome), names(nnt), value = TRUE)
  }))) %>%
  select(risk_group, contains("_nnt")) %>%
  inner_join(pred_nnt2, by="risk_group")

nnt_overall <- cohort %>% summarise(median_benefit = median(survdiff_ckd_egfr40)) %>% mutate(nnt = 1/median_benefit) %>% select(nnt) %>% as.numeric()
nnt_guideline_recommended <- cohort %>% filter(risk_group != "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol") %>% summarise(median_benefit = median(survdiff_ckd_egfr40)) %>% mutate(nnt = 1/median_benefit) %>% select(nnt) %>% as.numeric()

############################3 CUMULATIVE INCIDENCE CURVES################################################################
#set default colours (colour-blind accessible from the Okabe-Ito palette) for different drug classes
#this is the colour allocation that will be used across the department
cols <- c("SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00")
#in further analyses, the dpp4/su group will be combined, and we will use the dpp4 colour for this (strongest contrast)
cols <- c(cols, "DPP4i/SU" = "#0072B2")

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
zoom_y <- c(0,0.20)
zoom_x <- c(0,3)
break_y <- 0.025
legend_in <- c(0.35, 0.55)
legend_out <- c(100,-100)

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
) %>%
  mutate(risk_group=
           ifelse(risk_group %in% c("eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", 
                                    "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"),
                  "uACR ≥30mg/mmol", ifelse(
                    risk_group %in% c("eGFR <60mL/min/1.73m2, uACR <3mg/mmol",
                                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol"),
                    "eGFR <60mL/min/1.73m2", as.character(risk_group)
                  )),
         risk_group=factor(risk_group, levels=c(
           "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", 
           paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat),
           paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat),
           "eGFR <60mL/min/1.73m2",
           "uACR ≥30mg/mmol"
         )))

fit_presegfr_noalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                              data = cohort[cohort$risk_group == "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol" &
                                              cohort$.imp == n.imp,],
                              conf.type = "log", conf.int = 0.95)

plots_bottom <- list()

plots_bottom[[1]] <- ggsurvplot(
  fit = fit_presegfr_noalb,
  fun = "cumhaz",
  data = cohort[cohort$risk_group == "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol" &
                  cohort$.imp == n.imp,],
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
  cumevents = F,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "40% decline in eGFR / ESKD",
  title = "",
  subtitle = "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol"
)

plots_bottom[[1]]$plot <- plots_bottom[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
plots_bottom[[1]]$table <- plots_bottom[[1]]$table + theme(plot.title = element_text(size = font_risktable_title))
plots_bottom[[1]]$cumevents <- plots_bottom[[1]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))

plots_top <- list()

fit_presegfr_microalb_lowrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                         data = cohort[cohort$risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat) &
                                                         cohort$.imp == n.imp,],
                                         conf.type = "log", conf.int = 0.95)


plots_top[[1]] <- ggsurvplot(
  fit = fit_presegfr_microalb_lowrisk,
  fun = "cumhaz",
  data = cohort[cohort$risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", low_risk_cat) &
                  cohort$.imp == n.imp,],
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
  cumevents = F,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "",
  ylab = "40% decline in eGFR / ESKD",
  title = low_risk_cat,
  subtitle = "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol"
)

plots_top[[1]]$plot <- plots_top[[1]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
plots_top[[1]]$table <- plots_top[[1]]$table + theme(plot.title = element_text(size = font_risktable_title))
plots_top[[1]]$cumevents <- plots_top[[1]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))



fit_presegfr_microalb_highrisk <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                          data = cohort[cohort$risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat) &
                                                          cohort$.imp == n.imp,],
                                          conf.type = "log", conf.int = 0.95)

plots_top[[2]] <- ggsurvplot(
  fit = fit_presegfr_microalb_highrisk,
  fun = "cumhaz",
  data = cohort[cohort$risk_group == paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, ", high_risk_cat) &
                  cohort$.imp == n.imp,],
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
  cumevents = F,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "",
  ylab = "",
  title = high_risk_cat,
  subtitle = "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol"
) 

plots_top[[2]]$plot <- plots_top[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
plots_top[[2]]$table <- plots_top[[2]]$table + theme(plot.title = element_text(size = font_risktable_title))
plots_top[[2]]$cumevents <- plots_top[[2]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))


fit_macroalb <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                        data = cohort[cohort$risk_group %in% c("uACR ≥30mg/mmol") &
                                        cohort$.imp == n.imp,],
                        conf.type = "log", conf.int = 0.95)

plots_bottom[[3]] <- ggsurvplot(
  fit = fit_macroalb,
  fun = "cumhaz",
  data = cohort[cohort$risk_group %in% c(", uACR ≥30mg/mmol") &
                  cohort$.imp == n.imp,],
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
  cumevents = F,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "",
  title = "",
  subtitle = "uACR ≥30mg/mmol"
) 

plots_bottom[[3]]$plot <- plots_bottom[[3]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
plots_bottom[[3]]$table <- plots_bottom[[3]]$table + theme(plot.title = element_text(size = font_risktable_title))
plots_bottom[[3]]$cumevents <- plots_bottom[[3]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))



fit_redegfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                       data = cohort[cohort$risk_group %in% c("eGFR <60mL/min/1.73m2") &
                                       cohort$.imp == n.imp,],
                       conf.type = "log", conf.int = 0.95)

plots_bottom[[2]] <- ggsurvplot(
  fit = fit_redegfr,
  fun = "cumhaz",
  data = cohort[cohort$risk_group %in% c("eGFR <60mL/min/1.73m2") &
                  cohort$.imp == n.imp,],
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
  cumevents = F,
  tables.height = risktable_height,
  fontsize = font_risktable,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Years",
  ylab = "",
  title = "",
  subtitle = "eGFR <60mL/min/1.73m2"
) 

plots_bottom[[2]]$plot <- plots_bottom[[2]]$plot + coord_cartesian(xlim = zoom_x, ylim = zoom_y)
plots_bottom[[2]]$table <- plots_bottom[[2]]$table + theme(plot.title = element_text(size = font_risktable_title))
plots_bottom[[2]]$cumevents <- plots_bottom[[2]]$cumevents + theme(plot.title = element_text(size = font_risktable_title))



arrange_ggsurvplots(plots_top, print = T, nrow = 1, ncol = 2)
arrange_ggsurvplots(plots_bottom, print = T, nrow = 1, ncol = 3)
