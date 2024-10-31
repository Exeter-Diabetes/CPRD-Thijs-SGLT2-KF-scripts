
############################0 SETUP################################################################

# Setup
library(tidyverse)
library(survival)
library(survminer)
library(rms)
library(broom)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())

# set seed
set.seed(123)
n.imp <- 10
n.quantiles <- 10

#today <- as.character(Sys.Date(), format="%Y%m%d")
today <- "2024-10-29"

# covariates for multivariable adjustment
covariates <- c("dstartdate_age", "malesex", "imd2015_10", "ethnicity_4cat", "initiation_year", "prebmi", "prehba1c",
                "pretotalcholesterol", "preegfr", "uacr", "presbp", "ckdpc_50egfr_score_cal", "ncurrtx", "statin", "INS", 
                "ACEi_or_ARB", "smoking_status", "dstartdate_dm_dur_all", "predrug_hypertension", "predrug_af", "hosp_admission_prev_year")

# define outcomes to be analysed
outcomes <- c("ckd_egfr40", "ckd_egfr50", "macroalb", "dka", "side_effect", "ckd_egfr50_5y"#, "death", "amputation"
)
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
load("2024-10-29_t2d_ckdpc_recalibrated.Rda")

############################1 PREPARE DATASET################################################################


cohort$studydrug2 <- as.factor(cohort$studydrug2)

# calculate predicted sglt2 benefit (absolute risk reduction = ARR):
# ARR = S0(t)^HR - S0(t)
trial_hr_kf_sglt2i <- 0.62

cohort <- cohort %>% 
  mutate(ckdpc_50egfr_survival_cal=(100-ckdpc_50egfr_score_cal)/100,
         ckdpc_50egfr_survival_cal_sglt2i=ckdpc_50egfr_survival_cal^trial_hr_kf_sglt2i,
         ckdpc_50egfr_sglt2i_benefit=ckdpc_50egfr_survival_cal_sglt2i - ckdpc_50egfr_survival_cal)

cohort$benefit_decile <- ntile(cohort$ckdpc_50egfr_sglt2i_benefit, n.quantiles)

cohort <- cohort %>% mutate(
  ethnicity_4cat = ifelse(!ethnicity_4cat %in% c("White", "South Asian", "Black"), "Other", as.character(ethnicity_4cat)),
  initiation_year = ifelse(initiation_year %in% c("2019", "2020"), "2019-2020", as.character(initiation_year)),
  ncurrtx = ifelse(ncurrtx %in% c("3.", "4+"), "3+", as.character(ncurrtx))
)

cohort <- cohort %>% filter(!.imp > n.imp)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated_with_riskgroup.Rda"))

## prep data for estimating absolute benefit:
# centre predictors and save dataset (performing analysis in entire dataset exceeds memory limit)

centre_and_reference <- function(df, covariates) {
  df %>%
    mutate(across(all_of(covariates), ~ if(is.numeric(.)) . - mean(., na.rm = TRUE) else . ),  # Center numeric variables
           across(all_of(covariates), ~ if(is.logical(.)) . == names(sort(table(.), decreasing = TRUE))[1] else . ) ,  # Set most frequent level as reference for logical
           across(all_of(covariates), ~ if(is.factor(.)) relevel(., ref = names(sort(table(.), decreasing = TRUE))[1]) else . ))  # Set most frequent level as reference for factor
}

# create regex pattern of censoring variables to select
outcome_variables <- paste0("(", paste(outcomes, collapse = "|"), ")(?!.*(5y|pp)).*(_censtime_yrs|_censvar)$")
outcome_variables_5y <- paste0("(", paste("ckd_egfr50", collapse = "|"), ")(?!.*(pp)).*(_censtime_yrs|_censvar)$")

cohort_5y <- cohort %>%
  mutate(across(contains("predrug_"), as.logical),
         hosp_admission_prev_year=as.logical(hosp_admission_prev_year),
         INS=as.logical(INS),
         MFN=as.logical(MFN),
         malesex=as.factor(malesex),
         initiation_year=as.factor(initiation_year),
         benefit_decile=as.factor(benefit_decile)) %>%
  select(patid, .imp, risk_group, studydrug2, benefit_decile,
         matches(outcome_variables_5y, perl = TRUE), 
         all_of(covariates)) %>% 
  centre_and_reference(covariates)

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
save(cohort_5y, file=paste0(today, "_recalibrated_data_centred_predictors_5y.Rda"))

############################2 CALCULATE ADJUSTED OBSERVED BENEFITS################################################################

rm(list = setdiff(ls(), c("n.imp", "covariates", "today", "outcomes")))

## estimate absolute benefit
for (k in outcomes) {
  
  for (i in 1:n.imp) {
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
    load("2024-10-29_recalibrated_data_centred_predictors.Rda")
    
    if (k == "macroalb") {
      # remove subjects with established macroalbuminuria from these analyses
      cohort <- cohort %>% filter(.imp == i & !risk_group %in% c("uACR ≥30mg/mmol") & macroalb_censtime_yrs >= 0)
    } else if (k == "ckd_egfr50_5y") {
      load("2024-10-29_recalibrated_data_centred_predictors_5y.Rda")
      cohort <- cohort_5y %>% filter(.imp == i)
      rm(cohort_5y)
    } else {
      cohort <- cohort %>% filter(.imp == i)
    }
    
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")
    
    # fit multivariable-adjusted model
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", paste(covariates, collapse=" + "))) 
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE)
    
    # create dataframe with similar covariate distribution as our cohort but with everyone treated with SGLT2i
    obs_SGLT2 <- cohort %>%
      mutate(studydrug_original=studydrug2,
             studydrug2="SGLT2i",
             rowno=row_number())
    print(paste0("Survival estimates for SGLT2i in imputation ", i, "  (outcome ", k, ")"))
    
    # get multivariable-adjusted survival estimates
    observed_sglt2 <- survfit(model, newdata=as.data.frame(obs_SGLT2)) %>%
      tidy() %>%
      filter(time == if (k == "ckd_egfr50_5y") 5 else 3) %>%
      pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
      select(group, estimate_sglt2=estimate, se_sglt2=std.error) %>%
      mutate(group=as.numeric(group)) %>%
      inner_join(obs_SGLT2, by=c("group"="rowno"))
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
    # save estimates
    save(observed_sglt2, file=paste0(today, "_adjusted_surv_",k,"_SGLT2i_imp.", i, ".Rda"))
    
    rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "outcomes")))
  }
  
  # similar for DPP4i/SU:
  
  for (i in 1:n.imp) {
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
    load("2024-10-29_recalibrated_data_centred_predictors.Rda")
    
    if (k == "macroalb") {
      cohort <- cohort %>% filter(.imp == i & !risk_group %in% c("uACR ≥30mg/mmol") & macroalb_censtime_yrs >= 0)
    } else if (k == "ckd_egfr50_5y") {
      load("2024-10-29_recalibrated_data_centred_predictors_5y.Rda")
      cohort <- cohort_5y %>% filter(.imp == i)
      rm(cohort_5y)
    } else {
      cohort <- cohort %>% filter(.imp == i)
    }
    
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")  
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", paste(covariates, collapse=" + ")))
    
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE)
    
    obs_DPP4SU <- cohort %>%
      mutate(studydrug_original=studydrug2,
             studydrug2="DPP4i/SU",
             rowno=row_number())
    
    print(paste0("Survival estimates for DPP4i/SU in imputation ", i, "  (outcome ", k, ")"))
    
    observed_dpp4su <- survfit(model, newdata=as.data.frame(obs_DPP4SU)) %>%
      tidy() %>%
      filter(time == if (k == "ckd_egfr50_5y") 5 else 3) %>%
      pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
      select(group, estimate_dpp4su=estimate, se_dpp4su=std.error) %>%
      mutate(group=as.numeric(group)) %>%
      inner_join(obs_DPP4SU, by=c("group"="rowno"))
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
    save(observed_dpp4su, file=paste0(today, "_adjusted_surv_",k,"_DPP4iSU_imp.", i, ".Rda"))
    
    rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "outcomes")))
  }
  
  temp_sglt2 <- temp_dpp4su <- data.frame()
  
  #for every outcome, join survival estimates from each imputation in one dataframe
  for (i in 1:n.imp) {
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
    load(paste0("2024-10-29_adjusted_surv_",k,"_SGLT2i_imp.", i, ".Rda"))
    load(paste0("2024-10-29_adjusted_surv_",k,"_DPP4iSU_imp.", i, ".Rda"))
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

### add adjusted (observed) survival estimates to main dataset
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-10-29_t2d_ckdpc_recalibrated_with_riskgroup.Rda")

for (k in outcomes) {
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
  load(paste0("2024-10-29_adjusted_surv_", k, ".Rda"))
  cohort <- cohort %>% left_join(benefits %>%
                                   select(.imp, patid, studydrug2, contains("survdiff")), 
                                 by=c(".imp", "patid", "studydrug2"))
  rm(benefits)
}

cohort$studydrug2 <- as.factor(cohort$studydrug2)

# save dataset with adjusted survival probabilities

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_recalibrated_with_adjsurv.Rda"))


