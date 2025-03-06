# in this script we aim to derive the counterfactual absolute risk differences between SGLT2i and DPP4i/SU.
# these are essentially the observed absolute risk differences that we want to validate the predictions against.
# in order to do this, we will fit a marginal structural model - a weighted cox model. 
# we will then apply this marginal structural model to two counterfactual datasets - two copies of the dataset with treatment for all individuals set as "SGLT2i" and "DPP4i/SU" respectively.
# the counterfactual absolute risk differences will then be stored for use later.
# given the computational power required, the script is designed to minimise the use of working memory by only loading one dataset at a time with minimal columns.


############################0 SETUP################################################################

# Setup
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/")
source("00 Setup.R")

############################1 PREPARE DATASET################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
load(paste0(today, "_t2d_ckdpc_imputed_data_withweights.Rda"))

cohort$studydrug2 <- as.factor(cohort$studydrug2)

## prep data:
# centre predictors and retain relevant columns only (performing analysis in entire dataset exceeds memory limit)

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
         initiation_year=as.factor(initiation_year)) %>%
  select(patid, .imp, albuminuria, studydrug2, overlap2,
         matches(outcome_variables_5y, perl = TRUE), 
         all_of(covariates)) %>% 
  centre_and_reference(covariates)

cohort <- cohort %>%
  mutate(across(contains("predrug_"), as.logical),
         hosp_admission_prev_year=as.logical(hosp_admission_prev_year),
         INS=as.logical(INS),
         MFN=as.logical(MFN),
         malesex=as.factor(malesex),
         initiation_year=as.factor(initiation_year)) %>%
  select(patid, .imp, albuminuria, studydrug2, overlap2,
         matches(outcome_variables, perl = TRUE), 
         all_of(covariates)) %>% 
  centre_and_reference(covariates)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
save(cohort, file=paste0(today, "_data_centred_predictors.Rda"))
save(cohort_5y, file=paste0(today, "_data_centred_predictors_5y.Rda"))

############################2 FIT WEIGHTED COX MODEL AND ESTIMATE COUNTERFACTUAL OBSERVED SURVIVAL################################################################

# clear R memory to ensure memory limit not exceeded
rm(list = setdiff(ls(), c("n.imp", "covariates", "today")))

# if evaluating other risk scores, other relevant outcomes may be added
outcomes_msm <- "ckd_egfr50"

for (k in outcomes_msm) {
  
  for (i in 1:n.imp) {
    # load minimal dataset
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
    load(paste0(today, "_data_centred_predictors.Rda"))
    
    cohort <- cohort %>% filter(.imp == i)

    # clear R memory
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")
    
    # fit overlap-weighted cox model
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", paste(covariates, collapse=" + "))) 
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE, weights=overlap2)
    
    # create counterfactual dataset where everyone is treated with SGLT2i
    obs_SGLT2 <- cohort %>%
      # create variable for factual study arm for reference later
      mutate(studydrug_original=studydrug2,
             studydrug2="SGLT2i",
             rowno=row_number())
    print(paste0("Survival estimates for SGLT2i in imputation ", i, "  (outcome ", k, ")"))
    
    # get counterfactual (observed weighted) survival
    observed_sglt2 <- survfit(model, newdata=as.data.frame(obs_SGLT2)) %>%
      tidy() %>%
      filter(time == if (k == "ckd_egfr50_5y") 5 else 3) %>%
      pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
      select(group, estimate_sglt2=estimate, se_sglt2=std.error) %>%
      mutate(group=as.numeric(group)) %>%
      inner_join(obs_SGLT2, by=c("group"="rowno"))
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
    # save results
    save(observed_sglt2, file=paste0(today, "_adjusted_surv_",k,"_SGLT2i_imp.", i, ".Rda"))
    
    # clear environment
    rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "outcomes_msm")))
  }
  
  # repeat for DPP4i/SU:
  
  for (i in 1:n.imp) {
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
    load(paste0(today, "_data_centred_predictors.Rda"))
    
    cohort <- cohort %>% filter(.imp == i)
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")  
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", paste(covariates, collapse=" + ")))
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE, weights = overlap2)
    
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
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
    save(observed_dpp4su, file=paste0(today, "_adjusted_surv_",k,"_DPP4iSU_imp.", i, ".Rda"))
    
    rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "outcomes_msm")))
  }
  
  # after separately estimating counterfactual survival probabilities in each imputation, we will now combine these in the main dataset
  # create empty dataframes to append each imputed dataset to
  temp_sglt2 <- temp_dpp4su <- data.frame()
  
  #for every outcome, join estimates from each imputation in one dataframe
  for (i in 1:n.imp) {
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
    load(paste0(today, "_adjusted_surv_",k,"_SGLT2i_imp.", i, ".Rda"))
    load(paste0(today, "_adjusted_surv_",k,"_DPP4iSU_imp.", i, ".Rda"))
    temp_sglt2 <- temp_sglt2 %>% rbind(observed_sglt2)
    temp_dpp4su <- temp_dpp4su %>% rbind(observed_dpp4su)
    rm(observed_sglt2)
    rm(observed_dpp4su)
  }
  
  # combine counterfactual observed survival probabilities
  benefits <- temp_sglt2 %>%
    select(group, .imp, estimate_sglt2, se_sglt2) %>%
    inner_join(temp_dpp4su, by=c("group", ".imp")) %>%
    # remove counterfactual studydrug annotations and replace with original study arm
    select(-studydrug2) %>%
    # calculate counterfactual observed absolute risk difference
    mutate(survdiff=estimate_sglt2-estimate_dpp4su,
           se_survdiff=sqrt(se_sglt2^2 + se_dpp4su^2),
           studydrug2=studydrug_original) %>%
    select(-studydrug_original) 
  
  #rename variables by outcome name
  benefits <- benefits %>% rename_with(
    ~ paste0(.x, "_", k),
    contains("survdiff")
  )
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
  save(benefits, file=paste0(today, "_adjusted_surv_",k,".Rda"))
  rm(benefits)
}

rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "outcomes_msm")))

### add counterfactual observed absolute risk reductions to main dataset
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
load(paste0(today, "_t2d_ckdpc_imputed_data_withweights.Rda"))

for (k in outcomes_msm) {
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
  load(paste0(today, "_adjusted_surv_", k, ".Rda"))
  cohort <- cohort %>% left_join(benefits %>%
                                   select(.imp, patid, studydrug2, contains("survdiff")), 
                                 by=c(".imp", "patid", "studydrug2"))
  rm(benefits)
}

cohort$studydrug2 <- as.factor(cohort$studydrug2)

# save dataset

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_data_with_adjsurv.Rda"))


