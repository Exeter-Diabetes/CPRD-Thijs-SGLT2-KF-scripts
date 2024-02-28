## In this script, we will select a cohort from CPRD and make it suitable for analyses.
## Contents:
# 0 setup
# 1 cohort selection
# 2 multiple imputation of missing values
# 3 calculate risk scores (for CKD and QRISK2)
# 4 select subjects with valid CKD risk scores (ie. input values within range for calculation of risk scores)
# 5 baseline table and store dataset for further analyses

########################0 SETUP####################################################################

# 0 Setup
library(tidyverse)
library(gt)
library(gtsummary)
library(mice)
library(EHRBiomarkr)
library(tableone)

options(dplyr.summarise.inform = FALSE)

rm(list=ls())

set.seed(123)

########################1 COHORT SELECTION####################################################################

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition_kf function for details)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-02-26_t2d_1stinstance_a.Rda")
load("2024-02-26_t2d_1stinstance_b.Rda")
load("2024-02-26_t2d_all_drug_periods.Rda")

t2d_1stinstance <- rbind(t2d_1stinstance_a, t2d_1stinstance_b)
rm(t2d_1stinstance_a)
rm(t2d_1stinstance_b)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-CKD-scripts/Functions/")
source("cohort_definition_kf.R")
cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

table(cohort$studydrug)
# DPP4 SGLT2    SU 
# 69132 50637 43434

## B Make variables for survival analysis of all endpoints (see survival_variables_kf function for details)

source("survival_variables_kf.R")

cohort <- add_surv_vars(cohort, main_only=FALSE) # add per-protocol survival variables as well


## C Just keep variables of interest

# create acr variable for ckdpc risk scores that uses further source of acr if acr not available
cohort <- cohort %>% 
  mutate(uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
         uacr=ifelse(uacr<0.6, 0.6, uacr),
         
         #and create variable to code whether someone is on oral hyperglycaemic agents
         oha=ifelse(Acarbose+MFN+DPP4+Glinide+GLP1+SGLT2+SU+TZD>add, 1L, 0L))

# add variables that we would want to use later for weights
cohort$statin <- !is.na(cohort$predrug_latest_statins)
cohort$ACEi <- !is.na(cohort$predrug_latest_ace_inhibitors)
cohort$ARB <- !is.na(cohort$predrug_latest_arb)
cohort$BB <- !is.na(cohort$predrug_latest_beta_blockers)
cohort$CCB <- !is.na(cohort$predrug_latest_calcium_channel_blockers)
cohort$ThZD <- !is.na(cohort$predrug_latest_thiazide_diuretics)
cohort$loopD <- !is.na(cohort$predrug_latest_loop_diuretics)
cohort$MRA <- !is.na(cohort$predrug_latest_ksparing_diuretics)
cohort$steroids <- !is.na(cohort$predrug_latest_oralsteroids)
cohort$immunosuppr <- !is.na(cohort$predrug_latest_immunosuppressants)
cohort$osteoporosis <- !is.na(cohort$predrug_latest_osteoporosis)
cohort$later_glp1 <- !is.na(cohort$next_glp1_start)

cohort <- cohort %>%
  #for hba1c take closest result to index date within window of 2 years prior and 7 days post, similar as other biomarkers
  mutate(prehba1c = prehba1c2yrs) %>%
  select(patid, malesex, ethnicity_5cat, ethnicity_qrisk2, imd2015_10, regstartdate, gp_record_end, death_date, 
         drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, 
         MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, height, prehba1c, prebmi, 
         prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, predbp, preegfr, preckdstage, 
         preacr, uacr, qrisk2_smoking_cat, contains("cens"), last_sglt2_stop, oha, next_glp1_start, later_glp1,
         #add variables necessary to calculate qrisk2/qhdf and ckdpc scores
         predrug_fh_premature_cvd, predrug_af, predrug_rheumatoidarthritis, tds_2011,
         predrug_angina, predrug_myocardialinfarction, predrug_stroke, predrug_revasc,
         predrug_heartfailure, predrug_hypertension, predrug_acutepancreatitis,
         predrug_earliest_ace_inhibitors, predrug_earliest_arb,
         predrug_earliest_beta_blockers, predrug_earliest_calcium_channel_blockers,
         predrug_earliest_thiazide_diuretics, predrug_latest_ace_inhibitors, 
         predrug_latest_arb, predrug_latest_beta_blockers, 
         predrug_latest_calcium_channel_blockers, predrug_latest_thiazide_diuretics,
         #add variables necessary for PS weights later on
         predrug_dka, predrug_falls, predrug_urinary_frequency, predrug_volume_depletion, 
         predrug_micturition_control, predrug_dementia, hosp_admission_prev_year, 
         statin, ACEi, ARB, BB, CCB, ThZD, loopD, MRA, steroids, immunosuppr, osteoporosis
  )

rm(list=setdiff(ls(), "cohort"))

# set SU as reference group

cohort$qrisk2_smoking_cat <- as.factor(cohort$qrisk2_smoking_cat)

cohort$studydrug <- relevel(as.factor(cohort$studydrug), ref = "SU")

# create variable for year of treatment initiation
cohort$initiation_year <- substring(as.character(cohort$dstartdate), 1, 4)

# ethnicity cannot be calculated in the imputation model due to it being a constant variable
# for the sake of imputation, we will class missing as a separate category "missing" (5-cat ethnicity: 5; QRISK2: 10)
cohort <- cohort %>%
  mutate(ethnicity_qrisk2=ifelse(is.na(ethnicity_qrisk2), "10", ethnicity_qrisk2),
         ethnicity_5cat=ifelse(is.na(ethnicity_5cat), "5", ethnicity_5cat),
         ethnicity_5cat=factor(ethnicity_5cat,
                               levels = c(0, 1, 2, 3, 4, 5),
                               labels = c("White", "South Asian", "Black", "Other", "Mixed", "Not stated/Unknown"))) %>% 
  relocate(ethnicity_5cat, .after = last_col())



########################2 MULTIPLE IMPUTATION####################################################################
# 2 Impute missing data

# inspect missing data

# #table with n missing and % missing per variable
# miss_data <- cohort %>% 
#   gather(key, value) %>% 
#   group_by(key) %>% 
#   count(na = is.na(value)) %>% 
#   pivot_wider(names_from = na, values_from = n, values_fill = 0) %>% 
#   mutate(pct_missing = (`TRUE`/sum(`TRUE`, `FALSE`))*100) %>% 
#   ungroup()
# 
# miss_data %>% 
#   gt()
# 
# #graph with proportion missing per variable
# miss_data %>% 
#   mutate(Present = 100 - pct_missing) %>% 
#   gather(Key, value, 4:5) %>% 
#   mutate(Key = recode(Key, pct_missing = "Missing")) %>% 
#   ggplot(aes(x = reorder(key, `TRUE`), y = value, fill = Key)) +
#   geom_col(alpha = 0.85) +
#   scale_fill_manual(name = "", 
#                     values = c('tomato3', 'steelblue'), 
#                     labels = c("Missing", "Present")) +
#   coord_flip() +
#   labs(x = NULL, y = "Missing (%)")


#dry run
ini <- mice(cohort, seed = 123, maxit = 0)

# there are a couple of variables that we do not need to impute, we can tell mice not to impute these
meth <- ini$meth

meth[c("death_date", "preacr", "last_sglt2_stop", "preckdstage", "predrug_earliest_ace_inhibitors", 
       "predrug_earliest_arb",
       "predrug_earliest_beta_blockers", "predrug_earliest_calcium_channel_blockers",
       "predrug_earliest_thiazide_diuretics", "predrug_latest_ace_inhibitors", 
       "predrug_latest_arb",
       "predrug_latest_beta_blockers", "predrug_latest_calcium_channel_blockers",
       "ethnicity_qrisk2", 
       "predrug_latest_thiazide_diuretics", "next_glp1_start" )] <- ""

meth["qrisk2_smoking_cat"] <- "polyreg"

meth["prebmi"] <- "~ I( preweight / (height/100)^2)"

# use quickpred function to build predictor matrix
# we can specify which variables to definitely include (inlist) and which ones to leave out (outlist)

inlist <- c("malesex", "imd2015_10", "dstartdate_age",                # main sociodemographic factors
            "studydrug",                                              # treatment variable
            "dstartdate_dm_dur_all", "prebmi", "pretotalcholesterol", # laboratory and vital sign measurements
            "presbp", "preegfr", "uacr", "qrisk2_smoking_cat", 
            "ckd_egfr40_censvar", "death_censvar"     # outcome variables
)


#list variables that are 100% complete and are not interesting for the imputation model
complete_vars <- names(ini$nmis[ini$nmis == 0])
#inspect complete_vars by printing it > print(complete_vars) then choose variables that we want to omit
outlist1 <- c("patid", "gp_record_end", "drugclass", "drugline_all", "ncurrtx", 
              "DPP4", "GLP1", "SGLT2", "SU", "INS", 
              "cens_itt", "cens_pp", "cens_itt_3_yrs", "cens_pp_3_yrs", 
              "ckd_345_censdate", "ckd_345_censtime_yrs", 
              "ckd_egfr40_censdate", "ckd_egfr40_censtime_yrs", 
              "death_censdate", "death_censtime_yrs", 
              "ckd_345_pp_censdate", "ckd_345_pp_censvar", "ckd_345_pp_censtime_yrs", 
              "ckd_egfr40_pp_censdate", "ckd_egfr40_pp_censvar", "ckd_egfr40_pp_censtime_yrs", 
              "death_pp_censdate", "death_pp_censvar", "death_pp_censtime_yrs", 
              "oha", "predrug_angina", "predrug_myocardialinfarction", "predrug_stroke", 
              "predrug_revasc", "predrug_heartfailure", "initiation_year", "ethnicity_5cat")

#list variables with outflux <0.5 
#outflux is an indicator of the potential usefulness for imputing other variables - 
#outflux depends on the proportion of missing data of the variable: 
#outflux of a completely observed variable is equal to 1, 
#whereas outflux of a completely missing variable is equal to 0
#for two variables having the same proportion of missing data, 
#the variable with higher outflux is better connected to the missing data, 
#and thus potentially more useful for imputing other variables)
fx <- flux(cohort)
outlist2 <- row.names(fx)[fx$outflux < 0.5]

#identify problematic variables from initial run (constant/collinear variables)
outlist3 <- as.character(ini$loggedEvents[, "out"])

#combine above variables in one list
outlist <- unique(c(outlist1, outlist2, outlist3))

pred <- quickpred(cohort, include = inlist, exclude = outlist)


# limit imputations to plausible range
post <- ini$post
post["preegfr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(60, 120))"
post["prebmi"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(20, 40))"
post["prehba1c"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(42, 97))"
post["uacr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0.6, 56.5))"
post["presbp"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(80, 180))"

n.imp <- 10

imp <- mice(data = cohort, 
            meth = meth, 
            pred = pred, 
            post = post, 
            m=n.imp, 
            seed = 123)

# check imputed vs original values (disabled as time-consuming)
#densityplot(x = imp, data = ~ imd2015_10 + dstartdate_dm_dur_all + preweight + height + prehba1c + prebmi + 
#              prehdl + preldl + pretriglyceride + pretotalcholesterol + prealt + presbp + predbp + preegfr + 
#              uacr + qrisk2_smoking_cat + tds_2011)

#extract imputations so we can add variables
temp <- complete(imp, action = "long", include = T)

#the post-processing of prebmi imputations does not work as the method is passive imputation
#I will set all imputed prebmi values that are <20 at 20 which is what the post-processing procedure would otherwise do
z <- cohort[is.na(cohort$prebmi),]$patid
temp <- temp %>% mutate(prebmi = ifelse(
  patid %in% z & prebmi < 20, 20, prebmi
))
rm(z)

########################3 CALCULATE RISK SCORES####################################################################

# 3 using imputed data, calculate risk scores (qrisk2, qhdf, ckdpc risk of CKD, kidney disease progression)


# for qrisk2/qhdf scores
temp <- temp %>% mutate(precholhdl=pretotalcholesterol/prehdl,
                        ckd45=ifelse(is.na(preckdstage), ifelse(preegfr < 30, T, F), 
                                     ifelse(preckdstage == "stage_4" | preckdstage == "stage_5", T, F)),
                        cvd=predrug_myocardialinfarction==1 | predrug_angina==1 | predrug_stroke==1,
                        sex=ifelse(malesex==1, "male", "female"),
                        dm_duration_cat=ifelse(dstartdate_dm_dur_all<=1, 0L,
                                               ifelse(dstartdate_dm_dur_all<4, 1L,
                                                      ifelse(dstartdate_dm_dur_all<7, 2L,
                                                             ifelse(dstartdate_dm_dur_all<11, 3L, 4L)))),
                        
                        earliest_bp_med=pmin(
                          ifelse(is.na(predrug_earliest_ace_inhibitors),as.Date("2050-01-01"),predrug_earliest_ace_inhibitors),
                          ifelse(is.na(predrug_earliest_arb),as.Date("2050-01-01"),predrug_earliest_arb),
                          ifelse(is.na(predrug_earliest_beta_blockers),as.Date("2050-01-01"),predrug_earliest_beta_blockers),
                          ifelse(is.na(predrug_earliest_calcium_channel_blockers),as.Date("2050-01-01"),predrug_earliest_calcium_channel_blockers),
                          ifelse(is.na(predrug_earliest_thiazide_diuretics),as.Date("2050-01-01"),predrug_earliest_thiazide_diuretics),
                          na.rm=TRUE
                        ),
                        latest_bp_med=pmax(
                          ifelse(is.na(predrug_latest_ace_inhibitors),as.Date("1900-01-01"),predrug_latest_ace_inhibitors),
                          ifelse(is.na(predrug_latest_arb),as.Date("1900-01-01"),predrug_latest_arb),
                          ifelse(is.na(predrug_latest_beta_blockers),as.Date("1900-01-01"),predrug_latest_beta_blockers),
                          ifelse(is.na(predrug_latest_calcium_channel_blockers),as.Date("1900-01-01"),predrug_latest_calcium_channel_blockers),
                          ifelse(is.na(predrug_latest_thiazide_diuretics),as.Date("1900-01-01"),predrug_latest_thiazide_diuretics),
                          na.rm=TRUE
                        ),
                        bp_meds_qrisk2=ifelse(earliest_bp_med!=as.Date("2050-01-01") & latest_bp_med!=as.Date("1900-01-01") & difftime(dstartdate, latest_bp_med, units="days")<=28 & earliest_bp_med!=latest_bp_med, 1L, 0L),
                        
                        type1=0L,
                        type2=1L,
                        surv_5yr=5L,
                        surv_10yr=10L)
# for ckcpc risk scores:
temp <- temp %>%
  
  mutate(preckdstage=ifelse(is.na(preckdstage), ifelse(preegfr<15, "stage_5", 
                                                       ifelse(preegfr<30, "stage_4", 
                                                              ifelse(preegfr<45, "stage_3b", 
                                                                     ifelse(preegfr<60, "stage_3a", 
                                                                            ifelse(preegfr<90, "stage_2", 
                                                                                   "stage_1"))))), preckdstage),
         
         black_ethnicity=ifelse(!is.na(ethnicity_qrisk2) & (ethnicity_qrisk2 == 6 | ethnicity_qrisk2 == 7), 
                                1L, 
                                ifelse(is.na(ethnicity_qrisk2), NA, 0L)),
         
         cvd=predrug_myocardialinfarction==1 | predrug_revasc==1 | predrug_heartfailure==1 | predrug_stroke==1,
         
         ever_smoker=ifelse(!qrisk2_smoking_cat == 0, 1L, 0L),
         
         bp_meds_ckdpc=ifelse(latest_bp_med!=as.Date("1900-01-01") & difftime(dstartdate, latest_bp_med, units="days")<=183, 1L, 0L),
         
         hypertension=ifelse((!is.na(presbp) & presbp>=140) | (!is.na(predbp) & predbp>=90) | bp_meds_ckdpc==1, 1L,0L),
         
         chd=predrug_myocardialinfarction==1 | predrug_revasc==1,
         
         current_smoker=ifelse(qrisk2_smoking_cat==2 | qrisk2_smoking_cat == 3 | qrisk2_smoking_cat == 4, 1L, 0L),
         
         ex_smoker=ifelse(qrisk2_smoking_cat==1, 1L, 0L))


# recalculate risk scores (qrisk2, qhdf, ckdpc risk scores)

#qrisk + qhdf
temp <- temp %>%
  
  mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
  
  calculate_qdiabeteshf(sex=sex2, age=dstartdate_age, ethrisk=ethnicity_qrisk2, smoking=qrisk2_smoking_cat, duration=dm_duration_cat, type1=type1, cvd=cvd, renal=ckd45, af=predrug_af, hba1c=prehba1c, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, town=tds_2011, surv=surv_5yr) %>%
  
  mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
  
  calculate_qrisk2(sex=sex2, age=dstartdate_age, ethrisk=ethnicity_qrisk2, smoking=qrisk2_smoking_cat, type1=type1, type2=type2, fh_cvd=predrug_fh_premature_cvd, renal=ckd45, af=predrug_af, rheumatoid_arth=predrug_rheumatoidarthritis, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, bp_med=bp_meds_qrisk2, town=tds_2011, surv=surv_5yr) %>%
  
  rename(qrisk2_score_5yr=qrisk2_score) %>%
  
  select(-qrisk2_lin_predictor) %>%
  
  calculate_qrisk2(sex=sex2, age=dstartdate_age, ethrisk=ethnicity_qrisk2, smoking=qrisk2_smoking_cat, type1=type1, type2=type2, fh_cvd=predrug_fh_premature_cvd, renal=ckd45, af=predrug_af, rheumatoid_arth=predrug_rheumatoidarthritis, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, bp_med=bp_meds_qrisk2, town=tds_2011, surv=surv_10yr) %>%
  
  rename(qrisk2_score_10yr=qrisk2_score) 

#ckdpc incident ckd and ckd progression risk scores
temp <- temp %>% 
  
  mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
  
  calculate_ckdpc_egfr60_risk(age=dstartdate_age, sex=sex2, black_eth=black_ethnicity, egfr=preegfr, cvd=cvd, hba1c=prehba1c, insulin=INS, oha=oha, ever_smoker=ever_smoker, hypertension=hypertension, bmi=prebmi, acr=uacr, remote=FALSE) %>%
  
  calculate_ckdpc_40egfr_risk(age=dstartdate_age, sex=sex2, egfr=preegfr, acr=uacr, sbp=presbp, bp_meds=bp_meds_ckdpc, hf=predrug_heartfailure, chd=chd, af=predrug_af, current_smoker=current_smoker, ex_smoker=ex_smoker, bmi=prebmi, hba1c=prehba1c, oha=oha, insulin=INS, remote=FALSE) 


########################4 REMOVE RISK SCORE VALUES OUTSIDE OF RANGE####################################################################

# 4 remove risk scores for subjects whose characteristics are outside of the reference range for the respective equations

# Look at counts of people whose characteristics are outside the reference range for the scores

ckdpc_outofrange <- temp[temp$.imp > 0,] %>%
  mutate(uacr_val=ifelse(uacr < 56.5, "in range",
                         ifelse(uacr < 113, "in range for 40eGFR but above range for CKD60",
                                "above range for both scores")),
         egfr_under_60=ifelse(!is.na(preegfr) & preegfr<60, 1, 0),
         hba1c_val=ifelse(prehba1c<42, "under",
                          ifelse(prehba1c>97, "over", "in range")),
         bmi_val=ifelse(prebmi<20, "under",
                        ifelse(prebmi>40, "over", "in range")))

prop.table(table(ckdpc_outofrange$uacr_val))
prop.table(table(ckdpc_outofrange$egfr_under_60))
prop.table(table(ckdpc_outofrange$hba1c_val))
prop.table(table(ckdpc_outofrange$bmi_val))

ckdpc_outofrange %>% filter(bmi_val=="over") %>% summarise(median=median(prebmi))
ckdpc_outofrange %>% filter(hba1c_val=="over") %>% summarise(median=median(prehba1c))



## Remove QDiabetes-HF score for those with biomarker values outside of range:
### CholHDL: missing or 1-11 (NOT 12)
### HbA1c: 40-150
### SBP: missing or 70-210
### Age: 25-84
### Also exclude if BMI<20 as v. different from development cohort

## Remove QRISK2 score for those with biomarker values outside of range:
### CholHDL: missing or 1-12
### SBP: missing or 70-210
### Age: 25-84
### Also exclude if BMI<20 as v. different from development cohort
temp <- temp %>% mutate(qdiabeteshf_5yr_score=ifelse((is.na(precholhdl) | (precholhdl>=1 & precholhdl<=11)) &
                                                       prehba1c>=40 & prehba1c<=150 &
                                                       (is.na(presbp) | (presbp>=70 & presbp<=210)) &
                                                       dstartdate_age>=25 & dstartdate_age<=84 &
                                                       prebmi>=20, qdiabeteshf_score, NA),
                        
                        qdiabeteshf_lin_predictor=ifelse((is.na(precholhdl) | (precholhdl>=1 & precholhdl<=11)) &
                                                           prehba1c>=40 & prehba1c<=150 &
                                                           (is.na(presbp) | (presbp>=70 & presbp<=210)) &
                                                           dstartdate_age>=25 & dstartdate_age<=84 &
                                                           prebmi>=20, qdiabeteshf_lin_predictor, NA),
                        
                        qrisk2_5yr_score=ifelse((is.na(precholhdl) | (precholhdl>=1 & precholhdl<=12)) &
                                                  (is.na(presbp) | (presbp>=70 & presbp<=210)) &
                                                  dstartdate_age>=25 & dstartdate_age<=84 &
                                                  prebmi>=20, qrisk2_score_5yr, NA),
                        
                        qrisk2_10yr_score=ifelse((is.na(precholhdl) | (precholhdl>=1 & precholhdl<=12)) &
                                                   (is.na(presbp) | (presbp>=70 & presbp<=210)) &
                                                   dstartdate_age>=25 & dstartdate_age<=84 &
                                                   prebmi>=20, qrisk2_score_10yr, NA),
                        
                        qrisk2_lin_predictor=ifelse((is.na(precholhdl) | (precholhdl>=1 & precholhdl<=12)) &
                                                      (is.na(presbp) | (presbp>=70 & presbp<=210)) &
                                                      dstartdate_age>=25 & dstartdate_age<=84 &
                                                      prebmi>=20, qrisk2_lin_predictor, NA))

## Also remove eGFR<60 score for those with sociodemographic/vital/laboratory measurements outside of range:
### Age: 20-80
### UACR: 0.6-56.5 (5-500 in mg/g)
### BMI: <20
### HbA1c 42-97 (6-11 in %)

## Also remove 40% decline in eGFR score for those with sociodemographic/vital/laboratory measurements outside of range:
### Age: 20-80
### UACR: 0.6-113 (5-1000 in mg/g)
### SBP: 80-180
### BMI: <20
### HbA1c 42-97 (6-11 in %)

temp <- temp %>%
  
  mutate(across(starts_with("ckdpc_egfr60"),
                ~ifelse(dstartdate_age>=20 & dstartdate_age<=80 &
                          prebmi>=20, .x, NA))) %>%
  
  mutate(across(starts_with("ckdpc_40egfr"),
                ~ifelse(dstartdate_age>=20 & dstartdate_age<=80 &
                          prebmi>=20, .x, NA)))

# retain those with available risk scores only
temp <- temp %>% filter(!is.na(ckdpc_40egfr_score))

# those left out:
ckdpc_outofrange <- ckdpc_outofrange %>% anti_join(temp, by = c("patid", ".imp"))

q <- ckdpc_outofrange %>% .$patid %>% unique() %>% length()
print(paste0("Number of subjects excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q))

q <- ckdpc_outofrange %>% nrow()
print(paste0("Number of drug episodes excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q/n.imp))

q <- temp %>% .$patid %>% unique() %>% length()
print(paste0("Number of subjects in study population ", q))
q <- temp %>% nrow()
print(paste0("Number of drug episodes in study population ", q/n.imp))

########################5 SAVE DATASET####################################################################

# 5 tabulate and save imputed dataset for further analyses


# add variable for "high CV risk" according to ADA

temp <- temp %>% mutate(
  obesity = ifelse(prebmi < 30, F, T),
  smoking_hx = ifelse(qrisk2_smoking_cat == 0, F, T),
  albuminuria = ifelse(uacr < 3, F, T),
  dyslipidaemia = ifelse(pretotalcholesterol < 5 &
                           preldl < 4 &
                           pretriglyceride < 2.3, F, T),
  cv_high_risk = ifelse(dstartdate_age >= 55 & 
                          obesity + predrug_hypertension + 
                          smoking_hx + dyslipidaemia +
                          albuminuria > 1, T, F),
  qrisk2_above_10_pct = ifelse(qrisk2_10yr_score >= 10, T, F),
  egfr_below_60 = ifelse(preckdstage %in% c("stage_1", "stage_2") & preegfr >=60, F, T),
  ACEi_or_ARB = ifelse(temp$ACEi + temp$ARB > 0, T, F),
  macroalbuminuria = ifelse(uacr < 30, F, T),
  microalbuminuria = ifelse(uacr <3, F, ifelse(macroalbuminuria == T, F, T))
)

# create table one: this will be an average of the imputed datasets (n to be divided by n.imp)

#variables to be shown
vars <- c("dstartdate_age", "malesex", "ethnicity_5cat", "imd2015_10",             # sociodemographic variables
          "prebmi", "preldl", "prehba1c", "presbp", "predbp", "preegfr",           # vital signs and laboratory measurements
          "preckdstage", "uacr", "albuminuria", "microalbuminuria", "macroalbuminuria",
          "dstartdate_dm_dur_all", "qrisk2_smoking_cat", "predrug_hypertension",   # comorbidities
          "predrug_af", "predrug_dka", "osteoporosis", 
          "predrug_acutepancreatitis", "predrug_falls", 
          "predrug_urinary_frequency", "predrug_volume_depletion", 
          "predrug_micturition_control", "predrug_dementia", "hosp_admission_prev_year",
          "initiation_year",
          "ncurrtx", "MFN", "TZD", "INS", "ACEi_or_ARB",                                   # medications
          "cv_high_risk", "qrisk2_above_10_pct"                                     # CV risk
          
)

#categorical variables
factors <- c("malesex", "ethnicity_qrisk2", "qrisk2_smoking_cat", "predrug_hypertension", 
             "predrug_af", "predrug_dka", "osteoporosis", "predrug_acutepancreatitis", 
             "predrug_falls", "predrug_urinary_frequency", "predrug_volume_depletion", 
             "predrug_micturition_control", "predrug_dementia", "hosp_admission_prev_year",
             "initiation_year",
             "ncurrtx", "MFN", "TZD", "INS", "ACEi_or_ARB",
             "cv_high_risk", "qrisk2_above_10_pct", 
             "preckdstage", "albuminuria", "microalbuminuria", "macroalbuminuria")

nonnormal <- c("imd2015_10", "uacr", "dstartdate_dm_dur_all")

table <- CreateTableOne(vars = vars, strata = "studydrug", data = temp[temp$.imp > 0,], 
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)
## Save to a CSV file
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv
today <- as.character(Sys.Date(), format="%Y%m%d")
write.csv2(tabforprint, file = paste0(today, "_baseline_table_incl_egfr_below_60.csv"))

# save imputed dataset so this can be used in the subsequent scripts
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(temp, file=paste0(today, "_t2d_ckdpc_imputed_data_incl_egfr_below_60.Rda"))
