## select cohort, impute missing values, and save imputed dataset for further analyses
 
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
load("2023-10-31_t2d_1stinstance_a.Rda")
load("2023-10-31_t2d_1stinstance_b.Rda")
load("2023-10-31_t2d_all_drug_periods.Rda")

t2d_1stinstance <- rbind(t2d_1stinstance_a, t2d_1stinstance_b)
rm(t2d_1stinstance_a)
rm(t2d_1stinstance_b)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Functions/")
source("cohort_definition_kf.R")
cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

table(cohort$studydrug)
#DPP4 SGLT2    SU 
#48387 21163 49055

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
  
  select(patid, malesex, ethnicity_qrisk2, imd2015_10, regstartdate, gp_record_end, death_date, 
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

# ethnicity cannot be calculated in the imputation model due to it being a constant variable
# for the sake of imputation, we will class missing as "missing" (10)
cohort <- cohort %>%
  mutate(ethnicity_qrisk2=ifelse(is.na(ethnicity_qrisk2), "10", ethnicity_qrisk2))

# select cohort with at least 1 year of data prior to start date and who did not have SGLT2 in year prior to start date

#one_year_reg_cohort <- cohort %>%
#  mutate(reg_before_drug_start=as.numeric(difftime(dstartdate, regstartdate, units="days")),
#         time_since_sglt2_stop=as.numeric(difftime(dstartdate, last_sglt2_stop, units="days"))) %>%
#  filter(reg_before_drug_start>=365 & (is.na(last_sglt2_stop) | last_sglt2_stop>=365))

#table(one_year_reg_cohort$studydrug)

########################2 MULTIPLE IMPUTATION####################################################################
# 2 Impute missing data and recalculate risk scores

# inspect missing data

#table with n missing and % missing per variable
miss_data <- cohort %>% 
  gather(key, value) %>% 
  group_by(key) %>% 
  count(na = is.na(value)) %>% 
  pivot_wider(names_from = na, values_from = n, values_fill = 0) %>% 
  mutate(pct_missing = (`TRUE`/sum(`TRUE`, `FALSE`))*100) %>% 
  ungroup()

miss_data %>% 
  gt()

#graph with proportion missing per variable
miss_data %>% 
  mutate(Present = 100 - pct_missing) %>% 
  gather(Key, value, 4:5) %>% 
  mutate(Key = recode(Key, pct_missing = "Missing")) %>% 
  ggplot(aes(x = reorder(key, `TRUE`), y = value, fill = Key)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(name = "", 
                    values = c('tomato3', 'steelblue'), 
                    labels = c("Missing", "Present")) +
  coord_flip() +
  labs(x = NULL, y = "Missing (%)")


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
       "predrug_latest_thiazide_diuretics", "preacr", "next_glp1_start" )] <- ""

meth["qrisk2_smoking_cat"] <- "polyreg"

meth["prebmi"] <- "~ I( preweight / (height/100)^2)"

# use quickpred function to build predictor matrix
# we can specify which variables to definitely include (inlist) and which ones to leave out (outlist)

inlist <- c("malesex", "imd2015_10", "dstartdate_age",                # main sociodemographic factors
            "studydrug",                                              # treatment variable
            "dstartdate_dm_dur_all", "prebmi", "pretotalcholesterol", # laboratory and vital sign measurements
            "presbp", "preegfr", "uacr", "qrisk2_smoking_cat", 
            "mace_censvar", "ckd_egfr40_censvar", "death_censvar"     # outcome variables
            )


#list variables that are 100% complete and are not interesting for the imputation model
complete_vars <- names(ini$nmis[ini$nmis == 0])
outlist1 <- complete_vars[c(1, 5:6, 10:13, 15:16, 18, 20:24, 26:30, 32:33, 35:36, 38:39, 41:42, 44:72, 75, 77:81)]

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

n.imp <- 10

imp <- mice(data = cohort, meth = meth, pred = pred, m=n.imp, seed = 123)

# check imputed vs original values
#densityplot(x = imp, data = ~ imd2015_10 + dstartdate_dm_dur_all + preweight + height + prehba1c + prebmi + 
#              prehdl + preldl + pretriglyceride + pretotalcholesterol + prealt + presbp + predbp + preegfr + 
#              uacr + qrisk2_smoking_cat + tds_2011)

#extract imputations so we can add variables
temp <- complete(imp, action = "long", include = T)

########################3 CALCULATE RISK SCORES####################################################################

# 3 using imputed data, calculate risk scores (qrisk2, qhdf, ckdpc risk of CKD, kidney disease progression)

# for the sake of risk score calculation, class missing ethnicity as other (9)
temp$ethnicity_qrisk2_2 <- NA
temp[temp$.imp >0,] <- temp[temp$.imp >0,] %>%
  mutate(ethnicity_qrisk2_2=ifelse(ethnicity_qrisk2 == "10", "9", ethnicity_qrisk2))

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
         
         black_ethnicity=ifelse(!is.na(ethnicity_qrisk2_2) & (ethnicity_qrisk2_2 == 6 | ethnicity_qrisk2_2 == 7), 
                                1L, 
                                ifelse(is.na(ethnicity_qrisk2_2), NA, 0L)),
         
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
  
  calculate_qdiabeteshf(sex=sex2, age=dstartdate_age, ethrisk=ethnicity_qrisk2_2, smoking=qrisk2_smoking_cat, duration=dm_duration_cat, type1=type1, cvd=cvd, renal=ckd45, af=predrug_af, hba1c=prehba1c, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, town=tds_2011, surv=surv_5yr) %>%

  mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
  
  calculate_qrisk2(sex=sex2, age=dstartdate_age, ethrisk=ethnicity_qrisk2_2, smoking=qrisk2_smoking_cat, type1=type1, type2=type2, fh_cvd=predrug_fh_premature_cvd, renal=ckd45, af=predrug_af, rheumatoid_arth=predrug_rheumatoidarthritis, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, bp_med=bp_meds_qrisk2, town=tds_2011, surv=surv_5yr) %>%
  
  rename(qrisk2_score_5yr=qrisk2_score) %>%
  
  select(-qrisk2_lin_predictor) %>%
  
  calculate_qrisk2(sex=sex2, age=dstartdate_age, ethrisk=ethnicity_qrisk2_2, smoking=qrisk2_smoking_cat, type1=type1, type2=type2, fh_cvd=predrug_fh_premature_cvd, renal=ckd45, af=predrug_af, rheumatoid_arth=predrug_rheumatoidarthritis, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, bp_med=bp_meds_qrisk2, town=tds_2011, surv=surv_10yr) %>%
  
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
  ### BMI: 20-40
  ### HbA1c 42-97 (6-11 in %)
  
  ## Also remove 40% decline in eGFR score for those with sociodemographic/vital/laboratory measurements outside of range:
  ### Age: 20-80
  ### UACR: 0.6-113 (5-1000 in mg/g)
  ### SBP: 80-180
  ### BMI: 20-40
  ### HbA1c 42-97 (6-11 in %)
  
temp <- temp %>%

  mutate(across(starts_with("ckdpc_egfr60"),
                ~ifelse((is.na(preckdstage) | preckdstage=="stage_1" | preckdstage=="stage_2") &
                          (is.na(preegfr) | preegfr>=60) &
                          dstartdate_age>=20 & dstartdate_age<=80 &
                          prebmi>=20, .x, NA))) %>%
  
  mutate(across(starts_with("ckdpc_40egfr"),
                ~ifelse((is.na(preckdstage) | preckdstage=="stage_1" | preckdstage=="stage_2") &
                          (is.na(preegfr) | preegfr>=60) &
                          dstartdate_age>=20 & dstartdate_age<=80 &
                          prebmi>=20, .x, NA)))

# retain those with available risk scores only
temp <- temp %>% filter(!is.na(ckdpc_egfr60_confirmed_score) & !is.na(ckdpc_40egfr_score))

# those left out:
ckdpc_outofrange <- ckdpc_outofrange[!ckdpc_outofrange$patid %in% temp$patid,]

q <- ckdpc_outofrange %>% summarise(n=n()/n.imp)

print(paste0("Number of subjects excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q))


########################5 SAVE DATASET####################################################################

# 5 tabulate and save imputed dataset for further analyses


# create table one: this will be an average of the imputed datasets (n to be divided by n.imp)

#variables to be shown
vars <- c("dstartdate_age", "malesex", "ethnicity_qrisk2", "imd2015_10", 
          "prebmi", "preegfr", "uacr", "preldl", "prehba1c", "presbp", "predbp",
          "dstartdate_dm_dur_all", "qrisk2_smoking_cat", "predrug_hypertension", 
          "predrug_dka", "osteoporosis", "predrug_acutepancreatitis", "predrug_falls", 
          "predrug_urinary_frequency", "predrug_volume_depletion", 
          "predrug_micturition_control", "predrug_dementia", "hosp_admission_prev_year",
          "ncurrtx", "MFN", "TZD", "ACEi", "ARB")

#categorical variables
factors <- c("malesex", "ethnicity_qrisk2", "qrisk2_smoking_cat", "predrug_hypertension", 
             "predrug_dka", "osteoporosis", "predrug_acutepancreatitis", "predrug_falls", 
             "predrug_urinary_frequency", "predrug_volume_depletion", 
             "predrug_micturition_control", "predrug_dementia", "hosp_admission_prev_year",
             "ncurrtx", "MFN", "TZD", "ACEi", "ARB")

table <- CreateTableOne(vars = vars, strata = "studydrug", data = temp[temp$.imp > 0,], factorVars = factors, test = F)
print(table)



# save imputed dataset so this can be used in the subsequent scripts
today <- as.character(Sys.Date(), format="%Y%m%d")
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
save(temp, file=paste0(today, "_t2d_ckdpc_imputed_data.Rda"))