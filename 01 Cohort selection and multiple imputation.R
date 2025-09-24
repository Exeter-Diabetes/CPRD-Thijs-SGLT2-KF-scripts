## In this script, we will select a cohort from CPRD and make it suitable for analyses.
## Contents:
# 0 setup
# 1 cohort selection
# 2 multiple imputation of missing values
# 3 calculate risk scores (for CKD and QRISK2)
# 4 select subjects with valid CKD risk scores (ie. input values within range for calculation of risk scores)
# 5 baseline table and store dataset for further analyses

########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/")
source("00 Setup.R")
########################1 COHORT SELECTION####################################################################

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition_kf function for details)

setwd("/slade/CPRD_data/Thijs/SGLT2/Raw data/")
load(paste0(today, "_t2d_1stinstance_a.Rda"))
load(paste0(today, "_t2d_1stinstance_b.Rda"))
load(paste0(today, "_t2d_all_drug_periods.Rda"))

t2d_1stinstance <- rbind(t2d_1stinstance_a, t2d_1stinstance_b)
rm(t2d_1stinstance_a)
rm(t2d_1stinstance_b)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/Functions/")
source("cohort_definition_kf.R")
cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

table(cohort$studydrug)

## for some reason the dataset contains dstopdate.x and dstopdate.x which are identical
# remove these if this is present
if ("dstopdate.x" %in% names(cohort)) {
  cohort <- cohort %>% mutate(dstopdate = dstopdate.x) %>% select(
    -c(dstopdate.x, dstopdate.y)
  )
}


## B Make variables for survival analysis of all endpoints (see survival_variables_kf function for details)

source("survival_variables_kf.R")

cohort <- add_surv_vars(cohort, main_only=FALSE) # add per-protocol survival variables as well

rm(list=setdiff(ls(), c("cohort", "today", "vars", "factors", "nonnormal")))

## C Just keep variables of interest

# create acr variable for ckdpc risk scores that uses further source of acr if acr not available
cohort <- cohort %>% 
  mutate(uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
         uacr=ifelse(uacr<0.6, 0.6, uacr),
         #and create variable to code whether someone is on oral hyperglycaemic agents
         oha=ifelse(Acarbose+MFN+DPP4+Glinide+GLP1+SGLT2+SU+TZD>add, 1L, 0L),
         statin=!is.na(predrug_latest_statins),
         ACEi=!is.na(predrug_latest_ace_inhibitors),
         ARB=!is.na(predrug_latest_arb),
         BB=!is.na(predrug_latest_beta_blockers),
         CCB=!is.na(predrug_latest_calcium_channel_blockers),
         ThZD=!is.na(predrug_latest_thiazide_diuretics),
         loopD=!is.na(predrug_latest_loop_diuretics),
         MRA=!is.na(predrug_latest_ksparing_diuretics),
         steroids=!is.na(predrug_latest_oralsteroids),
         immunosuppr=!is.na(predrug_latest_immunosuppressants),
         osteoporosis=!is.na(predrug_latest_osteoporosis),
         genital_infection=as.logical(predrug_medspecific_gi)
         )


cohort <- cohort %>%
  #default hba1c variable is from previous 6 months to index date
  #take hba1c within window of 2 years prior and 7 days post, similar as other biomarkers
  mutate(prehba1c = prehba1c2yrs) %>%
  select(patid, malesex, ethnicity_5cat, ethnicity_qrisk2, imd2015_10, regstartdate, gp_record_end, death_date, 
         drugclass, studydrug, dstartdate, dstopdate, drugline_all, drugsubstances, ncurrtx, DPP4, GLP1, 
         MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, height, prehba1c, prebmi, 
         prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, predbp, preegfr, preckdstage, 
         preacr, uacr, qrisk2_smoking_cat, contains("cens"), last_sglt2_stop, oha,
         #add variables necessary to calculate qrisk2/qhdf and ckdpc scores
         predrug_fh_premature_cvd, predrug_af, predrug_rheumatoidarthritis, tds_2011,
         predrug_angina, predrug_myocardialinfarction, predrug_stroke, predrug_revasc,
         predrug_heartfailure, predrug_hypertension, predrug_acutepancreatitis,
         predrug_earliest_ace_inhibitors, predrug_earliest_arb,
         predrug_earliest_beta_blockers, predrug_earliest_calcium_channel_blockers,
         predrug_earliest_thiazide_diuretics, predrug_latest_ace_inhibitors, 
         predrug_latest_arb, predrug_latest_beta_blockers, 
         predrug_latest_calcium_channel_blockers, predrug_latest_thiazide_diuretics,
         predrug_dka, predrug_falls, predrug_urinary_frequency, predrug_volume_depletion, 
         predrug_micturition_control, predrug_dementia, hosp_admission_prev_year,
         statin, ACEi, ARB, BB, CCB, ThZD, loopD, MRA, steroids, immunosuppr, 
         osteoporosis, genital_infection,
         ckd_egfr50_outcome_type, preacr_confirmed, preacr_previous, preacr_previous_date, preacr_next, preacr_next_date
  )

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

# imd_2015 is not a continuous variable - we will categorise this as quantiles
cohort <- cohort %>% mutate(
  imd2015_10 = ifelse(imd2015_10 %in% c(1,2), "1/2",
                      ifelse(imd2015_10 %in% c(3,4), "3/4",
                             ifelse(imd2015_10 %in% c(5,6), "5/6",
                                    ifelse(imd2015_10 %in% c(7,8), "7/8",
                                           ifelse(imd2015_10 %in% c(9,10), "9/10", NA))))),
  imd2015_10 = factor(imd2015_10)
)

# variable preacr_confirmed indicates whether a person had their presence of albuminuria (3mg/mmol) confirmed on 2 readings
# this shows as NA if no second reading available to confirm - replace with NA
cohort <- cohort %>% mutate(preacr_confirmed = ifelse(is.na(preacr_confirmed), F, preacr_confirmed),
                            preacr_confirmed = ifelse(uacr<3, F, preacr_confirmed))

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
       "predrug_latest_thiazide_diuretics",
       "ckd_egfr50_outcome_type", "preacr_confirmed", 
       "preacr_previous", "preacr_previous_date", "preacr_next", "preacr_next_date")] <- ""

meth[c("qrisk2_smoking_cat", "imd2015_10")] <- "polyreg"

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
              names(cohort)[grep("cens", names(cohort))], 
              "oha", "predrug_angina", "predrug_myocardialinfarction", "predrug_stroke", 
              "predrug_revasc", "predrug_heartfailure", "initiation_year", "ethnicity_5cat",
              "ckd_egfr50_outcome_type", "preacr_confirmed", "preacr_previous", "preacr_previous_date", "preacr_next", "preacr_next_date")

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
                          if_else(is.na(predrug_earliest_ace_inhibitors),as.Date("2050-01-01"),predrug_earliest_ace_inhibitors),
                          if_else(is.na(predrug_earliest_arb),as.Date("2050-01-01"),predrug_earliest_arb),
                          if_else(is.na(predrug_earliest_beta_blockers),as.Date("2050-01-01"),predrug_earliest_beta_blockers),
                          if_else(is.na(predrug_earliest_calcium_channel_blockers),as.Date("2050-01-01"),predrug_earliest_calcium_channel_blockers),
                          if_else(is.na(predrug_earliest_thiazide_diuretics),as.Date("2050-01-01"),predrug_earliest_thiazide_diuretics),
                          na.rm=TRUE
                        ) %>% as.Date(),
                        latest_bp_med=pmax(
                          if_else(is.na(predrug_latest_ace_inhibitors),as.Date("1900-01-01"),predrug_latest_ace_inhibitors),
                          if_else(is.na(predrug_latest_arb),as.Date("1900-01-01"),predrug_latest_arb),
                          if_else(is.na(predrug_latest_beta_blockers),as.Date("1900-01-01"),predrug_latest_beta_blockers),
                          if_else(is.na(predrug_latest_calcium_channel_blockers),as.Date("1900-01-01"),predrug_latest_calcium_channel_blockers),
                          if_else(is.na(predrug_latest_thiazide_diuretics),as.Date("1900-01-01"),predrug_latest_thiazide_diuretics),
                          na.rm=TRUE
                        ) %>% as.Date(),
                        bp_meds_qrisk2=if_else(earliest_bp_med!=as.Date("2050-01-01") & latest_bp_med!=as.Date("1900-01-01") & difftime(dstartdate, latest_bp_med, units="days")<=28 & earliest_bp_med!=latest_bp_med, 1L, 0L),
                        
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

#ckdpc kidney disease progression risk score

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/Functions/")
source("calculate_ckdpc_50egfr_risk.R")

temp <- temp %>% 
  
  mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
  calculate_ckdpc_50egfr_risk(age=dstartdate_age, sex=sex2, egfr=preegfr, acr=uacr, sbp=presbp, bp_meds=bp_meds_ckdpc, hf=predrug_heartfailure, chd=chd, af=predrug_af, current_smoker=current_smoker, ex_smoker=ex_smoker, bmi=prebmi, hba1c=prehba1c, oha=oha, insulin=INS) 


########################4 REMOVE RISK SCORE VALUES OUTSIDE OF RANGE####################################################################

# 4 remove risk scores for subjects whose characteristics are outside of the reference range for the respective equations

# Look at counts of people whose characteristics are outside the reference range for the scores

ckdpc_outofrange <- temp[temp$.imp > 0,] %>%
  mutate(uacr_val=ifelse(uacr < 113, "in range for 50eGFR score",
                                "above range for both scores"),
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


## Also remove 40% decline in eGFR score for those with sociodemographic/vital/laboratory measurements outside of range:
### Age: 20-80
### UACR: 0.6-113 (5-1000 in mg/g)
### SBP: 80-180
### BMI: <20
### HbA1c 42-97 (6-11 in %)

temp <- temp %>%
  
  mutate(across(starts_with("ckdpc_50egfr"),
                ~ifelse(dstartdate_age>=20 & dstartdate_age<=80 &
                          prebmi>=20, .x, NA)))

q <- temp %>% filter(.imp !=0 & (dstartdate_age<20 | dstartdate_age>80 | prebmi < 20))

# those left out:
q1 <- q %>% .$patid %>% unique() %>% length()
print(paste0("Number of subjects excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q1))

q2 <- q %>% nrow()
print(paste0("Number of drug episodes excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q2/n.imp))


# retain those with available risk scores only
temp <- temp %>% filter(!is.na(ckdpc_50egfr_score))

# study cohort:

# add variables
temp <- temp %>% mutate(
  obesity = ifelse(prebmi < 30, F, T),
  smoking_hx = ifelse(qrisk2_smoking_cat == 0, F, T),
  smoking_status = ifelse(qrisk2_smoking_cat == 0, "never", ifelse(qrisk2_smoking_cat == 1, "ex", "current")),
  albuminuria_unconfirmed = ifelse(uacr < 3, F, T),
  albuminuria = preacr_confirmed,        # 
  dyslipidaemia = ifelse(pretotalcholesterol < 5 &
                           preldl < 4 &
                           pretriglyceride < 2.3, F, T),
  cv_high_risk = ifelse(dstartdate_age >= 55 & 
                          obesity + predrug_hypertension + 
                          smoking_hx + dyslipidaemia +
                          albuminuria > 1, T, F),
  ACEi_or_ARB = ifelse(temp$ACEi + temp$ARB > 0, T, F),
  studydrug = ifelse(studydrug == "SGLT2", "SGLT2i", ifelse(studydrug == "DPP4", "DPP4i", "SU")),
  studydrug2 = ifelse(!studydrug == "SGLT2i", "DPP4i/SU", "SGLT2i"),
  ncurrtx = ifelse(ncurrtx==1, "1.", ifelse(ncurrtx==2, "2.", "3+"))
)


q <- temp %>% filter(!.imp == 0) %>% nrow()
p <- temp %>% filter(!.imp == 0) %>%  ## dataset at present contains separate drug episodes if a subject started a DPP4i and later a sulfonylurea
  group_by(.imp, patid) %>% filter(!duplicated(studydrug2)) %>% ungroup() %>% nrow()
print(paste0("Number of duplicate drug episodes removed ", (q-p)/n.imp))
print(paste0("Number of drug episodes in study population ", p/n.imp))
rm(p)
q <- temp %>% .$patid %>% unique() %>% length()
print(paste0("Number of subjects in study population ", q))


########################5 SAVE DATASET####################################################################

# 5 tabulate and save imputed dataset for further analyses

temp <- temp %>% group_by(.imp, patid) %>% filter(
   !duplicated(studydrug)
   ) %>% ungroup()

# save imputed dataset so this can be used in the subsequent scripts
setwd("/slade/CPRD_data/Thijs/SGLT2/Processed data/")
save(temp, file=paste0(today, "_t2d_ckdpc_imputed_data.Rda"))


# create table one: this will be an average of the imputed datasets (n to be divided by n.imp)

table <- CreateTableOne(vars = vars, strata = "studydrug2", data = temp %>% filter(!.imp == 0) %>%  ## dataset at present contains separate drug episodes if a subject started a DPP4i and later a sulfonylurea
                          group_by(.imp, patid) %>% filter(!duplicated(studydrug2)) %>% ungroup(),  ## these "duplicate" episodes will be removed after we have done the drug-specific analyses
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)
## Save to a CSV file
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv

write.csv2(tabforprint, file = paste0(today, "_baseline_table.csv"))

# baseline table with DPP4i/SU split
table <- CreateTableOne(vars = vars, strata = "studydrug", data = temp %>% filter(!.imp == 0) %>%  ## dataset at present contains separate drug episodes if a subject started a DPP4i and later a sulfonylurea
                          group_by(.imp, patid) %>% filter(!duplicated(studydrug)) %>% ungroup(),  ## these "duplicate" episodes will be removed after we have done the drug-specific analyses
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)

write.csv2(tabforprint, file = paste0(today, "_baseline_table_dpp4isu_split.csv"))


# get baseline table by albuminuria status:
vars <- c(vars, "studydrug2")
factors <- c(factors, "studydrug2")
table2 <- CreateTableOne(vars = vars, strata = "albuminuria", data = temp %>% filter(!.imp == 0) %>%
                          group_by(.imp, patid) %>% filter(!duplicated(studydrug2)) %>% ungroup(),
                        factorVars = factors, test = F)

tabforprint2 <- print(table2, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)

write.csv2(tabforprint2, file = paste0(today, "_baseline_table_by_albuminuria.csv"))


# events rates (sum of events divided by sum of person-years) by studydrug
outcomes <- c("ckd_egfr40", "ckd_egfr50", "ckd_egfr50_5y", "macroalb", "dka", "side_effect", "death", "amputation")

for (k in outcomes) {
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")  
  
  for (m in levels(as.factor(temp$studydrug2))) {
    events <- temp %>% filter(.imp !=0 & studydrug2 == m) %>% select(censvar_var) %>% sum()
    pyears <- temp %>% filter(.imp !=0 & studydrug2 == m) %>% select(censtime_var) %>% sum()
    print(paste0(m, " event rate for ", k, ": ", round(events/pyears*1000,1), " per 1000 patient-years"))
    rm(events)
    rm(pyears)
  }
  rm(censvar_var)
  rm(censtime_var)
}


## get weighted baseline table
#write functions to summarise weighted table:
cont <- function(x, var_name) {
  if (var_name %in% c("uacr", "dstartdate_dm_dur_all")) {
    # For specific variables, calculate Median (IQR)
    with(stats.apply.rounding(stats.default(x)), 
         c("Median (IQR)" = sprintf("%s (%s-%s)", 
                                    round_pad(as.numeric(MEDIAN), 1), 
                                    round_pad(as.numeric(Q1), 1), 
                                    round_pad(as.numeric(Q3), 1))))
  } else {
    # For other continuous variables, calculate Mean (SD)
    with(stats.apply.rounding(stats.default(x)), 
         c("Mean (SD)" = sprintf("%s (%s)", 
                                 round_pad(as.numeric(MEAN), 1), 
                                 round_pad(as.numeric(SD), 1))))
  }
}
missing <- function(x, ...) {
  with(stats.apply.rounding(stats.default(x)), c("Missing"=sprintf("%s", prettyNum(NMISS, big.mark=","))))
}


rndr <- function(x, name, ...) {
  if (is.logical(x)) {
    y <- render.default(x, name, ...)
    y[2]
  } else if (is.numeric(x)) {
    cont(x, name)  # pass both x and variable name to your cont() function
  } else {
    render.default(x, name, ...)
  }
}


strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s</span>", 
          label, prettyNum(n, big.mark=","))
}

cat <- function(x, ...) {
  vals <- stats.default(x)  # get raw stats without rounding
  c("", sapply(vals, function(y) {
    # assume y$PCT is numeric; convert if necessary
    sprintf("%.4f%%", as.numeric(y$PCT))
  }))
}

gc()

setwd("/slade/CPRD_data/Thijs/SGLT2/Processed data/")

library(readr)

# Define your full list of covariates
covariates_list <- c("dstartdate_age", "malesex", "ethnicity_5cat", "imd2015_10", 
                     "prebmi", "presbp", "predbp", "pretotalcholesterol", "prehdl", 
                     "preldl", "pretriglyceride", "prehba1c", "preegfr", "uacr", 
                     "albuminuria", "dstartdate_dm_dur_all", 
                     "smoking_status", "predrug_hypertension", "predrug_af", "predrug_dka", 
                     "genital_infection", "hosp_admission_prev_year", "initiation_year", 
                     "ncurrtx", "statin", "INS", "ACEi_or_ARB")

# Create chunks of 5 covariates
covariates_chunks <- split(covariates_list, ceiling(seq_along(covariates_list) / 4))

# Define full set of known columns in correct order
all_colnames <- c("studydrug2", covariates_list)

for (covariates_chunk in covariates_chunks) {
  
  columns_to_load <- c("studydrug2", covariates_chunk)
  
  load(paste0(today, "_weighted_imputed_data_for_table.Rda"))
  
  # Read just the needed columns
  data_chunk <- weighted_cohort %>% select(all_of(columns_to_load))
  
  rm(weighted_cohort)
  gc()
  
  # Build formula for table1
  formula <- as.formula(paste("~ ", paste(covariates_chunk, collapse = " + "), "| studydrug2"))
  
  # Run the table1 function on the selected chunk
  table1(formula, data=data_chunk, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat)
  
  # Optionally, save the output if needed
  
  # Clear the data from memory and run garbage collection to free memory
  rm(data_chunk)
  gc()
  
}
