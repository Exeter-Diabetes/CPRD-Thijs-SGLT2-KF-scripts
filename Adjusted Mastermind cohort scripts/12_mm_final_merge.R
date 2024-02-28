
# Extract dataset of all first instance drug periods (i.e. the first time patient has taken this particular drug class) for T2Ds WITH HES LINKAGE
## Exclude drug periods starting within 91 days of registration
## Set drugline to missing where diagnosed before registration

## Do not exclude where first line
## Do not exclude where patient is on insulin at drug initiation
## Do not exclude where only 1 prescription (dstartdate=dstopdate)

## Set hosp_admission_prev_year to 0/1 rather than NA/1

# Also extract all T2D drug start and stop dates for T2Ds so that you can see if people later initiate SGLT2is/GLP1s etc.

############################################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")


############################################################################################

# Today's date for table names

today <- as.character(Sys.Date(), format="%Y%m%d")


############################################################################################

# Get handles to pre-existing data tables

## Cohort and patient characteristics including Townsend scores
analysis = cprd$analysis("all")
diabetes_cohort <- diabetes_cohort %>% analysis$cached("diabetes_cohort")
townsend_score <- townsend_score %>% analysis$cached("patid_townsend_score")

## Drug info
analysis = cprd$analysis("mm")
drug_start_stop <- drug_start_stop %>% analysis$cached("drug_start_stop")
combo_start_stop <- combo_start_stop %>% analysis$cached("combo_start_stop")

## Biomarkers inc. CKD
#baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("baseline_biomarkers")
response_biomarkers <- response_biomarkers %>% analysis$cached("response_biomarkers") #includes baseline biomarker values for first instance drug periods so no need to use baseline_biomakers table
ckd_stages <- ckd_stages %>% analysis$cached("ckd_stages")

## Comorbidities
analysis = cprd$analysis("Thijs_ckd") # nieuwe versie daarom onder mijn naam opgeslagen ipv onder mm
comorbidities <- comorbidities %>% analysis$cached("comorbidities")

## Non-diabetes meds
analysis = cprd$analysis("Thijs_ckd") # nieuwe versie daarom onder mijn naam opgeslagen ipv onder mm
non_diabetes_meds <- non_diabetes_meds %>% analysis$cached("nondiabetes_meds")

## Smoking status at drug start
analysis = cprd$analysis("mm")
smoking <- smoking %>% analysis$cached("smoking")

## Discontinuation
discontinuation <- discontinuation %>% analysis$cached("discontinuation")

## Glycaemic failure
glycaemic_failure <- glycaemic_failure %>% analysis$cached("glycaemic_failure")

## Death causes
analysis = cprd$analysis("Thijs_test") #omdat ik kf_death_cause ergens anders moest toevoegen
death_causes <- death_causes %>% analysis$cached("death_causes")

############################################################################################

# Make first instance drug period dataset
analysis = cprd$analysis("Thijs_ckd")

## Define T2D cohort (1 line per patient) with HES linkage
## Add in Townsend Deprivation Scores
t2ds <- diabetes_cohort %>%
  left_join((townsend_score %>% select(patid, tds_2011)), by="patid") %>%
  relocate(tds_2011, .after=imd2015_10) %>%
  filter(diabetes_type=="type 2" & with_hes==1)


## Get info for first instance drug periods for cohort (1 line per patid-drugclass period)
### Make new drugline variable which is missing where diagnosed before registration or within 90 days following

t2d_drug_periods <- t2ds %>%
  inner_join(drug_start_stop, by="patid") %>%
  inner_join(combo_start_stop, by=c("patid", c("dstartdate"="dcstartdate"))) %>%
  inner_join((glycaemic_failure %>% select(-c(dstopdate, timetochange, timetoaddrem, nextdrugchange, nextdcdate, prehba1c, prehba1cdate, threshold_7.5, threshold_8.5, threshold_baseline, threshold_baseline_0.5))), by=c("patid", "dstartdate", "drugclass")) %>%
  mutate(drugline=ifelse(dm_diag_date_all<regstartdate | is.na(dm_diag_date), NA, drugline_all)) %>%
  relocate(drugline, .after=drugline_all) %>%
  analysis$cached(paste0(today, "_t2d_1stinstance_interim_1"), indexes=c("patid", "dstartdate", "drugclass"))

t2d_drug_periods %>% distinct(patid) %>% count()
# 865,054


### Keep first instance only
t2d_1stinstance <- t2d_drug_periods %>%
  filter(druginstance==1)

t2d_1stinstance %>% distinct(patid) %>% count()
# 865,054 as above


### Exclude drug periods starting within 90 days of registration
t2d_1stinstance <- t2d_1stinstance %>%
  filter(datediff(dstartdate, regstartdate)>90)

t2d_1stinstance %>% count()
# 1,663,398

t2d_1stinstance %>% distinct(patid) %>% count()
# 769,394



## Merge in biomarkers, comorbidities, non-diabetes meds, smoking status
### Could merge on druginstance too, but quicker not to
### Remove some variables to avoid duplicates
### Make new variables: age at drug start, diabetes duration at drug start, CV risk scores
### Now in two stages to speed it up

t2d_1stinstance <- t2d_1stinstance %>%
  inner_join((response_biomarkers %>% select(-c(druginstance, timetochange, timetoaddrem, multi_drug_start, timeprevcombo))), by=c("patid", "dstartdate", "drugclass")) %>%
  inner_join((ckd_stages %>% select(-druginstance)), by=c("patid", "dstartdate", "drugclass")) %>%
  inner_join((comorbidities %>% select(-druginstance)), by=c("patid", "dstartdate", "drugclass")) %>%
  analysis$cached(paste0(today, "_t2d_1stinstance_interim_2"), indexes=c("patid", "dstartdate", "drugclass"))

t2d_1stinstance <- t2d_1stinstance %>%
  inner_join((non_diabetes_meds %>% select(-druginstance)), by=c("patid", "dstartdate", "drugclass")) %>%
  inner_join((smoking %>% select(-druginstance)), by=c("patid", "dstartdate", "drugclass")) %>%
  inner_join((discontinuation %>% select(-c(druginstance, timeondrug, nextremdrug, timetolastpx))), by=c("patid", "dstartdate", "drugclass")) %>%
  left_join(death_causes, by="patid") %>%
  mutate(dstartdate_age=datediff(dstartdate, dob)/365.25,
         dstartdate_dm_dur_all=datediff(dstartdate, dm_diag_date_all)/365.25,
         dstartdate_dm_dur=datediff(dstartdate, dm_diag_date)/365.25,
         hosp_admission_prev_year=ifelse(is.na(hosp_admission_prev_year) & with_hes==1, 0L,
                                         ifelse(hosp_admission_prev_year==1, 1L, NA))) %>%
  analysis$cached(paste0(today, "_t2d_1stinstance_interim_3"), indexes=c("patid", "dstartdate", "drugclass"))


# Check counts

t2d_1stinstance %>% count()
# 1,663,398

t2d_1stinstance %>% distinct(patid) %>% count()
# 769,394

############################################################################################

# Export to R data object
## Convert integer64 datatypes to double

t2d_1stinstance_a <- collect(t2d_1stinstance %>% filter(patid<2000000000000) %>% mutate(patid=as.character(patid)))

is.integer64 <- function(x){
  class(x)=="integer64"
}

t2d_1stinstance_a <- t2d_1stinstance_a %>%
  mutate_if(is.integer64, as.integer)

save(t2d_1stinstance_a, file=paste0(today, "_t2d_1stinstance_a.Rda"))

rm(t2d_1stinstance_a)


t2d_1stinstance_b <- collect(t2d_1stinstance %>% filter(patid>=2000000000000) %>% mutate(patid=as.character(patid)))

t2d_1stinstance_b <- t2d_1stinstance_b %>%
  mutate_if(is.integer64, as.integer)

save(t2d_1stinstance_b, file=paste0(today, "_t2d_1stinstance_b.Rda"))

rm(t2d_1stinstance_b)


############################################################################################

# Make dataset of all T2D drug starts so that can see whether people later initiate SGLT2i/GLP1 etc.

t2d_all_drug_periods <- t2ds %>%
  inner_join(drug_start_stop, by="patid") %>%
  select(patid, drugclass, dstartdate, dstopdate) %>%
  analysis$cached(paste0(today, "_t2d_all_drug_periods"))


## Export to R data object
### No integer64 datatypes

t2d_all_drug_periods <- collect(t2d_all_drug_periods %>% mutate(patid=as.character(patid)))

save(t2d_all_drug_periods, file=paste0(today, "_t2d_all_drug_periods.Rda"))
