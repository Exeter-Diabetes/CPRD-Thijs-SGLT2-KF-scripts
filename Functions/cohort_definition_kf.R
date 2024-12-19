
# Inclusion/exclusion criteria:
## a) Subjects with T2D
## b) With HES linkage
## c) 1st instance
## d) Exclude if start drug within 91 days of registration

## e) Aged 18+
## f) SGLT2/DPP4/SU
## g) Initiated between 01/01/2013 and end of data (31/10/2020)
## h) No CVD (broad definition: angina, IHD, MI, PAD, revasc, stroke, TIA (as per NICE but with TIA))
## i) No HF before index date
## j) No missing eGFR/uACR
## k) No advanced CKD (egfr < 20 mL/min/1.73m2 or CKD stage 5) before index date, or eGFR <60 or albuminuria ≥30mg/mmol
## l) Exclude if also on GLP1 agonist
## m) Remove further episodes of starting DPP4/SU if already taking SGLT2i in previous episode (these episodes would overlap)


# Use "t2d_1stinstance" cohort_dataset which already has a)-d) applied
# all_drug_periods_dataset is used to define later start dates of meds for censoring


define_cohort <- function(cohort_dataset, all_drug_periods_dataset) {
  
  
  # e)-g) Keep those aged >=18 and within study period
  cohort <- cohort_dataset %>%
    filter(dstartdate_age>=18 &
             (drugclass=="SGLT2" | drugclass=="DPP4" | drugclass=="SU") &
             dstartdate>=as.Date("2013-01-01")
    ) %>%
    mutate(studydrug=ifelse(drugclass=="SGLT2", "SGLT2", ifelse(drugclass=="DPP4", "DPP4", "SU")))
  
  # Remove if on SGLT2 (except SGLT2 arm)
  cohort <- cohort %>%
    filter(!(drugclass!="SGLT2" & SGLT2==1))
  
  
  
  # Remove if on SGLT2 (except SGLT2 arm)
  cohort <- cohort %>%
    filter((drugclass=="SGLT2" | SGLT2==0))
  
  
  
  # create variables that denote whether this is first/only/last episode
  cohort <- cohort %>%
    group_by(patid) %>%
    mutate(earliest_start = min(dstartdate, na.rm=TRUE),
           latest_start   = max(dstartdate, na.rm=TRUE),
           episode_order  = ifelse(earliest_start == latest_start, "only", 
                                   ifelse(earliest_start == dstartdate, "first", 
                                          ifelse(latest_start == dstartdate, "last", "other")))) 
  
  
  q <- cohort %>% .$patid %>% unique() %>% length()
  print(paste0("Number of subjects of starting an SGLT2/DPP4/SU between 2013-2020: ", q))
  
  q <- cohort %>% nrow()
  print(paste0("Number of drug episodes of starting an SGLT2/DPP4/SU between 2013-2020: ", q))
  

  # h) Remove if CVD before index date
  
  cohort <- cohort %>%
    mutate(predrug_cvd=ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1 | predrug_tia==1, 1, 0)) 
  
  q <- cohort %>% filter(predrug_cvd == 1) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with established CVD: ", q))
  
  cohort <- cohort %>%
    filter(predrug_cvd==0)
  
  
  # i) Remove if HF before index date
  q <- cohort %>% filter(predrug_heartfailure == 1) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with established HF: ", q))
  
  cohort <- cohort %>%
    filter(predrug_heartfailure==0)
  
  # j) Remove if missing CKD status
  
   # create acr variable for ckdpc risk scores that uses further source of acr if acr not available
  cohort <- cohort %>% 
    mutate(uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
           uacr=ifelse(uacr<0.6, 0.6, uacr))
             
  q <- cohort %>% filter(is.na(preckdstage) | is.na(preegfr) | is.na(uacr)) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with unknown CKD status: ", q))
  
  cohort <- cohort %>%
    filter(
      !(is.na(preckdstage) | is.na(preegfr) | is.na(uacr))
    )

  # k) Remove if CKD before index date
  
  q <- cohort %>% filter(preckdstage=="stage_5" | predrug_ckd5_code == 1 | preegfr < 20) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with established eGFR <20 mL/min/1.73m2 or ESKD: ", q))
  
  cohort <- cohort %>%
    filter(!(preegfr < 20 | preckdstage=="stage_5" | predrug_ckd5_code == 1) )

  q <- cohort %>% filter(preckdstage=="stage_3a" | preckdstage=="stage_3b" | preckdstage=="stage_4" | preegfr < 60) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with established eGFR <60 mL/min/1.73m2: ", q))
  
  cohort <- cohort %>%
    filter(!(preckdstage=="stage_3a" | preckdstage=="stage_3b" | preckdstage=="stage_4" | preegfr < 60) )
  
  q <- cohort %>% filter(uacr >= 30) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with uACR ≥30mg/mmol: ", q))
  
  cohort <- cohort %>%
    filter(uacr < 30)
  
  # l) Remove if on GLP1 agonist at start 
  q <- cohort %>% filter(GLP1==1) %>% nrow()
  
  print(paste0("Number of drug episodes with concurrent treatment with GLP1 receptor agonist: ", q))
  
  cohort <- cohort %>%
    filter(GLP1==0)
  
  
  # m) Remove further episodes of starting DPP4/SU if already taking SGLT2i in previous episode (these episodes would overlap)
  #    or episodes of taking DPP4 following episode of SU and vice versa
  q <- cohort %>% filter(episode_order %in% c("last", "other") & (
    ((drugclass == "DPP4" | drugclass == "SU") & SGLT2 == 1) |
      (drugclass == "DPP4" & SU == 1) |
      (drugclass == "SU" & DPP4 == 1))
  ) %>% nrow()
  
  print(paste0("Number of drug episodes removed (e.g. subsequent episode of starting DPP4/SU after episode of SGLT2): ", q))
  
  cohort <- cohort %>%
    filter(!(episode_order %in% c("last", "other") & (
      ((drugclass == "DPP4" | drugclass == "SU") & SGLT2 == 1) |
        (drugclass == "DPP4" & SU == 1) |
        (drugclass == "SU" & DPP4 == 1))))
  
  
  
  q <- cohort %>% .$patid %>% unique() %>% length()
  print(paste0("Number of subjects included: ", q))
  
  q <- cohort %>% nrow()
  print(paste0("Number of drug episodes included: ", q))
  
  rm(q)
  
  ## Use all SGLT2, GLP1 starts to code up later censoring
  
  #
  ### Also get latest GLP1 and SGLT2 stop dates before drug start for DPP4/SU arms 
  ### we need this for sensitivity analysis where we exclude people who tried SGLT2 in the year before starting a DPP4/SU.
  
  later_sglt2 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="SGLT2") %>%
                  select(patid, next_sglt2=dstartdate)), by="patid") %>%
    filter(next_sglt2>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE)) %>%
    ungroup()
  
  later_dpp4 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="DPP4") %>%
                  select(patid, next_dpp4=dstartdate)), by="patid") %>%
    filter(next_dpp4>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_dpp4_start=min(next_dpp4, na.rm=TRUE)) %>%
    ungroup()
  
  later_su <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="SU") %>%
                  select(patid, next_su=dstartdate)), by="patid") %>%
    filter(next_su>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_su_start=min(next_su, na.rm=TRUE)) %>%
    ungroup()
  
  later_glp1 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="GLP1") %>%
                  select(patid, next_glp1=dstartdate)), by="patid") %>%
    filter(next_glp1>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_glp1_start=min(next_glp1, na.rm=TRUE)) %>%
    ungroup()
  
  
  later_tzd <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="TZD") %>%
                  select(patid, next_tzd=dstartdate)), by="patid") %>%
    filter(next_tzd>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_tzd_start=min(next_tzd, na.rm=TRUE)) %>%
    ungroup()
  
  
  last_sglt2_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="SGLT2") %>%
                  select(patid, last_sglt2=dstopdate)), by="patid") %>%
    filter(last_sglt2<dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(last_sglt2_stop=min(last_sglt2, na.rm=TRUE)) %>%
    ungroup()
  
  last_dpp4_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="DPP4") %>%
                  select(patid, last_dpp4=dstopdate)), by="patid") %>%
    filter(last_dpp4<dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(last_dpp4_stop=min(last_dpp4, na.rm=TRUE)) %>%
    ungroup()
  
  last_su_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="SU") %>%
                  select(patid, last_su=dstopdate)), by="patid") %>%
    filter(last_su<dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(last_su_stop=min(last_su, na.rm=TRUE)) %>%
    ungroup()
  
  last_glp1_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drugclass=="GLP1") %>%
                  select(patid, last_glp1=dstopdate)), by="patid") %>%
    filter(last_glp1<dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(last_glp1_stop=min(last_glp1, na.rm=TRUE)) %>%
    ungroup()
  
  
  cohort <- cohort %>%
    left_join(later_sglt2, by=c("patid", "dstartdate")) %>%
    left_join(later_dpp4, by=c("patid", "dstartdate")) %>%
    left_join(later_su, by=c("patid", "dstartdate")) %>%
    left_join(later_glp1, by=c("patid", "dstartdate")) %>%
    left_join(later_tzd, by=c("patid", "dstartdate")) %>%
    left_join(last_sglt2_stop, by=c("patid", "dstartdate")) %>%
    left_join(last_dpp4_stop, by=c("patid", "dstartdate")) %>%
    left_join(last_su_stop, by=c("patid", "dstartdate")) %>%
    left_join(last_glp1_stop, by=c("patid", "dstartdate"))
  
  
  
  # Tidy up gender, ncurrtx, drugline and ethnicity variables
  ## Also code up death cause variables
  
  cohort <- cohort %>%
    
    mutate(malesex=ifelse(gender==1, 1, 0),
           
           ncurrtx=DPP4+GLP1+MFN+SU+SGLT2+TZD+INS,          #INS, GLP1 and TZD should be 0 but include anyway; ignore Acarbose and Glinide
           
           drugline_all=as.factor(ifelse(drugline_all>=5, 5, drugline_all)),
           
           drugsubstances=ifelse(grepl("&", drugsubstances), NA, drugsubstances),
           
           ethnicity_5cat_decoded=case_when(ethnicity_5cat==0 ~"White",
                                            ethnicity_5cat==1 ~"South Asian",
                                            ethnicity_5cat==2 ~"Black",
                                            ethnicity_5cat==3 ~"Other",
                                            ethnicity_5cat==4 ~"Mixed"),
           
           cv_death_date_any_cause=if_else(!is.na(death_date) & !is.na(cv_death_any_cause) & cv_death_any_cause==1, death_date, as.Date(NA)),
           cv_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(cv_death_primary_cause) & cv_death_primary_cause==1, death_date, as.Date(NA)),
           kf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(kf_death_any_cause) & kf_death_any_cause==1, death_date, as.Date(NA)),
           kf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(kf_death_primary_cause) & kf_death_primary_cause==1, death_date, as.Date(NA)),
           hf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(hf_death_any_cause) & hf_death_any_cause==1, death_date, as.Date(NA)),
           hf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(hf_death_primary_cause) & hf_death_primary_cause==1, death_date, as.Date(NA)))
  
  return(cohort)
  
}