
# Produce survival variables for all endpoints (including for sensitivity analysis)
## All censored at 5 years post drug start (3 years for 'ckd_egfr40') / end of GP records / death / starting a different diabetes med which affects CV risk (TZD/GLP1/SGLT2), and also drug stop date + 6 months for per-protocol analysis

# Main analysis:
## 'mace': stroke, MI, CV death
## 'expanded_mace': stroke, MI, CV death, revasc, HES unstable angina
## 'hf'
## 'ckd_345': onset of CKD stage 3a-5
## 'ckd_egfr40': fall in eGFR of <=40% from baseline or onset of CKD stage 5 OR death from renal causes
## 'hosp': all-cause hospitalisation
## 'death': all-cause mortality

# Sensitivity analysis:
## 'narrow_mace': hospitalisation for incident MI (subset of HES codes), incident stroke (subset of HES codes, includes ischaemic only), CV death - all as primary cause for hospitalisation/death only
## 'narrow_hf': hospitalisation or death with HF as primary cause
## '{outcome}_pp': all of main analysis but per-protocol rather than intention to treat
## intention to treat: censoring if starting an SGLT2 inhibitor (if in DPP4 or SU arm) or GLP1 agonist.
## per-protocol: censoring if starting any other treatment arm or GLP1 agonist.
## not using per-protocol analyses other than to compare DPP4i with SU arm - therefore removing censoring 185 days after starting


add_surv_vars <- function(cohort_dataset, main_only=FALSE) {
  
  # Add survival variables for outcomes for main analysis
  main_outcomes <- c("ckd_egfr40", "death", "macroalb", "dka", "amputation", "side_effect", "any_ae")
  
  cohort <- cohort_dataset %>%
    
    mutate(cens_itt=pmin(dstartdate+(365.25*5),
                         gp_record_end,
                         death_date,
                         if_else(studydrug!="GLP1", next_glp1_start, as.Date("2050-01-01")),
                         if_else(studydrug!="SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                         na.rm=TRUE),
           
           cens_pp=pmin(dstartdate+(365.25*5),
                        gp_record_end,
                        death_date,
                        if_else(studydrug!="GLP1", next_glp1_start, as.Date("2050-01-01")),
                        if_else(studydrug!="SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                        if_else(studydrug!="SU", next_su_start, as.Date("2050-01-01")),
                        if_else(studydrug!="DPP4", next_dpp4_start, as.Date("2050-01-01")),
                        #   dstopdate+183,
                        na.rm=TRUE),
           
           cens_itt_3_yrs=pmin(dstartdate+(365.25*3),
                               gp_record_end,
                               death_date,
                               if_else(studydrug!="GLP1", next_glp1_start, as.Date("2050-01-01")),
                               if_else(studydrug!="SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                               na.rm=TRUE),
           
           cens_pp_3_yrs=pmin(dstartdate+(365.25*3),
                              gp_record_end,
                              death_date,
                              if_else(studydrug!="GLP1", next_glp1_start, as.Date("2050-01-01")),
                              if_else(studydrug!="SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                              if_else(studydrug!="SU", next_su_start, as.Date("2050-01-01")),
                              if_else(studydrug!="DPP4", next_dpp4_start, as.Date("2050-01-01")),
                              #    dstopdate+183,
                              na.rm=TRUE),
           
           
           #     mace_outcome=pmin(postdrug_first_myocardialinfarction,
           #                       postdrug_first_stroke,
           #                       cv_death_date_any_cause,
           #                       na.rm=TRUE),
           #    
           #     expanded_mace_outcome=pmin(postdrug_first_myocardialinfarction,
           #                                postdrug_first_stroke,
           #                                cv_death_date_any_cause,
           #                                postdrug_first_revasc,
           #                                postdrug_first_unstableangina,
           #                                na.rm=TRUE),
           #       
           #     hf_outcome=pmin(postdrug_first_heartfailure,
           #                     hf_death_date_any_cause,
           #                     na.rm=TRUE),
           #                      
           #     
           #     ckd_345_outcome=postckdstage345date,
           
           ckd_egfr40_outcome=pmin(egfr_40_decline_date,
                                   postckdstage5date,
                                   kf_death_date_any_cause,
                                   na.rm=TRUE),
           
           macroalb_outcome=macroalb_date,
           
           dka_outcome=postdrug_first_dka,
           
           amputation_outcome=postdrug_first_amputation,
           
           side_effect_outcome=pmin(postdrug_first_urinary_frequency,
                                    postdrug_first_micturition_control,
                                    postdrug_first_volume_depletion,
                                    postdrug_first_medspecific_gi,
                                    na.rm=TRUE),
           
           any_ae_outcome=pmin(macroalb_date,
                               postdrug_first_dka,
                               postdrug_first_amputation,
                               postdrug_first_urinary_frequency,
                               postdrug_first_micturition_control,
                               postdrug_first_volume_depletion,
                               postdrug_first_medspecific_gi,
                               na.rm=TRUE),
           
           #     hosp_outcome=postdrug_first_all_cause_hosp,
           #     
           death_outcome=death_date)
  
  
  for (i in main_outcomes) {
    
    outcome_var=paste0(i, "_outcome")
    censdate_var=paste0(i, "_censdate")
    censvar_var=paste0(i, "_censvar")
    censtime_var=paste0(i, "_censtime_yrs")
    censdate_var_5y=paste0(i, "_5y_censdate")
    censvar_var_5y=paste0(i, "_5y_censvar")
    censtime_var_5y=paste0(i, "_5y_censtime_yrs")
    
    if (i=="ckd_egfr40") {
      cohort <- cohort %>%
        mutate({{censdate_var}}:=pmin(!!sym(outcome_var), cens_itt_3_yrs, na.rm=TRUE),
               {{censvar_var}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var)==!!sym(outcome_var), 1, 0),
               {{censtime_var}}:=as.numeric(difftime(!!sym(censdate_var), dstartdate, unit="days"))/365.25,
               {{censdate_var_5y}}:=pmin(!!sym(outcome_var), cens_itt, na.rm=TRUE),
               {{censvar_var_5y}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var_5y)==!!sym(outcome_var), 1, 0),
               {{censtime_var_5y}}:=as.numeric(difftime(!!sym(censdate_var_5y), dstartdate, unit="days"))/365.25)
    } else {
      cohort <- cohort %>%
        mutate({{censdate_var}}:=pmin(!!sym(outcome_var), cens_itt, na.rm=TRUE),
               {{censvar_var}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var)==!!sym(outcome_var), 1, 0),
               {{censtime_var}}:=as.numeric(difftime(!!sym(censdate_var), dstartdate, unit="days"))/365.25)
    }
    
  }
  
  if (main_only==TRUE) {
    message(paste("survival variables for", paste(main_outcomes, collapse=", "), "added"))
  }
  
  
  # Add survival variables for outcomes for sensitivity analyses
  
  else {
    
    # Split by whether ITT or PP
    sensitivity_outcomes <- c("ckd_egfr40_pp", "death_pp")
    
    # cohort <- cohort %>%
    #    
    #  mutate(narrow_mace_outcome=pmin(postdrug_first_incident_mi,
    #                                    postdrug_first_incident_stroke,
    #                                   cv_death_date_primary_cause,
    #                                  na.rm=TRUE),
    #        
    #         narrow_hf_outcome=pmin(postdrug_first_primary_hhf,
    #                                hf_death_date_primary_cause,
    #                                na.rm=TRUE))
    
    
    for (i in sensitivity_outcomes) {
      
      censdate_var=paste0(i, "_censdate")
      censvar_var=paste0(i, "_censvar")
      censtime_var=paste0(i, "_censtime_yrs")
      
      
      outcome_var=paste0(substr(i, 1,  nchar(i)-3), "_outcome")
      
      if (i=="ckd_egfr40_pp") {
        cohort <- cohort %>%
          mutate({{censdate_var}}:=pmin(!!sym(outcome_var), cens_pp_3_yrs, na.rm=TRUE))
      } else {
        cohort <- cohort %>%
          mutate({{censdate_var}}:=pmin(!!sym(outcome_var), cens_pp, na.rm=TRUE))
      }
      
      
      
      cohort <- cohort %>%
        mutate({{censvar_var}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var)==!!sym(outcome_var), 1, 0),
               {{censtime_var}}:=as.numeric(difftime(!!sym(censdate_var), dstartdate, unit="days"))/365.25)
      
    }
    
    if (main_only==FALSE) {
      message(paste("survival variables for", paste(main_outcomes, collapse=", "), ",", paste(unlist(sensitivity_outcomes), collapse=", "), "added"))
    }
    
  }
  
  return(cohort) 
  
}