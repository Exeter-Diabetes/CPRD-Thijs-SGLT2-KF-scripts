
# Produce survival variables for all endpoints 
## All censored at 3 years post drug start / end of GP records / death / starting a SGLT2/GLP1 if not taking one
## Additional variables where maximal censoring time is lengthened to 5 years

# Main analysis:
## 'ckd_egfr40': decline in eGFR of >=40% from baseline or onset of CKD stage 5 OR death from renal causes
## 'ckd_egfr50': decline in eGFR of >=50% from baseline or onset of CKD stage 5 OR death from renal causes
## 'death': all-cause mortality
## 'macroalb': progression to albuminuria >=30mg/mmol
## 'dka': occurrence of diabetic ketoacidosis
## 'amputation': major limb amputation
## 'side_effect': mycotic genital infection

# Sensitivity analysis:
## '{outcome}_pp': all of main analysis but per-protocol rather than intention to treat (see below)
## intention to treat: censoring if starting an SGLT2 inhibitor (if in DPP4 or SU arm) or GLP1 agonist.
## per-protocol: censoring if starting any other treatment arm or GLP1 agonist.
## the per-protocol analyses are only used to compare treatment effects of DPP4i vs SU.
## no censoring is undertaken if patients do not continue their prescription.


add_surv_vars <- function(cohort_dataset, main_only=FALSE) {
  
  # Add survival variables for outcomes for main analysis
  main_outcomes <- c("ckd_egfr40", "ckd_egfr50", "death", "macroalb", "dka", "amputation", "side_effect")
  
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
           
           ckd_egfr40_outcome=pmin(egfr_40_decline_date,
                                   postckdstage5date,
                                   kf_death_date_any_cause,
                                   na.rm=TRUE),
           
           ckd_egfr50_outcome=pmin(egfr_50_decline_date,
                                   postckdstage5date,
                                   kf_death_date_any_cause,
                                   na.rm=TRUE),
           
           ckd_egfr50_outcome_type=ifelse(
             is.na(ckd_egfr50_outcome), 
             NA,
             ifelse(ckd_egfr50_outcome == egfr_50_decline_date, 
                    "50% eGFR decline",
                    ifelse(ckd_egfr50_outcome == postckdstage5date, 
                           "ESKD", 
                           "Kidney-related death")
             )
           ),
           
           macroalb_outcome=macroalb_date,
           
           dka_outcome=postdrug_first_dka,
           
           amputation_outcome=postdrug_first_amputation,
           
           side_effect_outcome=pmin(postdrug_first_medspecific_gi,
                                    na.rm=TRUE),
    
           death_outcome=death_date)
  
  
  for (i in main_outcomes) {
    
    outcome_var=paste0(i, "_outcome")
    censdate_var=paste0(i, "_censdate")
    censvar_var=paste0(i, "_censvar")
    censtime_var=paste0(i, "_censtime_yrs")
    censdate_var_5y=paste0(i, "_5y_censdate")
    censvar_var_5y=paste0(i, "_5y_censvar")
    censtime_var_5y=paste0(i, "_5y_censtime_yrs")
    

    
    cohort <- cohort %>%
        mutate({{censdate_var}}:=pmin(!!sym(outcome_var), cens_itt_3_yrs, na.rm=TRUE),
               {{censvar_var}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var)==!!sym(outcome_var), 1, 0),
               {{censtime_var}}:=as.numeric(difftime(!!sym(censdate_var), dstartdate, unit="days"))/365.25,
               {{censdate_var_5y}}:=pmin(!!sym(outcome_var), cens_itt, na.rm=TRUE),
               {{censvar_var_5y}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var_5y)==!!sym(outcome_var), 1, 0),
               {{censtime_var_5y}}:=as.numeric(difftime(!!sym(censdate_var_5y), dstartdate, unit="days"))/365.25)

    
  }
  
  if (main_only==TRUE) {
    message(paste("survival variables for", paste(main_outcomes, collapse=", "), "added"))
  }
  
  
  # Add survival variables for outcomes for sensitivity analyses
  
  else {
    
    # Split by whether ITT or PP
    sensitivity_outcomes <- c("ckd_egfr40_pp", "ckd_egfr50_pp")
    
    
    for (i in sensitivity_outcomes) {
      
      censdate_var=paste0(i, "_censdate")
      censvar_var=paste0(i, "_censvar")
      censtime_var=paste0(i, "_censtime_yrs")
      
      
      outcome_var=paste0(substr(i, 1,  nchar(i)-3), "_outcome")
      
        cohort <- cohort %>%
          mutate({{censdate_var}}:=pmin(!!sym(outcome_var), cens_pp_3_yrs, na.rm=TRUE))
      
      
      
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