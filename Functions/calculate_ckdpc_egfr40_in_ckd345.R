calculate_ckdpc_egfr40_in_ckd345 <- function (dataframe, age, sex, egfr, acr, sbp, bp_meds, hf, 
          chd, af, current_smoker, ex_smoker, bmi, hba1c, oha, insulin, 
          remote) 
  
  ####NEED TO ADJUST CONSTANTS BELOW ACCORDING TO PAPER
  
{
  message("Note that values may be incorrect if 'remote' is not specified correctly (TRUE = on SQL server; FALSE = local in R)")
  age_col <- as.symbol(deparse(substitute(age)))
  sex_col <- as.symbol(deparse(substitute(sex)))
  egfr_col <- as.symbol(deparse(substitute(egfr)))
  acr_col <- as.symbol(deparse(substitute(acr)))
  sbp_col <- as.symbol(deparse(substitute(sbp)))
  bp_meds_col <- as.symbol(deparse(substitute(bp_meds)))
  hf_col <- as.symbol(deparse(substitute(hf)))
  chd_col <- as.symbol(deparse(substitute(chd)))
  af_col <- as.symbol(deparse(substitute(af)))
  current_smoker_col <- as.symbol(deparse(substitute(current_smoker)))
  ex_smoker_col <- as.symbol(deparse(substitute(ex_smoker)))
  bmi_col <- as.symbol(deparse(substitute(bmi)))
  hba1c_col <- as.symbol(deparse(substitute(hba1c)))
  oha_col <- as.symbol(deparse(substitute(oha)))
  insulin_col <- as.symbol(deparse(substitute(insulin)))
  dataframe <- dataframe %>% mutate(id_col = row_number())
  new_dataframe <- dataframe %>% mutate(join_col = 1L)
  ckdpc_40egfr_risk_vars <- data.frame(unlist(lapply(EHRBiomarkr::ckdpc40EgfrRiskConstants, 
                                                     function(y) lapply(y, as.numeric)), recursive = "FALSE")) %>% 
    mutate(join_col = 1L)
  new_dataframe <- new_dataframe %>% inner_join(ckdpc_40egfr_risk_vars, 
                                                by = "join_col", copy = TRUE)
  new_dataframe <- new_dataframe %>% mutate(male_sex = ifelse(!!sex_col == 
                                                                "male", 1L, 0L), hba1c_percent = (!!hba1c_col * 0.09148) + 
                                              2.152, oha_var = ifelse(!!insulin_col == 1, 1L, !!oha_col), 
                                            ex_smoker_var = ifelse(!!current_smoker_col == 1, 0L, 
                                                                   !!ex_smoker_col), acr_mgg = !!acr_col * 8.8402, 
                                            log_acr_var = if (remote == TRUE) 
                                              sql("LN(acr_mgg/10)")
                                            else log(acr_mgg/10), ckdpc_40egfr_lin_predictor = ifelse(is.na(!!acr_col) | 
                                                                                                        !!acr_col == 0, NA, (age_cons * ((!!age_col - 60)/10)) + 
                                                                                                        (-(male_cons * (male_sex - 0.5))) + (-(egfr_cons * 
                                                                                                                                                 ((!!egfr_col - 85)/5))) + (acr_cons * log_acr_var) + 
                                                                                                        (sbp_cons * ((!!sbp_col - 130)/20)) + (bp_med_cons * 
                                                                                                                                                 !!bp_meds_col) + (-(sbp_bp_med_cons * ((!!sbp_col - 
                                                                                                                                                                                           130)/20) * !!bp_meds_col)) + (hf_cons * (!!hf_col - 
                                                                                                                                                                                                                                      0.05)) + (chd_cons * (!!chd_col - 0.15)) + (af_cons * 
                                                                                                                                                                                                                                                                                    !!af_col) + (current_smoker_cons * !!current_smoker_col) + 
                                                                                                        (ex_smoker_cons * ex_smoker_var) + (bmi_cons * ((!!bmi_col - 
                                                                                                                                                           30)/5)) + (hba1c_cons * (hba1c_percent - 7)) + (-(oha_cons * 
                                                                                                                                                                                                               oha_var)) + (insulin_cons * !!insulin_col)), ckdpc_40egfr_score = 100 * 
                                              (exp(ckdpc_40egfr_lin_predictor + intercept)/(1 + 
                                                                                              exp(ckdpc_40egfr_lin_predictor + intercept))))
  new_dataframe <- new_dataframe %>% select(id_col, ckdpc_40egfr_score, 
                                            ckdpc_40egfr_lin_predictor)
  dataframe <- dataframe %>% inner_join(new_dataframe, by = "id_col") %>% 
    select(-id_col)
  message("New columns 'ckdpc_40egfr_score' and 'ckdpc_40egfr_lin_predictor' added")
  return(dataframe)
}