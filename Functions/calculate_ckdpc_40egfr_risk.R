calculate_ckdpc_40egfr_risk <- function(dataframe, age, sex, egfr, acr, sbp, bp_meds, hf, 
                                        chd, af, current_smoker, ex_smoker, bmi, hba1c, oha, insulin, 
                                        remote) {
  message("Note that values may be incorrect if 'remote' is not specified correctly (TRUE = on SQL server; FALSE = local in R)")
  
  # Convert column names to symbols
  age_col <- deparse(substitute(age))
  sex_col <- deparse(substitute(sex))
  egfr_col <- deparse(substitute(egfr))
  acr_col <- deparse(substitute(acr))
  sbp_col <- deparse(substitute(sbp))
  bp_meds_col <- deparse(substitute(bp_meds))
  hf_col <- deparse(substitute(hf))
  chd_col <- deparse(substitute(chd))
  af_col <- deparse(substitute(af))
  current_smoker_col <- deparse(substitute(current_smoker))
  ex_smoker_col <- deparse(substitute(ex_smoker))
  bmi_col <- deparse(substitute(bmi))
  hba1c_col <- deparse(substitute(hba1c))
  oha_col <- deparse(substitute(oha))
  insulin_col <- deparse(substitute(insulin))
  
  dataframe <- dataframe %>% mutate(id_col = row_number())
  new_dataframe <- data.frame()
  
  # Filter and calculate for eGFR < 60
  if (any(dataframe[[egfr_col]] < 60)) {
    new_dataframe_lt60 <- dataframe %>% mutate(join_col = 1L) %>% filter(!!sym(egfr_col) < 60)
    
    ckdpc_40egfr_risk_vars_lt60 <- data.frame(
      intercept = -3.2874,
      age_cons = -0.1725,
      male_cons = -0.1452,
      egfr_cons = -0.0769,
      acr_cons = 0.4658,
      sbp_cons = 0.2068,
      bp_med_cons = 0.1665,
      sbp_bp_med_cons = -0.0549,
      hf_cons = 0.4216,
      chd_cons = 0.2173,
      af_cons = 0.0456,
      current_smoker_cons = -0.0338,
      ex_smoker_cons = 0.1386,
      bmi_cons = 0.0253,
      hba1c_cons = -0.0031,
      oha_cons = -0.1232,
      insulin_cons = 0.0980,
      join_col = 1L
    )
    
    new_dataframe_lt60 <- new_dataframe_lt60 %>% inner_join(ckdpc_40egfr_risk_vars_lt60, by = "join_col")
    
    new_dataframe_lt60 <- new_dataframe_lt60 %>% mutate(
      male_sex = ifelse(!!sym(sex_col) == "male", 1L, 0L),
      hba1c_percent = (!!sym(hba1c_col) * 0.09148) + 2.152,
      oha_var = ifelse(!!sym(insulin_col) == 1, 1L, !!sym(oha_col)),
      ex_smoker_var = ifelse(!!sym(current_smoker_col) == 1, 0L, !!sym(ex_smoker_col)),
      acr_mgg = !!sym(acr_col) * 8.8402,
      log_acr_var = ifelse(remote, sql("LN(acr_mgg/10)"), log(acr_mgg/10)),
      ckdpc_40egfr_lin_predictor = ifelse(
        is.na(!!sym(acr_col)) | !!sym(acr_col) == 0, NA,
        (age_cons * ((!!sym(age_col) - 60) / 10)) + 
          ((male_cons * (male_sex - 0.5))) + 
          ((egfr_cons * ((!!sym(egfr_col) - 45) / 5))) + 
          (acr_cons * log_acr_var) + 
          (sbp_cons * ((!!sym(sbp_col) - 130) / 20)) + 
          (bp_med_cons * !!sym(bp_meds_col)) + 
          ((sbp_bp_med_cons * ((!!sym(sbp_col) - 130) / 20) * !!sym(bp_meds_col))) + 
          (hf_cons * (!!sym(hf_col) - 0.05)) + 
          (chd_cons * (!!sym(chd_col) - 0.15)) + 
          (af_cons * !!sym(af_col)) + 
          (current_smoker_cons * !!sym(current_smoker_col)) + 
          (ex_smoker_cons * ex_smoker_var) + 
          (bmi_cons * ((!!sym(bmi_col) - 30) / 5)) + 
          (hba1c_cons * (hba1c_percent - 7)) + 
          ((oha_cons * oha_var)) + 
          (insulin_cons * !!sym(insulin_col))
      ),
      ckdpc_40egfr_score = 100 * (exp(ckdpc_40egfr_lin_predictor + intercept) / (1 + exp(ckdpc_40egfr_lin_predictor + intercept)))
    )
    
    new_dataframe_lt60 <- new_dataframe_lt60 %>% select(id_col, ckdpc_40egfr_score, ckdpc_40egfr_lin_predictor)
    new_dataframe <- bind_rows(new_dataframe, new_dataframe_lt60)
  }
  
  # Filter and calculate for eGFR >= 60
  if (any(dataframe[[egfr_col]] >= 60)) {
    new_dataframe_ge60 <- dataframe %>% mutate(join_col = 1L) %>% filter(!!sym(egfr_col) >= 60)
    
    ckdpc_40egfr_risk_vars_ge60 <- data.frame(
      intercept = -4.2125,
      age_cons = 0.1465,
      male_cons = -0.2481,
      egfr_cons = -0.0562,
      acr_cons = 0.4098,
      sbp_cons = 0.1518,
      bp_med_cons = 0.2846,
      sbp_bp_med_cons = -0.0321,
      hf_cons = 0.9234,
      chd_cons = 0.2168,
      af_cons = 0.3087,
      current_smoker_cons = 0.1186,
      ex_smoker_cons = 0.0777,
      bmi_cons = 0.0308,
      hba1c_cons = 0.0987,
      oha_cons = -0.0638,
      insulin_cons = 0.2359,
      join_col = 1L
    )
    
    new_dataframe_ge60 <- new_dataframe_ge60 %>% inner_join(ckdpc_40egfr_risk_vars_ge60, by = "join_col")
    
    new_dataframe_ge60 <- new_dataframe_ge60 %>% mutate(
      male_sex = ifelse(!!sym(sex_col) == "male", 1L, 0L),
      hba1c_percent = (!!sym(hba1c_col) * 0.09148) + 2.152,
      oha_var = ifelse(!!sym(insulin_col) == 1, 1L, !!sym(oha_col)),
      ex_smoker_var = ifelse(!!sym(current_smoker_col) == 1, 0L, !!sym(ex_smoker_col)),
      acr_mgg = !!sym(acr_col) * 8.8402,
      log_acr_var = ifelse(remote, sql("LN(acr_mgg/10)"), log(acr_mgg/10)),
      ckdpc_40egfr_lin_predictor = ifelse(
        is.na(!!sym(acr_col)) | !!sym(acr_col) == 0, NA,
        (age_cons * ((!!sym(age_col) - 60) / 10)) + 
          ((male_cons * (male_sex - 0.5))) + 
          ((egfr_cons * ((!!sym(egfr_col) - 85) / 5))) + 
          (acr_cons * log_acr_var) + 
          (sbp_cons * ((!!sym(sbp_col) - 130) / 20)) + 
          (bp_med_cons * !!sym(bp_meds_col)) + 
          ((sbp_bp_med_cons * ((!!sym(sbp_col) - 130) / 20) * !!sym(bp_meds_col))) + 
          (hf_cons * (!!sym(hf_col) - 0.05)) + 
          (chd_cons * (!!sym(chd_col) - 0.15)) + 
          (af_cons * !!sym(af_col)) + 
          (current_smoker_cons * !!sym(current_smoker_col)) + 
          (ex_smoker_cons * ex_smoker_var) + 
          (bmi_cons * ((!!sym(bmi_col) - 30) / 5)) + 
          (hba1c_cons * (hba1c_percent - 7)) + 
          ((oha_cons * oha_var)) + 
          (insulin_cons * !!sym(insulin_col))
      ),
      ckdpc_40egfr_score = 100 * (exp(ckdpc_40egfr_lin_predictor + intercept) / (1 + exp(ckdpc_40egfr_lin_predictor + intercept)))
    )
    
    new_dataframe_ge60 <- new_dataframe_ge60 %>% select(id_col, ckdpc_40egfr_score, ckdpc_40egfr_lin_predictor)
    new_dataframe <- bind_rows(new_dataframe, new_dataframe_ge60)
  }
  
  dataframe <- dataframe %>% inner_join(new_dataframe, by = "id_col") %>% 
    select(-id_col)
  message("New columns 'ckdpc_40egfr_score' and 'ckdpc_40egfr_lin_predictor' added")
  return(dataframe)
}
