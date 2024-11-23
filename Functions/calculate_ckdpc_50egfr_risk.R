calculate_ckdpc_50egfr_risk <- function(dataframe, age, sex, egfr, acr, sbp, bp_meds, hf, 
                                            chd, af, current_smoker, ex_smoker, bmi, hba1c, oha, insulin) {

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
    
    ckdpc_50egfr_risk_vars_lt60 <- data.frame(
      intercept = -3.896094,
      age_cons = -0.284654,
      male_cons = -0.0596991,
      egfr_cons = -0.1079417,
      acr_cons = 0.5004966,
      sbp_cons = 0.2658852,
      bp_med_cons = 0.2755529,
      sbp_bp_med_cons = -0.1184013,
      hf_cons = 0.3461398,
      chd_cons = 0.1899315,
      af_cons = -0.0792505,
      current_smoker_cons = 0.0506796,
      ex_smoker_cons = 0.0884935,
      bmi_cons = 0.0068514,
      hba1c_cons = 0.001295,
      oha_cons = -0.1892067,
      insulin_cons = 0.0665215,
      join_col = 1L
    )
    
    new_dataframe_lt60 <- new_dataframe_lt60 %>% inner_join(ckdpc_50egfr_risk_vars_lt60, by = "join_col")
    
    new_dataframe_lt60 <- new_dataframe_lt60 %>% mutate(
      male_sex = ifelse(!!sym(sex_col) == "male", 1L, 0L),
      hba1c_percent = (!!sym(hba1c_col) * 0.09148) + 2.152,
      oha_var = ifelse(!!sym(insulin_col) == 1, 1L, !!sym(oha_col)),
      ex_smoker_var = ifelse(!!sym(current_smoker_col) == 1, 0L, !!sym(ex_smoker_col)),
      acr_mgg = !!sym(acr_col) * 8.8402,
      log_acr_var = log(acr_mgg/10),
      ckdpc_50egfr_lin_predictor = ifelse(
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
      ckdpc_50egfr_score = 100 * (exp(ckdpc_50egfr_lin_predictor + intercept) / (1 + exp(ckdpc_50egfr_lin_predictor + intercept)))
    )
    
    new_dataframe_lt60 <- new_dataframe_lt60 %>% select(id_col, ckdpc_50egfr_score, ckdpc_50egfr_lin_predictor)
    new_dataframe <- bind_rows(new_dataframe, new_dataframe_lt60)
  }
  
  # Filter and calculate for eGFR >= 60
  if (any(dataframe[[egfr_col]] >= 60)) {
    new_dataframe_ge60 <- dataframe %>% mutate(join_col = 1L) %>% filter(!!sym(egfr_col) >= 60)
    
    ckdpc_50egfr_risk_vars_ge60 <- data.frame(
      intercept = -4.941086,
      age_cons = 0.0794321,
      male_cons = -0.1141482,
      egfr_cons = -0.0459708,
      acr_cons = 0.4326073,
      sbp_cons = 0.1917602,
      bp_med_cons = 0.306044,
      sbp_bp_med_cons = -0.0670274,
      hf_cons = 0.9185919,
      chd_cons = 0.2386131,
      af_cons = 0.2157927,
      current_smoker_cons = 0.1122236,
      ex_smoker_cons = 0.1327815,
      bmi_cons = 0.0154271,
      hba1c_cons = 0.1103419,
      oha_cons = -0.0460205,
      insulin_cons = 0.2571656,
      join_col = 1L
    )
    
    new_dataframe_ge60 <- new_dataframe_ge60 %>% inner_join(ckdpc_50egfr_risk_vars_ge60, by = "join_col")
    
    new_dataframe_ge60 <- new_dataframe_ge60 %>% mutate(
      male_sex = ifelse(!!sym(sex_col) == "male", 1L, 0L),
      hba1c_percent = (!!sym(hba1c_col) * 0.09148) + 2.152,
      oha_var = ifelse(!!sym(insulin_col) == 1, 1L, !!sym(oha_col)),
      ex_smoker_var = ifelse(!!sym(current_smoker_col) == 1, 0L, !!sym(ex_smoker_col)),
      acr_mgg = !!sym(acr_col) * 8.8402,
      log_acr_var = log(acr_mgg/10),
      ckdpc_50egfr_lin_predictor = ifelse(
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
      ckdpc_50egfr_score = 100 * (exp(ckdpc_50egfr_lin_predictor + intercept) / (1 + exp(ckdpc_50egfr_lin_predictor + intercept)))
    )
    
    new_dataframe_ge60 <- new_dataframe_ge60 %>% select(id_col, ckdpc_50egfr_score, ckdpc_50egfr_lin_predictor)
    new_dataframe <- bind_rows(new_dataframe, new_dataframe_ge60)
  }
  
  dataframe <- dataframe %>% inner_join(new_dataframe, by = "id_col") %>% 
    select(-id_col)
  message("New columns 'ckdpc_50egfr_score' and 'ckdpc_50egfr_lin_predictor' added")
  return(dataframe)
}
