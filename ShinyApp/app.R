



library(shiny)
library(tidyverse)

# Function to calculate pARR
calculate_parr <- function(dataframe, age, sex, egfr, acr, sbp, bp_meds, hf, chd, af, smoking_status, diabetes_med, bmi, hba1c) {
  
  trial_hr_kf_sglt2i <- 0.62
  
  # constants for risk calculation
  ckdpc_50egfr_risk_vars <- list(
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
    insulin_cons = 0.2571656
  )
  # coefficients for risk calculation
  dataframe <- dataframe %>%
    mutate(
      male_sex = ifelse(sex == "Male", 1L, 0L),
      current_smoker = ifelse(smoking_status == "Current Smoker", 1, 0),
      ex_smoker = ifelse(smoking_status == "Ex Smoker", 1, 0),
      oha = ifelse(diabetes_med == "Oral Medications Only", 1, 0),
      insulin = ifelse(diabetes_med == "Insulin", 1, 0),
      hba1c_percent = (as.numeric(hba1c) * 0.09148) + 2.152,
      oha_var = ifelse(insulin == 1, 1L, oha),
      ex_smoker_var = ifelse(current_smoker == 1, 0L, ex_smoker),
      acr_mgg = as.numeric(acr) * 8.8402,
      log_acr_var = log(acr_mgg / 10),
      ckdpc_50egfr_lin_predictor = (ckdpc_50egfr_risk_vars$age_cons * ((as.numeric(age) - 60) / 10)) +
        (ckdpc_50egfr_risk_vars$male_cons * (male_sex - 0.5)) +
        (ckdpc_50egfr_risk_vars$egfr_cons * ((as.numeric(egfr) - 85) / 5)) +
        (ckdpc_50egfr_risk_vars$acr_cons * log_acr_var) +
        (ckdpc_50egfr_risk_vars$sbp_cons * ((as.numeric(sbp) - 130) / 20)) +
        (ckdpc_50egfr_risk_vars$bp_med_cons * bp_meds) +
        (ckdpc_50egfr_risk_vars$sbp_bp_med_cons * ((as.numeric(sbp) - 130) / 20) * bp_meds) +
        (ckdpc_50egfr_risk_vars$hf_cons * (hf - 0.05)) +
        (ckdpc_50egfr_risk_vars$chd_cons * (chd - 0.15)) +
        (ckdpc_50egfr_risk_vars$af_cons * af) +
        (ckdpc_50egfr_risk_vars$current_smoker_cons * current_smoker) +
        (ckdpc_50egfr_risk_vars$ex_smoker_cons * ex_smoker_var) +
        (ckdpc_50egfr_risk_vars$bmi_cons * ((as.numeric(bmi) - 30) / 5)) +
        (ckdpc_50egfr_risk_vars$hba1c_cons * (hba1c_percent - 7)) +
        (ckdpc_50egfr_risk_vars$oha_cons * oha_var) +
        (ckdpc_50egfr_risk_vars$insulin_cons * insulin),
      
      ckdpc_50egfr_score = 100 * (exp(ckdpc_50egfr_lin_predictor + ckdpc_50egfr_risk_vars$intercept) / (1 + exp(ckdpc_50egfr_lin_predictor + ckdpc_50egfr_risk_vars$intercept))),
      
      ckdpc_50egfr_survival=(100-ckdpc_50egfr_score)/100,
      ckdpc_50egfr_survival_sglt2i=ckdpc_50egfr_survival^trial_hr_kf_sglt2i,
      parr = (ckdpc_50egfr_survival_sglt2i - ckdpc_50egfr_survival) * 100  # Multiply by 100 to get percentage
    )
  
  return(dataframe)
}


# UI for the Shiny app
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      .title {
      text-align: center;
      margin-top: 20px;
      max-width: 700px; 
      min-width: 650px
      }
    .beta-notice {
      text-align: center;
      font-size: 12px;
      font-style: italic;
      color: red;
      margin-top: 20px;
    }
    .explanation {
      text-align: center;
      font-size: 14px;
      margin-top: 10px;
      max-width: 700px; 
      min-width: 650px
    }    
    .subtitle {
      text-align: center;
      font-size: 18px;
      font-weight: bold;
      margin-top: 10px;
      max-width: 700px; 
      min-width: 650px
    }
      .result-text {
        font-size: 48px;
        font-weight: bold;
        text-align: center;
        margin-bottom: 0px;
      }
      .blue { color: #0072B2; text-shadow: 1px 1px 2px black; }
      .red { color: #FF5733; text-shadow: 1px 1px 2px black; }
      .green { color: #D55E00; text-shadow: 1px 1px 2px black; }
      .validation-message {
        font-size: 16px;
        font-weight: normal;
        color: black;
        text-align: center;
        margin-top: 0px;
        margin-bottom: 0px;
      }
      .references-title {
        font-weight: bold;
      max-width: 700px; 
      min-width: 650px
      }
    "))
  ),
  column(width = 12, align = "center", div(class="title", titlePanel("Predicted kidney protection benefit from SGLT2-inhibitors in people with type 2 diabetes"))),
  column(width = 12, align = "center", div(class="beta-notice", 
                                           HTML("This calculator is in beta and is for research purposes only. Please leave any feedback <a href='https://forms.gle/ewYM9jmQfNb2fD2U8' target='_blank'>here</a>."))),
  column(width = 12, align = "center", div(class="explanation", uiOutput("dynamic_explanation"))),
  column(width = 12, align = "center", div(class="subtitle", h4("3-year risk of kidney disease progression (≥50% decline in eGFR or kidney failure):"))),
  column(width = 12, align="center",
         div(style = "max-width: 700px; min-width: 650px", # these values set the min and max width for the elements inside
             uiOutput("result_text"),  # Result displayed here
             plotOutput("risk_plot", height = "120px"),  # Add the bar chart
             uiOutput("traffic_light"),
         )
  ),
  fluidRow(
    column(12, align = "center",
           div(style = "max-width: 650px; display: flex; justify-content: space-between; text-align: center;",
               column(6, 
                      selectInput("age", "Age (years):",
                                  choices = as.list(seq(20, 80, 1)), selected = 50),
                      selectInput("hba1c", "HbA1c (mmol/mol):", 
                                  choices = as.list(c(seq(42, 96, 1), "≥97 (model will use 97 as input)")), selected = 58),
                      selectInput("bmi", "BMI (kg/m2):", 
                                  choices = as.list(c(seq(20, 40, 1), "≥40 (model will use 40 as input)")), selected = 30),
                      selectInput("egfr", "eGFR (ml/min per 1.73m2):", 
                                  choices = as.list(c("<60", seq(60, 140, 1))), selected = 90),
                      selectInput(
                        inputId = "acr",
                        label = "ACR (mg/mmol)",
                        choices = as.list(c("≤0.6 or incalculable", seq(1, 29.5, 0.5), "≥30")),  # Default choices
                        selected = 1
                      ),
                      actionLink(
                        inputId = "toggle_units",
                        label = "Change to conventional units:",  # Dynamic label
                        style = "display: block; margin-top: 1px; margin-bottom: 15px; color: #007BFF; text-decoration: underline; cursor: pointer;"  # Styling for hyperlink
                      )
               ),
               column(6, 
                      selectInput("sex", "Sex:", 
                                  choices = list(Male = "Male", Female = "Female")),
                      selectInput("diabetes_med", "Diabetes Medication:", 
                                  choices = list("No diabetes medications" = "No Diabetes Medication", 
                                                 "Oral medications only" = "Oral Medications Only", 
                                                 "Insulin (with or without oral medications)" = "Insulin")),
                      selectInput("smoking_status", "Smoking Status:", 
                                  choices = list("Never smoker" = "Never Smoker", 
                                                 "Ex smoker" = "Ex Smoker", 
                                                 "Current smoker" = "Current Smoker")),
                      #      checkboxInput("type2_dm", "Type 2 diabetes", value = TRUE),
                      selectInput("sbp", "SBP (mmHg):", 
                                  choices = as.list(c(seq(90, 179, 1), "≥180 (model will use 180 as input)")), selected = 130),
                      checkboxInput("bp_meds", "On antihypertensive medications", value = FALSE),
                      checkboxInput("hf", "History of heart failure", value = FALSE),
                      checkboxInput("chd", "History of ischaemic heart disease, stroke, or peripheral vascular disease", value = FALSE),
                      checkboxInput("af", "History of atrial fibrillation", value = FALSE)
               )
           )
    ),
  ),
  div(style = "text-align:center; margin-top:20px;", 
      actionButton("predict", "Predict", class = "btn-primary btn-lg")
  ),
  column(width = 12, align = "center",
         div(style = "text-align:center; margin-top:50px; max-width: 700px; min-width: 650px",
             p("This prediction model is based on the CKD Prognosis Consortium risk score for ≥50% decline in eGFR or kidney failure over 3 years integrated with the relative treatment effect from SGLT2-inhibitor trial meta-analysis. The model was independently validated using UK routine general practice data of 141,500 participants with type 2 diabetes, preserved eGFR (≥60mL/min/1.73m2), normal or low-level albuminuria (<30mg/mmol), and without a history of atherosclerotic vascular disease or heart failure. In this population, a predicted benefit threshold of 0.65% (corresponding to an NNT 154) would prevent over 10% more kidney disease progression events than using an albuminuria threshold ≥3mg/mmol, while targeting a comparable proportion of the population."),
             p(class = "references-title", "References:"),
             p("Nuffield Department of Population Health Renal Studies Group, SGLT2 inhibitor Meta-Analysis Cardio-Renal Trialists' Consortium. Impact of diabetes on the effects of sodium glucose co-transporter-2 inhibitors on kidney outcomes: collaborative meta-analysis of large placebo-controlled trials. Lancet 2022; 400(10365): 1788-801."),
             p("Grams ME, Brunskill NJ, Ballew SH, et al. Development and Validation of Prediction Models of Adverse Kidney Outcomes in the Population With and Without Diabetes. Diabetes Care 2022; 45(9): 2055-63.")
         )
  ),

)

# Server logic for the Shiny app
server <- function(input, output, session) {
  # Ensure decimal separator uses a period
  Sys.setlocale("LC_NUMERIC", "C")
  
  # Track the current units (default: SI)
  unit <- reactiveVal("SI")
  
  # Toggle unit state when the action link is clicked
  observeEvent(input$toggle_units, {
    if (unit() == "SI") {
      unit("Conventional")
    } else {
      unit("SI")
    }
  })
  
  output$dynamic_explanation <- renderUI({
    acr_limit_value <- if (unit() == "SI") "<30mg/mmol" else "<265mg/g"
    HTML(paste0("This calculator predicts kidney protection benefit from SGLT2-inhibitor treatment for people with type 2 diabetes, preserved eGFR (<60mL/min/1.73m2), normal or low-level albuminuria (uACR", acr_limit_value, "), and without a history of atherosclerotic vascular disease or heart failure."))
  })
  
  # Dynamically update ACR choices based on the unit
  observe({
    if (unit() == "SI") {
      updateSelectInput(session, "acr", label = "ACR (mg/mmol):",
                        choices = as.list(c("≤0.6 or incalculable", seq(1, 29.5, 0.5), "≥30")), selected = 1)
      updateSelectInput(session, "hba1c", label = "HbA1c (mmol/mol):", 
                        choices = as.list(c(seq(42, 96, 1), "≥97 (model will use 97 as input)")), selected = 58)
      updateActionButton(session, "toggle_units", label = "Change to conventional units")
    } else {
      updateSelectInput(session, "acr", label = "ACR (mg/g):",
                        choices = as.list(c("≤5.3 or incalculable", seq(10, 260, 5), "≥265")), selected = 10)
      updateSelectInput(session, "hba1c", label = "HbA1c (%):", 
                        choices = as.list(c(seq(6, 10.9, 0.1), "≥11 (model will use 11 as input)")), selected = 7.5)
      updateActionButton(session, "toggle_units", label = "Change to SI units")
      
    }
  })
  
  # EventReactive for input processing
  patient <- eventReactive(input$predict, {
    # Convert ACR back to mg/mmol if needed
    acr_value <- if (unit() == "Conventional") {
      if (input$acr == "≤5.3 or incalculable") {
        0.6  # Corresponding value in mg/mmol
      } else if (input$acr == "≥265") {
        30
    } else {
        as.numeric(input$acr) / 8.84  # Convert mg/g to mg/mmol
      }
    } else {
      if (input$acr == "≤0.6 or incalculable") {
        0.6
      } else if (input$acr == "≥30") {
        30
      } else {
        as.numeric(input$acr)
      }
    }
    
    hba1c_value <- if (unit() == "Conventional") {
      if (input$hba1c == "≥11 (model will use 11 as input)") {
        97
      } else 
        {
      (as.numeric(input$hba1c) - 2.15) * 10.929 # Convert % to mmol/mol
      }
    } else 
      {
      if (input$hba1c == "≥97 (model will use 97 as input)") {
        97
      } else 
        {
      as.numeric(input$hba1c)
      }
    }
    
    egfr_value <- if (input$egfr == "<60") {
      59 } else { as.numeric(input$egfr) } 
    
    sbp_value <- if (input$sbp == "≥180 (model will use 180 as input)") {
      180
    } else { as.numeric(input$sbp) }
    
    bmi_value <- if (input$bmi == "≥40 (model will use 40 as input)") {
      40
    } else { as.numeric(input$bmi) }
    
    # Combine everything into a data.frame
    data.frame(
      age = input$age,
      sex = input$sex,
      egfr = egfr_value,
      acr = acr_value,  # Use the converted ACR value
      sbp = sbp_value,
      bp_meds = as.numeric(input$bp_meds),
      hf = as.numeric(input$hf),
      chd = as.numeric(input$chd),
      af = as.numeric(input$af),
      smoking_status = input$smoking_status,
      diabetes_med = input$diabetes_med,
      bmi = bmi_value,
      hba1c = hba1c_value
    )
  })
  
  result <- eventReactive(input$predict, {
    patient <- patient()
    if (patient$hf | patient$chd | patient$egfr < 60 | patient$acr >= 30
    ) {
      return(NULL)
    }
    calculate_parr(
      dataframe = data.frame(
        age = patient$age,
        sex = patient$sex,
        egfr = patient$egfr,
        acr = patient$acr,
        sbp = patient$sbp,
        bp_meds = as.numeric(patient$bp_meds),
        hf = as.numeric(patient$hf),
        chd = as.numeric(patient$chd),
        af = as.numeric(patient$af),
        smoking_status = patient$smoking_status,
        diabetes_med = patient$diabetes_med,
        bmi = patient$bmi,
        hba1c = patient$hba1c
      ),
      age = patient$age,
      sex = patient$sex,
      egfr = patient$egfr,
      acr = patient$acr,
      sbp = patient$sbp,
      bp_meds = as.numeric(patient$bp_meds),
      hf = as.numeric(patient$hf),
      chd = as.numeric(patient$chd),
      af = as.numeric(patient$af),
      smoking_status = patient$smoking_status,
      diabetes_med = patient$diabetes_med,
      bmi = patient$bmi,
      hba1c = patient$hba1c
    )
  })
  
  
  parr_threshold <- 0.65059
  
  output$result_text <- renderUI({

    
    patient <- patient()
    if (input$predict == 0) {
      NULL
      # } else if (!patient$type2_dm) {
      #   div(class = "validation-message", 
      #       "This prediction model is validated for individuals with type 2 diabetes only.")
    } else if (patient$hf || patient$chd) {
      div(class = "validation-message", 
          "This patient has an indication for SGLT2-inhibitor treatment due to a history of atherosclerotic vascular disease or heart failure.")
    } else if (patient$acr >= 30) {
      div(class = "validation-message",
          "This patient has an indication for SGLT2-inhibitor treatment due to significant albuminuria.")
    } else if (patient$egfr < 60) {
      div(class = "validation-message",
          "This patient has an indication for SGLT2-inhibitor treatment due to an eGFR <60mL/min/1.73m2.")
    } else {
      df <- result()
      if (is.null(df)) {
        div(class = "validation-message", 
            "This prediction model cannot calculate a valid result with the current inputs.")
      } else {
        current_risk <- df$ckdpc_50egfr_score[1]
        risk_with_treatment <- df$ckdpc_50egfr_score[1] - df$parr[1]
        parr <- df$parr[1]
        nnt <- ifelse(df$parr[1] > 0, round(1 / (df$parr[1]/100), 0), NA)
        
        
        result_colour <- ifelse(parr < parr_threshold, "#FF5733", "#00C853") 
        
        div(
          div(style = "display: flex; justify-content: space-between; text-align: center;",
              div(class = "result-text", style = "width: 33%;",
                  h4("Risk without SGLT2i"),
              ),
              div(class = "result-text", style = "width: 33%;",
                  h4("Risk with SGLT2i"),
              ),
              div(class = "result-text", style = "width: 33%; margin-bottom: 0px; padding: 0; ",
                  h4("3-year NNT"),
              ),
          ),
          div(style = "display: flex; justify-content: space-between; text-align: center; margin-top: 0px;",
              div(class = "result-text", style = "width: 33%;",
                  span(style = "color: #0072B2; text-shadow: 1px 1px 2px black;",
                       sprintf("%.1f%%", current_risk))
              ),
              div(class = "result-text", style = "width: 33%;",
                  span(style = paste("color: ", result_colour, "; text-shadow: 1px 1px 2px black;"),
                       sprintf("%.1f%%", risk_with_treatment))
              ),
              div(class = "result-text", style = "width: 33%;",
                  span(style = paste("color: ", result_colour, "; text-shadow: 1px 1px 2px black;"),
                       ifelse(!is.na(nnt), sprintf("%d", nnt), "N/A")),
              ),
          ),
          
        )
        
        
        
      }
    }
  })
  
  output$risk_plot <- renderPlot({
    patient <- patient()
    if (input$predict == 0 #|| !patient$type2_dm 
        || patient$hf || patient$chd || patient$egfr < 60 || patient$acr >= 30) {
      return(NULL)
    }
    
    if (is.null(result())) return()
    
    # Extract the necessary values
    df <- result()
    current_risk <- df$ckdpc_50egfr_score[1]
    treated_risk <- 100 * (1 - df$ckdpc_50egfr_survival_sglt2i[1])
    parr <- df$parr[1]
    
    bar_colour <- ifelse(parr < parr_threshold, "#FF5733", "#00C853")
    
    # Create a data frame for the two bars
    plot_data <- data.frame(
      xmin = c(0, 0),  # Add a small gap for yellow bar inset
      xmax = c(treated_risk, current_risk),
      ymin = c(0.2, 0.45),  # 
      ymax = c(0.4, 0.65),  # 
      fill = c(bar_colour, "#0072B2")  # Blue for untreated, yellow for treated
    )
    
    # Generate the plot
    ggplot() +
      
      # Full-width blue bar (untreated risk)
      geom_rect(data = plot_data[1, ], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
                color = "black", 
                size = 0.2) +  # No border for blue bar
      
      
      # Narrower yellow bar (treated risk)
      geom_rect(data = plot_data[2, ], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
                color = "black", 
                size = 0.2) +  # 
      
      # Use specified colors for the bars
      scale_fill_identity() +
      
      # X-axis configuration
      scale_x_continuous(limits = c(0, max(current_risk, 5)), 
                         breaks = seq(0, max(current_risk, 5), 1), 
                         labels = paste0(seq(0, max(current_risk, 5), 1), "%")) +
      
      # Minimal theme with adjusted aesthetics
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
      ) # +
      
      # Chart title
      #labs(title = "3-year risk of kidney disease progression")
  })
  
  output$traffic_light <- renderUI({
    if (input$predict > 0 && !is.null(result())) {
      
      df <- result()
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      parr <- df$parr[1]
      parr_threshold <- 0.65059
      is_above_threshold <- parr < parr_threshold
      
      # Render the traffic light and text
      div(
        style = "display: flex; align-items: flex-start; margin-top: 40px; margin-bottom: 40px",  # Align top edges of elements
        
        # Traffic Light
        div(
          style = "width: 90px; height: 120px; background-color: #333; border: 2px solid black; border-radius: 5px; padding: 10px; margin-right: 20px; position: relative;",
          
          # Top Light (Red or Grey)
          div(style = sprintf("width: 30px; height: 30px; background-color: %s; border-radius: 50%%; position: absolute; top: 15px; left: 50%%; transform: translateX(-50%%);",
                              ifelse(!is_above_threshold, "#808080", "#FF5733")),
              
          ),
          
          # Bottom Light (Green or Grey)
          div(style = sprintf("width: 30px; height: 30px; background-color: %s; border-radius: 50%%; position: absolute; bottom: 15px; left: 50%%; transform: translateX(-50%%);",
                              ifelse(!is_above_threshold, "#00C853", "#808080")),
              
          )
        ),
        
        # Text and Arrows
        div(
          style = "display: flex; flex-direction: column; justify-content: space-between; height: 120px;",
          
          # Text aligned with Red Light (if below threshold)
          if (is_above_threshold) {
            div(style = "display: flex; align-items: center; margin-bottom: 15px;",
                div(style = "font-size: 20px; font-weight: bold; margin-right: 5px; color: black;", "→"),  # Unicode arrow
                div(style = "font-size: 14px;", "SGLT2-inhibitor treatment not suggested for kidney protection in this patient. This patient's predicted benefit is below the threshold of 0.65%. Targeting treatment using this threshold would result in the same number of individuals treated, but over 10% more kidney disease progression events prevented over 3 years than when using a ≥3mg/mmol albuminuria threshold (currently recommended in NICE, KDIGO, ADA, and EASD guidelines).")
            )
          },
          
          # Spacer to align Green Light text
          div(style = "flex-grow: 1;"),
          
          # Text aligned with Green Light (if above threshold)
          if (!is_above_threshold) {
            div(style = "display: flex; align-items: center; margin-top: 15px;",
                div(style = "font-size: 20px; font-weight: bold;margin-right: 5px; color: black;", "→"),  # Unicode arrow
                div(style = "font-size: 14px;", "Consider SGLT2-inhibitor treatment for kidney protection in this patient. This patient's predicted benefit is above the threshold of 0.65%. Targeting treatment using this threshold would result in the same number of individuals treated, but over 10% more kidney disease progression events prevented over 3 years than when using a ≥3mg/mmol albuminuria threshold (currently recommended in NICE, KDIGO, ADA, and EASD guidelines).")
            )
          }
        )
      )
    
  
    # 
    #   
    # div(style = "display: flex; align-items: flex-start; margin-top: 40px; margin-bottom: 40px; ",
    #     
    #     div(
    #       style = "width: 90px; height: 120px; background-color: #333; border: 2px solid black; border-radius: 5px; padding: 10px; margin-right: 20px; position: relative;",
    #       
    #       # Red Light
    #       div(style = "width: 30px; height: 30px; background-color: %s; border-radius: 50%; position: absolute; top: 15px; left: 50%; transform: translateX(-50%);",
    #           ifelse(!is_above_threshold, "#808080", "#FF5733")),
    #       
    #       # Green Light
    #       div(style = "width: 30px; height: 30px; background-color: %s; border-radius: 50%; position: absolute; bottom: 15px; left: 50%; transform: translateX(-50%);",
    #           ifelse(!is_above_threshold, "#00C853", "#808080")),
    #     ),
    #     
    # 
    #     
    #     # Text and Arrows
    #     div(
    #       style = "display: flex; flex-direction: column; justify-content: space-between; height: 120px;",
    #       
    #       # Text aligned with Red Light
    #       div(style = "display: flex; align-items: center; margin-top: 15px;",
    #           div(style = "font-size: 20px; font-weight: bold; margin-right: 5px; color: black;", "→"),
    #           div(style = "font-size: 14px;", "SGLT2-inhibitor treatment not suggested for kidney protection. This indicates a predicted benefit below the threshold (0.65%) where a comparable proportion of the population would be treated as currently recommended by NICE and international guidelines (KDIGO/ADA/EASD).")
    #       ),
    #       
    #       # Spacer to align Green Light text
    #       div(style = "flex-grow: 1;"),
    #       
    #       # Text aligned with Green Light
    #       div(style = "display: flex; align-items: center; margin-top: 15px;",
    #           div(style = "font-size: 20px; font-weight: bold; margin-right: 5px; color: black;", "→"),
    #           div(style = "font-size: 14px;", "Consider SGLT2-inhibitor treatment for kidney protection. This indicates a predicted benefit above the threshold (0.65%) where a comparable proportion of the population would be treated as currently recommended by NICE and international guidelines (KDIGO/ADA/EASD).")
    #       )
    #     )
    # )
        
    #     # Text aligned with the lights
    #     div(
    #       style = "display: flex; flex-direction: column; justify-content: space-between; height: 120px;",
    #       
    #       # Text aligned with Red Light
    #       div(style = "margin-top: 15px; font-size: 14px;", "SGLT2-inhibitor treatment not suggested for kidney protection."),
    #       
    #       # Spacer to align Green Light text
    #       div(style = "flex-grow: 1;"),
    #       
    #       # Text aligned with Green Light
    #       div(style = "margin-top: 15px; font-size: 14px;", "Consider SGLT2-inhibitor treatment for kidney protection. This indicates a predicted benefit above the threshold (0.65%) where a comparable proportion of the population would be treated as currently recommended by NICE and international guidelines (KDIGO/ADA/EASD).")
    #     )
    # )
        
    #     # Traffic Light
    #     div(style = "width: 90px; height: 100px; background-color: #333; border: 2px solid black; border-radius: 10px; display: flex; flex-direction: column; justify-content: space-around; align-items: center; margin-right: 0px;",
    #         div(style = "width: 40px; height: 40px; background-color: #FF5733; border: 0.5px solid black; border-radius: 50%;"),
    #         div(style = "width: 40px; height: 40px; background-color: #00C853; border: 0.5px solid black; border-radius: 50%;")
    #     ),
    #     
    #     # Explanatory Text
    #     div(
    #       style = "display: flex; flex-direction: column; justify-content: space-between; height: 120px;",
    #       div(style = "color: red; font-weight: bold;", " SGLT2-inhibitor treatment not suggested for kidney protection."),
    #       div(style = "color: green; font-weight: bold;", " Consider SGLT2-inhibitor treatment for kidney protection. This indicates a predicted benefit above the threshold (0.65%) where a comparable proportion of the population would be treated as currently recommended by NICE and international guidelines (KDIGO/ADA/EASD).")
    #     )
    # )
    #     
        
        # div(style = "flex-grow: 1; margin-left: 00px; font-size: 14px; color: black;",
        #     
        #     tags$ul(
        #       tags$li(tags$b("Red Light:"), " SGLT2-inhibitor treatment not suggested for kidney protection."),
        #       tags$li(tags$b("Green Light:"), " Consider SGLT2-inhibitor treatment for kidney protection. This indicates a predicted benefit above the threshold (0.65%) where a comparable proportion of the population would be treated as currently recommended by NICE and international guidelines (KDIGO/ADA/EASD).")
        #     )
        # )
        # )
  
  } else {
      NULL  # Do not display anything
    }
  })

  
}

shinyApp(ui = ui, server = server)
