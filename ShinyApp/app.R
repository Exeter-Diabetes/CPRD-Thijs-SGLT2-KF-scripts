library(shiny)
library(tidyverse)
library(shinyscroll)

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
      .toggle-pictogram {
      background: none;
      border: none;
      padding: 0;
      color: #007BFF;
      text-decoration: underline;
      font-style: italic;
      font-size: 16px;
      display: block;
      text-align: center;
      margin: 20px auto;
      }
          .toggle-traffic-light {
      background: none;
      border: none;
      padding: 0;
      color: #007BFF;
      text-decoration: underline;
      font-style: italic;
      font-size: 16px;
      display: block;
      text-align: center;
      margin: 20px auto;
    }
      .references-title {
        font-weight: bold;
      max-width: 700px; 
      min-width: 650px
      }
    "))
  ),
  column(width = 12, align = "center", div(class="title", titlePanel("Predicted kidney protection benefit with SGLT2-inhibitors in people with type 2 diabetes"))),
  column(width = 12, align = "center", div(class="beta-notice", 
                                           HTML("This calculator is in beta and is for research purposes only. Please leave any feedback <a href='https://forms.gle/ewYM9jmQfNb2fD2U8' target='_blank'>here</a>."))),
  column(width = 12, align = "center", div(class="explanation", uiOutput("dynamic_explanation"))),

  column(width = 12, align = "center", div(class="subtitle", h3(id = "result_title",
                                                                style = "font-weight: bold;",
                                                                       HTML("3-year risk of kidney disease progression<br>(≥50% decline in eGFR or kidney failure)")))),
  
  
    column(width = 12, align = "center",
         div(style = "max-width: 700px; min-width: 650px", 
             uiOutput("result_text"),  # Result displayed here
             
             # Wrap plotOutput in a div to allow styling
             div(plotOutput("risk_plot", height = "120px"), 
                 style = "margin-bottom: 20px;"),  # Add margin-bottom to the plot
             
             # Add space below the pictogram
             uiOutput("risk_pictogram", style = "margin-bottom: 10px;"),  # Add margin-bottom to the pictogram
             
             uiOutput("toggle_pictogram"),
             # actionButton("toggle_pictogram", "Click here to show the risk of possible side effects", class = "toggle-pictogram"),
             
             uiOutput("harm_pictogram", style = "margin_bottom: 10px;"),
             
             uiOutput("toggle_traffic_light"),
             
             # Traffic light with less added space
             uiOutput("traffic_light", style = "margin-bottom: 20px;")
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
                      selectInput("bmi", "BMI (kg/m²):", 
                                  choices = as.list(c(seq(20, 40, 1), "≥40 (model will use 40 as input)")), selected = 30),
                      selectInput("egfr", "eGFR (ml/min per 1.73m²):", 
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
             p("This prediction model is based on the CKD Prognosis Consortium risk score² for ≥50% decline in eGFR or kidney failure over 3 years integrated with the relative treatment effect from SGLT2-inhibitor trial meta-analysis¹. The model was independently validated using UK routine general practice data of 141,500 participants with type 2 diabetes, preserved eGFR (≥60mL/min/1.73m²), normal or low-level albuminuria (<30mg/mmol), and without a history of atherosclerotic vascular disease or heart failure³."),
             p(class = "references-title", "References:"),
             p("1. Nuffield Department of Population Health Renal Studies Group, SGLT2 inhibitor Meta-Analysis Cardio-Renal Trialists' Consortium. Impact of diabetes on the effects of sodium glucose co-transporter-2 inhibitors on kidney outcomes: collaborative meta-analysis of large placebo-controlled trials. Lancet 2022; 400(10365): 1788-801."),
             p("2. Grams ME, Brunskill NJ, Ballew SH, et al. Development and Validation of Prediction Models of Adverse Kidney Outcomes in the Population With and Without Diabetes. Diabetes Care 2022; 45(9): 2055-63."),
             p("3. Jansz TT, Young KG, Hopkins R, et al. Precision Medicine in Type 2 Diabetes: Targeting SGLT2-inhibitor Treatment for Kidney Protection. medRxiv 2024.09.01.24312905")
         )
  ),
  
  shinyscroll::use_shinyscroll(),
  
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
    HTML(paste0("This calculator predicts kidney protection benefit with SGLT2-inhibitor treatment for people with type 2 diabetes, preserved eGFR (≥60mL/min/1.73m²), normal or low-level albuminuria (uACR", acr_limit_value, "), and without a history of atherosclerotic vascular disease or heart failure."))
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
    
    # scroll to 
    shinyscroll::scroll("result_title")
    
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
          "This patient has an indication for SGLT2-inhibitor treatment due to severe albuminuria.")
    } else if (patient$egfr < 60) {
      div(class = "validation-message",
          "This patient has an indication for SGLT2-inhibitor treatment due to an eGFR <60mL/min/1.73m².")
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
        
        
        # result_colour <- ifelse(parr < parr_threshold, "#FF5733", "#00C853") 

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
                  span(style = "color: grey25; text-shadow: 1px 1px 2px black;",
                       sprintf("%.1f%%", current_risk))
              ),
              div(class = "result-text", style = "width: 33%;",
                  span(style = paste("color: #0072B2; text-shadow: 1px 1px 2px black;"),
                       sprintf("%.1f%%", risk_with_treatment))
              ),
              div(class = "result-text", style = "width: 33%;",
                  span(style = paste("color: #E69F00; text-shadow: 1px 1px 2px black;"),
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
    
    # bar_colour <- ifelse(parr < parr_threshold, "#FF5733", "#00C853")
    bar_colour <- "#0072B2"
    
    # Create a data frame for the two bars
    plot_data <- data.frame(
      xmin = c(0, 0),  # Add a small gap for yellow bar inset
      xmax = c(treated_risk, current_risk),
      ymin = c(0.2, 0.45),  # 
      ymax = c(0.4, 0.65),  # 
      fill = c(bar_colour, "grey25")  # Blue for untreated, yellow for treated
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
  
  output$risk_pictogram <- renderUI({
    patient <- patient()
    if (input$predict == 0 #|| !patient$type2_dm 
        || patient$hf || patient$chd || patient$egfr < 60 || patient$acr >= 30) {
      return(NULL)
    }
    
    if (is.null(result())) return()   
    
    # Extract the necessary values
    df <- result()
    current_risk <- df$ckdpc_50egfr_score[1]
    risk_with_treatment <- df$ckdpc_50egfr_score[1] - df$parr[1]
    
    total <- 100
    
    no_event <- round(total - current_risk) # white
    
    prevented <- round(max(current_risk - risk_with_treatment, 0))    # yellow
    still_event <- round(ifelse(prevented == 0, current_risk, risk_with_treatment))  # blue
    no_event <- total - still_event - prevented                       # white
    
    
    
    # Create the circle divs
    circles <- c(
      rep("<div class='circle white'></div>", no_event),
      rep("<div class='circle yellow'></div>", prevented),
      rep("<div class='circle blue'></div>", still_event)
      
    )
    
    HTML(sprintf("
<style>
  .circle-grid {
    display: grid;
    grid-template-columns: repeat(10, 20px);
    grid-gap: 4px;
    width: max-content;
    margin-top: 20px;
    margin-bottom: 20px;
  }
  .circle {
    width: 20px;
    height: 20px;
    border-radius: 50%%;
    border: 1px solid black;
  }
  .yellow {
    background-color: #E69F00;
  }
  .blue {
    background-color: #0072B2;
  }
  .white {
    background-color: white;
  }

  .legend-wrapper {
    display: flex;
    flex-direction: row;
    align-items: flex-start;
  }

  .legend {
    display: flex;
    flex-direction: column;
    gap: 10px;
    margin-top: 10px;
    margin-left: 40px;
    font-size: 14px;
    max-width: 400px;
  }

  .legend-item {
    display: flex;
    align-items: center;
    text-align: left;
    gap: 10px;
  }

  .legend-item .circle {
    flex-shrink: 0;
    display: inline-block;
  }
</style>

  <!-- Title for the circle diagram -->
  <h3 style='text-align: center;font-weight: bold;'>Making a decision</h3>
  <h4 style='text-align: center;font-style: italic;'>Predicted benefit</h3>

<div class='legend-wrapper'>
  <div class='circle-grid'>
    %s
  </div>
  <div class='legend'>
    <div><strong>If 100 people with this predicted risk take an SGLT2-inhibitor, over 3 years on average:</strong></div>
    <div class='legend-item'>
      <div class='circle white'></div>
      <div>about %d will not get kidney failure or a halving of their kidney function, but would not even if they had not taken an SGLT2-inhibitor</div>
    </div>
    <div class='legend-item'>
      <div class='circle yellow'></div>
      <div>about %d will not get kidney failure or a halving of their kidney function, because they take an SGLT2-inhibitor</div>
    </div>
    <div class='legend-item'>
      <div class='circle blue'></div>
      <div>about %d will get kidney failure or a halving of their kidney function, even though they take an SGLT2-inhibitor</div>
    </div>
  </div>
</div>
", paste(circles, collapse = "\n"), no_event, prevented, still_event))
    
    
    
  })
  
  output$toggle_pictogram <- renderUI({
    # Require result to be present and predict clicked
    if (is.null(result()) || input$predict == 0) return(NULL)
    
    # Require patient to meet conditions
    patient <- patient()
    if (patient$hf || patient$chd || patient$egfr < 60 || patient$acr >= 30) return(NULL)
    
    # If all conditions pass, render the button
    actionButton(
      inputId = "toggle_pictogram",
      label = "Possible side effects",
      class = "toggle-pictogram"
    )
  })
  
  
  
  
  
  # Create toggle state
  show_harm_pictogram <- reactiveVal("hide")
  
  # Toggle button logic
  observeEvent(input$toggle_pictogram, {
    if (show_harm_pictogram() == "hide") {
      show_harm_pictogram("show")
    } else {
      show_harm_pictogram("hide")
    }
  })
  
  # Update the button label based on state
  observe({
    # Only update label if the button is visible
    if (!is.null(result()) && input$predict != 0 && !(patient()$hf || patient()$chd || patient()$egfr < 60 || patient()$acr >= 30)) {
      label <- if (show_harm_pictogram() == "hide") {
        "Possible side effects"
      } else {
        "Possible side effects"
      }
      updateActionButton(session, "toggle_pictogram", label = label)
    }
  })
  
  # Render the harm pictogram
  output$harm_pictogram <- renderUI({
    if (show_harm_pictogram() == "hide") return(NULL)
    
    patient <- patient()
    if (input$predict == 0 || patient$hf || patient$chd || patient$egfr < 60 || patient$acr >= 30) return(NULL)
    if (is.null(result())) return(NULL)
    
    # Example values — replace with your actual logic
    current_risk <- 0.7
    risk_with_treatment <- 1.8
    total <- 100
    
    added <- round(max(risk_with_treatment - current_risk, 0))
    event_without_treatment <- round(current_risk)
    no_event <- total - event_without_treatment - added
    
    circles <- c(
      rep("<div class='circle white'></div>", no_event),
      rep("<div class='circle yellow'></div>", added),
      rep("<div class='circle blue'></div>", event_without_treatment)
    )
    
    HTML(sprintf("
  <style>
    .circle-grid {
      display: grid;
      grid-template-columns: repeat(10, 20px);
      grid-gap: 4px;
      width: max-content;
      margin: 20px auto;
    }
    .circle {
      width: 20px;
      height: 20px;
      border-radius: 50%%;
      border: 1px solid black;
    }
    .yellow { background-color: #E69F00; }
    .blue { background-color: #0072B2; }
    .white { background-color: white; }

    .legend-wrapper {
      display: flex;
      flex-direction: row;
      align-items: flex-start;
      justify-content: center;
      gap: 40px;
      max-width: 900px;
      margin: auto;
    }

    .legend {
      display: flex;
      flex-direction: column;
      gap: 10px;
      font-size: 14px;
      max-width: 400px;
    }

    .legend-item {
      display: flex;
      align-items: center;
      text-align: left;
      gap: 10px;
    }

    .legend-item .circle {
      flex-shrink: 0;
      display: inline-block;
    }

    .description {
      font-size: 14px;
      text-align: center;
      margin-bottom: 00px;
    }

    .footnote {
      font-size: 14px;
      text-align: center;
      margin-top: 00px;
    }
  </style>
  
  <div class='description'>
    SGLT2-inhibitors can cause genital thrush, but some people get this from time to time whether they take an SGLT2-inhibitor or not. The information below is a summary of results from many large studies¹.
  </div>

  <div class='legend-wrapper'>
    <div class='circle-grid'>
      %s
    </div>
    <div class='legend'>
      <div><strong>On average, for every 100 people who took an SGLT2-inhibitor over 3 years:</strong></div>
      <div class='legend-item'>
        <div class='circle white'></div>
        <div>about %d did not get genital thrush</div>
      </div>
      <div class='legend-item'>
        <div class='circle yellow'></div>
        <div>about %d got genital thrush, because they took an SGLT2-inhibitor</div>
      </div>
      <div class='legend-item'>
        <div class='circle blue'></div>
        <div>about %d got genital thrush, but would have done if they had not taken an SGLT2-inhibitor</div>
      </div>
    </div>
  </div>

  <div class='footnote'>
    More rarely, people can get diabetic keto-acidosis. This happens anyway in about 6 in 10,000 people (so 9,994 do not get this). If all 10,000 people took an SGLT2-inhibitor, an extra 9 might get diabetic keto-acidosis and 9,985 do not get diabetic keto-acidosis. Even more rarely, people can get a severe genital infection, called Fournier's gangrene. This might happen slightly more often in people taking an SGLT2-inhibitor, but this is extremely rare (about 1 in 10,000 people).
  </div>
  ", paste(circles, collapse = "\n"), no_event, added, event_without_treatment))
  })
  
  output$toggle_pictogram <- renderUI({
    # Require result to be present and predict clicked
    if (is.null(result()) || input$predict == 0) return(NULL)
    
    # Require patient to meet conditions
    patient <- patient()
    if (patient$hf || patient$chd || patient$egfr < 60 || patient$acr >= 30) return(NULL)
    
    # If all conditions pass, render the button
    actionButton(
      inputId = "toggle_pictogram",
      label = "Possible side effects",
      class = "toggle-pictogram"
    )
  })
  
  
  
  
  
  # Create toggle state
  show_traffic_light <- reactiveVal("hide")
  
  # Toggle button logic
  observeEvent(input$toggle_traffic_light, {
    if (show_traffic_light() == "hide") {
      show_traffic_light("show")
    } else {
      show_traffic_light("hide")
    }
  })
  
  # Update the button label based on state
  observe({
    # Only update label if the button is visible
    if (!is.null(result()) && input$predict != 0 && !(patient()$hf || patient()$chd || patient()$egfr < 60 || patient()$acr >= 30)) {
      label <- if (show_traffic_light() == "hide") {
        "Suggested treatment threshold"
      } else {
        "Suggested treatment threshold"
      }
      updateActionButton(session, "toggle_traffic_light", label = label)
    }
  })
  
  # Render the traffic light
  
  output$traffic_light <- renderUI({
    
    if (show_traffic_light() == "hide") return(NULL)
    
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
        style = "display: flex; flex-direction: column; align-items: center; margin-top: 20px; margin-bottom: 40px;",  # Align elements vertically and center them
        
        # # Title for the traffic light
        # div(
        #   style = "margin-bottom: 00px;",
        #   h4(style = "font-style: italic;", "Suggested treatment threshold")  # Title text
        # ),
        
        # Traffic Light and Text
        div(
          style = "display: flex; align-items: flex-start; width: 100%;",
          
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
              div(style = "display: flex; align-items: center; margin-bottom: 25px;",
                  div(style = "font-size: 20px; font-weight: bold; margin-right: 5px; color: black;", "→"),  # Unicode arrow
                  div(style = "font-size: 14px;", HTML("SGLT2-inhibitor treatment not suggested for kidney protection in this patient.<br>This patient's predicted benefit is below the suggested threshold (NNT of 154 or less), a threshold that would target the same number of patients as a ≥3mg/mmol albuminuria threshold (currently recommended in NICE, KDIGO, ADA, and EASD guidelines), but could prevent over 10% more kidney disease progression events over 3 years³. For context, the NNT seen in kidney outcome trials is 30¹."))
              )
            },
            
            # Spacer to align Green Light text
            div(style = "flex-grow: 1;"),
            
            # Text aligned with Green Light (if above threshold)
            if (!is_above_threshold) {
              div(style = "display: flex; align-items: center; margin-top: 25px;",
                  div(style = "font-size: 20px; font-weight: bold;margin-right: 5px; color: black;", "→"),  # Unicode arrow
                  div(style = "font-size: 14px;", HTML("SGLT2-inhibitor treatment suggested for kidney protection in this patient.<br>This patient's predicted benefit is above the suggested threshold (NNT of 154 or less), a threshold that would target the same number of patients as a ≥3mg/mmol albuminuria threshold (currently recommended in NICE, KDIGO, ADA, and EASD guidelines), but could prevent over 10% more kidney disease progression events over 3 years³. For context, the NNT seen in kidney outcome trials was about 30¹."))
              )
            }
          )
        )
      )
    } else {
      NULL  # Do not display anything
    }
  })
  
  
  output$toggle_traffic_light <- renderUI({
    # Require result to be present and predict clicked
    if (is.null(result()) || input$predict == 0) return(NULL)
    
    # Require patient to meet conditions
    patient <- patient()
    if (patient$hf || patient$chd || patient$egfr < 60 || patient$acr >= 30) return(NULL)
    
    # If all conditions pass, render the button
    actionButton(
      inputId = "toggle_traffic_light",
      label = "Suggested treatment threshold",
      class = "toggle-traffic-light"
    )
  })
  
}

shinyApp(ui = ui, server = server)
