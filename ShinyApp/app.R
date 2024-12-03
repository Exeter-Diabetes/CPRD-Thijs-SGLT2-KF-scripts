



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
    }
    .subtitle {
      text-align: center;
      font-size: 18px;
      margin-top: 10px;
    }
      .result-text {
        font-size: 48px;
        font-weight: bold;
        text-align: center;
        margin-bottom: 20px;
      }
      .blue { color: #0072B2; text-shadow: 1px 1px 2px black; }
      .yellow { color: #E69F00; text-shadow: 1px 1px 2px black; }
      .orange { color: #D55E00; text-shadow: 1px 1px 2px black; }
      .validation-message {
        font-size: 16px;
        font-weight: normal;
        color: black;
        text-align: center;
        margin-top: 20px;
        margin-bottom: 10px;
      }
      .references-title {
        font-weight: bold;
      }
    "))
  ),
  div(class="title", titlePanel("Predicted kidney protection benefit from SGLT2-inhibitors")),
  div(class="subtitle", h4("3-year absolute risk reduction in kidney disease progression (≥50% decline in eGFR or kidney failure)")),
  uiOutput("result_text"),  # Result displayed here
  plotOutput("risk_plot", height = "120px"),  # Add the bar chart
  fluidRow(
    column(6,
           selectInput("age", "Age (years):", 
                       choices = as.list(seq(20, 80, 1)), selected = 50),
           selectInput("sex", "Sex:", 
                       choices = list(Male = "Male", Female = "Female")),
           selectInput("egfr", "eGFR (ml/min/1.73m2):", 
                       choices = as.list(seq(60, 140, 1)), selected = 90),
           selectInput("acr", "ACR (mg/mmol):", 
                       choices = as.list(seq(0.6, 29.9, 0.1)), selected = 1),
           selectInput("bmi", "BMI (kg/m2):", 
                       choices = as.list(seq(20, 40, 1)), selected = 30),
           selectInput("sbp", "SBP (mmHg):", 
                       choices = as.list(seq(90, 180, 1)), selected = 130),
           checkboxInput("bp_meds", "On BP Medications", value = FALSE)
    ),
    column(6,
           checkboxInput("type2_dm", "Type 2 Diabetes", value = TRUE),
           selectInput("hba1c", "HbA1c (mmol/mol):", 
                       choices = as.list(seq(42, 97, 1)), selected = 58),
           selectInput("diabetes_med", "Diabetes Medication:", 
                       choices = list("No Diabetes Medication" = "No Diabetes Medication", 
                                      "Oral Medications Only" = "Oral Medications Only", 
                                      Insulin = "Insulin")),
           checkboxInput("hf", "History of Heart Failure (HF)", value = FALSE),
           checkboxInput("chd", "History of Coronary Heart Disease (CHD)", value = FALSE),
           checkboxInput("af", "History of Atrial Fibrillation (AF)", value = FALSE),
           selectInput("smoking_status", "Smoking Status:", 
                       choices = list("Never Smoker" = "Never Smoker", 
                                      "Current Smoker" = "Current Smoker", 
                                      "Ex Smoker" = "Ex Smoker"))
    )
  ),
  div(style = "text-align:center; margin-top:20px;", 
      actionButton("predict", "Predict", class = "btn-primary btn-lg")
  ),
  div(style = "text-align:center; margin-top:50px;",
      p("This prediction model integrates the relative risk reduction estimate from SGLT2-inhibitor trial meta-analysis with an existing prediction model for risk of a composite of ≥50% decline in eGFR or kidney failure over 3 years. It was validated using UK routine general practice data of 120,315 participants with type 2 diabetes, preserved eGFR (≥60mL/min/1.73m2), normal or low-level albuminuria (<30mg/mmol), and without a history of atherosclerotic vascular disease or heart failure."),
      p(class = "references-title", "References:"),
      p("Nuffield Department of Population Health Renal Studies Group, SGLT2 inhibitor Meta-Analysis Cardio-Renal Trialists' Consortium. Impact of diabetes on the effects of sodium glucose co-transporter-2 inhibitors on kidney outcomes: collaborative meta-analysis of large placebo-controlled trials. Lancet 2022; 400(10365): 1788-801."),
      p("Grams ME, Brunskill NJ, Ballew SH, et al. Development and Validation of Prediction Models of Adverse Kidney Outcomes in the Population With and Without Diabetes. Diabetes Care 2022; 45(9): 2055-63.")
  )
)

# Server logic for the Shiny app
server <- function(input, output, session) {
  # Ensure decimal separator uses a period
  Sys.setlocale("LC_NUMERIC", "C")
  
  result <- eventReactive(input$predict, {
    if (input$hf || input$chd || !input$type2_dm) {
      return(NULL)
    }
    calculate_parr(
      dataframe = data.frame(
        age = input$age,
        sex = input$sex,
        egfr = input$egfr,
        acr = input$acr,
        sbp = input$sbp,
        bp_meds = as.numeric(input$bp_meds),
        hf = as.numeric(input$hf),
        chd = as.numeric(input$chd),
        af = as.numeric(input$af),
        smoking_status = input$smoking_status,
        diabetes_med = input$diabetes_med,
        bmi = input$bmi,
        hba1c = input$hba1c
      ),
      age = input$age,
      sex = input$sex,
      egfr = input$egfr,
      acr = input$acr,
      sbp = input$sbp,
      bp_meds = as.numeric(input$bp_meds),
      hf = as.numeric(input$hf),
      chd = as.numeric(input$chd),
      af = as.numeric(input$af),
      smoking_status = input$smoking_status,
      diabetes_med = input$diabetes_med,
      bmi = input$bmi,
      hba1c = input$hba1c
    )
  })
  
  output$result_text <- renderUI({
    if (input$predict == 0) {
      NULL
    } else if (!input$type2_dm) {
      div(class = "validation-message", 
          "This prediction model is validated for individuals with type 2 diabetes only.")
    } else if (input$hf || input$chd) {
      div(class = "validation-message", 
          "This prediction model is validated for individuals without a history of atherosclerotic vascular disease or heart failure.")
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
        
        div(style = "display: flex; justify-content: space-between; text-align: center;",
            div(class = "result-text", style = "width: 33%;", 
                h4("Current risk"), 
                span(style = "color: #0072B2; text-shadow: 1px 1px 2px black;", 
                     sprintf("%.2f%%", current_risk))),
            # div(class = "result-text", style = "width: 33%;", 
            #     h4("Residual risk with treatment"), 
            #     span(style = "color: #0072B2; text-shadow: 1px 1px 2px black;", 
            #          sprintf("%.2f%%", risk_with_treatment))),
            div(class = "result-text", style = "width: 33%;", 
                h4("Reduction with treatment"), 
                span(style = "color: #E69F00; text-shadow: 1px 1px 2px black;", 
                     sprintf("%.2f%%", parr))),
            div(class = "result-text", style = "width: 33%;", 
                h4("3-year NNT"), 
                span(style = "color: #E69F00; text-shadow: 1px 1px 2px black;", 
                     ifelse(!is.na(nnt), sprintf("%d", nnt), "N/A")))
        )
      }
    }
  })
  
  output$risk_plot <- renderPlot({
    if (input$predict == 0 || !input$type2_dm || input$hf || input$chd) {
      return(NULL)
    }
    
    if (is.null(result())) return()
    
    # Extract the necessary values
    df <- result()
    current_risk <- df$ckdpc_50egfr_score[1]
    treated_risk <- 100 * (1 - df$ckdpc_50egfr_survival_sglt2i[1])
    
    # Create a data frame for the two bars
    plot_data <- data.frame(
      xmin = c(0, treated_risk),  # Add a small gap for yellow bar inset
      xmax = c(current_risk, current_risk),
      ymin = c(0.4, 0.4),  # Blue bar slightly taller
      ymax = c(0.6, 0.6),  # Yellow bar slightly shorter
      fill = c("#0072B2", "#E69F00")  # Blue for untreated, yellow for treated
    )
    
    # Generate the plot
    ggplot() +
      # Full-width blue bar (untreated risk)
      geom_rect(data = plot_data[1, ], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
                color = "#0072B2", size = 1.2) +  # No border for blue bar
      
      # Narrower yellow bar (treated risk) with blue border
      geom_rect(data = plot_data[2, ], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
                color = "#0072B2", size = 1.2) +  # Blue border
      
      # Use specified colors for the bars
      scale_fill_identity() +
      
      # X-axis configuration
      scale_x_continuous(limits = c(0, max(current_risk, 10)), 
                         breaks = seq(0, 10, 1), 
                         labels = paste0(seq(0, 10, 1), "%")) +
      
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
      ) +
      
      # Chart title
      labs(title = "3-year risk of kidney disease progression")
  })
  
  
  # output$risk_plot <- renderPlot({
  #   if (input$predict == 0 || !input$type2_dm || input$hf || input$chd) {
  #     return(NULL)
  #   }
  #   
  #   if (is.null(result())) return()
  #   
  #   # Extract the necessary values
  #   df <- result()
  #   current_risk <- df$ckdpc_50egfr_score[1]
  #   treated_risk <- 100 * (1 - df$ckdpc_50egfr_survival_sglt2i[1])
  #   
  #   # Create a data frame for plotting
  #   plot_data <- data.frame(
  #     xmin = c(0, treated_risk),
  #     xmax = c(treated_risk, current_risk),
  #     ymin = c(0.4, 0.4),  # Bar thickness
  #     ymax = c(0.6, 0.6),  # Bar thickness
  #     fill = c("#0072B2", "#E69F00")  # Blue for treated, yellow for untreated
  #   )
  #   
  #   # Generate the plot
  #   ggplot() +
  #     geom_rect(data = plot_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
  #               color = NA) +  # Remove borders
  #     scale_fill_identity() +  # Use specified colors
  #     scale_x_continuous(limits = c(0, max(current_risk, 10)), breaks = seq(0, 10, 1), labels = paste0(seq(0, 10, 1), "%")) +
  #     theme_minimal() +
  #     theme(
  #       axis.title.y = element_blank(),
  #       axis.text.y = element_blank(),
  #       axis.ticks.y = element_blank(),
  #       axis.title.x = element_blank(),
  #       axis.text.x = element_text(size = 14),
  #       plot.title = element_text(size = 16, hjust = 0.5),
  #       panel.grid.major.y = element_blank(),
  #       panel.grid.minor.y = element_blank(),
  #       panel.grid.major.x = element_line(color = "grey80"),
  #       panel.grid.minor.x = element_blank(),
  #       plot.margin = margin(10, 10, 10, 10) 
  #     ) +
  #     labs(title = "Current 3-year risk of kidney disease progression")
  # })
}

shinyApp(ui = ui, server = server)
