######################################0 SET UP#########################################
# install.packages("remotes")                         # to install ggsankey package
# remotes::install_github("davidsjoberg/ggsankey")    # to install ggsankey package
library(tidyverse)
library(ggsankey)
library(survival)
library(survminer)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/2023-12-05_t2d_csection_2020-07.Rda")
cohort <- prev_2020
rm(prev_2020)

######################################1 DEFINE COHORT AND INDICATIONS#########################################

cohort <- cohort %>%
  mutate(pre_index_date_cvd=ifelse(pre_index_date_angina==1 | pre_index_date_ihd==1 | pre_index_date_myocardialinfarction==1 | pre_index_date_pad==1 | pre_index_date_revasc==1 | pre_index_date_stroke==1 | pre_index_date_tia==1, 1, 0),
         uacr = preacr) 

cohort <- cohort %>%
  mutate(
    macroalbuminuria = ifelse(uacr < 30, F, T),
    microalbuminuria = ifelse(uacr <3, F, ifelse(macroalbuminuria == T, F, T)),
    egfr_below_60=ifelse( # preckdstage %in% c("stage_1", "stage_2") & 
      preegfr >=60, F, T)
  )

cohort$uacr %>% is.na() %>% summary()
# only 60% have a uACR

exclude <- cohort[is.na(cohort$uacr) | is.na(cohort$preegfr),]

cohort <- cohort[!is.na(cohort$uacr) & !is.na(cohort$preegfr),]

cohort <- cohort %>% mutate(
  indication_NICE = ifelse(pre_index_date_ckd5_code == 1, "No renal indication",
                           ifelse(pre_index_date_heartfailure == 1 | pre_index_date_cvd == 1,
                                  "Offer SGLT2-inhibitor (cardiovascular indication)",
                                  ifelse(macroalbuminuria == T & preegfr >=20,
                                         "Offer SGLT2-inhibitor (renal indication)",

                                         ifelse(
                                           microalbuminuria == T & preegfr >=20 | (preegfr < 90 & preegfr >= 20),
                                           "Consider SGLT2-inhibitor / option",
                                           "No renal indication"
                                         )))),
  indication_NICE = ifelse(is.na(indication_NICE), "No renal indication", indication_NICE) %>%
    factor(levels = c("No renal indication", "Consider SGLT2-inhibitor / option", "Offer SGLT2-inhibitor (renal indication)", "Offer SGLT2-inhibitor (cardiovascular indication)")),
  indication_KDIGO = ifelse(pre_index_date_ckd5_code == 1, "No renal indication", 
                            ifelse(pre_index_date_heartfailure == 1 | pre_index_date_cvd == 1,
                                   "Offer SGLT2-inhibitor (cardiovascular indication)", 
                                   ifelse(macroalbuminuria == T & preegfr >= 20,
                                          "Offer SGLT2-inhibitor (renal indication)", 
                                          ifelse(
                                            (microalbuminuria == T & preegfr >=20 | (preegfr < 60 & preegfr >= 20)),
                                            "Offer SGLT2-inhibitor (renal indication)",
                                            "No renal indication"
                                          )))),
  indication_KDIGO = ifelse(is.na(indication_KDIGO), "No renal indication", indication_KDIGO) %>% 
    factor(levels = c("No renal indication", "Offer SGLT2-inhibitor (renal indication)", "Offer SGLT2-inhibitor (cardiovascular indication)")),
  Population = paste0(c("Sample from CPRD (n=", nrow(cohort), ")"), collapse=""),
  Percentage = 0,
  SGLT2i = ifelse(is.na(pre_index_date_earliest_sglt2), F, T)
)

cohort$indication_NICE %>% summary()
cohort$indication_KDIGO %>% summary()

exclude2 <- cohort %>% filter(!indication_KDIGO == "Offer SGLT2-inhibitor (renal indication)")
cohort <- cohort %>% filter(indication_KDIGO == "Offer SGLT2-inhibitor (renal indication)")

## as "new indication" we will set everyone with microalbuminuria, preserved eGFR, and high risk (>90 percentile) as consider
## and everyone with reduced eGFR and normoalbuminuria as consider
## everyone with reduced eGFR and microalbuminuria as offer
## and everyone with preserved eGFR otherwise as no renal indication

threshold_1 <- as.numeric(quantile(cohort[cohort$egfr_below_60 == F & cohort$microalbuminuria == T,]$uacr, 
                                   0.90))

cohort <- cohort %>% mutate(
  newindication = ifelse(macroalbuminuria == T, 
                      ifelse(egfr_below_60 == F, "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol", "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"), 
                      ifelse(egfr_below_60 == T, 
                             ifelse(microalbuminuria == T, "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol", "eGFR <60mL/min/1.73m2, uACR <3mg/mmol"),
                             ifelse(microalbuminuria == F, "eGFR ≥60mL/min/1.73m2, uACR <3mg/mmol", "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol"))),
  newindication = ifelse(newindication == "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol", 
                         ifelse(uacr < threshold_1, 
                                paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score 0-90th percentile"),
                                paste0("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score ≥90th percentile")), 
                         as.character(newindication)),
  newindication = newindication %>%
    factor(levels = c("eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score 0-90th percentile",
                      "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score ≥90th percentile", 
                      "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                      "eGFR <60mL/min/1.73m2, uACR <3mg/mmol", 
                      "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol", 
                      "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol"))
)

cohort$newindication %>% summary()

######################################2 PLOTS#########################################

## sankey plot: KDIGO/EASD/ADA recommendations

# make long data frame
test2 <- cohort %>%
  make_long(indication_KDIGO, newindication) 

test2 <- test2 %>% mutate(
  node = ifelse(node == "Offer SGLT2-inhibitor (renal indication)", "Treatment recommended\nby KDIGO/ADA/EASD guidelines", as.character(node)),
  node = factor(node, levels = c("Treatment recommended\nby KDIGO/ADA/EASD guidelines", 
                                 "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score 0-90th percentile",
                                 "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score ≥90th percentile", 
                                 "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                 "eGFR <60mL/min/1.73m2, uACR <3mg/mmol", 
                                 "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol", 
                                 "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol")),
  next_node = factor(next_node, levels = c("Treatment recommended\nby KDIGO/ADA/EASD guidelines", 
                                           "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score 0-90th percentile",
                                           "eGFR ≥60mL/min/1.73m2, uACR 3-30mg/mmol, CKD-PC risk score ≥90th percentile", 
                                           "eGFR ≥60mL/min/1.73m2, uACR ≥30mg/mmol",
                                           "eGFR <60mL/min/1.73m2, uACR <3mg/mmol", 
                                           "eGFR <60mL/min/1.73m2, uACR 3-30mg/mmol", 
                                           "eGFR <60mL/min/1.73m2, uACR ≥30mg/mmol")),
  x = ifelse(x == "indication_KDIGO", "KDIGO/EASD/ADA recommendations", "Targeted SGLT2-inhibitor treatment")
) 

ggplot(test2, aes(x = x, 
                  next_x = next_x, 
                  node = fct_rev(node), 
                  next_node = fct_rev(next_node),
                  fill = node,
                  label = fct_rev(node))) +
  theme_sankey() +
  scale_fill_manual(values = c("#E69F00", "grey", "#D55E00", "#D55E00", "#D55E00", "#D55E00", "#D55E00")) + 
  ggtitle("People with type 2 diabetes that KDIGO/ADA/EASD guidelines recommend treating for kidney protection",
          subtitle = paste0(c("UK-representative cohort from CPRD with non-missing eGFR/uACR on 01/07/2020"), collapse = "")) +
  guides(fill = guide_legend(title="")) +
  geom_sankey(flow.alpha = 0.4, node.color = 1, flow.fill = "grey", flow.color = "grey85", width = .25, smooth = 4, space = 2800) +
  # geom_sankey_label(aes(group = next_node, label = fct_rev(node)), hjust = 0.5, vjust = 0.5, color = 1, fill = "white", size=rel(3)) +
  theme(plot.subtitle=element_text(size=rel(1))) +
  xlab("")

# percentages per indication_KDIGO:
for (k in levels(as.factor(test2$x))) {
  for (i in levels(test2$node)) {
    print(paste0(c(k, ": ", i, " ", round(nrow(test2[test2$x == k & test2$node == i,])/nrow(test2[test2$x == k,])*100,2), "%"), collapse = ""))
  }
}
