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

load("2024-03-06_t2d_ckdpc_recalibrated_incl_egfr_below_60.Rda")
noncal_cohort$studydrug2 <- as.factor(noncal_cohort$studydrug2)

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

cohort <- cohort[!is.na(cohort$uacr) & !is.na(cohort$preegfr),]

cohort <- cohort %>% mutate(
  indication_NICE = ifelse(pre_index_date_ckd5_code == 1, "No renal indication", 
                           ifelse(macroalbuminuria == T,
                                  "Offer SGLT2 inhibitor (renal indication)", 
                                  ifelse(pre_index_date_heartfailure == 1 | pre_index_date_cvd == 1,
                                         "Offer SGLT2 inhibitor (cardiovascular indication)", 
                                         ifelse(
                                           microalbuminuria == T | (preegfr < 90 & preegfr >= 20),
                                           "Consider SGLT2 inhibitor / option",
                                           "No renal indication"
                                         )))),
  indication_NICE = ifelse(is.na(indication_NICE), "No renal indication", indication_NICE) %>% 
    factor(levels = c("No renal indication", "Consider SGLT2 inhibitor / option", "Offer SGLT2 inhibitor (renal indication)", "Offer SGLT2 inhibitor (cardiovascular indication)")),
  indication_KDIGO = ifelse(pre_index_date_ckd5_code == 1, "No renal indication", 
                            ifelse(macroalbuminuria == T,
                                   "Offer SGLT2 inhibitor (renal indication)", 
                                   ifelse(pre_index_date_heartfailure == 1 | pre_index_date_cvd == 1,
                                          "Offer SGLT2 inhibitor (cardiovascular indication)", 
                                          ifelse(
                                            (microalbuminuria == T | (preegfr < 60 & preegfr >= 20)),
                                            "Consider SGLT2 inhibitor / option",
                                            "No renal indication"
                                          )))),
  indication_KDIGO = ifelse(is.na(indication_KDIGO), "No renal indication", indication_KDIGO) %>% 
    factor(levels = c("No renal indication", "Consider SGLT2 inhibitor / option", "Offer SGLT2 inhibitor (renal indication)", "Offer SGLT2 inhibitor (cardiovascular indication)")),
  Population = paste0(c("Sample from CPRD (n=", nrow(cohort), ")"), collapse=""),
  Percentage = 0,
  SGLT2i = ifelse(is.na(pre_index_date_earliest_sglt2), F, T)
)

cohort$indication_NICE %>% summary()
cohort$indication_KDIGO %>% summary()

## as "new indication" we will set everyone with microalbuminuria, preserved eGFR, and high risk (>90 percentile) as consider
## and everyone with reduced eGFR and normoalbuminuria as consider
## everyone with reduced eGFR and microalbuminuria as offer
## and everyone with preserved eGFR otherwise as no renal indication

threshold_1 <- as.numeric(quantile(cohort[cohort$indication_NICE == "Consider SGLT2 inhibitor / option" & 
                                            cohort$egfr_below_60 == F & cohort$microalbuminuria == T,]$uacr, 
                                   0.90))

threshold_2 <- as.numeric(quantile(cohort[cohort$indication_NICE == "Consider SGLT2 inhibitor / option" & 
                                            cohort$egfr_below_60 == T & cohort$microalbuminuria == F,]$uacr, 
                                   0))

cohort <- cohort %>% mutate(
  newindication = ifelse(
    indication_NICE == "Consider SGLT2 inhibitor / option" & egfr_below_60 == T & microalbuminuria == T,
    "Offer SGLT2 inhibitor (renal indication)", ifelse(
      indication_NICE == "Consider SGLT2 inhibitor / option" & egfr_below_60 == F & microalbuminuria == F,
      "No renal indication", ifelse(
        indication_NICE == "Consider SGLT2 inhibitor / option" & egfr_below_60 == F & microalbuminuria == T, 
        ifelse(uacr < threshold_1, "No renal indication (as per risk score)", 
               "Consider SGLT2 inhibitor / option (as per risk score)"), ifelse(
            indication_NICE == "Consider SGLT2 inhibitor / option" & egfr_below_60 == T & microalbuminuria == F & 
              uacr < threshold_2, "No renal indication (as per risk score)", as.character(indication_NICE))
      )))
   %>% 
    factor(levels = c("No renal indication", "No renal indication (as per risk score)", "Consider SGLT2 inhibitor / option (as per risk score)", "Consider SGLT2 inhibitor / option", "Offer SGLT2 inhibitor (renal indication)", "Offer SGLT2 inhibitor (cardiovascular indication)"))
)

cohort$newindication %>% as.factor() %>% summary()

######################################2 PLOTS#########################################

# # simple bar charts for current and targeted indications:


# cohort <- cohort %>%
#   add_count(indication_NICE) %>%
#   ungroup() %>%
#   mutate(percent = round(n/nrow(cohort)*100))
# 
# ggplot(cohort, aes(fill=indication_NICE, x = Population, y = percent)) + 
#   geom_bar(position="fill", stat="identity") + 
#   # geom_text(aes(label=paste0(sprintf("%1.1f", percent),"%")),
#   #                                                      position=position_stack(vjust=0.5), colour="white", size = 2) +
#   scale_fill_manual(values = c("grey", "#E69F00", "#D55E00")) +
#   labs(x ="", y = "Percentage") +
#   scale_y_continuous(labels = scales::percent) + 
#   ggtitle("Prevalent patients with type 2 diabetes (01/07/2020)") +
#   theme_bw() +
#   theme(legend.position = "right")

# cohort <- cohort %>% 
#   add_count(newindication) %>%
#   ungroup() %>%
#   mutate(percent2 = round(nn/nrow(cohort)*100))
# 
# ggplot(cohort, aes(fill=newindication, x = Population, y = percent)) + 
#   geom_bar(position="fill", stat="identity") + 
#   # geom_text(aes(label=paste0(sprintf("%1.1f", percent),"%")),
#   #                                                      position=position_stack(vjust=0.5), colour="white", size = 2) +
#   scale_fill_manual(values = c("grey", "grey", "#E69F00", "#E69F00", "#D55E00")) +
#   labs(x ="", y = "Percentage") +
#   scale_y_continuous(labels = scales::percent) + 
#   ggtitle("Prevalent patients with type 2 diabetes (01/07/2020)") +
#   theme_bw() +
#   theme(legend.position = "right")

## sankey plot: NICE recommendations

# make long data frame
test <- cohort %>%
  make_long(indication_NICE, newindication) 

test <- test %>% mutate(
  node = factor(node, levels = c("Offer SGLT2 inhibitor (cardiovascular indication)", "Offer SGLT2 inhibitor (renal indication)", "Consider SGLT2 inhibitor / option", "Consider SGLT2 inhibitor / option (as per risk score)", "No renal indication (as per risk score)", "No renal indication")),
  next_node = factor(next_node, levels = c("Offer SGLT2 inhibitor (cardiovascular indication)", "Offer SGLT2 inhibitor (renal indication)", "Consider SGLT2 inhibitor / option", "Consider SGLT2 inhibitor / option (as per risk score)", "No renal indication (as per risk score)", "No renal indication")),
  x = ifelse(x == "indication_NICE", "NICE recommendations", "Targeted SGLT2i treatment")
) 

ggplot(test, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node),
               label = node)) +
  theme_sankey() +
  scale_fill_manual(values = c("#D55E00", "#D55E00", "#E69F00", "#E69F00", "grey", "grey")) + 
  ggtitle("Population effect of targeted SGLT2i treatment", 
          subtitle = paste0(c("Type 2 diabetes and non-missing eGFR/uACR in CPRD on 01/07/2020 (n=", nrow(cohort), ")"), collapse = "")) +
  guides(fill = guide_legend(title = "")) +
  geom_sankey(flow.alpha = 0.4, node.color = 1, flow.fill = "grey", flow.color = "grey85", width = .35, smooth = 10, space = 6000) +
  # geom_sankey_label(color = 1, fill = "white") +
  theme(axis.text=element_text(size=rel(.85)),
        plot.subtitle=element_text(size=rel(1))) +
  xlab("")
  
# percentages per indication_NICE:
for (k in levels(as.factor(test$x))) {
  for (i in levels(test$node)) {
  print(paste0(c(k, ": ", i, " ", round(nrow(test[test$x == k & test$node == i,])/nrow(test[test$x == k,])*100,2), "%"), collapse = ""))
  }
}

# percentage of people under current recommendations who are not on SGLT2i and should consider:
old <- round(nrow(cohort[cohort$indication_NICE == "Consider SGLT2 inhibitor / option" & cohort$SGLT2i==F,])/nrow(cohort)*100,2)

# percentage of people under targeted treatment who are not on SGLT2i and should consider:
new_1 <- round(nrow(cohort[cohort$newindication == "Consider SGLT2 inhibitor / option" & cohort$SGLT2i==F,])/nrow(cohort)*100,2)

# percentage of people under targeted treatment who are not on SGLT2i and should now be offered:
new_2 <- round(nrow(cohort[cohort$indication_NICE == "Consider SGLT2 inhibitor / option" & 
                    cohort$newindication == "Offer SGLT2 inhibitor (renal indication)" &
                    cohort$SGLT2i==F,])/nrow(cohort)*100,2)
# percentage that would additionally get treated for renal indications (previously in consider group, now in consider/offer)
new_1 + new_2
old


## sankey plot: KDIGO/EASD/ADA recommendations

# make long data frame
test2 <- cohort %>%
  make_long(indication_KDIGO, newindication) 

test2 <- test2 %>% mutate(
  node = factor(node, levels = c("Offer SGLT2 inhibitor (cardiovascular indication)", "Offer SGLT2 inhibitor (renal indication)", "Consider SGLT2 inhibitor / option", "Consider SGLT2 inhibitor / option (as per risk score)", "No renal indication (as per risk score)", "No renal indication")),
  next_node = factor(next_node, levels = c("Offer SGLT2 inhibitor (cardiovascular indication)", "Offer SGLT2 inhibitor (renal indication)", "Consider SGLT2 inhibitor / option", "Consider SGLT2 inhibitor / option (as per risk score)", "No renal indication (as per risk score)", "No renal indication")),
  x = ifelse(x == "indication_KDIGO", "KDIGO/EASD/ADA recommendations", "Targeted SGLT2i treatment")
) 

ggplot(test2, aes(x = x, 
                 next_x = next_x, 
                 node = node, 
                 next_node = next_node,
                 fill = factor(node),
                 label = node)) +
  theme_sankey() +
  scale_fill_manual(values = c("#D55E00", "#D55E00", "#E69F00", "#E69F00", "grey", "grey")) + 
  ggtitle("Population effect of targeted SGLT2i treatment", 
          subtitle = paste0(c("Type 2 diabetes and non-missing eGFR/uACR in CPRD on 01/07/2020 (n=", nrow(cohort), ")"), collapse = "")) +
  guides(fill = guide_legend(title = "")) +
  geom_sankey(flow.alpha = 0.4, node.color = 1, flow.fill = "grey", flow.color = "grey85", width = .35, smooth = 10, space = 6000) +
  # geom_sankey_label(color = 1, fill = "white") +
  theme(axis.text=element_text(size=rel(.85)),
        plot.subtitle=element_text(size=rel(1))) +
  xlab("")

# percentages per indication_KDIGO:
for (k in levels(as.factor(test2$x))) {
  for (i in levels(test2$node)) {
    print(paste0(c(k, ": ", i, " ", round(nrow(test2[test2$x == k & test2$node == i,])/nrow(test2[test2$x == k,])*100,2), "%"), collapse = ""))
  }
}

# percentage of people under current recommendations who are not on SGLT2i and should consider:
old <- round(nrow(cohort[cohort$indication_KDIGO == "Consider SGLT2 inhibitor / option" & cohort$SGLT2i==F,])/nrow(cohort)*100,2)

# percentage of people under targeted treatment who are not on SGLT2i and should consider:
new_1 <- round(nrow(cohort[cohort$newindication == "Consider SGLT2 inhibitor / option" & cohort$SGLT2i==F,])/nrow(cohort)*100,2)

# percentage of people under targeted treatment who are not on SGLT2i and should now be offered:
new_2 <- round(nrow(cohort[cohort$indication_KDIGO == "Consider SGLT2 inhibitor / option" & 
                             cohort$newindication == "Offer SGLT2 inhibitor (renal indication)" &
                             cohort$SGLT2i==F,])/nrow(cohort)*100,2)
# percentage that would additionally get treated for renal indications (previously in consider group, now in consider/offer)
new_1 + new_2
old

## sankey plot: current use vs NICE and KDIGO recommendations

cohort <- cohort %>% mutate(
  prev_users_NICE = ifelse(SGLT2i == T, ifelse(
    indication_NICE == "No renal indication", "No renal indication but currently on SGLT2i", as.character(indication_NICE)),
    "No renal indication"),
  pot_users_NICE = ifelse(prev_users_NICE == "No renal indication but currently on SGLT2i", "No renal indication but currently on SGLT2i",
                     as.character(indication_NICE)),
  prev_users_KDIGO = ifelse(SGLT2i == T, ifelse(
    indication_KDIGO == "No renal indication", "No renal indication but currently on SGLT2i", as.character(indication_KDIGO)),
    "No renal indication"),
  pot_users_KDIGO = ifelse(prev_users_KDIGO == "No renal indication but currently on SGLT2i", "No renal indication but currently on SGLT2i",
                          as.character(indication_KDIGO))
)

test3 <- cohort %>%
  make_long(prev_users_NICE, pot_users_NICE) 

test3 <- test3 %>% mutate(
  node = factor(node, levels = c("No renal indication but currently on SGLT2i", "Offer SGLT2 inhibitor (cardiovascular indication)", "Offer SGLT2 inhibitor (renal indication)", "Consider SGLT2 inhibitor / option", "No renal indication")),
  next_node = factor(next_node, levels = c("No renal indication but currently on SGLT2i", "Offer SGLT2 inhibitor (cardiovascular indication)", "Offer SGLT2 inhibitor (renal indication)", "Consider SGLT2 inhibitor / option", "No renal indication"))#,
  # x = ifelse(x == "pot_users", "Eligible users as per NICE recommendations", "SGLT2i users on 01/07/2020")
) 

ggplot(test3, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node),
                  label = node)) +
  theme_sankey() +
  scale_fill_manual(values = c("white", "#D55E00", "#D55E00", "#E69F00", "grey")) + 
  ggtitle("Population effect of NICE recommendations for SGLT2i treatment", 
          subtitle = paste0(c("Type 2 diabetes and non-missing eGFR/uACR in CPRD on 01/07/2020 (n=", nrow(cohort), ")"), collapse = "")) +
  guides(fill = guide_legend(title = "")) +
  geom_sankey(flow.alpha = 0.4, node.color = 1, flow.fill = "grey", flow.color = "grey75", width = .35, smooth = 10, space =6000) +
  # geom_sankey_label(color = 1, fill = "white") +
  theme(axis.text=element_text(size=rel(.85)),
        plot.subtitle=element_text(size=rel(1))) +
  xlab("")

for (k in levels(as.factor(test3$x))) {
  for (i in levels(as.factor(test3$node))) {
    print(paste0(c(k, ": ", i, " ", round(nrow(test3[test3$x == k & test3$node == i,])/nrow(test3[test3$x == k,])*100,2), "%"), collapse = ""))
  }
}



test4 <- cohort %>%
  make_long(prev_users_KDIGO, pot_users_KDIGO) 

test4 <- test4 %>% mutate(
  node = factor(node, levels = c("No renal indication but currently on SGLT2i", "Offer SGLT2 inhibitor (cardiovascular indication)", "Offer SGLT2 inhibitor (renal indication)", "Consider SGLT2 inhibitor / option", "No renal indication")),
  next_node = factor(next_node, levels = c("No renal indication but currently on SGLT2i", "Offer SGLT2 inhibitor (cardiovascular indication)", "Offer SGLT2 inhibitor (renal indication)", "Consider SGLT2 inhibitor / option", "No renal indication"))#,
) 

ggplot(test4, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node),
                  label = node)) +
  theme_sankey() +
  scale_fill_manual(values = c("white", "#D55E00", "#D55E00", "#E69F00", "grey")) + 
  ggtitle("Population effect of KDIGO/EASD/ADA recommendations for SGLT2i treatment", 
          subtitle = paste0(c("Type 2 diabetes and non-missing eGFR/uACR in CPRD on 01/07/2020 (n=", nrow(cohort), ")"), collapse = "")) +
  guides(fill = guide_legend(title = "")) +
  geom_sankey(flow.alpha = 0.4, node.color = 1, flow.fill = "grey", flow.color = "grey75", width = .35, smooth = 10, space =6000) +
  # geom_sankey_label(color = 1, fill = "white") +
  theme(axis.text=element_text(size=rel(.85)),
        plot.subtitle=element_text(size=rel(1))) +
  xlab("")

for (k in levels(as.factor(test4$x))) {
  for (i in levels(as.factor(test4$node))) {
    print(paste0(c(k, ": ", i, " ", round(nrow(test4[test4$x == k & test4$node == i,])/nrow(test4[test4$x == k,])*100,2), "%"), collapse = ""))
  }
}

######################################3 NNT IN NICE/KDIGO RECOMMENDATIONS#########################################

noncal_cohort <- noncal_cohort %>% filter(macroalbuminuria == F)

temp <- noncal_cohort %>% filter(preegfr < 90 | albuminuria == T)

surv_indication_NICE <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                data = temp,
                                weights = overlap2,
                                conf.type = "log", conf.int = 0.95)

NNT_NICE <- round(1/(summary(surv_indication_NICE, times=3)$surv[2]-summary(surv_indication_NICE, times=3)$surv[1])
                  , 0)

temp2 <- noncal_cohort %>% filter(preegfr < 60 | albuminuria == T)

surv_indication_KDIGO <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                 data = temp2,
                                 weights = overlap2,
                                 conf.type = "log", conf.int = 0.95)

NNT_KDIGO <- round(1/(summary(surv_indication_KDIGO, times=3)$surv[2]-summary(surv_indication_KDIGO, times=3)$surv[1])
                   , 0)


temp3 <- noncal_cohort %>% filter(preegfr>90)
surv_noindication_NICE <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                  data = temp3,
                                  weights = overlap2,
                                  conf.type = "log", conf.int = 0.95)

NNT_noindication_NICE <- round(1/(summary(surv_noindication_NICE, times=3)$surv[2]-summary(surv_noindication_NICE, times=3)$surv[1])
                               , 0)

temp4 <- noncal_cohort %>% filter(preegfr>60 & albuminuria == F)
surv_noindication_KDIGO <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                   data = temp4,
                                   weights = overlap2,
                                   conf.type = "log", conf.int = 0.95)

NNT_noindication_KDIGO <- round(1/(summary(surv_noindication_KDIGO, times=3)$surv[2]-summary(surv_noindication_KDIGO, times=3)$surv[1])
                                , 0)

temp5 <- noncal_cohort %>% filter(preegfr>60 & albuminuria == F |
                                    preegfr>60 & albuminuria == T & 
                                    uacr < threshold_1)
surv_noindication_tt <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
                                data = temp5,
                                weights = overlap2,
                                conf.type = "log", conf.int = 0.95)

NNT_noindication_tt <- round(1/(summary(surv_noindication_tt, times=3)$surv[2]-summary(surv_noindication_tt, times=3)$surv[1])
                             , 0)
