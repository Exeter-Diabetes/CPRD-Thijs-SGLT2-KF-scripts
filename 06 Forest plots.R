## Forest plots:
# 1 per drug class
# 2 by analysis strategy
# 3 by albuminuria status
# 4 secondary outcomes by albuminuria status
# 5 secondary outcomes by pARR threshold

############################0 SETUP################################################################
# 0 Setup
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
load(paste0(today, "_all_hrs.Rda"))
load(paste0(today, "_all_SGLT2ivsDPP4iSU_hrs.Rda"))
load(paste0(today, "_SGLT2ivsDPP4iSU_hrs.Rda"))
load(paste0(today, "_subgroup_hrs.Rda"))
load(paste0(today, "_subgroup_SGLT2ivsDPP4iSU_hrs.Rda"))
load(paste0(today, "_subgroup_parr_hrs.Rda"))
load(paste0(today, "_subgroup_parr_SGLT2ivsDPP4iSU_hrs.Rda"))

############################1 FOREST PLOT FOR HR BY DRUG CLASS (SUPPLEMENTAL FIGURE)################################################################

# create labels
labels <- data.frame(matrix("", nrow = 1, ncol = length(all_hrs)))
names(labels) <- names(all_hrs)
labels <- labels %>% mutate(contrast = "Overall (by analytic approach)", 
                            analysis = "Overall (by analytic approach)",
                            string = "Hazard Ratio (95% CI)",
                            SU_nN = "Events/subjects (SU)",
                            DPP4i_nN = "Events/subjects (DPP4i)",
                            SGLT2i_nN = "Events/subjects (SGLT2i)")

labels_plot <- all_hrs %>% mutate(
  string = ifelse(outcome == "ckd_egfr50_pp" & contrast == "SGLT2i vs SU" & analysis == "IPTW", paste(string, " "), string),
  string = as.factor(string)
)

for (k in unique(all_hrs$contrast)) {
  for (m in unique(all_hrs$outcome)) {
    labels_temp <- labels %>% data.frame()
    labels_temp$outcome <- m
    labels_temp$contrast <- k
    labels_plot <- rbind(labels_temp, labels_plot)
  }
}

labels_plot$analysis <- factor(labels_plot$analysis, levels = rev(unique(labels_plot$analysis)))



# plot
p_hr_1 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.3, 0.5, 0.75, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp",]$analysis)) + 1), 
                  xlim=c(0.3, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.25,
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp",]$analysis)) + 1, 
           label = "Favours SU") +
  labs(x="", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_1 <-
  labels_plot %>%
  mutate(analysis = ifelse(analysis == "Overall (by analytic approach)",
                           "SGLT2-inhibitors vs sulfonylureas",
                           as.character(analysis)),
         analysis = factor(analysis, levels = c(
           "IPTW",
           "Overlap-weighted",
           "Adjusted",           
           "SGLT2-inhibitors vs sulfonylureas"
         ))) %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = (analysis))) + 
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs SU",]$
                        analysis == "Overall (by analytic approach)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_1 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50" & labels_plot$contrast == "SGLT2i vs SU",]$
                        string == "Hazard Ratio (95% CI)", "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_counts_1 <- labels_plot %>% filter(outcome=="ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs SU") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
            fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs SU",]$SGLT2i_nN == labels$SGLT2i_nN, "bold", "plain")) +
  geom_text(aes(x = 4, label = `SU_nN`), hjust = 1, fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs SU",]$SU_nN == labels$SU_nN, "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(-2, 5))

## same for SGLT2i vs DPP4

p_hr_2 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.3, 0.5, 0.75, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1), 
                  xlim=c(0.3, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.25,
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours DPP4i") +
  labs(x="", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 

p_left_2 <-
  labels_plot %>%
  mutate(analysis = ifelse(analysis == "Overall (by analytic approach)",
                           "SGLT2-inhibitors vs DPP4-inhibitors",
                           as.character(analysis)),
         analysis = factor(analysis, levels = c(
           "IPTW",
           "Overlap-weighted",
           "Adjusted",           
           "SGLT2-inhibitors vs DPP4-inhibitors"
         ))) %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = (analysis))) + 
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50" & labels_plot$contrast == "SGLT2i vs DPP4i",]$
                        analysis == "Overall (by analytic approach)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_2 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    colour = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50" & labels_plot$contrast == "SGLT2i vs DPP4i",]$
                      string == "Hazard Ratio (95% CI)", "white", "black")) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_counts_2 <- labels_plot %>% filter(outcome=="ckd_egfr50_pp") %>%
  filter(contrast == "SGLT2i vs DPP4i") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
            fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs DPP4i",]$SGLT2i_nN == labels$SGLT2i_nN, "bold", "plain")) +
  geom_text(aes(x = 4, label = `DPP4i_nN`), hjust = 1, fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "SGLT2i vs DPP4i",]$DPP4i_nN == labels$DPP4i_nN, "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(-2, 5))

## same for DPP4i vs SU

p_hr_3 <- 
  all_hrs %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.3, 0.5, 0.75, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1), 
                  xlim=c(0.3, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours DPP4i") +
  annotate("text", x = 1.25,
           y = length(unique(all_hrs[all_hrs$outcome == "ckd_egfr50_pp" & !all_hrs$analysis == "Meta-analysis of RCTs",]$analysis)) + 1, 
           label = "Favours SU") +
  labs(x="Hazard Ratio", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 

p_left_3 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  mutate(analysis = ifelse(analysis == "Overall (by analytic approach)",
                           "DPP4-inhibitors vs sulfonylureas",
                           as.character(analysis)),
         analysis = factor(analysis, levels = c(
           "IPTW",
           "Overlap-weighted",
           "Adjusted",           
           "DPP4-inhibitors vs sulfonylureas"
         ))) %>%
  ggplot(aes(y = (analysis))) + 
  # geom_text(aes(x = 0, label = analysis), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "DPP4i vs SU",]$
                        analysis == "Overall (by analytic approach)", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_3 <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    colour = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50" & labels_plot$contrast == "DPP4i vs SU",]$
                      string == "Hazard Ratio (95% CI)", "white", "black")) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_counts_3 <- labels_plot %>% filter(outcome=="ckd_egfr50_pp") %>%
  filter(contrast == "DPP4i vs SU") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(aes(x = 1, label = DPP4i_nN), hjust = 1, 
            fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "DPP4i vs SU",]$DPP4i_nN == labels$DPP4i_nN, "bold", "plain")) +
  geom_text(aes(x = 4, label = `SU_nN`), hjust = 1, fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50_pp" & labels_plot$contrast == "DPP4i vs SU",]$SU_nN == labels$SU_nN, "bold", "plain")) +
  theme_void() +
  coord_cartesian(xlim = c(-2, 5))

layout <- c(
  area(t = 0, l = 0, b = 5, r = 4), 
  area(t = 0, l = 5, b = 5, r = 9), 
  area(t = 0, l = 10, b = 5, r = 14),
  area(t = 0, l = 15, b = 5, r = 19),
  area(t = 6, l = 0, b = 11, r = 4), 
  area(t = 6, l = 5, b = 11, r = 9),
  area(t = 6, l = 10, b = 11, r = 14),
  area(t = 6, l = 15, b = 11, r = 19),
  area(t = 12, l = 0, b = 17, r = 4), 
  area(t = 12, l = 5, b = 17, r = 9),
  area(t = 12, l = 10, b = 17, r = 14),
  area(t = 12, l = 15, b = 17, r = 19)
)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_HR_by_drugclass.tiff"), width=18, height=5.5, units = "in", res=800) 
p_left_1 + p_counts_1 + p_hr_1 + p_right_1 + 
  p_left_2 + p_counts_2 + p_hr_2 + p_right_2 + 
  p_left_3 + p_counts_3 + p_hr_3 + p_right_3 + plot_layout(design = layout)
dev.off()
############################2 FOREST PLOT OF OVERALL HR BY ANALYSIS STRATEGY (SUPPLEMENTAL FIGURE)################################################################


# prep all_hrs dataframe for forest plot to show hazard ratios and add literature-reported HR
trial_hr <- data.frame(matrix("", nrow=1, ncol=length(all_SGLT2ivsDPP4iSU_hrs))) 
names(trial_hr) <- names(all_SGLT2ivsDPP4iSU_hrs)
trial_hr <- trial_hr %>% mutate(
  outcome = "ckd_egfr50", contrast = "SGLT2i vs SU", analysis = "Nuffield Group, 2022", 
  HR = 0.62, LB = 0.56, UB = 0.68, string = paste0(HR, " (", LB, ", ", UB, ")"), 
  model = paste0(string, " [", analysis, "]"))

labels1 <- data.frame(matrix("", nrow = 1, ncol = length(trial_hr)))
names(labels1) <- names(trial_hr)
labels1 <- labels1 %>% mutate(contrast = "Meta-analysis of RCTs", 
                              analysis = "Meta-analysis of RCTs",
                              string = "Hazard Ratio (95% CI)")
labels_plot <- trial_hr

for (k in unique(trial_hr$contrast)) {
  for (m in unique(trial_hr$outcome)) {
    labels_temp <- labels1
    labels_temp$outcome <- m
    labels_temp$contrast <- k
    labels_plot <- rbind(labels_temp, labels_plot)
  }
}

# create labels for forest plot with studydrug2 (SGLT2i vs DPP4i/SU combined)
labels_plot2 <- all_SGLT2ivsDPP4iSU_hrs

labels2 <- data.frame(matrix("", nrow = 1, ncol = length(all_SGLT2ivsDPP4iSU_hrs)))
names(labels2) <- names(all_SGLT2ivsDPP4iSU_hrs)
labels2 <- labels2 %>% mutate(contrast = "Study cohort (by analysis approach)", 
                              analysis = "Study cohort (by analysis approach)",
                              string = "Hazard Ratio (95% CI)",
)



for (k in unique(all_SGLT2ivsDPP4iSU_hrs$contrast)) {
  for (m in unique(all_SGLT2ivsDPP4iSU_hrs$outcome)) {
    labels_temp <- labels2 
    labels_temp$outcome <- m
    labels_temp$contrast <- k
    labels_plot2 <- rbind(labels_temp, labels_plot2)
  }
}

labels_plot2$analysis <- factor(labels_plot2$analysis, levels = rev(unique(labels_plot2$analysis)))


class(trial_hr$HR) <- class(trial_hr$LB) <- class(trial_hr$UB) <-
  class(all_SGLT2ivsDPP4iSU_hrs$HR) <- class(all_SGLT2ivsDPP4iSU_hrs$LB) <- class(all_SGLT2ivsDPP4iSU_hrs$UB) <- "numeric"

# plot
p_hr_all <- 
  all_SGLT2ivsDPP4iSU_hrs %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.45, 0.6, 0.80, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$outcome == "ckd_egfr50",]$analysis)) + 1), 
                  xlim=c(0.45, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  labs(x="Hazard Ratio", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5)) 


p_left_all <-
  labels_plot2 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = (analysis))) + 
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot2[labels_plot2$outcome == "ckd_egfr50",]$
                        analysis == labels2$analysis, "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_all <-
  labels_plot2 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    fontface = "plain",
    colour = ifelse(labels_plot2[labels_plot2$outcome == "ckd_egfr50",]$string == labels2$string, "white", "black")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

# trial hr
p_hr_trial <- 
  trial_hr %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(analysis, levels = rev(unique(analysis))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.45, 0.60, 0.80, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(trial_hr[trial_hr$outcome == "ckd_egfr50",]$analysis)) + 1), 
                  xlim=c(0.45, 1.5)) +
  theme_classic() +
  geom_point(aes(x=HR), shape=15, size=3, colour = "#0072B2") +
  geom_linerange(aes(xmin=LB, xmax=UB), colour = "#0072B2") +
  geom_vline(xintercept = 1, linetype="dashed") +
  annotate("text", x = .65, 
           y = length(unique(trial_hr[trial_hr$outcome == "ckd_egfr50",]$analysis)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.25,
           y = length(unique(trial_hr[trial_hr$outcome == "ckd_egfr50",]$analysis)) + 1, 
           label = "Favours\nDPP4i/SU") +
  labs(x="", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 


p_left_trial <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = rev(analysis))) + 
  geom_text(
    aes(x = 1, label = analysis),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50",]$
                        analysis == labels1$analysis, "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_trial <-
  labels_plot %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  # geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    fontface = ifelse(labels_plot[labels_plot$outcome == "ckd_egfr50",]$
                        string == labels1$string, "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

# layout for plots below
layout <- c(
  area(t = 21, l = 0, b = 55, r = 4), 
  area(t = 21, l = 5, b = 55, r = 9),
  area(t = 21, l = 10, b = 55, r = 14),
  area(t = 0, l = 0, b = 23, r = 4), 
  area(t = 0, l = 5, b = 23, r = 9), 
  area(t = 0, l = 10, b = 23, r = 14)
  
)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_HR_by_analysis_strategy.tiff"), width=15, height=4, units = "in", res=800) 
# Final plot arrangement
p_left_all + p_hr_all + p_right_all + 
  p_left_trial + p_hr_trial + p_right_trial + plot_layout(design = layout)
dev.off()

############################3 FOREST PLOT OF HR BY ALBUMINURIA STATUS (FIGURE 1A)################################################################

# if analyses in 02 Calculate weights and hazard ratios.R not run, then object p_value_interaction will not be defined
# if that is the case, define p-value at 0.99 (as per analyses dd 17/12/2024)
if (!exists("p_value_interaction")) {
  p_value_interaction <- 0.99
}

# prep data frames with row for overall
overall <- all_SGLT2ivsDPP4iSU_hrs[all_SGLT2ivsDPP4iSU_hrs$analysis == "Overlap-weighted",]
overall$contrast <- " "

labels_plot5 <- overall

labels5 <- data.frame(matrix("", nrow = 1, ncol = length(overall)))
names(labels5) <- names(overall)
labels5 <- labels5 %>% mutate(contrast = " All subjects", 
                              analysis = " ", 
                              string = "Hazard Ratio (95% CI)",
                              SGLT2i_nN = "Events/subjects (SGLT2i)",
                              `DPP4i/SU_nN` = "(DPP4i/SU)")

for (m in unique(overall$outcome)) {
  labels_temp <- labels5 
  labels_temp$outcome <- m
  labels_plot5 <- rbind(labels_temp, labels_plot5)
}

# have to coerce HR and CI to class numeric as they sometimes default to character
class(subgroup_hrs$HR) <- class(subgroup_hrs$LB) <- class(subgroup_hrs$UB) <- "numeric"

labels3 <- data.frame(matrix("", nrow = 1, ncol = length(subgroup_hrs)))
names(labels3) <- names(subgroup_hrs)
labels3 <- labels3 %>% mutate(contrast = paste0("By albuminuria status (p = ", sprintf("%.2f", round(p_value_interaction,2)), " for trend)"), 
                              analysis = "By albuminuria status", 
                              string = "Hazard Ratio (95% CI)",
                              SGLT2i_nN = "Events/subjects",
                              `DPP4i/SU_nN` = "")

labels_plot3 <- subgroup_hrs %>% filter(analysis == "Overlap-weighted") %>%
  mutate(contrast = case_when(
  contrast == "uACR <3mg/mmol" ~ "Albuminuria <3mg/mmol",
  contrast == "uACR 3-30mg/mmol" ~ "Albuminuria 3-30mg/mmol"))

for (k in unique(subgroup_hrs$outcome)) {
  labels_temp <- labels3
  labels_temp$outcome <- k
  labels_plot3 <- rbind(labels_temp, labels_plot3)
}

labels_plot3$contrast <- factor(labels_plot3$contrast, levels = c(labels3$contrast,
                                                                  "Albuminuria <3mg/mmol", 
                                                                  "Albuminuria 3-30mg/mmol"))

# plot by risk group
p_counts_subgroup <- labels_plot3 %>% filter(outcome=="ckd_egfr50") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
            colour = ifelse(labels_plot3[labels_plot3$outcome == "ckd_egfr50",]$SGLT2i_nN == labels3$SGLT2i_nN, "white", "black")) +
  geom_text(aes(x = 3, label = `DPP4i/SU_nN`), hjust = 1, fontface = "plain") +
  theme_void(base_size = 14) +
  coord_cartesian(xlim = c(-2, 5))

p_hr_subgroup <- 
  subgroup_hrs %>% filter(analysis == "Overlap-weighted") %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0)) +
  coord_cartesian(ylim=c(1,length(unique(subgroup_hrs[subgroup_hrs$outcome == "ckd_egfr50",]$contrast)) + 1), 
                  xlim=c(0.25, 2)) +
  theme_classic(base_size = 14) +
  geom_point(aes(x=HR), shape=15, size=3) +
  geom_linerange(aes(xmin=LB, xmax=UB)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  geom_vline(xintercept = 0.62, linetype="twodash", size = 1, colour = "#D55E00") +
  geom_vline(xintercept = 0.56, linetype="twodash", size = 0.5, colour = "#D55E00") +
  geom_vline(xintercept = 0.68, linetype="twodash", size = 0.5, colour = "#D55E00") +
  annotate("text", x = .65, 
           y = length(unique(subgroup_hrs[subgroup_hrs$outcome == "ckd_egfr50",]$contrast)) + 1, 
           label = "") +
  annotate("text", x = 1.25,
           y = length(unique(subgroup_hrs[subgroup_hrs$outcome == "ckd_egfr50",]$contrast)) + 1, 
           label = "") +
  labs(x="Hazard ratio", y="") + 
  geom_point(aes(x = ifelse(UB > 2, 2.2, NA),
                 y=1.01), 
             shape = 62, size = 5, color = "black") + # add symbol to indicate CI extends beyond axis limits
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_subgroup <-
  labels_plot3 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = rev(unique((contrast))))) + 
  geom_text(
    aes(x = 1, label = contrast),
    hjust = 0,
    fontface = ifelse(labels_plot3[labels_plot3$outcome == "ckd_egfr50",]$
                        contrast == labels3$contrast, "bold", "plain")
  ) +
  theme_void(base_size = 14) +
  coord_cartesian(xlim = c(0, 4), clip = "off")

p_right_subgroup <-
  labels_plot3 %>%
  filter(outcome == "ckd_egfr50") %>%
  mutate(string = 
           ifelse(coalesce(string == lead(string), FALSE), paste0(string, " "), string)) %>% # if hazard ratios are the same, add a space to the end of one of them so they do not get taken as one
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    colour = ifelse(labels_plot3[labels_plot3$outcome == "ckd_egfr50",]$string == "Hazard Ratio (95% CI)", "white", "black")) +
  theme_void(base_size = 14) +
  coord_cartesian(xlim = c(0, 4))


# plot for overall HR
p_counts_overall <- labels_plot5 %>% filter(outcome=="ckd_egfr50") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
            fontface = ifelse(labels_plot5[labels_plot5$outcome == "ckd_egfr50",]$SGLT2i_nN == labels5$SGLT2i_nN, "bold", "plain")) +
  geom_text(aes(x = 3, label = `DPP4i/SU_nN`), hjust = 1, 
            fontface = ifelse(labels_plot5[labels_plot5$outcome == "ckd_egfr50",]$`DPP4i/SU_nN` == labels5$`DPP4i/SU_nN`, "bold", "plain")) +
  theme_void(base_size = 14) +
  coord_cartesian(xlim = c(-2, 5), clip="off")

p_hr_overall <- 
  overall %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
  scale_x_continuous(trans = "log10", breaks = c(0.25, 0.5, 0.75, 1.0, 1.5)) +
  coord_cartesian(ylim=c(1,length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast)) + 1), 
                  xlim=c(0.25, 2)) +
  theme_classic(base_size = 14) +
  geom_point(aes(x=HR), shape=15, size=3, colour = "black") +
  geom_linerange(aes(xmin=LB, xmax=UB), colour = "black") +
  geom_vline(xintercept = 1, linetype="dashed") +
  geom_segment(aes(x = 0.62, xend = 0.62, y = 0, yend = length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast))+0.75), 
                                                             linetype = "twodash", size = 1, colour = "#D55E00") +
  geom_segment(aes(x = 0.56, xend = 0.56, y = 0, yend = length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast))+0.75), 
               linetype = "twodash", size = 0.5, colour = "#D55E00") +
  geom_segment(aes(x = 0.68, xend = 0.68, y = 0, yend = length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast))+0.75), 
               linetype = "twodash", size = 0.5, colour = "#D55E00") +
  annotate("text", x = .6, 
           y = length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast)) + 1, 
           label = "Favours SGLT2i") +
  annotate("text", x = 1.5,
           y = length(unique(overall[overall$outcome == "ckd_egfr50",]$contrast)) + 1, 
           label = "Favours\nDPP4i/SU") +
  labs(x="Hazard Ratio", y="") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) 


p_left_overall <-
  labels_plot5 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = rev(unique((contrast))))) + 
  geom_text(
    aes(x = 1, label = contrast),
    hjust = 0,
    fontface = ifelse(labels_plot5[labels_plot5$outcome == "ckd_egfr50",]$
                        contrast == overall$contrast[1], "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_right_overall <-
  labels_plot5 %>%
  filter(outcome == "ckd_egfr50") %>%
  ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
  geom_text(
    aes(x = 0, label = string),
    hjust = 0,
    fontface = ifelse(labels_plot5[labels_plot5$outcome == "ckd_egfr50",]$string == "Hazard Ratio (95% CI)", "bold", "plain")) +
  theme_void(base_size = 14) +
  coord_cartesian(xlim = c(0, 4))


# # layout for plots below
# layout <- c(
#   area(t = 13, l = 7, b = 30, r = 13),
#   area(t = 13, l = 0, b = 30, r = 7), 
#   area(t = 13, l = 12, b = 30, r = 18),
#   area(t = 13, l = 19, b = 30, r = 24),
#   area(t = 0, l = 7, b = 13, r = 13), 
#   area(t = 0, l = 0, b = 13, r = 7),
#   area(t = 0, l = 12, b = 13, r = 18),
#   area(t = 0, l = 19, b = 13, r = 24)
#   
# )

layout <- c(
  area(t = 8, l = 2, b = 18, r = 5), 
  area(t = 8, l = 0, b = 18, r = 2), 
  area(t = 8, l = 5.5, b = 18, r = 6),
  area(t = 8, l = 7, b = 18, r = 8),
  area(t = 0, l = 2, b = 8, r = 5), 
  area(t = 0, l = 0, b = 8, r = 2),
  area(t = 0, l = 5.5, b = 8, r = 6),
  area(t = 0, l = 7, b = 8, r = 8)
)


# Final plot arrangement
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_HR_by_subgroup.tiff"), width=11, height=3.5, units = "in", res=800) 
p_counts_subgroup + p_left_subgroup + p_hr_subgroup + p_right_subgroup + 
  p_counts_overall + p_left_overall + p_hr_overall + p_right_overall + plot_layout(design = layout)
dev.off()

############################4 FOREST PLOT FOR HR OF SECONDARY OUTCOMES BY ALBUMINURIA STATUS (SUPPLEMENTAL FIGURE)################################################################

secondary <- subgroup_hrs %>% 
  filter(analysis=="Overlap-weighted") %>%
  filter(!outcome %in% c("ckd_egfr50", "death")) %>%
  mutate(albuminuria_status = case_when(
    contrast == "uACR <3mg/mmol" ~ "Albuminuria <3mg/mmol",
    contrast == "uACR 3-30mg/mmol" ~ "Albuminuria 3-30mg/mmol"
  ))

# Prepare labels for each outcome by albuminuria status
labels_plot6 <- secondary

labels6 <- data.frame(matrix("", nrow = 1, ncol = length(labels_plot6)))
names(labels6) <- names(labels_plot6)
labels6 <- labels6 %>% mutate(analysis = "Overall",
                              string = "Hazard Ratio (95% CI)",
                              SGLT2i_nN = "Events/subjects (SGLT2i)",
                              `DPP4i/SU_nN` = "(DPP4i/SU)")
 
for (m in unique(secondary$outcome)) {
  labels_temp <- labels6
  labels_temp$outcome <- m
  labels_plot6 <- rbind(labels_temp, labels_plot6)
}

labels_plot6 <- labels_plot6 %>%
  mutate(
    outcome_label = case_when(
      outcome == "ckd_egfr40" ~ "40% eGFR decline / ESKD",
      outcome == "ckd_egfr50_5y" ~ "50% eGFR decline / ESKD (5 years)",
      outcome == "macroalb" ~ "Progression to albuminuria ≥30mg/mmol",
      outcome == "dka" ~ "Diabetic keto-acidosis",
      outcome == "amputation" ~ "Amputation",
      outcome == "side_effect" ~ "Mycotic genital infection"
    )
  )

labels_plot6 <- labels_plot6 %>%
  mutate(
    contrast = ifelse(analysis == "Overall", paste0(" ", outcome_label), albuminuria_status)
    )


# have to coerce HR and CI to class numeric as they sometimes default to character
class(secondary$HR) <- class(secondary$LB) <- class(secondary$UB) <- "numeric"

plot_expression <- ""

for (m in rev(unique(secondary$outcome))) {
  
  p_counts <- labels_plot6 %>% filter(outcome==m) %>%
    ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
    geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
              colour = ifelse(!m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$SGLT2i_nN == labels6$SGLT2i_nN,
                              "white", "black"),
              fontface = ifelse(m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$SGLT2i_nN == labels6$SGLT2i_nN,
                                "bold", "plain")) +
    geom_text(aes(x = 3, label = `DPP4i/SU_nN`), hjust = 1, 
              colour = ifelse(!m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$`DPP4i/SU_nN` == labels6$`DPP4i/SU_nN`,
                              "white", "black"),
              fontface = ifelse(m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$`DPP4i/SU_nN` == labels6$`DPP4i/SU_nN`,
                                "bold", "plain")) +
    theme_void() +
    coord_cartesian(xlim = c(-2, 5))
  
  p_hr <- 
    secondary %>%
    filter(outcome == m) %>%
    ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
    scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.5, 2.25, 3.25)) +
    coord_cartesian(ylim=c(1,length(unique(labels_plot6[labels_plot6$outcome == m,]$contrast)) + 1), 
                    xlim=c(0.5, 3.25)) +
    theme_classic() +
    geom_point(aes(x=HR), shape=15, size=3) +
    geom_linerange(aes(xmin=LB, xmax=UB)) +
    geom_vline(xintercept = 1, linetype="dashed") +
    annotate("text", x = .65, 
             y = length(unique(secondary[secondary$outcome == m,]$contrast)) + 2, 
             label = ifelse(m==unique(secondary$outcome)[1], "Favours SGLT2i", "")) +
    annotate("text", x = 1.5,
             y = length(unique(secondary[secondary$outcome == m,]$contrast)) + 2, 
             label = ifelse(m==unique(secondary$outcome)[1], "Favours DPP4i/SU", "")) +
    labs(x=ifelse(m==unique(secondary$outcome)[nlevels(as.factor(secondary$outcome))], "Hazard ratio", ""), y="") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y= element_blank(),
          axis.title.y= element_blank(),
          axis.line.x = if (!m==unique(secondary$outcome)[nlevels(as.factor(secondary$outcome))]) {element_blank()},
          axis.text.x = if (!m==unique(secondary$outcome)[nlevels(as.factor(secondary$outcome))]) {element_blank()},
          axis.ticks.x = if (!m==unique(secondary$outcome)[nlevels(as.factor(secondary$outcome))]) {element_blank()},
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) 
  
  p_left <-
    labels_plot6 %>%
    filter(outcome == m) %>%
    ggplot(aes(y = rev(unique((contrast))))) + 
    geom_text(
      aes(x = 1, label = contrast),
      hjust = 0,
      fontface = ifelse(labels_plot6[labels_plot6$outcome == m,]$
                          contrast == paste0(" ", labels_plot6[labels_plot6$outcome == m,]$outcome_label), "bold", "plain"),
      colour = ifelse(m==unique(secondary$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$
                        contrast == "Outcome" | !labels_plot6[labels_plot6$outcome == m,]$
                        contrast == "Outcome", "black", "white")
    ) +
    theme_void() +
    coord_cartesian(xlim = c(0, 4))
  
  p_right <-
    labels_plot6 %>%
    filter(outcome == m) %>%
    ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
    geom_text(
      aes(x = 0, label = string),
      hjust = 0,
      fontface = ifelse(labels_plot6[labels_plot6$outcome == m,]$string == "Hazard Ratio (95% CI)", "bold", "plain"),
      colour = ifelse(labels_plot6[labels_plot6$outcome == m,]$string == "Hazard Ratio (95% CI)" & !m==unique(secondary$outcome)[1], 
                      "white", "black")) +
    theme_void() +
    coord_cartesian(xlim = c(0, 4))
  
  assign(paste0("p_counts_", m), p_counts)
  assign(paste0("p_hr_", m), p_hr)
  assign(paste0("p_left_", m), p_left)
  assign(paste0("p_right_", m), p_right)
  
  plot_expression <- paste0(plot_expression, "p_counts_", m, " + p_left_", m, " + p_hr_", m, " + p_right_", m, " + ")
}

n.plots <- nlevels(as.factor(secondary$outcome))

# layout for plots below

i <- 1
layout <- paste("area(t = ",(i-1)*24, ", l = 7, b = ",(i-1)*24+24,", r = 13), area(t = ",(i-1)*24, ", l = 0, b = ",(i-1)*24+24,", r = 7), area(t = ",(i-1)*24, ", l = 12, b = ",(i-1)*24+24,", r = 18), area(t = ",(i-1)*24, ", l = 19, b = ", (i-1)*24+24,", r = 24)")


for (i in 2:n.plots) {
  layout <- paste("area(t = ",(i-1)*24, ", l = 7, b = ",(i-1)*24+24,", r = 13), area(t = ",(i-1)*24, ", l = 0, b = ",(i-1)*24+24,", r = 7), area(t = ",(i-1)*24, ", l = 12, b = ",(i-1)*24+24,", r = 18), area(t = ",(i-1)*24, ", l = 19, b = ", (i-1)*24+24,", r = 24), ", layout)
}

layout <- paste0("c(", layout, ")")

layout <- eval(str2lang(layout))

# Final plot arrangement

plot_expression <- paste0(plot_expression, "plot_layout(design = layout)")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_HR_secondary_outcomes_by_albuminuria.tiff"), width=18, height=5.5, units = "in", res=800) 
eval(str2lang(plot_expression))
dev.off()

############################5 FOREST PLOT FOR HRs OF SECONDARY OUTCOMES BY pARR (SUPPLEMENTAL FIGURE)################################################################

secondary2 <- subgroup_parr_hrs %>% 
  filter(analysis=="Overlap-weighted") %>%
  filter(!outcome %in% c("ckd_egfr50", "death")) %>%
  mutate(parr_status = case_when(
    contrast == "pARR below threshold" ~ "pARR <0.65%",
    contrast == "pARR above threshold" ~ "pARR ≥0.65%"
  ))

# Prepare labels for each outcome by pARR status
labels_plot6 <- secondary2

labels6 <- data.frame(matrix("", nrow = 1, ncol = length(labels_plot6)))
names(labels6) <- names(labels_plot6)
labels6 <- labels6 %>% mutate(analysis = "Overall",
                              string = "Hazard Ratio (95% CI)",
                              SGLT2i_nN = "Events/subjects (SGLT2i)",
                              `DPP4i/SU_nN` = "(DPP4i/SU)")

for (m in unique(secondary2$outcome)) {
  labels_temp <- labels6
  labels_temp$outcome <- m
  labels_plot6 <- rbind(labels_temp, labels_plot6)
}

labels_plot6 <- labels_plot6 %>%
  mutate(
    outcome_label = case_when(
      outcome == "ckd_egfr40" ~ "40% eGFR decline / ESKD",
      outcome == "ckd_egfr50_5y" ~ "50% eGFR decline / ESKD (5 years)",
      outcome == "macroalb" ~ "Progression to albuminuria ≥30mg/mmol",
      outcome == "dka" ~ "Diabetic keto-acidosis",
      outcome == "amputation" ~ "Amputation",
      outcome == "side_effect" ~ "Mycotic genital infection"
    )
  )

labels_plot6 <- labels_plot6 %>%
  mutate(
    contrast = ifelse(analysis == "Overall", paste0(" ", outcome_label), parr_status),
    contrast = factor(contrast, levels = c(" 40% eGFR decline / ESKD",
                                           " 50% eGFR decline / ESKD (5 years)",
                                           " Progression to albuminuria ≥30mg/mmol",
                                           " Diabetic keto-acidosis",
                                           " Amputation",
                                           " Mycotic genital infection",
                                           "pARR <0.65%",
                                           "pARR ≥0.65%"))
  )


# have to coerce HR and CI to class numeric as they sometimes default to character
class(secondary2$HR) <- class(secondary2$LB) <- class(secondary2$UB) <- "numeric"

plot_expression <- ""

for (m in rev(unique(secondary2$outcome))) {
  
  p_counts <- labels_plot6 %>% filter(outcome==m) %>%
    ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
    geom_text(aes(x = 1, label = SGLT2i_nN), hjust = 1, 
              colour = ifelse(!m==unique(secondary2$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$SGLT2i_nN == labels6$SGLT2i_nN,
                              "white", "black"),
              fontface = ifelse(m==unique(secondary2$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$SGLT2i_nN == labels6$SGLT2i_nN,
                                "bold", "plain")) +
    geom_text(aes(x = 3, label = `DPP4i/SU_nN`), hjust = 1, 
              colour = ifelse(!m==unique(secondary2$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$`DPP4i/SU_nN` == labels6$`DPP4i/SU_nN`,
                              "white", "black"),
              fontface = ifelse(m==unique(secondary2$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$`DPP4i/SU_nN` == labels6$`DPP4i/SU_nN`,
                                "bold", "plain")) +
    theme_void() +
    coord_cartesian(xlim = c(-2, 5))
  
  p_hr <- 
    secondary2 %>%
    filter(outcome == m) %>%
    ggplot(aes(y = factor(contrast, levels = rev(unique(contrast))))) + 
    scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.5, 2.25, 3.25)) +
    coord_cartesian(ylim=c(1,length(unique(labels_plot6[labels_plot6$outcome == m,]$contrast)) + 1), 
                    xlim=c(0.5, 3.25)) +
    theme_classic() +
    geom_point(aes(x=HR), shape=15, size=3) +
    geom_linerange(aes(xmin=LB, xmax=UB)) +
    geom_vline(xintercept = 1, linetype="dashed") +
    annotate("text", x = .65, 
             y = length(unique(secondary2[secondary2$outcome == m,]$contrast)) + 2, 
             label = ifelse(m==unique(secondary2$outcome)[1], "Favours SGLT2i", "")) +
    annotate("text", x = 1.5,
             y = length(unique(secondary2[secondary2$outcome == m,]$contrast)) + 2, 
             label = ifelse(m==unique(secondary2$outcome)[1], "Favours DPP4i/SU", "")) +
    labs(x=ifelse(m==unique(secondary2$outcome)[nlevels(as.factor(secondary2$outcome))], "Hazard ratio", ""), y="") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y= element_blank(),
          axis.title.y= element_blank(),
          axis.line.x = if (!m==unique(secondary2$outcome)[nlevels(as.factor(secondary2$outcome))]) {element_blank()},
          axis.text.x = if (!m==unique(secondary2$outcome)[nlevels(as.factor(secondary2$outcome))]) {element_blank()},
          axis.ticks.x = if (!m==unique(secondary2$outcome)[nlevels(as.factor(secondary2$outcome))]) {element_blank()},
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) 
  
  p_left <-
    labels_plot6 %>%
    filter(outcome == m) %>%
    ggplot(aes(y = rev(contrast))) + 
    geom_text(
      aes(x = 1, label = contrast),
      hjust = 0,
      fontface = ifelse(labels_plot6[labels_plot6$outcome == m,]$
                          contrast == paste0(" ", labels_plot6[labels_plot6$outcome == m,]$outcome_label), "bold", "plain"),
      colour = ifelse(m==unique(secondary2$outcome)[1] & labels_plot6[labels_plot6$outcome == m,]$
                        contrast == "Outcome" | !labels_plot6[labels_plot6$outcome == m,]$
                        contrast == "Outcome", "black", "white")
    ) +
    theme_void() +
    coord_cartesian(xlim = c(0, 4))
  
  p_right <-
    labels_plot6 %>%
    filter(outcome == m) %>%
    ggplot(aes(y = factor(string, levels = rev(unique(string))))) + 
    geom_text(
      aes(x = 0, label = string),
      hjust = 0,
      fontface = ifelse(labels_plot6[labels_plot6$outcome == m,]$string == "Hazard Ratio (95% CI)", "bold", "plain"),
      colour = ifelse(labels_plot6[labels_plot6$outcome == m,]$string == "Hazard Ratio (95% CI)" & !m==unique(secondary2$outcome)[1], 
                      "white", "black")) +
    theme_void() +
    coord_cartesian(xlim = c(0, 4))
  
  assign(paste0("p_counts_", m), p_counts)
  assign(paste0("p_hr_", m), p_hr)
  assign(paste0("p_left_", m), p_left)
  assign(paste0("p_right_", m), p_right)
  
  plot_expression <- paste0(plot_expression, "p_counts_", m, " + p_left_", m, " + p_hr_", m, " + p_right_", m, " + ")
}

n.plots <- nlevels(as.factor(secondary2$outcome))

# layout for plots below

i <- 1
layout <- paste("area(t = ",(i-1)*24, ", l = 7, b = ",(i-1)*24+24,", r = 13), area(t = ",(i-1)*24, ", l = 0, b = ",(i-1)*24+24,", r = 7), area(t = ",(i-1)*24, ", l = 12, b = ",(i-1)*24+24,", r = 18), area(t = ",(i-1)*24, ", l = 19, b = ", (i-1)*24+24,", r = 24)")


for (i in 2:n.plots) {
  layout <- paste("area(t = ",(i-1)*24, ", l = 7, b = ",(i-1)*24+24,", r = 13), area(t = ",(i-1)*24, ", l = 0, b = ",(i-1)*24+24,", r = 7), area(t = ",(i-1)*24, ", l = 12, b = ",(i-1)*24+24,", r = 18), area(t = ",(i-1)*24, ", l = 19, b = ", (i-1)*24+24,", r = 24), ", layout)
}

layout <- paste0("c(", layout, ")")

layout <- eval(str2lang(layout))

# Final plot arrangement

plot_expression <- paste0(plot_expression, "plot_layout(design = layout)")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
tiff(paste0(today, "_HR_secondary_outcomes_by_parr.tiff"), width=18, height=5.5, units = "in", res=800) 
eval(str2lang(plot_expression))
dev.off()

