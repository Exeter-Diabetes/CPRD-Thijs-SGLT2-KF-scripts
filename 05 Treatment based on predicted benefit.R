
############################0 SETUP################################################################

# Setup
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-11-06_t2d_ckdpc_recalibrated_with_adjsurv.Rda")


############################1 DEFINE CUTOFFS################################################################

## consider following methods to define cut-off:
# 1 Strategy A: choose cut-off that matches treatment proportion of guidelines but improves outcomes through better-targeted treatment
# 2 Strategy B: take as cut-off the top decile  (10%) of predicted benefit

cutoff1 <- cohort %>% .$ckdpc_50egfr_sglt2i_benefit %>% quantile( 
  1 - (cohort %>% filter(albuminuria == T) %>% nrow() / cohort  %>% nrow())) %>% as.numeric()

cutoff1_equivalent_40egfr_score <- cohort  %>% .$ckdpc_40egfr_score %>% quantile( 
  1 - (cohort %>% filter(albuminuria == T) %>% nrow() / cohort  %>% nrow())) %>% as.numeric()

cutoff2 <- cohort %>% .$ckdpc_50egfr_sglt2i_benefit %>% quantile(0.90) %>% as.numeric()

cutoff2_equivalent_40egfr_score <- cohort %>% .$ckdpc_40egfr_score %>% quantile(0.90) %>% as.numeric()

cohort <- cohort %>% mutate(
  treat_guideline = ifelse(albuminuria == F, F, T),                          # everyone current guidelines recommend treating (uACR cutoff)
  treat_model1 = ifelse(ckdpc_50egfr_sglt2i_benefit > cutoff1, T, F),        # strategy A (matching uACR treatment proportion)
  treat_model2 = ifelse(ckdpc_50egfr_sglt2i_benefit > cutoff2, T, F),        # strategy B (top 10%)
)

############################2 EVALUATE OBSERVED BENEFIT DISTRIBUTION/CALIBRATION################################################################

## FIGURE 2A: histogram of predicted benefit
benefit_histogram <- ggplot(cohort %>%
                      mutate(predicted_benefit_percent = ckdpc_50egfr_sglt2i_benefit * 100), 
                    aes(x = predicted_benefit_percent)) +
  geom_histogram(aes(y = ..count.. / 10, fill = predicted_benefit_percent > cutoff1*100),
                 binwidth = 0.02, color = "black") +  # Adjust binwidth as necessary
  scale_fill_manual(values = c("TRUE" = "#E69F00", "FALSE" = "grey")) +
  geom_vline(xintercept = cutoff1*100, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = cutoff2*100, linetype = "dashed", color = "black", size = 1) +
  # geom_vline(xintercept = cutoff3*100, linetype = "dashed", color = "black", size = 1) +

  annotate("text", x = cutoff1*100, y = Inf, label = "pARR strategy A", vjust = -0.5, hjust = 1.1, angle = 90, size = 4, color = "black") +
  annotate("text", x = cutoff2*100, y = Inf, label = "pARR strategy B", vjust = -0.5, hjust = 1.1, angle = 90, size = 4, color = "black") +

  labs(x = "Predicted absolute risk reduction with SGLT2i (%)", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "none") + 
  coord_cartesian(xlim=c(0,2.5))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_predicted_benefit_histogram.tiff"), width=8, height=4.5, units = "in", res=800) 
benefit_histogram
dev.off()

##FIGURE 2B: calibration plot of predicted vs observed benefit
obs_v_pred_for_plot <- cohort %>%
  mutate(benefit_decile = ntile(ckdpc_50egfr_sglt2i_benefit, n.quantiles)) %>%
  group_by(benefit_decile) %>%
  summarise(median_predicted_benefit=median(ckdpc_50egfr_sglt2i_benefit, na.rm=T),
            mean_predicted_benefit=mean(ckdpc_50egfr_sglt2i_benefit, na.rm=T),
            mean_benefit=mean(survdiff_ckd_egfr50),
            se_benefit=mean(se_survdiff_ckd_egfr50),
            median_benefit=median(survdiff_ckd_egfr50),
            lq_benefit=quantile(survdiff_ckd_egfr50, prob=c(.25)),
            uq_benefit=quantile(survdiff_ckd_egfr50, prob=c(.75)),
            upper_ci=mean_benefit + 1.96*se_benefit,
            lower_ci=mean_benefit - 1.96*se_benefit)


empty_tick <- data.frame(matrix(NA, nrow = 1, ncol = length(obs_v_pred_for_plot)))
names(empty_tick) <- names(obs_v_pred_for_plot)
empty_tick <- empty_tick %>%
  mutate(benefit_decile=0)

## SGLT2i benefit predicted vs observed - mean
p_benefit_bydeciles_mean <- ggplot(data=bind_rows(empty_tick,obs_v_pred_for_plot), aes(x=mean_predicted_benefit*100)) +
  geom_errorbar(aes(ymax=upper_ci*100,ymin=lower_ci*100, color= "#E69F00"),width=0.1,size=1) +
  geom_point(aes(y = mean_benefit*100, color="#E69F00"), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Predicted SGLT2-inhibitor benefit (%)") + ylab("Observed benefit (%)")+
  scale_colour_manual(values = "#E69F00") +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  ggtitle("Mean predicted versus observed SGLT2-inhibitor benefit", subtitle = "By predicted benefit decile") +
  coord_cartesian(xlim = c(0,2.5), ylim = c(0,2.5))

p_benefit_bydeciles_mean

## SGLT2i benefit predicted vs observed - median
p_benefit_bydeciles_median <- ggplot(data=bind_rows(empty_tick,obs_v_pred_for_plot), aes(x=median_predicted_benefit*100)) +
  geom_errorbar(aes(ymax=uq_benefit*100,ymin=lq_benefit*100, color= "#E69F00"),width=0.1,size=1) +
  geom_point(aes(y = median_benefit*100, color="#E69F00"), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Predicted absolute risk reduction with SGLT2i (%)") + ylab("Observed absolute risk reduction (%)")+
  scale_colour_manual(values = "#E69F00") +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
#  ggtitle("Median predicted versus observed SGLT2-inhibitor benefit", subtitle = "By predicted benefit decile") +
  coord_cartesian(xlim = c(0,1.5), ylim = c(-.1,1.5))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_predicted_benefit_calibration.tiff"), width=6, height=5.5, units = "in", res=800) 
p_benefit_bydeciles_median
dev.off()

# mean / median pARR overall and by group
print(paste0("Overall mean pARR: ", round(100*mean(cohort$ckdpc_50egfr_sglt2i_benefit),2), "% ±", round(100*sd(cohort$ckdpc_50egfr_sglt2i_benefit),2)))
print(paste0("Overall median pARR: ", round(100*median(cohort$ckdpc_50egfr_sglt2i_benefit),2), "%, IQR ", round(100*quantile(cohort$ckdpc_50egfr_sglt2i_benefit, 0.25),2), "-",round(100*quantile(cohort$ckdpc_50egfr_sglt2i_benefit, 0.75),2)))

print(paste0("Normoalbuminuria mean pARR: ", round(100*mean(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit),2), "% ±", round(100*sd(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit),2)))
print(paste0("Normoalbuminuria median pARR: ", round(100*median(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit),2), "%, IQR ", round(100*quantile(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit, 0.25),2), "-",round(100*quantile(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit, 0.75),2)))

print(paste0("Low-level albuminuria mean pARR: ", round(100*mean(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit),2), "% ±", round(100*sd(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit),2)))
print(paste0("Low-level albuminuria median pARR: ", round(100*median(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit),2), "%, IQR ", round(100*quantile(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit, 0.25),2), "-",round(100*quantile(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit, 0.75),2)))

############################3 GUIDELINES VS MODEL TREAT/NOT TREAT################################################################

covariates_ps <- setdiff(covariates, "initiation_year")

ps.formula2 <- formula(paste("studydrug2 ~ ", paste(covariates_ps, collapse=" + ")))

#guideline-based strata
cohort_guideline_N <- cohort %>% filter(treat_guideline==F) %>% mutate(subgp="cohort_guideline_N")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_guideline_N), weight="overlap") # calculate overlap weights by subgroup
cohort_guideline_N$overlap_bygroup <- overlap$ps.weights$overlap

cohort_guideline_Y <- cohort %>% filter(treat_guideline==T) %>% mutate(subgp="cohort_guideline_Y")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_guideline_Y), weight="overlap")
cohort_guideline_Y$overlap_bygroup <- overlap$ps.weights$overlap

#Model-based strata

#Define the two strata:
#cutoff 1: benchmark scenario (same proportion treated as under guidelines but based on benefit rather than albuminuria)
cohort_model1_N <- cohort %>% filter(treat_model1 == F) %>% mutate(subgp="cohort_model1_N")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_model1_N), weight="overlap")
cohort_model1_N$overlap_bygroup <- overlap$ps.weights$overlap

cohort_model1_Y <- cohort %>% filter(treat_model1 == T) %>% mutate(subgp="cohort_model1_Y")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_model1_Y), weight="overlap")
cohort_model1_Y$overlap_bygroup <- overlap$ps.weights$overlap

#cutoff 2: top quintile of predicted benefit of those that guidelines recommend treating.
cohort_model2_N <- cohort %>% filter(treat_model2 == F) %>% mutate(subgp="cohort_model2_N")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_model2_N), weight="overlap")
cohort_model2_N$overlap_bygroup <- overlap$ps.weights$overlap

cohort_model2_Y <- cohort %>% filter(treat_model2 == T) %>% mutate(subgp="cohort_model2_Y")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_model2_Y), weight="overlap")
cohort_model2_Y$overlap_bygroup <- overlap$ps.weights$overlap

## observed data at 5 years:

#guidelines
ddist <- datadist(cohort_guideline_Y);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_Y, x=TRUE, y=TRUE, surv=TRUE)
obs_guideline_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_guideline_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv         
arr_guideline_Y <- (obs_guideline_Y_SGLT2-obs_guideline_Y_DPP4)*100
print(paste0(c("5-year ARR - guidelines recommend treating: ", sprintf("%.2f", arr_guideline_Y))))

ddist <- datadist(cohort_guideline_N);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_N, x=TRUE, y=TRUE, surv=TRUE)
obs_guideline_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_guideline_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv           
arr_guideline_N <- (obs_guideline_N_SGLT2-obs_guideline_N_DPP4)*100
print(paste0(c("5-year ARR - guidelines do not recommend treating: ", sprintf("%.2f", arr_guideline_N))))

#strategy A
ddist <- datadist(cohort_model1_Y);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model1_Y, x=TRUE, y=TRUE, surv=TRUE)
obs_model1_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_model1_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv
arr_model1_Y <- (obs_model1_Y_SGLT2-obs_model1_Y_DPP4)*100
print(paste0(c("5-year ARR - model recommends treating: ", sprintf("%.2f", arr_model1_Y))))


ddist <- datadist(cohort_model1_N);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model1_N, x=TRUE, y=TRUE, surv=TRUE)
obs_model1_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_model1_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv     
arr_model1_N <- (obs_model1_N_SGLT2-obs_model1_N_DPP4)*100
print(paste0(c("5-year ARR - model does not recommend treating: ", sprintf("%.2f", arr_model1_N))))

#strategy B
ddist <- datadist(cohort_model2_Y);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model2_Y, x=TRUE, y=TRUE, surv=TRUE)
obs_model2_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_model2_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv           
arr_model2_Y <- (obs_model2_Y_SGLT2-obs_model2_Y_DPP4)*100
print(paste0(c("5-year ARR - model recommends treating: ", sprintf("%.2f", arr_model2_Y))))


ddist <- datadist(cohort_model2_N);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model2_N, x=TRUE, y=TRUE, surv=TRUE)
obs_model2_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_model2_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv           
arr_model2_N <- (obs_model2_N_SGLT2-obs_model2_N_DPP4)*100
print(paste0(c("5-year ARR - model does not recommend treating: ", sprintf("%.2f", arr_model2_N))))


## Km plots

#guidelines vs model 1
cohort1 <- rbind(cohort_guideline_N, cohort_guideline_Y, cohort_model1_N, cohort_model1_Y)

survfit_list_1 <- lapply(split(cohort1, f=cohort1$subgp), function(x) survfit(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, data=x, weights=x$overlap_bygroup))

names(survfit_list_1[[1]][["strata"]]) <- c("DPP4i/SU   ", "SGLT2i")

# Create each plot separately
p1_1 <- ggsurvplot(
  survfit_list_1[[1]],
  data = cohort1,
  title = paste0("uACR <3mg/mmol (", round(100 * (nrow(cohort_guideline_N) / nrow(cohort)), 1), "%)"),
  subtitle = paste0("SGLT2i not recommended; 5-year ARR: ", sprintf("%.2f", arr_guideline_N), "%"),
  size = 1.5,
  fun = function(x) {100 - x * 100},
  conf.int = TRUE,
  risk.table = TRUE,
  ylim = c(0, 10),
  xlab = "Years",
  ylab = "Incidence of\nkidney disease progression (%)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  font.x = c(16), font.y = c(16), font.tickslab = c(14),
  axes.offset = FALSE,
  palette = c("#0072B2", "#E69F00")
)

p1_2 <- ggsurvplot(
  survfit_list_1[[2]],
  data = cohort1,
  title = paste0("uACR 3-30mg/mmol (", round(100 * (nrow(cohort_guideline_Y) / nrow(cohort)), 1), "%)"),
  subtitle = paste0("SGLT2i recommended; 5-year ARR: ", sprintf("%.2f", arr_guideline_Y), "%"),
  size = 1.5,
  fun = function(x) {100 - x * 100},
  conf.int = TRUE,
  risk.table = TRUE,
  ylim = c(0, 10),
  xlab = "Years",
  ylab = "Incidence of\nkidney disease progression (%)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  font.x = c(16), font.y = c(16), font.tickslab = c(14),
  axes.offset = FALSE,
  palette = c("#0072B2", "#E69F00")
)

p1_3 <- ggsurvplot(
  survfit_list_1[[3]],
  data = cohort1,
  title = paste0("pARR below threshold (", round(100 * (nrow(cohort_model1_N) / nrow(cohort)), 1), "%)"),
  subtitle = paste0("Threshold matched to uACR-recommended proportion\nSGLT2i not recommended; 5-year ARR: ", sprintf("%.2f", arr_model1_N), "%"),
  size = 1.5,
  fun = function(x) {100 - x * 100},
  conf.int = TRUE,
  risk.table = TRUE,
  ylim = c(0, 10),
  xlab = "Years",
  ylab = "Incidence of\nkidney disease progression (%)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  font.x = c(16), font.y = c(16), font.tickslab = c(14),
  axes.offset = FALSE,
  palette = c("#0072B2", "#E69F00")
)

p1_4 <- ggsurvplot(
  survfit_list_1[[4]],
  data = cohort1,
  title = paste0("pARR above threshold (", round(100 * (nrow(cohort_model1_Y) / nrow(cohort)), 1), "%)"),
  subtitle = paste0("Threshold matched to uACR-recommended proportion\nSGLT2i recommended; 5-year ARR: ", sprintf("%.2f", arr_model1_Y), "%"),
  size = 1.5,
  fun = function(x) {100 - x * 100},
  conf.int = TRUE,
  risk.table = TRUE,
  ylim = c(0, 10),
  xlab = "Years",
  ylab = "Incidence of\nkidney disease progression (%)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  font.x = c(16), font.y = c(16), font.tickslab = c(14),
  axes.offset = FALSE,
  palette = c("#0072B2", "#E69F00")
)

# Combine the plots
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_treat_asper_guidelines_vs_model1.tiff"), width=16, height=10, units = "in", res=800) 
grid.arrange(arrangeGrob(p1_1[["plot"]] + theme(legend.position=c(0.5, 0.3), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
                         p1_2[["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p1_3[["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p1_4[["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)), ncol=2, nrow=2, widths=c(1,1)))
dev.off()


table(cohort1$subgp)

#guidelines vs model 2
cohort2 <- rbind(cohort_guideline_N, cohort_guideline_Y, cohort_model2_N, cohort_model2_Y)

survfit_list_2 <- lapply(split(cohort2, f=cohort2$subgp), function(x) survfit(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, data=x, weights=x$overlap_bygroup))

names(survfit_list_2[[1]][["strata"]]) <- c("DPP4i/SU   ", "SGLT2i")

# Create each plot separately
p2_1 <- ggsurvplot(
  survfit_list_2[[1]],
  data = cohort2,
  title = paste0("uACR <3mg/mmol (", round(100 * (nrow(cohort_guideline_N) / nrow(cohort)), 1), "%)"),
  subtitle = paste0("SGLT2i not recommended; 5-year ARR: ", sprintf("%.2f", arr_guideline_N), "%"),
  size = 1.5,
  fun = function(x) {100 - x * 100},
  conf.int = TRUE,
  risk.table = TRUE,
  ylim = c(0, 10),
  xlab = "Years",
  ylab = "Incidence of\nkidney disease progression (%)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  font.x = c(16), font.y = c(16), font.tickslab = c(14),
  axes.offset = FALSE,
  palette = c("#0072B2", "#E69F00")
)

p2_2 <- ggsurvplot(
  survfit_list_2[[2]],
  data = cohort2,
  title = paste0("uACR 3-30mg/mmol (", round(100 * (nrow(cohort_guideline_Y) / nrow(cohort)), 1), "%)"),
  subtitle = paste0("SGLT2i recommended; 5-year ARR: ", sprintf("%.2f", arr_guideline_Y), "%"),
  size = 1.5,
  fun = function(x) {100 - x * 100},
  conf.int = TRUE,
  risk.table = TRUE,
  ylim = c(0, 10),
  xlab = "Years",
  ylab = "Incidence of\nkidney disease progression (%)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  font.x = c(16), font.y = c(16), font.tickslab = c(14),
  axes.offset = FALSE,
  palette = c("#0072B2", "#E69F00")
)

p2_3 <- ggsurvplot(
  survfit_list_2[[3]],
  data = cohort2,
  title = paste0("pARR below 90th percentile (", round(100 * (nrow(cohort_model2_N) / nrow(cohort)), 1), "%)"),
  subtitle = paste0("SGLT2i not recommended; 5-year ARR: ", sprintf("%.2f", arr_model2_N), "%"),
  size = 1.5,
  fun = function(x) {100 - x * 100},
  conf.int = TRUE,
  risk.table = TRUE,
  ylim = c(0, 10),
  xlab = "Years",
  ylab = "Incidence of\nkidney disease progression (%)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  font.x = c(16), font.y = c(16), font.tickslab = c(14),
  axes.offset = FALSE,
  palette = c("#0072B2", "#E69F00")
)

p2_4 <- ggsurvplot(
  survfit_list_2[[4]],
  data = cohort2,
  title = paste0("pARR above 90th percentile (", round(100 * (nrow(cohort_model2_Y) / nrow(cohort)), 1), "%)"),
  subtitle = paste0("SGLT2i recommended; 5-year ARR: ", sprintf("%.2f", arr_model2_Y), "%"),
  size = 1.5,
  fun = function(x) {100 - x * 100},
  conf.int = TRUE,
  risk.table = TRUE,
  ylim = c(0, 10),
  xlab = "Years",
  ylab = "Incidence of\nkidney disease progression (%)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  font.x = c(16), font.y = c(16), font.tickslab = c(14),
  axes.offset = FALSE,
  palette = c("#0072B2", "#E69F00")
)

# Combine the plots
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_treat_asper_guidelines_vs_model2.tiff"), width=16, height=10, units = "in", res=800) 
grid.arrange(arrangeGrob(p2_1[["plot"]] + theme(legend.position=c(0.5, 0.3), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
                         p2_2[["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p2_3[["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p2_4[["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)), ncol=2, nrow=2, widths=c(1,1)))
dev.off()


table(cohort2$subgp)

options(datadist=NULL)

############################4 NUMBERS TREATED AND EVENTS AVOIDED################################################################

### predicted data at 3 years:
#No treatment
describe(cohort$ckdpc_50egfr_score_cal) #
#estimated number of events with no one treated
ckd.notx <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score_cal/100)) 

#Treat all
describe(1-cohort$ckdpc_50egfr_survival_cal_sglt2i) #
#estimated number of events with everyone treated
ckd.tx_all <- round(nrow(cohort)*mean(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)) 


#Treat as per guideline recommendations
describe(cohort$treat_guideline)
cohort <- cohort %>% mutate(ckdpc_50egfr_score_cal.applied.guideline = 
                              ifelse(treat_guideline == F ,
                                     ckdpc_50egfr_score_cal, 100*(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)))
#estimated number of events with treatment as per guidelines
ckd.tx.guideline <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score_cal.applied.guideline/100)) 


#Treat as per strategy A
describe(cohort$treat_model1)
cohort <- cohort %>% mutate(ckdpc_50egfr_score_cal.applied.model1 = 
                              ifelse(treat_model1 == F, 
                                     ckdpc_50egfr_score_cal, 100*(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)))
#estimated number of events with treatment as per model 1
ckd.tx.model1 <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score_cal.applied.model1/100)) 


#Treat as per strategy B
describe(cohort$treat_model2)
cohort <- cohort %>% mutate(ckdpc_50egfr_score_cal.applied.model2 = 
                              ifelse(treat_model2 == F, 
                                     ckdpc_50egfr_score_cal, 100*(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)))
#estimated number of events with treatment as per model 2
ckd.tx.model2 <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score_cal.applied.model2/100)) 


print(paste0(c("Number of people treated if no one treated: 0 (0%)")))
print(paste0(c("Number of events if no one treated: ", round(ckd.notx/n.imp), " (", round(100*ckd.notx/nrow(cohort), 1), "%)"), collapse = ""))

print(paste0(c("Number of people treated with treat-everyone strategy: ", round(nrow(cohort)/n.imp), " (", round(nrow(cohort)/nrow(cohort)*100,1), "%)"), collapse = ""))
print(paste0(c("Number of events with treat-everyone strategy: ", round(ckd.tx_all/n.imp), " (", round(100*(ckd.tx_all/nrow(cohort)),1), "%)"), collapse = ""))
print(paste0(c("Number of events avoided with treat-everyone strategy: ", round((ckd.notx - ckd.tx_all)/n.imp), " (", round((ckd.tx_all-ckd.notx)/(ckd.tx_all-ckd.notx)*100,1), "%)"), collapse = ""))
print(paste0(c("NNT with treat-everyone strategy: ", round(1/(mean(cohort$ckdpc_50egfr_score_cal/100)-mean(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)))), collapse = ""))

print(paste0(c("Number of people treated with guideline treatment strategy: ", round(nrow(cohort[cohort$treat_guideline == T,])/n.imp), " (", round(nrow(cohort[cohort$treat_guideline == T,])/nrow(cohort)*100,1), "%)"), collapse = ""))
print(paste0(c("Number of events with guideline treatment strategy: ", round(ckd.tx.guideline/n.imp), " (", round(100*(ckd.tx.guideline/nrow(cohort)),1), "%)"), collapse = ""))
print(paste0(c("Number of events avoided with guideline treatment strategy: ", round(abs(ckd.tx.guideline-ckd.notx)/n.imp), " (", round((ckd.tx.guideline-ckd.notx)/(ckd.tx_all-ckd.notx)*100,1), "%)"), collapse = ""))
print(paste0(c("NNT with guideline treatment strategy: ", 
               round(1/(mean(cohort[cohort$treat_guideline == T,]$ckdpc_50egfr_score_cal/100) - 
                          mean(cohort[cohort$treat_guideline == T,]$ckdpc_50egfr_score_cal.applied.guideline/100)))), collapse = ""))

print(paste0(c("Number of people treated with treatment strategy A: ", round(nrow(cohort[cohort$treat_model1 == T,])/n.imp), " (", round(nrow(cohort[cohort$treat_model1 == T,])/nrow(cohort)*100,1), "%)"), collapse = ""))
print(paste0(c("Number of events with treatment strategy A: ", round(ckd.tx.model1/n.imp), " (", round(100*(ckd.tx.model1/nrow(cohort)),1), "%)"), collapse = ""))
print(paste0(c("Number of events avoided with treatment strategy A: ", round(abs(ckd.tx.model1-ckd.notx)/n.imp), " (", round((ckd.tx.model1-ckd.notx)/(ckd.tx_all-ckd.notx)*100,1), "%)"), collapse = ""))
print(paste0(c("NNT with treatment strategy A: ", 
               round(1/(mean(cohort[cohort$treat_model1 == T,]$ckdpc_50egfr_score_cal/100) - 
                          mean(cohort[cohort$treat_model1 == T,]$ckdpc_50egfr_score_cal.applied.model1/100)))), collapse = ""))

print(paste0(c("Number of people treated with treatment strategy B: ", round(nrow(cohort[cohort$treat_model2 == T,])/n.imp), " (", round(nrow(cohort[cohort$treat_model2 == T,])/nrow(cohort)*100,1), "%)"), collapse = ""))
print(paste0(c("Number of events with treatment strategy B: ", round(ckd.tx.model2/n.imp), " (", round(100*(ckd.tx.model2/nrow(cohort)),1), "%)"), collapse = ""))
print(paste0(c("Number of events avoided with treatment strategy B: ", round(abs(ckd.tx.model2-ckd.notx)/n.imp), " (", round((ckd.tx.model2-ckd.notx)/(ckd.tx_all-ckd.notx)*100,1), "%)"), collapse = ""))
print(paste0(c("NNT with treatment strategy B: ", 
               round(1/(mean(cohort[cohort$treat_model2 == T,]$ckdpc_50egfr_score_cal/100) - 
                          mean(cohort[cohort$treat_model2 == T,]$ckdpc_50egfr_score_cal.applied.model2/100)))), collapse = ""))


############################5 HR BY pARR THRESHOLD################################################################


#create empty data frame to which we can append the hazard ratios once calculated
subgroup_parr_SGLT2ivsDPP4iSU_hrs <- subgroup_parr_hrs <- data.frame()

#ensure treat_model2 recognised as logical variable:
cohort <- cohort %>% mutate(treat_model2 = as.logical(treat_model2))

p_value_interaction <- setNames(vector("numeric", length(outcomes)), outcomes)

## analyses stratified by subgroup
for (k in outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  
  
  # calculate number of subjects in each group
  count <- cohort %>%
    group_by(studydrug2,treat_model2) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort %>%
    group_by(studydrug2,treat_model2) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- cohort %>%
    group_by(studydrug2,treat_model2) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug2, treat_model2, events) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*treat_model2"))
  
  f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*treat_model2 + ", paste(covariates, collapse=" + ")))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.lowparr.unadj <- SE.lowparr.unadj <-
    COEFS.highparr.unadj <- SE.highparr.unadj <-
    # for the adjusted survival models
    COEFS.lowparr.adj <- SE.lowparr.adj <-
    COEFS.highparr.adj <- SE.highparr.adj <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f2, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.lowparr.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"]
    SE.lowparr.unadj[i] <- sqrt(fit.unadj$var[1,1])
    
    COEFS.highparr.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:treat_model2TRUE"]
    SE.highparr.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var),nrow(fit.unadj$var)]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted2, cohort[cohort$.imp == i,])
    
    COEFS.lowparr.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"]
    SE.lowparr.adj[i] <- sqrt(fit.adj$var[1,1])
    
    COEFS.highparr.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:treat_model2TRUE"]
    SE.highparr.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var),nrow(fit.adj$var)]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)])
    
    
    if (i == n.imp) {
      f_adjusted3 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + treat_model2 + ", paste(covariates, collapse=" + ")))
      fit.no_interaction <- coxph(f_adjusted3, cohort[cohort$.imp == i,])
      
      loglikelihood_test <- anova(fit.no_interaction, fit.adj, test = "Chisq")
      p_value_interaction[k] <- loglikelihood_test$`Pr(>|Chi|)`[2]      
    }
    
  }
  
  # pool hazard ratios
  unadjusted_lowparr <- pool.rubin.HR(COEFS.lowparr.unadj, SE.lowparr.unadj, n.imp)
  adjusted_lowparr <- pool.rubin.HR(COEFS.lowparr.adj, SE.lowparr.adj, n.imp)
  
  unadjusted_highparr <- pool.rubin.HR(COEFS.highparr.unadj, SE.highparr.unadj, n.imp)
  adjusted_highparr <- pool.rubin.HR(COEFS.highparr.adj, SE.highparr.adj, n.imp)
  
  # save pooled HR and 95% confidence interval
  unadjusted_lowparr_string <- paste0(sprintf("%.2f", round(unadjusted_lowparr[1], 2)), " (", sprintf("%.2f", round(unadjusted_lowparr[2], 2)), ", ", sprintf("%.2f", round(unadjusted_lowparr[3], 2)), ")")
  adjusted_lowparr_string <- paste0(sprintf("%.2f", round(adjusted_lowparr[1], 2)), " (", sprintf("%.2f", round(adjusted_lowparr[2], 2)), ", ", sprintf("%.2f", round(adjusted_lowparr[3], 2)), ")")
  
  unadjusted_highparr_string <- paste0(sprintf("%.2f", round(unadjusted_highparr[1], 2)), " (", sprintf("%.2f", round(unadjusted_highparr[2], 2)), ", ", sprintf("%.2f", round(unadjusted_highparr[3], 2)), ")")
  adjusted_highparr_string <- paste0(sprintf("%.2f", round(adjusted_highparr[1], 2)), " (", sprintf("%.2f", round(adjusted_highparr[2], 2)), ", ", sprintf("%.2f", round(adjusted_highparr[3], 2)), ")")
  
  # combine in dataframe that we can tabulate
  presegfr_lowparr_hr <- cbind(outcome=k, count[1,c(2:3)], followup[1,c(2:3)], events[1,c(2:3)],
                             unadjusted=unadjusted_lowparr_string, adjusted=adjusted_lowparr_string
  )
  presegfr_highparr_hr <- cbind(outcome=k, count[2,c(2:3)], followup[2,c(2:3)], events[2,c(2:3)],
                                unadjusted=unadjusted_highparr_string, adjusted=adjusted_highparr_string
  )
  
  outcome_subgroup_parr_SGLT2ivsDPP4iSU_hrs <- rbind(presegfr_lowparr_hr, presegfr_highparr_hr)
  
  temp <- rbind(
    cbind(outcome = k, contrast = "pARR below 90th percentile", analysis = "Adjusted",
          HR = adjusted_lowparr[1], LB = adjusted_lowparr[2], UB = adjusted_lowparr[3], string = adjusted_lowparr_string),
    cbind(outcome = k, contrast = "pARR above 90th percentile", analysis = "Adjusted",
          HR = adjusted_highparr[1], LB = adjusted_highparr[2], UB = adjusted_highparr[3], string = adjusted_highparr_string)
  )
  
  subgroup_parr_hrs <- rbind(subgroup_parr_hrs, temp)
  subgroup_parr_SGLT2ivsDPP4iSU_hrs <- rbind(subgroup_parr_SGLT2ivsDPP4iSU_hrs, outcome_subgroup_parr_SGLT2ivsDPP4iSU_hrs)
  
}


subgroup_parr_hrs <- subgroup_parr_hrs %>% cbind(subgroup_parr_SGLT2ivsDPP4iSU_hrs %>% select(-c(outcome, adjusted, unadjusted)))

subgroup_parr_hrs <- subgroup_parr_hrs %>%
  separate(`DPP4i/SU_events`, into = c("DPP4i/SU_events_number", "DPP4i/SU_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(SGLT2i_events, into = c("SGLT2i_events_number", "SGLT2i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `DPP4i/SU_events_percentage` = str_replace(`DPP4i/SU_events_percentage`, "\\)", ""),
    SGLT2i_events_percentage = str_replace(SGLT2i_events_percentage, "\\)", ""),
    `DPP4i/SU_nN` = paste0(`DPP4i/SU_events_number`, "/", `DPP4i/SU_count`),
    SGLT2i_nN = paste0(SGLT2i_events_number, "/", SGLT2i_count)
  )


# save all_hrs table and SGLT2i vs DPP4i/su table
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
save(subgroup_parr_hrs, file=paste0(today, "_subgroup_parr_hrs.Rda"))
save(subgroup_parr_SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_subgroup_parr_SGLT2ivsDPP4iSU_hrs.Rda"))

