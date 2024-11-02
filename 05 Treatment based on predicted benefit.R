
############################0 SETUP################################################################

# Setup
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Raw data/")
load("2024-10-29_t2d_ckdpc_recalibrated_with_adjsurv.Rda")


############################1 DEFINE CUTOFFS################################################################

## consider following methods to define cut-off:
# 1 take as cut-off the upper top decile value of predicted benefit in the group of people that guidelines recommend treating
# 2 take as cut-off the percentile corresponding to the percentage of people that would be treated based on current guidelines (ie same number treated but targeted based on our approach rather than purely uACR cut-off)

cutoff2 <- cohort %>% filter(albuminuria == T) %>% .$ckdpc_50egfr_sglt2i_benefit %>% quantile(0.90) %>% as.numeric()

cutoff2_equivalent_40egfr_score <- cohort %>% filter(albuminuria == T) %>% .$ckdpc_40egfr_score %>% quantile(0.90) %>% as.numeric()

cutoff1 <- cohort %>% .$ckdpc_50egfr_sglt2i_benefit %>% quantile( 
  1 - (cohort %>% filter(albuminuria == T) %>% nrow() / cohort  %>% nrow())) %>% as.numeric()

cutoff1_equivalent_40egfr_score <- cohort  %>% .$ckdpc_40egfr_score %>% quantile( 
  1 - (cohort %>% filter(albuminuria == T) %>% nrow() / cohort  %>% nrow())) %>% as.numeric()

cohort <- cohort %>% mutate(
  treat_guideline = ifelse(albuminuria == F, F, T),            # everyone current guidelines recommend treating (uACR cutoff)
  treat_model2 = ifelse(ckdpc_50egfr_sglt2i_benefit > cutoff2, T, F),        # anyone with benefit above upper quartile of those currently recommended treatment in 3-30mg/mmol
  treat_model1 = ifelse(ckdpc_50egfr_sglt2i_benefit > cutoff1, T, F)         # treating same proportion of people as guidelines would recommend but based on risk-score rather than uACR cutoff
  )

############################2 EVALUATE OBSERVED BENEFIT DISTRIBUTION/CALIBRATION################################################################

## FIGURE 2A: histogram of predicted benefit (people with albuminuria <30mg/mmol only)
benefit_histogram <- ggplot(cohort %>%
                      mutate(predicted_benefit_percent = ckdpc_50egfr_sglt2i_benefit * 100), 
                    aes(x = predicted_benefit_percent)) +
  geom_histogram(aes(y = ..count.. / 10, fill = predicted_benefit_percent > cutoff2*100),
                 binwidth = 0.02, color = "black") +  # Adjust binwidth as necessary
  scale_fill_manual(values = c("TRUE" = "#E69F00", "FALSE" = "grey")) +
  geom_vline(xintercept = cutoff2*100, linetype = "dashed", color = "black", size = 1) +
  labs(x = "Predicted absolute risk reduction with SGLT2i (%)", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "none") + 
  coord_cartesian(xlim=c(0,3))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_predicted_benefit_histogram.tiff"), width=5, height=2.5, units = "in", res=800) 
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
#cutoff 1: same number of people treated as per guidelines except based on risk rather than albuminuria 3 mg/mmol threshold
cohort_model1_N <- cohort %>% filter(treat_model1 == F) %>% mutate(subgp="cohort_model1_N")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_model1_N), weight="overlap")
cohort_model1_N$overlap_bygroup <- overlap$ps.weights$overlap

cohort_model1_Y <- cohort %>% filter(treat_model1 == T) %>% mutate(subgp="cohort_model1_Y")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_model1_Y), weight="overlap")
cohort_model1_Y$overlap_bygroup <- overlap$ps.weights$overlap

#cutoff 2: top decile of predicted benefit amongst those with uACR <30 that guidelines recommend treating.
cohort_model2_N <- cohort %>% filter(treat_model2 == F) %>% mutate(subgp="cohort_model2_N")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_model2_N), weight="overlap")
cohort_model2_N$overlap_bygroup <- overlap$ps.weights$overlap

cohort_model2_Y <- cohort %>% filter(treat_model2 == T) %>% mutate(subgp="cohort_model2_Y")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_model2_Y), weight="overlap")
cohort_model2_Y$overlap_bygroup <- overlap$ps.weights$overlap

## Km plots

#guidelines vs model 2
cohort2 <- rbind(cohort_guideline_N, cohort_guideline_Y, cohort_model2_N, cohort_model2_Y)

survfit_list_2 <- lapply(split(cohort2, f=cohort2$subgp), function(x) survfit(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, data=x, weights=x$overlap_bygroup))

names(survfit_list_2[[1]][["strata"]]) <- c("DPP4i/SU   ", "SGLT2i")


p2 <- ggsurvplot_list(survfit_list_2,
                      data=cohort2 %>% filter(.imp == n.imp),
                      title=c(
                        paste0("Guidelines do not recommend SGLT2i (uACR <3mg/mmol; ",round(100*(nrow(cohort_guideline_N)/nrow(cohort)),1),"%)"), 
                        paste0("Guidelines recommend SGLT2i (uACR ≥3mg/mmol; ",round(100*(nrow(cohort_guideline_Y)/nrow(cohort)),1),"%)"), 
                        paste0("Model: No SGLT2i (",round(100*(nrow(cohort_model2_N)/nrow(cohort)),1),"%)"),
                        paste0("Model: SGLT2i (",round(100*(nrow(cohort_model2_Y)/nrow(cohort)),1),"%)")),
                      size = 1.5,
                      fun = function(x) {100 - x*100},
                      conf.int = T,
                      risk.table=TRUE,
                      ylim = c(0, 10.2),
                      xlab = "Years",
                      ylab = "Incidence of\nkidney disease progression (%)",
                      break.time.by = 1,
                      ggtheme = theme_classic(),
                      font.x = c(16),font.y = c(16),font.tickslab = c(14),
                      axes.offset = FALSE, # start the axis at the origin
                      palette = c("#0072B2","#E69F00"))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_treat_asper_guidelines_vs_model2.tiff"), width=16, height=10, units = "in", res=800) 
grid.arrange(arrangeGrob(p2[[1]][["plot"]] + theme(legend.position=c(0.5, 0.7), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
                         p2[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p2[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p2[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)), 
                         ncol=2, nrow=2, widths=c(1,1)))
dev.off()

table(cohort2$subgp)

#guidelines vs model 2
cohort1 <- rbind(cohort_guideline_N, cohort_guideline_Y, cohort_model1_N, cohort_model1_Y)

survfit_list_1 <- lapply(split(cohort1, f=cohort1$subgp), function(x) survfit(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, data=x, weights=x$overlap_bygroup))

names(survfit_list_1[[1]][["strata"]]) <- c("DPP4i/SU   ", "SGLT2i")


p1 <- ggsurvplot_list(survfit_list_1,
                      data=cohort1,
                      title=c(
                        paste0("Guidelines: No SGLT2i (",round(100*(nrow(cohort_guideline_N)/nrow(cohort)),1),"%)"), 
                        paste0("Guidelines: SGLT2i (",round(100*(nrow(cohort_guideline_Y)/nrow(cohort)),1),"%)"), 
                        paste0("Model: No SGLT2i (",round(100*(nrow(cohort_model1_N)/nrow(cohort)),1),"%)"),
                        paste0("Model: SGLT2i (",round(100*(nrow(cohort_model1_Y)/nrow(cohort)),1),"%)")),
                      size = 1.5,
                      fun = function(x) {100 - x*100},
                      conf.int = T,
                      risk.table=TRUE,
                      ylim = c(0, 10),
                      xlab = "Years",
                      ylab = "Incidence of\nkidney disease progression (%)",
                      break.time.by = 1,
                      ggtheme = theme_classic(),
                      font.x = c(16),font.y = c(16),font.tickslab = c(14),
                      axes.offset = FALSE, # start the axis at the origin
                      palette = c("#0072B2","#E69F00"))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_treat_asper_guidelines_vs_model1.tiff"), width=16, height=10, units = "in", res=800) 
grid.arrange(arrangeGrob(p1[[1]][["plot"]] + theme(legend.position=c(0.5, 0.7), legend.title=element_blank(), legend.text=element_text(size=16, face="bold"), plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)) + guides(colour = guide_legend(nrow = 1)),
                         p1[[2]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p1[[3]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)),
                         p1[[4]][["plot"]] + theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5, margin=margin(b=5, t=12), size = 22), plot.subtitle=element_text(hjust=0.5, margin=margin(b=-25), size = 16), axis.title.y=element_text(vjust=-0.3), plot.margin=margin(l=5, r=5)), ncol=2, nrow=2, widths=c(1,1)))
dev.off()

table(cohort1$subgp)

options(datadist=NULL)

############################4 NUMBERS TREATED AND EVENTS AVOIDED################################################################

### predicted data at 3 years:
#No treatment
describe(cohort$ckdpc_50egfr_score_cal) #
#estimated number of events with no one treated
ckd.notx <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score_cal/100)) 
ckd.notx

#Treat all
describe(1-cohort$ckdpc_50egfr_survival_cal_sglt2i) #
#estimated number of events with everyone treated
ckd.tx_all <- round(nrow(cohort)*mean(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)) 
ckd.tx_all
#estimated number of events avoided (with everyone treated)
ckd.notx - ckd.tx_all
100*(ckd.tx/nrow(cohort))

1/(mean(cohort$ckdpc_50egfr_score_cal/100)-mean(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)) #NNT

#Treat as per guideline recommendations
describe(cohort$treat_guideline)
cohort <- cohort %>% mutate(ckdpc_50egfr_score_cal.applied = 
                              ifelse(treat_guideline == F ,
                                     ckdpc_50egfr_score_cal, 100*(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)))
#estimated number of events with treatment as per guidelines
ckd.tx <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score_cal.applied/100)) 
ckd.tx
#estimated number of events avoided (with guidelines)
ckd.notx - ckd.tx
(ckd.notx - ckd.tx)/ckd.tx_all
100*(ckd.tx/nrow(cohort))

1/(mean(cohort[cohort$treat_guideline == T,]$ckdpc_50egfr_score_cal/100)-mean(cohort[cohort$treat_guideline == T,]$ckdpc_50egfr_score_cal.applied/100)) #NNT

#Treat as per model 1
describe(cohort$treat_model2)
cohort <- cohort %>% mutate(ckdpc_50egfr_score_cal.applied = 
                              ifelse(treat_model2 == F, 
                                     ckdpc_50egfr_score_cal, 100*(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)))
#estimated number of events with treatment as per model 1
ckd.tx <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score_cal.applied/100)) 
ckd.tx
#estimated number of events avoided (with model 1)
ckd.notx - ckd.tx  
(ckd.notx - ckd.tx)/ckd.tx_all
100*(ckd.tx/nrow(cohort))

1/(mean(cohort[cohort$treat_model2 == T,]$ckdpc_50egfr_score_cal/100)-mean(cohort[cohort$treat_model2 == T,]$ckdpc_50egfr_score_cal.applied/100)) #NNT

#Treat as per model 2
describe(cohort$treat_model1)
cohort <- cohort %>% mutate(ckdpc_50egfr_score_cal.applied = 
                              ifelse(treat_model1 == F, 
                                     ckdpc_50egfr_score_cal, 100*(1-cohort$ckdpc_50egfr_survival_cal_sglt2i)))
#estimated number of events with treatment as per model 2
ckd.tx <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score_cal.applied/100)) 
ckd.tx
#estimated number of events avoided (with model 2)
ckd.notx - ckd.tx  
(ckd.notx - ckd.tx)/ckd.tx_all
100*(ckd.tx/nrow(cohort))

1/(mean(cohort[cohort$treat_model1 == T,]$ckdpc_50egfr_score_cal/100)-mean(cohort[cohort$treat_model1 == T,]$ckdpc_50egfr_score_cal.applied/100)) #NNT



## observed data for 3 years:

ddist <- datadist(cohort_guideline_Y);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_Y, x=TRUE, y=TRUE, surv=TRUE)
obs_guideline_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times = 3)$surv
obs_guideline_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times = 3)$surv           
print(paste0(c("3-year ARR - guidelines recommend treating: ", (obs_guideline_Y_SGLT2-obs_guideline_Y_DPP4)*100)))

ddist <- datadist(cohort_guideline_N);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_N, x=TRUE, y=TRUE, surv=TRUE)
obs_guideline_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times = 3)$surv
obs_guideline_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times = 3)$surv           
print(paste0(c("3-year ARR - guidelines do not recommend treating: ", (obs_guideline_N_SGLT2-obs_guideline_N_DPP4)*100)))

ddist <- datadist(cohort_model2_Y);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model2_Y, x=TRUE, y=TRUE, surv=TRUE)
obs_model2_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times = 3)$surv
obs_model2_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times = 3)$surv           
print(paste0(c("3-year ARR - model recommends treating: ", (obs_model2_Y_SGLT2-obs_model2_Y_DPP4)*100)))

ddist <- datadist(cohort_model2_N);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model2_N, x=TRUE, y=TRUE, surv=TRUE)
obs_model2_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times = 3)$surv
obs_model2_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times = 3)$surv           
print(paste0(c("3-year ARR - model does not recommend treating: ", (obs_model2_N_SGLT2-obs_model2_N_DPP4)*100)))

ddist <- datadist(cohort_model1_Y);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model1_Y, x=TRUE, y=TRUE, surv=TRUE)
obs_model1_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times = 3)$surv
obs_model1_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times = 3)$surv           
print(paste0(c("3-year ARR - model recommends treating: ", (obs_model1_Y_SGLT2-obs_model1_Y_DPP4)*100)))

ddist <- datadist(cohort_model1_N);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model1_N, x=TRUE, y=TRUE, surv=TRUE)
obs_model1_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times = 3)$surv
obs_model1_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times = 3)$surv           
print(paste0(c("3-year ARR - model does not recommend treating: ", (obs_model1_N_SGLT2-obs_model1_N_DPP4)*100)))

nnt_median_overall <- cohort %>% summarise(median_benefit = median(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/median_benefit) %>% select(nnt) %>% as.numeric()
nnt_median_guideline_recommended <- cohort %>% filter(risk_group != "uACR <3mg/mmol") %>% summarise(median_benefit = median(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/median_benefit) %>% select(nnt) %>% as.numeric()
nnt_median_guideline_not_recommended <- cohort %>% filter(albuminuria == F) %>% summarise(median_benefit = median(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/median_benefit) %>% select(nnt) %>% as.numeric()
nnt_median_model2_recommended <- cohort %>% filter(treat_model2 == T) %>% summarise(median_benefit = median(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/median_benefit) %>% select(nnt) %>% as.numeric()
nnt_median_model2_not_recommended <- cohort %>% filter(treat_model2 == F) %>% summarise(median_benefit = median(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/median_benefit) %>% select(nnt) %>% as.numeric()

nnt_mean_overall <- cohort %>% summarise(mean_benefit = mean(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/mean_benefit) %>% select(nnt) %>% as.numeric()
nnt_mean_guideline_recommended <- cohort %>% filter(risk_group != "uACR <3mg/mmol") %>% summarise(mean_benefit = mean(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/mean_benefit) %>% select(nnt) %>% as.numeric()
nnt_mean_guideline_not_recommended <- cohort %>% filter(albuminuria == F) %>% summarise(mean_benefit = mean(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/mean_benefit) %>% select(nnt) %>% as.numeric()
nnt_mean_model2_recommended <- cohort %>% filter(treat_model2 == T) %>% summarise(mean_benefit = mean(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/mean_benefit) %>% select(nnt) %>% as.numeric()
nnt_mean_model2_not_recommended <- cohort %>% filter(treat_model2 == F) %>% summarise(mean_benefit = mean(survdiff_ckd_egfr50)) %>% mutate(nnt = 1/mean_benefit) %>% select(nnt) %>% as.numeric()

# 
# #### counterfactual observed events if no one treated
# events_none_treated <- (1-obs_model2_Y_DPP4) * nrow(cohort_model2_Y) + (1-obs_model2_N_DPP4) * nrow(cohort_model2_N)
# events_none_treated / n.imp
# #### counterfactual observed events if everyone treated
# events_everyone_treated <- (1-obs_model2_Y_SGLT2) * nrow(cohort_model2_Y) + (1-obs_model2_N_SGLT2) * nrow(cohort_model2_N)
# events_everyone_treated / n.imp
# # number of events avoided if everyone treated
# events_avoided_max <- events_none_treated - events_everyone_treated
# events_avoided_max / n.imp
# nnt_everyone_treated <- nrow(cohort) / events_avoided_max
# nnt_everyone_treated
# 
# #### counterfactual observed events if guideline treated
# events_guideline_treated <- (1-obs_guideline_Y_SGLT2) * nrow(cohort_guideline_Y) + (1-obs_guideline_N_DPP4) * nrow(cohort_guideline_N)
# events_guideline_treated / n.imp
# events_avoided_guideline <- events_none_treated - events_guideline_treated
# events_avoided_guideline / n.imp
# events_avoided_guideline / events_avoided_max
# nnt_guideline_treated <- nrow(cohort_guideline_Y) / events_avoided_guideline
# nnt_guideline_treated
# 
# #### counterfactual observed events if treated per model 1
# events_model2_treated <- (1-obs_model2_Y_SGLT2) * nrow(cohort_model2_Y) + (1-obs_model2_N_DPP4) * nrow(cohort_model2_N)
# events_model2_treated / n.imp
# events_avoided_model2 <- events_none_treated - events_model2_treated
# events_avoided_model2 / n.imp
# events_avoided_model2 / events_avoided_max
# nnt_model2_treated <- nrow(cohort_model2_Y) / events_avoided_model2
# nnt_model2_treated
# 
# #### counterfactual observed events if treated per model 2
# events_model1_treated <- (1-obs_model1_Y_SGLT2) * nrow(cohort_model1_Y) + (1-obs_model1_N_DPP4) * nrow(cohort_model1_N)
# events_model1_treated / n.imp
# events_avoided_model1 <- events_none_treated - events_model1_treated
# events_avoided_model1 / n.imp
# events_avoided_model1 / events_avoided_max
# nnt_model1_treated <- nrow(cohort_model1_Y) / events_avoided_model1
# nnt_model1_treated
# 
