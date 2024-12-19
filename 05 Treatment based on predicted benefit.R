
############################0 SETUP################################################################

# Setup
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
load(paste0(today, "_t2d_ckdpc_recalibrated_with_adjsurv.Rda"))


############################1 DEFINE CUTOFFS################################################################

# calculate predicted sglt2 benefit (absolute risk reduction = ARR):
# ARR = S0(t)^HR - S0(t)
trial_hr_kf_sglt2i <- 0.62

cohort <- cohort %>% 
  mutate(ckdpc_50egfr_survival=(100-ckdpc_50egfr_score)/100,
         ckdpc_50egfr_survival_sglt2i=ckdpc_50egfr_survival^trial_hr_kf_sglt2i,
         ckdpc_50egfr_sglt2i_benefit=ckdpc_50egfr_survival_sglt2i - ckdpc_50egfr_survival)

## consider following methods to define cut-off:
# 1 Strategy A: choose cut-off that matches treatment proportion of guidelines but improves outcomes through better-targeted treatment
# 2 Strategy B: take as cut-off the top decile  (10%) of predicted benefit

cutoff1 <- cohort %>% .$ckdpc_50egfr_sglt2i_benefit %>% quantile( 
  1 - (cohort %>% filter(albuminuria == T) %>% nrow() / cohort  %>% nrow())) %>% as.numeric()

cutoff1_equivalent_50egfr_score <- cohort  %>% .$ckdpc_50egfr_score %>% quantile( 
  1 - (cohort %>% filter(albuminuria == T) %>% nrow() / cohort  %>% nrow())) %>% as.numeric()

cutoff2 <- cohort %>% .$ckdpc_50egfr_sglt2i_benefit %>% quantile(0.90) %>% as.numeric()

cutoff2_equivalent_50egfr_score <- cohort %>% .$ckdpc_50egfr_score %>% quantile(0.90) %>% as.numeric()

cohort <- cohort %>% mutate(
  treat_guideline = ifelse(albuminuria == F, F, T),                          # everyone current guidelines recommend treating (uACR cutoff)
  treat_model1 = ifelse(ckdpc_50egfr_sglt2i_benefit > cutoff1, T, F),        # strategy A (matching uACR treatment proportion)
  treat_model2 = ifelse(ckdpc_50egfr_sglt2i_benefit > cutoff2, T, F),        # strategy B (top 10%)
)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv

vars <- c(vars, "studydrug2")
factors <- c(factors, "studydrug2")
table3 <- CreateTableOne(vars = vars, strata = "treat_model1", data = cohort,
                        factorVars = factors, test = F)

tabforprint3 <- print(table3, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)

write.csv2(tabforprint3, file = paste0(today, "_baseline_table_by_parr.csv"))
############################2 HR BY RISK SCORE################################################################
## check whether there is evidence of treatment heterogeneity by baseline risk (figure 1B)

# fit model using interaction term of treatment with risk score, 
# modelled with restricted cubic splines [rcs()] with 3-5 knots

ddist <- cohort %>% datadist()
options(datadist='ddist')


# Define the range of knots to test
k_range <- 3:5

# Initialize empty vectors to store results
aic_values <- numeric(length(k_range))
bic_values <- numeric(length(k_range))

# Loop over each value of k
for (i in seq_along(k_range)) {
  k <- k_range[i]
  
  # Fit the model with k knots
  model <- cph(
    as.formula(paste0(
      "Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ studydrug2*rcs(ckdpc_50egfr_score,", k, ") + ",
      paste(setdiff(covariates, "ckdpc_50egfr_score"), collapse=" + ") # need to remove risk score from covariate list as already specified in interaction term
    )),
    data = cohort %>% filter(.imp == n.imp), x = TRUE, y = TRUE
  )
  
  # Store the AIC and BIC values
  aic_values[i] <- AIC(model)
  bic_values[i] <- BIC(model)
}

# Find the optimal k based on minimum AIC or BIC
optimal_k_aic <- k_range[which.min(aic_values)]
optimal_k_bic <- k_range[which.min(bic_values)]

# Fit the final model using the optimal number of knots based on AIC
final_model <- cph(
  as.formula(paste0(
    "Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ studydrug2*rcs(ckdpc_50egfr_score,", optimal_k_bic, ") + ",
    paste(setdiff(covariates, "ckdpc_50egfr_score"), collapse=" + ")
  )),
  data = cohort %>% filter(.imp == n.imp), x = TRUE, y = TRUE
)

# Print optimal k values
cat("Optimal number of knots based on AIC:", optimal_k_aic, "\n")
cat("Optimal number of knots based on BIC:", optimal_k_bic, "\n")

anova(final_model)
p_value_non_linear <- anova(final_model)[2,3] # p value for non-linear interaction term

# create data frame with range of scores by study drug
contrast_spline <- contrast(final_model, 
                            list(studydrug2 = "SGLT2i", ckdpc_50egfr_score = seq(quantile(cohort$ckdpc_50egfr_score, .01, na.rm=TRUE), quantile(cohort$ckdpc_50egfr_score, .99, na.rm=TRUE), by=0.05)), 
                            list(studydrug2 = "DPP4i/SU", ckdpc_50egfr_score = seq(quantile(cohort$ckdpc_50egfr_score, .01, na.rm=TRUE), quantile(cohort$ckdpc_50egfr_score, .99, na.rm=TRUE), by=0.05))
)

contrast_spline_df <- as.data.frame(contrast_spline[c('ckdpc_50egfr_score','Contrast','Lower','Upper')])
# plot
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_HR_by_ckd_egfr50_risk.tiff"), width=10, height=4, units = "in", res=800) 
ggplot(data=contrast_spline_df, aes(x=ckdpc_50egfr_score, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=ckdpc_50egfr_score, y=exp(Contrast)), size=1) +
  xlab(expression(paste("Predicted 3-year risk of kidney disease progression"))) +
  ylab("Hazard ratio") +
  coord_trans(y = "log10") +
  scale_x_continuous(breaks = seq(0,20,.5)) +
  scale_y_continuous(breaks = c(seq(0, 0.8, 0.1), seq(0.8, 1.6, 0.2))) +
  geom_ribbon(data=contrast_spline_df, aes(x=ckdpc_50egfr_score, ymin=exp(Lower), ymax=exp(Upper)), alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  geom_hline(aes(yintercept = 0.62, linetype = "hr", size="hr"), color="#D55E00")  +
  geom_hline(aes(yintercept = 0.68, linetype = "hr_95", size="hr_95"), color="#D55E00")  +
  geom_hline(aes(yintercept = 0.57, linetype = "hr_95", size="hr_95"), color="#D55E00")  +
  annotate("text", x = mean(range(contrast_spline_df$ckdpc_50egfr_score)), y = 0.35, 
           label = paste0("italic(p)~'='~italic(", sprintf("%.2f", p_value_non_linear), "~'for non-linear interaction term')"), 
           size = 5, hjust = 0.5, parse = T) +
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="bottom",
        legend.title = element_text(size=14, face = "italic"),
        legend.text = element_text(face="italic")) +
  scale_linetype_manual(values = c(hr = "twodash", hr_95 = "twodash"), labels = c(hr = "0.62", hr_95 = "95% CI 0.56-0.68"), name="Trial meta-analysis hazard ratio") +
  scale_size_manual(values = c(hr = 1, hr_95 = 0.5), labels = c(hr = "0.62", hr_95 = "95% CI 0.56-0.68"), name="Trial meta-analysis hazard ratio")
dev.off()

options(datadist = NULL)
## there is no significant treatment heterogeneity by baseline risk score

############################3 EVALUATE OBSERVED BENEFIT DISTRIBUTION/CALIBRATION################################################################

## FIGURE 2A: histogram of predicted benefit
benefit_histogram <- ggplot(cohort %>%
                      mutate(predicted_benefit_percent = ckdpc_50egfr_sglt2i_benefit * 100), 
                    aes(x = predicted_benefit_percent)) +
  geom_histogram(aes(y = ..count.. / n.imp, fill = predicted_benefit_percent > cutoff1*100),
                 binwidth = 0.02, color = "black") +  # Adjust binwidth as necessary
  scale_fill_manual(values = c("TRUE" = "#E69F00", "FALSE" = "grey")) +
  geom_vline(xintercept = cutoff1*100, linetype = "dashed", color = "black", size = 1) +
#  geom_vline(xintercept = cutoff2*100, linetype = "dashed", color = "black", size = 1) +

  annotate("text", x = cutoff1*100, y = Inf, label = "pARR threshold", vjust = -0.5, hjust = 1.1, angle = 90, size = 4, color = "black") +
  #  annotate("text", x = cutoff1*100, y = Inf, label = "pARR threshold A", vjust = -0.5, hjust = 1.1, angle = 90, size = 4, color = "black") +
  #  annotate("text", x = cutoff2*100, y = Inf, label = "pARR threshold B", vjust = -0.5, hjust = 1.1, angle = 90, size = 4, color = "black") +

  labs(x = "Predicted absolute risk reduction with SGLT2i (%)", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "none") + 
  coord_cartesian(xlim=c(0,2.5))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_predicted_benefit_histogram.tiff"), width=6, height=3.75, units = "in", res=800) 
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
  coord_cartesian(xlim = c(0,1.5), ylim = c(-.1,1.56))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_predicted_benefit_calibration.tiff"), width=6, height=5.5, units = "in", res=800) 
p_benefit_bydeciles_median
dev.off()

# calibration slope
x <- lm(mean_benefit ~ mean_predicted_benefit, data = obs_v_pred_for_plot)
slope <- coef(x)
slope <- slope[2]
x <- confint(x)
slope_lc <- x[2,1]
slope_uc <- x[2,2]
rm(x)
print(paste0("Calibration slope for mean pARR and observed ARR by decile: ", round(slope, 3), ", 95% CI ", round(slope_lc, 3), "-", round(slope_uc,3)))

# mean / median pARR overall and by group
print(paste0("Overall mean pARR: ", round(100*mean(cohort$ckdpc_50egfr_sglt2i_benefit),2), "% ±", round(100*sd(cohort$ckdpc_50egfr_sglt2i_benefit),2)))
print(paste0("Overall median pARR: ", round(100*median(cohort$ckdpc_50egfr_sglt2i_benefit),2), "%, IQR ", round(100*quantile(cohort$ckdpc_50egfr_sglt2i_benefit, 0.25),2), "-",round(100*quantile(cohort$ckdpc_50egfr_sglt2i_benefit, 0.75),2)))

print(paste0("Normoalbuminuria mean pARR: ", round(100*mean(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit),2), "% ±", round(100*sd(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit),2)))
print(paste0("Normoalbuminuria median pARR: ", round(100*median(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit),2), "%, IQR ", round(100*quantile(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit, 0.25),2), "-",round(100*quantile(cohort[cohort$albuminuria == F,]$ckdpc_50egfr_sglt2i_benefit, 0.75),2)))

print(paste0("Low-level albuminuria mean pARR: ", round(100*mean(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit),2), "% ±", round(100*sd(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit),2)))
print(paste0("Low-level albuminuria median pARR: ", round(100*median(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit),2), "%, IQR ", round(100*quantile(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit, 0.25),2), "-",round(100*quantile(cohort[cohort$albuminuria == T,]$ckdpc_50egfr_sglt2i_benefit, 0.75),2)))

############################4 GUIDELINES VS MODEL TREAT/NOT TREAT################################################################

ps.formula2 <- formula(paste("studydrug2 ~ ", paste(covariates, collapse=" + ")))


## strata based on whether people are recommended treatment my albuminuria (guidelines) and/or the model

# not recommended by model or based on albuminuria
cohort_guideline_N_model_N <- cohort %>% filter(treat_guideline == F & treat_model1 == F) %>% mutate(subgp="guideline_N_model_N")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_guideline_N_model_N), weight="overlap") # calculate overlap weights by subgroup
cohort_guideline_N_model_N$overlap_bygroup <- overlap$ps.weights$overlap

ddist <- datadist(cohort_guideline_N_model_N);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_N_model_N, x=TRUE, y=TRUE, surv=TRUE)
obs_guideline_N_model_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_guideline_N_model_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv         
arr_guideline_N_model_N <- (obs_guideline_N_model_N_SGLT2-obs_guideline_N_model_N_DPP4)*100

se_guideline_N_model_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$std.err
se_guideline_N_model_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$std.err
se_guideline_N_model_N_ARR <- sqrt(se_guideline_N_model_N_SGLT2^2 + se_guideline_N_model_N_DPP4^2)

print(paste0(c("5-year ARR - Albuminuria <3mg/mmol, model does not recommend: ", sprintf("%.2f", arr_guideline_N_model_N))))


# recommended based on albuminuria only
cohort_guideline_Y_model_N <- cohort %>% filter(treat_guideline == T & treat_model1 == F) %>% mutate(subgp="guideline_Y_model_N")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_guideline_Y_model_N), weight="overlap")
cohort_guideline_Y_model_N$overlap_bygroup <- overlap$ps.weights$overlap

ddist <- datadist(cohort_guideline_Y_model_N);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_Y_model_N, x=TRUE, y=TRUE, surv=TRUE)
obs_guideline_Y_model_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_guideline_Y_model_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv         
arr_guideline_Y_model_N <- (obs_guideline_Y_model_N_SGLT2-obs_guideline_Y_model_N_DPP4)*100

se_guideline_Y_model_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$std.err
se_guideline_Y_model_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$std.err
se_guideline_Y_model_N_ARR <- sqrt(se_guideline_Y_model_N_SGLT2^2 + se_guideline_Y_model_N_DPP4^2)

print(paste0(c("5-year ARR - Albuminuria ≥3mg/mmol, model does not recommend: ", sprintf("%.2f", arr_guideline_Y_model_N))))


# recommended based on model only
cohort_guideline_N_model_Y <- cohort %>% filter(treat_guideline == F & treat_model1 == T) %>% mutate(subgp="guideline_N_model_Y")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_guideline_N_model_Y), weight="overlap") # calculate overlap weights by subgroup
cohort_guideline_N_model_Y$overlap_bygroup <- overlap$ps.weights$overlap

ddist <- datadist(cohort_guideline_N_model_Y);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_N_model_Y, x=TRUE, y=TRUE, surv=TRUE)
obs_guideline_N_model_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_guideline_N_model_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv         
arr_guideline_N_model_Y <- (obs_guideline_N_model_Y_SGLT2-obs_guideline_N_model_Y_DPP4)*100

se_guideline_N_model_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$std.err
se_guideline_N_model_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$std.err
se_guideline_N_model_Y_ARR <- sqrt(se_guideline_N_model_Y_SGLT2^2 + se_guideline_N_model_Y_DPP4^2)

print(paste0(c("5-year ARR - Model recommends SGLT2i, albuminuria <3mg/mmol: ", sprintf("%.2f", arr_guideline_N_model_Y))))


#recommended by both model and based on albuminuria
cohort_guideline_Y_model_Y <- cohort %>% filter(treat_guideline == T & treat_model1 == T) %>% mutate(subgp="guideline_Y_model_Y")
overlap <- SumStat(ps.formula=ps.formula2, data=as.data.frame(cohort_guideline_Y_model_Y), weight="overlap")
cohort_guideline_Y_model_Y$overlap_bygroup <- overlap$ps.weights$overlap

ddist <- datadist(cohort_guideline_Y_model_Y);options(datadist='ddist')
model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_Y_model_Y, x=TRUE, y=TRUE, surv=TRUE)
obs_guideline_Y_model_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
obs_guideline_Y_model_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv         
arr_guideline_Y_model_Y <- (obs_guideline_Y_model_Y_SGLT2-obs_guideline_Y_model_Y_DPP4)*100

se_guideline_Y_model_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$std.err
se_guideline_Y_model_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$std.err
se_guideline_Y_model_Y_ARR <- sqrt(se_guideline_Y_model_Y_SGLT2^2 + se_guideline_Y_model_Y_DPP4^2)

print(paste0(c("5-year ARR - SGLT2i recommended by both model and based on albuminuria: ", sprintf("%.2f", arr_guideline_Y_model_Y))))

options(datadist=NULL)

cohort1 <- rbind(cohort_guideline_N_model_N, cohort_guideline_Y_model_N, cohort_guideline_N_model_Y, cohort_guideline_Y_model_Y)

for (i in levels(as.factor(cohort1$subgp))) {
  
  print(i)
  
  for (m in levels(as.factor(cohort1$studydrug2))) {
    
    df_name <- paste0("df_", m)
    
    fit <- survfit(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ 1, 
                   data = cohort1 %>% filter(subgp == i & studydrug2 == m), 
                   weights = cohort1 %>% filter(subgp == i & studydrug2 == m) %>% .$overlap_bygroup)
    
    df <- data.frame(time = fit$time, surv = fit$surv, std.err = fit$std.err)
    
    assign(df_name, df)
  }
  
  df_diff_name <- paste0("df_diff_", i)
  
  df_diff <- merge(df_SGLT2i, `df_DPP4i/SU`, by = "time", suffixes = c("_SGLT2i", "_DPP4i/SU"))
  
  df_diff <- df_diff %>% mutate(
    difference = surv_SGLT2i - `surv_DPP4i/SU`,
    se_difference = sqrt(std.err_SGLT2i^2 + `std.err_DPP4i/SU`^2),
    lower_ci = difference - 1.96 * se_difference,
    upper_ci = difference + 1.96 * se_difference,
    subgp = i,
    smoothed_diff = loess(difference ~ time, data = df_diff, span = 1)$fitted,
    smoothed_lower_ci = loess(lower_ci ~ time, data = df_diff, span = 1)$fitted,
    smoothed_upper_ci = loess(upper_ci ~ time, data = df_diff, span = 1)$fitted
  )
  
  assign(df_diff_name, df_diff)
}

df_combined <- rbind(df_diff_guideline_N_model_N, df_diff_guideline_Y_model_N, df_diff_guideline_N_model_Y, df_diff_guideline_Y_model_Y) 
df_combined <- df_combined %>% mutate(
  subgp_label = ifelse(subgp == "guideline_N_model_N", paste0("Albuminuria <3mg/mmol, model does not recommend (n=", round(nrow(cohort_guideline_N_model_N)/n.imp), ")"), 
                       ifelse(subgp == "guideline_Y_model_N", paste0("Albuminuria ≥3mg/mmol, model does not recommend (n=", round(nrow(cohort_guideline_Y_model_N)/n.imp), ")"),
                              ifelse(subgp == "guideline_Y_model_Y", paste0("Model recommends SGLT2i, albuminuria ≥3mg/mmol (n=", round(nrow(cohort_guideline_Y_model_Y)/n.imp), ")"),
                                     paste0("Model recommends SGLT2i, albuminuria <3mg/mmol (n=", round(nrow(cohort_guideline_N_model_Y)/n.imp), ")")
                              ))
  ),
  subgp_label = factor(subgp_label, levels = c(
    paste0("Model recommends SGLT2i, albuminuria ≥3mg/mmol (n=", round(nrow(cohort_guideline_Y_model_Y)/n.imp), ")"),
    paste0("Model recommends SGLT2i, albuminuria <3mg/mmol (n=", round(nrow(cohort_guideline_N_model_Y)/n.imp), ")"),
    paste0("Albuminuria ≥3mg/mmol, model does not recommend (n=", round(nrow(cohort_guideline_Y_model_N)/n.imp), ")"),    
    paste0("Albuminuria <3mg/mmol, model does not recommend (n=", round(nrow(cohort_guideline_N_model_N)/n.imp), ")")
  )))



# plot
ci_plot <- df_combined %>% ggplot(aes(x = time, y = smoothed_diff*100, color = subgp_label)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = smoothed_lower_ci*100, ymax = smoothed_upper_ci*100, fill = subgp_label), alpha = 0.2, color = NA) +
  labs(#title = "Kidney protection benefit by SGLT2i treatment recommendation",
       x = "Time (years)",
       y = "Absolute risk reduction in kidney disease progression") +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("#D55E00", "#E69F00", "#56B4E9", "#0072B2")) +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#56B4E9", "#0072B2")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(legend.position = c(0.44, 0.875), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) +
  coord_cartesian(ylim = c(0,3.5), xlim = c(0, 4.8))


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_arr_by_treatment_recommendation.tiff"), width=7.5, height=6, units = "in", res=800)
ci_plot
dev.off()

arr_data <- data.frame(
  stratum = c(    
    "Albuminuria <3mg/mmol, model does not recommend",
    "Albuminuria ≥3mg/mmol, model does not recommend",
    "Model recommends SGLT2i, albuminuria <3mg/mmol",
    "Model recommends SGLT2i, albuminuria ≥3mg/mmol"
  ),
  ARR = c(
    arr_guideline_N_model_N,
    arr_guideline_Y_model_N,
    arr_guideline_N_model_Y,
    arr_guideline_Y_model_Y
  ),
  SE = c(
    se_guideline_N_model_N_ARR,
    se_guideline_Y_model_N_ARR,
    se_guideline_N_model_Y_ARR,
    se_guideline_Y_model_Y_ARR
  )
)

arr_data <- arr_data %>% mutate(
  stratum = factor(stratum, levels = c(    
    "Albuminuria <3mg/mmol, model does not recommend",
    "Albuminuria ≥3mg/mmol, model does not recommend",
    "Model recommends SGLT2i, albuminuria <3mg/mmol",
    "Model recommends SGLT2i, albuminuria ≥3mg/mmol"
  )),
  lower_bound = ARR - SE,
  upper_bound = ARR + SE)


bar_plot <- ggplot(arr_data, aes(x = stratum, y = ARR, fill = stratum)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +  
  labs(
    #   title = "5-year observed kidney protection benefit",
    x = "",
    y = "5-year absolute risk reduction in kidney disease progression         "
  ) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2, color = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_bw(base_size = 16) +  # Adjust font size
  theme(
    axis.text.x = element_text(hjust = 1), 
    legend.position = "none",
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
  ) +
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#E69F00", "#D55E00"))  + 
  coord_flip()

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_arr_barplot_by_treatment_recommendation.tiff"), width=7.5, height=6, units = "in", res=800)
bar_plot
dev.off()

# 
# ## observed data at 5 years:
# 
# #guidelines
# ddist <- datadist(cohort_guideline_Y);options(datadist='ddist')
# model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_Y, x=TRUE, y=TRUE, surv=TRUE)
# obs_guideline_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
# obs_guideline_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv         
# arr_guideline_Y <- (obs_guideline_Y_SGLT2-obs_guideline_Y_DPP4)*100
# print(paste0(c("5-year ARR - guidelines recommend treating: ", sprintf("%.2f", arr_guideline_Y))))
# 
# ddist <- datadist(cohort_guideline_N);options(datadist='ddist')
# model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_guideline_N, x=TRUE, y=TRUE, surv=TRUE)
# obs_guideline_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
# obs_guideline_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv           
# arr_guideline_N <- (obs_guideline_N_SGLT2-obs_guideline_N_DPP4)*100
# print(paste0(c("5-year ARR - guidelines do not recommend treating: ", sprintf("%.2f", arr_guideline_N))))
# 
# #strategy A
# ddist <- datadist(cohort_model1_Y);options(datadist='ddist')
# model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model1_Y, x=TRUE, y=TRUE, surv=TRUE)
# obs_model1_Y_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
# obs_model1_Y_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv
# arr_model1_Y <- (obs_model1_Y_SGLT2-obs_model1_Y_DPP4)*100
# print(paste0(c("5-year ARR - model recommends treating: ", sprintf("%.2f", arr_model1_Y))))
# 
# 
# ddist <- datadist(cohort_model1_N);options(datadist='ddist')
# model <- cph(Surv(ckd_egfr50_5y_censtime_yrs, ckd_egfr50_5y_censvar) ~ studydrug2, weights=overlap_bygroup, data=cohort_model1_N, x=TRUE, y=TRUE, surv=TRUE)
# obs_model1_N_SGLT2 <- survest(model, newdata=expand.grid(studydrug2="SGLT2i"), times=5)$surv
# obs_model1_N_DPP4 <- survest(model, newdata=expand.grid(studydrug2="DPP4i/SU"), times=5)$surv     
# arr_model1_N <- (obs_model1_N_SGLT2-obs_model1_N_DPP4)*100
# print(paste0(c("5-year ARR - model does not recommend treating: ", sprintf("%.2f", arr_model1_N))))
# 
# 
# 
# options(datadist=NULL)

############################5 NUMBERS TREATED AND EVENTS AVOIDED################################################################

### predicted data at 3 years:
#No treatment
describe(cohort$ckdpc_50egfr_score) #
#estimated number of events with no one treated
ckd.notx <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score/100)) 

#Treat all
describe(1-cohort$ckdpc_50egfr_survival_sglt2i) #
#estimated number of events with everyone treated
ckd.tx_all <- round(nrow(cohort)*mean(1-cohort$ckdpc_50egfr_survival_sglt2i)) 


#Treat as per guideline recommendations
describe(cohort$treat_guideline)
cohort <- cohort %>% mutate(ckdpc_50egfr_score.applied.guideline = 
                              ifelse(treat_guideline == F ,
                                     ckdpc_50egfr_score, 100*(1-cohort$ckdpc_50egfr_survival_sglt2i)))
#estimated number of events with treatment as per guidelines
ckd.tx.guideline <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score.applied.guideline/100)) 


#Treat as per strategy A
describe(cohort$treat_model1)
cohort <- cohort %>% mutate(ckdpc_50egfr_score.applied.model1 = 
                              ifelse(treat_model1 == F, 
                                     ckdpc_50egfr_score, 100*(1-cohort$ckdpc_50egfr_survival_sglt2i)))
#estimated number of events with treatment as per model 1
ckd.tx.model1 <- round(nrow(cohort)*mean(cohort$ckdpc_50egfr_score.applied.model1/100)) 



print(paste0(c("Number of people treated if no one treated: 0 (0%)")))
print(paste0(c("Number of events if no one treated: ", round(ckd.notx/n.imp), " (", round(100*ckd.notx/nrow(cohort), 1), "%)"), collapse = ""))

print(paste0(c("Number of people treated with treat-everyone strategy: ", round(nrow(cohort)/n.imp), " (", round(nrow(cohort)/nrow(cohort)*100,1), "%)"), collapse = ""))
print(paste0(c("Number of events with treat-everyone strategy: ", round(ckd.tx_all/n.imp), " (", round(100*(ckd.tx_all/nrow(cohort)),1), "%)"), collapse = ""))
print(paste0(c("Number of events avoided with treat-everyone strategy: ", round((ckd.notx - ckd.tx_all)/n.imp), " (", round((ckd.tx_all-ckd.notx)/(ckd.tx_all-ckd.notx)*100,1), "%)"), collapse = ""))
print(paste0(c("NNT with treat-everyone strategy: ", round(1/(mean(cohort$ckdpc_50egfr_score/100)-mean(1-cohort$ckdpc_50egfr_survival_sglt2i)))), collapse = ""))

print(paste0(c("Number of people treated with guideline treatment strategy: ", round(nrow(cohort[cohort$treat_guideline == T,])/n.imp), " (", round(nrow(cohort[cohort$treat_guideline == T,])/nrow(cohort)*100,1), "%)"), collapse = ""))
print(paste0(c("Number of events with guideline treatment strategy: ", round(ckd.tx.guideline/n.imp), " (", round(100*(ckd.tx.guideline/nrow(cohort)),1), "%)"), collapse = ""))
print(paste0(c("Number of events avoided with guideline treatment strategy: ", round(abs(ckd.tx.guideline-ckd.notx)/n.imp), " (", round((ckd.tx.guideline-ckd.notx)/(ckd.tx_all-ckd.notx)*100,1), "%)"), collapse = ""))
print(paste0(c("NNT with guideline treatment strategy: ", 
               round(1/(mean(cohort[cohort$treat_guideline == T,]$ckdpc_50egfr_score/100) - 
                          mean(cohort[cohort$treat_guideline == T,]$ckdpc_50egfr_score.applied.guideline/100)))), collapse = ""))

print(paste0(c("Number of people treated with pARR strategy: ", round(nrow(cohort[cohort$treat_model1 == T,])/n.imp), " (", round(nrow(cohort[cohort$treat_model1 == T,])/nrow(cohort)*100,1), "%)"), collapse = ""))
print(paste0(c("Number of events with pARR strategy: ", round(ckd.tx.model1/n.imp), " (", round(100*(ckd.tx.model1/nrow(cohort)),1), "%)"), collapse = ""))
print(paste0(c("Number of events avoided with pARR strategy: ", round(abs(ckd.tx.model1-ckd.notx)/n.imp), " (", round((ckd.tx.model1-ckd.notx)/(ckd.tx_all-ckd.notx)*100,1), "%)"), collapse = ""))
print(paste0(c("NNT with pARR strategy: ", 
               round(1/(mean(cohort[cohort$treat_model1 == T,]$ckdpc_50egfr_score/100) - 
                          mean(cohort[cohort$treat_model1 == T,]$ckdpc_50egfr_score.applied.model1/100)))), collapse = ""))



############################6 DECISION CURVE ANALYSIS################################################################
cohort <- cohort %>% mutate(ckdpc_50egfr_score_risk = ckdpc_50egfr_score/100) 
dca_data <- dca(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ treat_guideline + ckdpc_50egfr_score_risk,
      thresholds = seq(0, 0.20, by = 0.005),
      time = 3,
      data=cohort)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Output/")
tiff(paste0(today, "_decision_curve_analysis.tiff"), width=6, height=5, units = "in", res=800) 
as_tibble(dca_data) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label, linetype = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", 
              span = 0.2, size = 1.25) +
  coord_cartesian(ylim = c(-0.00105984276971715, 0.0105984276971715
  )) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Risk Tolerance", y = "Net Utility", color = "") +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "black", "grey40", "#0072B2", "#E69F00"
    ),
    labels = c("Treat all", 
               "Treat none", 
               "Treat if uACR ≥3mg/mmol", 
               "Treat according to pARR")) +
  scale_linetype_manual(
    values = c(
      "solid", 
      "solid", 
      "longdash", 
      "solid"),
    labels = c("Treat all", 
               "Treat none", 
               "Treat if uACR ≥3mg/mmol", 
               "Treat according to pARR")) + 
  guides(color = guide_legend("Treatment strategy"), 
         linetype = guide_legend("Treatment strategy")) +
  theme(legend.position = c(0.75, 0.8)) +
  coord_cartesian(xlim = c(0.0012,0.03), ylim = c(0,0.01))
dev.off()

# #DCA with static thresholds (pARR thresholds A + B used in table)
# dca_data2 <- dca(Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ treat_guideline + treat_model1 + treat_model2,
#                  data = cohort,
#                  thresholds = seq(0, 0.20, by = 0.01),
#                  time = 3)
# 
# #plot
# as_tibble(dca_data2) %>%
#   dplyr::filter(!is.na(net_benefit)) %>%
#   ggplot(aes(x = threshold, y = net_benefit, color = label, linetype = label)) +
#   stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", 
#               span = 0.2, size = 1.25) +
#   coord_cartesian(ylim = c(-0.00105984276971715, 0.0105984276971715
#   )) +
#   scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
#   labs(x = "Threshold Probability", y = "Net Utility", color = "") +
#   theme_bw() +
#   scale_color_manual(
#     values = c(
#       "black", "grey40", "#0072B2", "#E69F00", "#E69F00" 
#     ),
#     labels = c("Treat all", 
#                "Treat none", 
#                "Treat if uACR ≥3mg/mmol", 
#                "Treat if pARR above threshold A", 
#                "Treat if pARR above threshold B")) +
#   scale_linetype_manual(
#     values = c(
#       "solid", 
#       "solid", 
#       "longdash", 
#       "solid", 
#       "dashed"),
#     labels = c("Treat all", 
#                "Treat none", 
#                "Treat if uACR ≥3mg/mmol", 
#                "Treat if pARR above threshold A", 
#                "Treat if pARR above threshold B")) + 
#   guides(color = guide_legend("Treatment strategy"), 
#          linetype = guide_legend("Treatment strategy")) +
#   theme(legend.position = c(0.75, 0.8)) +
#   coord_cartesian(xlim = c(0.0012,0.03), ylim = c(0,0.01))
############################7 HR BY pARR THRESHOLD################################################################


#create empty data frame to which we can append the hazard ratios once calculated
subgroup_parr_SGLT2ivsDPP4iSU_hrs <- subgroup_parr_hrs <- data.frame()

#ensure treat_model1 recognised as logical variable:
cohort <- cohort %>% mutate(treat_model1 = as.logical(treat_model1))

p_value_interaction <- setNames(vector("numeric", length(outcomes)), outcomes)

## analyses stratified by subgroup
for (k in outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  
  
  # calculate number of subjects in each group
  count <- cohort %>%
    group_by(studydrug2,treat_model1) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort %>%
    group_by(studydrug2,treat_model1) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- cohort %>%
    group_by(studydrug2,treat_model1) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug2, treat_model1, events) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*treat_model1"))
  
  f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*treat_model1 + ", paste(covariates, collapse=" + ")))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.lowparr.unadj <- SE.lowparr.unadj <-
    COEFS.highparr.unadj <- SE.highparr.unadj <-
    # for the adjusted survival models
    COEFS.lowparr.adj <- SE.lowparr.adj <-
    COEFS.highparr.adj <- SE.highparr.adj <-
    # for the overlap-weighted survival models
    COEFS.lowparr.ow <- SE.lowparr.ow <-
    COEFS.highparr.ow <- SE.highparr.ow <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f2, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.lowparr.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"]
    SE.lowparr.unadj[i] <- sqrt(fit.unadj$var[1,1])
    
    COEFS.highparr.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:treat_model1TRUE"]
    SE.highparr.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var),nrow(fit.unadj$var)]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted2, cohort[cohort$.imp == i,])
    
    COEFS.lowparr.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"]
    SE.lowparr.adj[i] <- sqrt(fit.adj$var[1,1])
    
    COEFS.highparr.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:treat_model1TRUE"]
    SE.highparr.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var),nrow(fit.adj$var)]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)])
    
    #overlap-weighted analyses
    fit.ow <- coxph(f_adjusted2, cohort[cohort$.imp == i,], weights = overlap2)
    
    COEFS.lowparr.ow[i] <- fit.ow$coefficients["studydrug2SGLT2i"]
    SE.lowparr.ow[i] <- sqrt(fit.ow$var[1,1])
    
    COEFS.highparr.ow[i] <- fit.ow$coefficients["studydrug2SGLT2i"] + fit.ow$coefficients["studydrug2SGLT2i:treat_model1TRUE"]
    SE.highparr.ow[i] <- sqrt(abs(fit.ow$var[1]) + abs(fit.ow$var[nrow(fit.ow$var),nrow(fit.ow$var)]) + 2 * vcov(fit.ow)[1,nrow(fit.ow$var)])
    
    
    if (i == n.imp) {
      f_adjusted3 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + treat_model1 + ", paste(covariates, collapse=" + ")))
      fit.no_interaction <- coxph(f_adjusted3, cohort[cohort$.imp == i,])
      
      loglikelihood_test <- anova(fit.no_interaction, fit.adj, test = "Chisq")
      p_value_interaction[k] <- loglikelihood_test$`Pr(>|Chi|)`[2]      
    }
    
  }
  
  # pool hazard ratios
  unadjusted_lowparr <- pool.rubin.HR(COEFS.lowparr.unadj, SE.lowparr.unadj, n.imp)
  adjusted_lowparr <- pool.rubin.HR(COEFS.lowparr.adj, SE.lowparr.adj, n.imp)
  ow_lowparr <- pool.rubin.HR(COEFS.lowparr.ow, SE.lowparr.ow, n.imp)
  
  unadjusted_highparr <- pool.rubin.HR(COEFS.highparr.unadj, SE.highparr.unadj, n.imp)
  adjusted_highparr <- pool.rubin.HR(COEFS.highparr.adj, SE.highparr.adj, n.imp)
  ow_highparr <- pool.rubin.HR(COEFS.highparr.ow, SE.highparr.ow, n.imp)
  
  # save pooled HR and 95% confidence interval
  unadjusted_lowparr_string <- paste0(sprintf("%.2f", round(unadjusted_lowparr[1], 2)), " (", sprintf("%.2f", round(unadjusted_lowparr[2], 2)), ", ", sprintf("%.2f", round(unadjusted_lowparr[3], 2)), ")")
  adjusted_lowparr_string <- paste0(sprintf("%.2f", round(adjusted_lowparr[1], 2)), " (", sprintf("%.2f", round(adjusted_lowparr[2], 2)), ", ", sprintf("%.2f", round(adjusted_lowparr[3], 2)), ")")
  ow_lowparr_string <- paste0(sprintf("%.2f", round(ow_lowparr[1], 2)), " (", sprintf("%.2f", round(ow_lowparr[2], 2)), ", ", sprintf("%.2f", round(ow_lowparr[3], 2)), ")")
  
  unadjusted_highparr_string <- paste0(sprintf("%.2f", round(unadjusted_highparr[1], 2)), " (", sprintf("%.2f", round(unadjusted_highparr[2], 2)), ", ", sprintf("%.2f", round(unadjusted_highparr[3], 2)), ")")
  adjusted_highparr_string <- paste0(sprintf("%.2f", round(adjusted_highparr[1], 2)), " (", sprintf("%.2f", round(adjusted_highparr[2], 2)), ", ", sprintf("%.2f", round(adjusted_highparr[3], 2)), ")")
  ow_highparr_string <- paste0(sprintf("%.2f", round(ow_highparr[1], 2)), " (", sprintf("%.2f", round(ow_highparr[2], 2)), ", ", sprintf("%.2f", round(ow_highparr[3], 2)), ")")
  
  # combine in dataframe that we can tabulate
  presegfr_lowparr_hr <- cbind(outcome=k, count[1,c(2:3)], followup[1,c(2:3)], events[1,c(2:3)], contrast = "pARR below threshold",
                             unadjusted=unadjusted_lowparr_string, adjusted=adjusted_lowparr_string, ow=ow_lowparr_string
  )
  presegfr_highparr_hr <- cbind(outcome=k, count[2,c(2:3)], followup[2,c(2:3)], events[2,c(2:3)], contrast = "pARR above threshold",
                                unadjusted=unadjusted_highparr_string, adjusted=adjusted_highparr_string, ow=ow_highparr_string
  )
  
  outcome_subgroup_parr_SGLT2ivsDPP4iSU_hrs <- rbind(presegfr_lowparr_hr, presegfr_highparr_hr)
  
  temp <- rbind(
    cbind(outcome = k, contrast = "pARR below threshold", analysis = "Adjusted",
          HR = adjusted_lowparr[1], LB = adjusted_lowparr[2], UB = adjusted_lowparr[3], string = adjusted_lowparr_string),
    cbind(outcome = k, contrast = "pARR above threshold", analysis = "Adjusted",
          HR = adjusted_highparr[1], LB = adjusted_highparr[2], UB = adjusted_highparr[3], string = adjusted_highparr_string),
    cbind(outcome = k, contrast = "pARR below threshold", analysis = "Overlap-weighted",
          HR = ow_lowparr[1], LB = ow_lowparr[2], UB = ow_lowparr[3], string = ow_lowparr_string),
    cbind(outcome = k, contrast = "pARR above threshold", analysis = "Overlap-weighted",
          HR = ow_highparr[1], LB = ow_highparr[2], UB = ow_highparr[3], string = ow_highparr_string)
  )
  
  subgroup_parr_hrs <- rbind(subgroup_parr_hrs, temp)
  subgroup_parr_SGLT2ivsDPP4iSU_hrs <- rbind(subgroup_parr_SGLT2ivsDPP4iSU_hrs, outcome_subgroup_parr_SGLT2ivsDPP4iSU_hrs)
  
}


subgroup_parr_hrs <- subgroup_parr_hrs %>% left_join(subgroup_parr_SGLT2ivsDPP4iSU_hrs %>% select(-c(adjusted, unadjusted)), by = c("outcome", "contrast"))

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
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
save(subgroup_parr_hrs, file=paste0(today, "_subgroup_parr_hrs.Rda"))
save(subgroup_parr_SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_subgroup_parr_SGLT2ivsDPP4iSU_hrs.Rda"))

