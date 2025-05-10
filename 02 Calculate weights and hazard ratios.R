## In this script we aim to see whether the relative treatment effect of SGLT2i vs control is similar
## in the literature and in the CPRD cohort
## we will compare hazard ratios for SGLT2i vs DPP4i/SU in CPRD to those from trials.

## for reference, the trial meta-analysis RR reported in Lancet (Nuffield Group 2022) is:
# HR 0.62 (0.56 - 0.68)

## Here we will be calculating HRs for SGLT2i and DPP4i vs SU for:
## primary outcome kidney disease progression: eGFR decline of ≥50% from baseline or onset of CKD stage 5
## secondary outcomes: progression to macroalbuminuria, all-cause mortality, dka, amputation, and mycotic genital infections

## Contents:
# 0 setup
# 1 calculate weights (IPTW, overlap weights)
# 2 calculate hazard ratios
# 3 subgroup analyses (by albuminuria status)
# 4 store hazard ratios

############################0 SETUP################################################################
# 0 Setup
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
load(paste0(today, "_t2d_ckdpc_imputed_data.Rda"))

############################1 CALCULATE WEIGHTS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in them, then stack them again at the end
temp$IPTW <- temp$overlap <- 
temp$IPTW2 <- temp$overlap2 <- NA

# take non-imputed dataset to append the imputed datasets to later on
temp2 <- temp[temp$.imp == 0,]

# propensity score formula
ps.formula <- formula(paste("studydrug ~ ", paste(covariates, collapse=" + ")))

ps.formula2 <- formula(paste("studydrug2 ~ ", paste(covariates, collapse=" + ")))


#calculate weights in each imputed dataset
for (i in 1:n.imp) {
  
  print(paste0("Calculating weights for imputed dataset number ", i))
  
  imp.x <- temp[temp$.imp == i,] %>% as.data.frame()
  
  w.overlap <- SumStat(ps.formula=ps.formula, 
                       data = imp.x, 
                       weight = c("IPW", "overlap"))
  
  w.overlap2 <- SumStat(ps.formula=ps.formula2, 
                       data = imp.x, 
                       weight = c("IPW", "overlap"))
  
  # truncate IPTW at 2nd and 98% percentile
  w.overlap$ps.weights$IPW <- ifelse(w.overlap$ps.weights$IPW < quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1], 
                                     quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                     ifelse(
                                       w.overlap$ps.weights$IPW > quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                       quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                       w.overlap$ps.weights$IPW
                                     ))
  
  
  w.overlap2$ps.weights$IPW <- ifelse(w.overlap2$ps.weights$IPW < quantile(w.overlap2$ps.weights$IPW, probs = c(0.02, 0.98))[1], 
                                      quantile(w.overlap2$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                      ifelse(
                                        w.overlap2$ps.weights$IPW > quantile(w.overlap2$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                        quantile(w.overlap2$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                        w.overlap2$ps.weights$IPW
                                      ))
  
  ## Add weights to data frame
  weights <- w.overlap$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
  imp.x$IPTW <- weights$IPW
  imp.x$overlap <- weights$overlap
  
  weights2 <- w.overlap2$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
  imp.x$IPTW2 <- weights2$IPW
  imp.x$overlap2 <- weights2$overlap
  
  # append imputed dataset with weights to dataframe with the original dataset + appended imputed datasets
  temp2 <- rbind(temp2, imp.x)
  
  rm(imp.x)
  rm(weights)
  
}

temp <- temp2
rm(temp2)

#plot SMD in weighted / unweighted populations (love plot - shows balance by variable)
plot(w.overlap2)
#histogram of propensity scores by treatment group (I have amended the code from the package to fit with our colours)
plot.hist.ps.2drugs<-function(x, weighted.var=TRUE, threshold=0.1, metric="ASD", breaks=50,...){
#get object info
m<-length(names(x$ps.weights)[-c(1,2)])+1
zname<-names(x$ps.weights)[1]
Z<-unlist(x$ps.weights[zname])
ncate<-length(unique(Z))
metric<-toupper(metric)

#extract original treatment labels
dic0<-rep(NA,ncate)
for (i in 1:ncate){
  dic0[i]<-as.character(x$ps.weights[which(x$ps.weights$zindex==i)[1],1])
}
e<-x$propensity[,2]
z<-x$ps.weights$zindex
#enough room for 8 character-strings in legend
par(mar=c(5,4,4,10.1),xpd=TRUE)
he1<-hist(e[z==1],breaks=breaks,plot = FALSE)
ylim1<-max(he1$counts)+0.5
he2<-hist(e[z==2],breaks=breaks,plot = FALSE)
ylim2<-max(he2$counts+0.5)
ylims<-max(c(ylim1,ylim2))

hist(e[z==1],breaks=breaks,col="#0072B2",border="black",bty='L',
     main=NULL,ylab=NULL,xlim=c(max(min(e)-0.2,0),min(max(e)+0.2,1)),xlab="Estimated propensity score",
     freq=T,cex.lab = 1.5, cex.axis = 1.5 ,cex.main = 2, cex = 2,bty='L',ylim = c(0,ylims))
hist(e[z==2],breaks=breaks,add=TRUE,col=rgb(230/255, 159/255, 0, alpha = 0.9),freq=T)

legend("right",inset=c(-0.2, 0), title="",legend=dic0,col=c("#0072B2","#E69F00"),
       lty=1,lwd=1.5,bty='n',cex=1.5,seg.len = 0.65,xjust = 1)
}
plot.hist.ps.2drugs(w.overlap2)

summary(w.overlap2, weighted.var = TRUE, metric = "ASD")

# save dataset with weights so this can be used in subsequent scripts
cohort <- temp %>% filter(!.imp == 0)
rm(temp)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/Processed data/")
save(cohort, file=paste0(today, "_t2d_ckdpc_imputed_data_withweights.Rda"))
#load(paste0(today, "_t2d_ckdpc_imputed_data_withweights.Rda"))
############################2 CALCULATE HAZARD RATIOS################################################################

## 2 calculate hazard ratios (unadjusted, adjusted, weighted) and n events per study drug

#create empty data frame to which we can append the hazard ratios once calculated
all_sglt2i_hrs <- 
  all_dpp4i_hrs <-
  all_SGLT2ivsDPP4i_hrs <- 
  all_hrs <- 
  data.frame()

for (k in outcomes_per_drugclass) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  # calculate number of subjects in each group
  count <- cohort %>%
    group_by(studydrug) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort %>%
    group_by(studydrug) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- cohort %>%
    group_by(studydrug) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug, events) %>%
    pivot_wider(names_from=studydrug,
                names_glue="{studydrug}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug"))
  
  f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug + ", paste(covariates, collapse=" + ")))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.SGLT2i.unadj <- SE.SGLT2i.unadj <-
    COEFS.DPP4i.unadj <- SE.DPP4i.unadj <-
    COEFS.SGLT2ivsDPP4i.unadj <- SE.SGLT2ivsDPP4i.unadj <-
  # for the adjusted survival models
    COEFS.SGLT2i.adj <- SE.SGLT2i.adj <-
    COEFS.DPP4i.adj <- SE.DPP4i.adj <-
    COEFS.SGLT2ivsDPP4i.adj <- SE.SGLT2ivsDPP4i.adj <-
  # for the overlap weighted survival models
    COEFS.SGLT2i.ow <- SE.SGLT2i.ow <-
    COEFS.DPP4i.ow <- SE.DPP4i.ow <-
    COEFS.SGLT2ivsDPP4i.ow <- SE.SGLT2ivsDPP4i.ow <-
  # for the inverse probability of treatment weighted survival models
    COEFS.SGLT2i.iptw <- SE.SGLT2i.iptw <-
    COEFS.DPP4i.iptw <- SE.DPP4i.iptw <-
    COEFS.SGLT2ivsDPP4i.iptw <- SE.SGLT2ivsDPP4i.iptw <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.DPP4i.unadj[i] <- fit.unadj$coefficients[1]
    COEFS.SGLT2i.unadj[i] <- fit.unadj$coefficients[2]
    SE.DPP4i.unadj[i] <- sqrt(fit.unadj$var[1,1])
    SE.SGLT2i.unadj[i] <- sqrt(fit.unadj$var[2,2])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted, cohort[cohort$.imp == i,])
    
    COEFS.DPP4i.adj[i] <- fit.adj$coefficients[1]
    COEFS.SGLT2i.adj[i] <- fit.adj$coefficients[2]
    SE.DPP4i.adj[i] <- sqrt(fit.adj$var[1,1])
    SE.SGLT2i.adj[i] <- sqrt(fit.adj$var[2,2])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = overlap)
    
    COEFS.DPP4i.ow[i] <- fit.ow$coefficients[1]
    COEFS.SGLT2i.ow[i] <- fit.ow$coefficients[2]
    SE.DPP4i.ow[i] <- sqrt(fit.ow$var[1,1])
    SE.SGLT2i.ow[i] <- sqrt(fit.ow$var[2,2])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = IPTW)
    
    COEFS.DPP4i.iptw[i] <- fit.iptw$coefficients[1]
    COEFS.SGLT2i.iptw[i] <- fit.iptw$coefficients[2]
    SE.DPP4i.iptw[i] <- sqrt(fit.iptw$var[1,1])
    SE.SGLT2i.iptw[i] <- sqrt(fit.iptw$var[2,2])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
    
    #calculate HRs for SGLT2i vs DPP4
    cohort$studydrug <- relevel(factor(cohort$studydrug), ref = "DPP4i")
    
    #unadjusted analyses first
    fit.unadj <- coxph(f, cohort[cohort$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.SGLT2ivsDPP4i.unadj[i] <- fit.unadj$coefficients[2]
    SE.SGLT2ivsDPP4i.unadj[i] <- sqrt(fit.unadj$var[2,2])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted, cohort[cohort$.imp == i,])
    
    COEFS.SGLT2ivsDPP4i.adj[i] <- fit.adj$coefficients[2]
    SE.SGLT2ivsDPP4i.adj[i] <- sqrt(fit.adj$var[2,2])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = overlap)
    
    COEFS.SGLT2ivsDPP4i.ow[i] <- fit.ow$coefficients[2]
    SE.SGLT2ivsDPP4i.ow[i] <- sqrt(fit.ow$var[2,2])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted, cohort[cohort$.imp ==i,], weights = IPTW)
    
    COEFS.SGLT2ivsDPP4i.iptw[i] <- fit.iptw$coefficients[2]
    SE.SGLT2ivsDPP4i.iptw[i] <- sqrt(fit.iptw$var[2,2])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
    
    cohort$studydrug <- relevel(factor(cohort$studydrug), ref = "SU")
    
  }
  
  # pool hazard ratios
  unadjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.unadj, SE.SGLT2i.unadj, n.imp)
  unadjusted_DPP4i <- pool.rubin.HR(COEFS.DPP4i.unadj, SE.DPP4i.unadj, n.imp)
  unadjusted_SGLT2ivsDPP4i <- pool.rubin.HR(COEFS.SGLT2ivsDPP4i.unadj, SE.SGLT2ivsDPP4i.unadj, n.imp)
  
  adjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.adj, SE.SGLT2i.adj, n.imp)
  adjusted_DPP4i <- pool.rubin.HR(COEFS.DPP4i.adj, SE.DPP4i.adj, n.imp)  
  adjusted_SGLT2ivsDPP4i <- pool.rubin.HR(COEFS.SGLT2ivsDPP4i.adj, SE.SGLT2ivsDPP4i.adj, n.imp)
  
  overlapweighted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.ow, SE.SGLT2i.ow, n.imp)
  overlapweighted_DPP4i <- pool.rubin.HR(COEFS.DPP4i.ow, SE.DPP4i.ow, n.imp)
  overlapweighted_SGLT2ivsDPP4i <- pool.rubin.HR(COEFS.SGLT2ivsDPP4i.ow, SE.SGLT2ivsDPP4i.ow, n.imp)
  
  iptw_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.iptw, SE.SGLT2i.iptw, n.imp)
  iptw_DPP4i <- pool.rubin.HR(COEFS.DPP4i.iptw, SE.DPP4i.iptw, n.imp)
  iptw_SGLT2ivsDPP4i <- pool.rubin.HR(COEFS.SGLT2ivsDPP4i.iptw, SE.SGLT2ivsDPP4i.iptw, n.imp)
  
  
  # save pooled HR and 95% confidence interval
  unadjusted_SGLT2i_string <- paste0(sprintf("%.2f", round(unadjusted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(unadjusted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(unadjusted_SGLT2i[3], 2)), ")")
  unadjusted_DPP4i_string <- paste0(sprintf("%.2f", round(unadjusted_DPP4i[1], 2)), " (", sprintf("%.2f", round(unadjusted_DPP4i[2], 2)), ", ", sprintf("%.2f", round(unadjusted_DPP4i[3], 2)), ")")
  unadjusted_SGLT2ivsDPP4i_string <- paste0(sprintf("%.2f", round(unadjusted_SGLT2ivsDPP4i[1], 2)), " (", sprintf("%.2f", round(unadjusted_SGLT2ivsDPP4i[2], 2)), ", ", sprintf("%.2f", round(unadjusted_SGLT2ivsDPP4i[3], 2)), ")")

  adjusted_SGLT2i_string <- paste0(sprintf("%.2f", round(adjusted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(adjusted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(adjusted_SGLT2i[3], 2)), ")")
  adjusted_DPP4i_string <- paste0(sprintf("%.2f", round(adjusted_DPP4i[1], 2)), " (", sprintf("%.2f", round(adjusted_DPP4i[2], 2)), ", ", sprintf("%.2f", round(adjusted_DPP4i[3], 2)), ")")
  adjusted_SGLT2ivsDPP4i_string <- paste0(sprintf("%.2f", round(adjusted_SGLT2ivsDPP4i[1], 2)), " (", sprintf("%.2f", round(adjusted_SGLT2ivsDPP4i[2], 2)), ", ", sprintf("%.2f", round(adjusted_SGLT2ivsDPP4i[3], 2)), ")")
  
  overlapweighted_SGLT2i_string <- paste0(sprintf("%.2f", round(overlapweighted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(overlapweighted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(overlapweighted_SGLT2i[3], 2)), ")")
  overlapweighted_DPP4i_string <- paste0(sprintf("%.2f", round(overlapweighted_DPP4i[1], 2)), " (", sprintf("%.2f", round(overlapweighted_DPP4i[2], 2)), ", ", sprintf("%.2f", round(overlapweighted_DPP4i[3], 2)), ")")
  overlapweighted_SGLT2ivsDPP4i_string <- paste0(sprintf("%.2f", round(overlapweighted_SGLT2ivsDPP4i[1], 2)), " (", sprintf("%.2f", round(overlapweighted_SGLT2ivsDPP4i[2], 2)), ", ", sprintf("%.2f", round(overlapweighted_SGLT2ivsDPP4i[3], 2)), ")")

  iptw_SGLT2i_string <- paste0(sprintf("%.2f", round(iptw_SGLT2i[1], 2)), " (", sprintf("%.2f", round(iptw_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(iptw_SGLT2i[3], 2)), ")")
  iptw_DPP4i_string <- paste0(sprintf("%.2f", round(iptw_DPP4i[1], 2)), " (", sprintf("%.2f", round(iptw_DPP4i[2], 2)), ", ", sprintf("%.2f", round(iptw_DPP4i[3], 2)), ")")
  iptw_SGLT2ivsDPP4i_string <- paste0(sprintf("%.2f", round(iptw_SGLT2ivsDPP4i[1], 2)), " (", sprintf("%.2f", round(iptw_SGLT2ivsDPP4i[2], 2)), ", ", sprintf("%.2f", round(iptw_SGLT2ivsDPP4i[3], 2)), ")")
  
  
  # combine in dataframe that we can tabulate
  outcome_sglt2i_hr <- cbind(outcome=k, 
                             count[c("SU_count", "SGLT2i_count")], 
                             followup[c("SU_followup", "SGLT2i_followup")], 
                             events[c("SU_events", "SGLT2i_events")], 
                            unadjusted_SGLT2i_string, 
                            adjusted_SGLT2i_string, 
                            overlapweighted_SGLT2i_string, 
                            iptw_SGLT2i_string)
  
  outcome_dpp4i_hr <- cbind(outcome=k, 
                            count[c("SU_count", "DPP4i_count")], 
                            followup[c("SU_followup", "DPP4i_followup")], 
                            events[c("SU_events", "DPP4i_events")], 
                            unadjusted_DPP4i_string, 
                            adjusted_DPP4i_string, 
                            overlapweighted_DPP4i_string, 
                            iptw_DPP4i_string)
  
  outcome_SGLT2ivsDPP4i_hr <- cbind(outcome=k, 
                                    count[c("DPP4i_count", "SGLT2i_count")], 
                                    followup[c("DPP4i_followup", "SGLT2i_followup")], 
                                    events[c("DPP4i_events", "SGLT2i_events")],  
                                    unadjusted_SGLT2ivsDPP4i_string, 
                                    adjusted_SGLT2ivsDPP4i_string, 
                                    overlapweighted_SGLT2ivsDPP4i_string, 
                                    iptw_SGLT2ivsDPP4i_string)
  
  all_sglt2i_hrs <- rbind(all_sglt2i_hrs, outcome_sglt2i_hr)
  all_dpp4i_hrs <- rbind(all_dpp4i_hrs, outcome_dpp4i_hr)
  all_SGLT2ivsDPP4i_hrs <- rbind(all_SGLT2ivsDPP4i_hrs, outcome_SGLT2ivsDPP4i_hr)
  
  outcome_hr <- rbind(
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "Unadjusted", 
          HR = unadjusted_SGLT2i[1], LB = unadjusted_SGLT2i[2], UB = unadjusted_SGLT2i[3], 
          string = unadjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "Adjusted", 
          HR = adjusted_SGLT2i[1], LB = adjusted_SGLT2i[2], UB = adjusted_SGLT2i[3], 
          string = adjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "Overlap-weighted", 
          HR = overlapweighted_SGLT2i[1], LB = overlapweighted_SGLT2i[2], UB = overlapweighted_SGLT2i[3], 
          string = overlapweighted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs SU", analysis = "IPTW", 
          HR = iptw_SGLT2i[1], LB = iptw_SGLT2i[2], UB = iptw_SGLT2i[3], 
          string = iptw_SGLT2i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "Unadjusted", 
          HR = unadjusted_DPP4i[1], LB = unadjusted_DPP4i[2], UB = unadjusted_DPP4i[3], 
          string = unadjusted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "Adjusted", 
          HR = adjusted_DPP4i[1], LB = adjusted_DPP4i[2], UB = adjusted_DPP4i[3], 
          string = adjusted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "Overlap-weighted", 
          HR = overlapweighted_DPP4i[1], LB = overlapweighted_DPP4i[2], UB = overlapweighted_DPP4i[3], 
          string = overlapweighted_DPP4i_string),
    cbind(outcome = k, contrast = "DPP4i vs SU", analysis = "IPTW", 
          HR = iptw_DPP4i[1], LB = iptw_DPP4i[2], UB = iptw_DPP4i[3], 
          string = iptw_DPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "Unadjusted", 
          HR = unadjusted_SGLT2ivsDPP4i[1], LB = unadjusted_SGLT2ivsDPP4i[2], UB = unadjusted_SGLT2ivsDPP4i[3], 
          string = unadjusted_SGLT2ivsDPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "Adjusted", 
          HR = adjusted_SGLT2ivsDPP4i[1], LB = adjusted_SGLT2ivsDPP4i[2], UB = adjusted_SGLT2ivsDPP4i[3], 
          string = adjusted_SGLT2ivsDPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "Overlap-weighted", 
          HR = overlapweighted_SGLT2ivsDPP4i[1], LB = overlapweighted_SGLT2ivsDPP4i[2], UB = overlapweighted_SGLT2ivsDPP4i[3], 
          string = overlapweighted_SGLT2ivsDPP4i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i", analysis = "IPTW", 
          HR = iptw_SGLT2ivsDPP4i[1], LB = iptw_SGLT2ivsDPP4i[2], UB = iptw_SGLT2ivsDPP4i[3], 
          string = iptw_SGLT2ivsDPP4i_string)
  )
  
  all_hrs <- rbind(all_hrs, outcome_hr)
  
}

# review HRs for different treatment contrasts
flextable(all_sglt2i_hrs)
flextable(all_dpp4i_hrs)
flextable(all_SGLT2ivsDPP4i_hrs) # need to compare pp outcomes as these censor DPP4i when starting SU and vice versa
# we can see that the DPP4i vs SU HRs are non-significant, and we can therefore pool these 2 groups together.

# store n/N for every drug class in all_hrs
all_hrs <- all_hrs %>% left_join (all_dpp4i_hrs %>% select(1:7)) %>% left_join(all_sglt2i_hrs %>% select(1,3,5,7))



## compare SGLT2i vs combined group of DPP4i/SU

#create empty data frame to which we can append the hazard ratios once calculated
all_SGLT2ivsDPP4iSU_hrs <- SGLT2ivsDPP4iSU_hrs <- data.frame()

## remove double overlapping entries for DPP4i and SU that overlap (take one only)
cohort <- cohort %>% group_by(.imp, patid) %>% filter(
  !duplicated(studydrug2)
) %>% ungroup()


for (k in outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  if (k == "macroalb") {
    temp <- cohort[cohort$uacr < 30,]
  } else {
    temp <- cohort
  }
  
  # calculate number of subjects in each group
  count <- temp[temp$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- temp[temp$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- temp[temp$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug2, events) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2"))
  
  f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", paste(covariates, collapse=" + ")))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.SGLT2i.unadj <- SE.SGLT2i.unadj <-
    # for the adjusted survival models
    COEFS.SGLT2i.adj <- SE.SGLT2i.adj <-
    # for the overlap weighted survival models
    COEFS.SGLT2i.ow <- SE.SGLT2i.ow <-
    # for the inverse probability of treatment weighted survival models
    COEFS.SGLT2i.iptw <- SE.SGLT2i.iptw <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f2, temp[temp$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.SGLT2i.unadj[i] <- fit.unadj$coefficients[1]
    SE.SGLT2i.unadj[i] <- sqrt(fit.unadj$var[1,1])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted2, temp[temp$.imp == i,])
    
    COEFS.SGLT2i.adj[i] <- fit.adj$coefficients[1]
    SE.SGLT2i.adj[i] <- sqrt(fit.adj$var[1,1])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted2, temp[temp$.imp ==i,], weights = overlap2)
    
    COEFS.SGLT2i.ow[i] <- fit.ow$coefficients[1]
    SE.SGLT2i.ow[i] <- sqrt(fit.ow$var[1,1])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted2, temp[temp$.imp ==i,], weights = IPTW2)
    
    COEFS.SGLT2i.iptw[i] <- fit.iptw$coefficients[1]
    SE.SGLT2i.iptw[i] <- sqrt(fit.iptw$var[1,1])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
    }
  
  # pool hazard ratios
  unadjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.unadj, SE.SGLT2i.unadj, n.imp)
  adjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.adj, SE.SGLT2i.adj, n.imp)
  overlapweighted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.ow, SE.SGLT2i.ow, n.imp)
  iptw_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.iptw, SE.SGLT2i.iptw, n.imp)
  
  # save pooled HR and 95% confidence interval
  unadjusted_SGLT2i_string <- paste0(sprintf("%.2f", round(unadjusted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(unadjusted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(unadjusted_SGLT2i[3], 2)), ")")
  adjusted_SGLT2i_string <- paste0(sprintf("%.2f", round(adjusted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(adjusted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(adjusted_SGLT2i[3], 2)), ")")
  overlapweighted_SGLT2i_string <- paste0(sprintf("%.2f", round(overlapweighted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(overlapweighted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(overlapweighted_SGLT2i[3], 2)), ")")
  iptw_SGLT2i_string <- paste0(sprintf("%.2f", round(iptw_SGLT2i[1], 2)), " (", sprintf("%.2f", round(iptw_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(iptw_SGLT2i[3], 2)), ")")
 
  # combine in dataframe that we can tabulate
  outcome_SGLT2ivsDPP4iSU_hr <- cbind(outcome=k, count[c(1:2)], followup[c(1:2)], events[c(1:2)],
                                      unadjusted_SGLT2i_string, adjusted_SGLT2i_string, overlapweighted_SGLT2i_string, iptw_SGLT2i_string)
  
  SGLT2ivsDPP4iSU_hrs <- rbind(SGLT2ivsDPP4iSU_hrs, outcome_SGLT2ivsDPP4iSU_hr)
  
  SGLT2ivsDPP4iSU_hr <- rbind(
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Unadjusted", 
          HR = unadjusted_SGLT2i[1], LB = unadjusted_SGLT2i[2], UB = unadjusted_SGLT2i[3], string = unadjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Adjusted", 
          HR = adjusted_SGLT2i[1], LB = adjusted_SGLT2i[2], UB = adjusted_SGLT2i[3], string = adjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Overlap-weighted", 
          HR = overlapweighted_SGLT2i[1], LB = overlapweighted_SGLT2i[2], UB = overlapweighted_SGLT2i[3], string = overlapweighted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "IPTW", 
          HR = iptw_SGLT2i[1], LB = iptw_SGLT2i[2], UB = iptw_SGLT2i[3], string = iptw_SGLT2i_string))
    
    all_SGLT2ivsDPP4iSU_hrs <- rbind(all_SGLT2ivsDPP4iSU_hrs, SGLT2ivsDPP4iSU_hr)
  
    rm(temp)
}

all_SGLT2ivsDPP4iSU_hrs2 <- SGLT2ivsDPP4iSU_hrs2 <- data.frame()

# sensitivity analyses - keep 1 drug episode per patid only
for (k in "ckd_egfr50") {
  
  # keep one observation per patid keeping the SGLT2i episode preferentially
  temp <- cohort %>%
    mutate(studydrug2 = factor(studydrug2, levels = c("SGLT2i", "DPP4i/SU"))) %>%
    arrange(.imp, patid, studydrug2) %>%
    distinct(.imp, patid, .keep_all = TRUE) %>%
    ungroup() %>% 
    mutate(studydrug2 = factor(studydrug2, levels = c("DPP4i/SU", "SGLT2i")))
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  # calculate number of subjects in each group
  count <- temp[temp$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_count",
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- temp[temp$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_followup",
                values_from=time)
  
  # summarise number of events per group
  events <- temp[temp$.imp > 0,] %>%
    group_by(studydrug2) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(studydrug2, events) %>%
    pivot_wider(names_from=studydrug2,
                names_glue="{studydrug2}_events",
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2"))
  
  f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", paste(covariates, collapse=" + ")))
  
  # create empty vectors to store the hazard ratios from every imputed dataset
  # for the unadjusted survival models
  COEFS.SGLT2i.unadj <- SE.SGLT2i.unadj <-
    # for the adjusted survival models
    COEFS.SGLT2i.adj <- SE.SGLT2i.adj <-
    # for the overlap weighted survival models
    COEFS.SGLT2i.ow <- SE.SGLT2i.ow <-
    # for the inverse probability of treatment weighted survival models
    COEFS.SGLT2i.iptw <- SE.SGLT2i.iptw <-
    rep(NA,n.imp)
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
    #unadjusted analyses first
    fit.unadj <- coxph(f2, temp[temp$.imp == i,])
    
    #store coefficients and standard errors from this model
    COEFS.SGLT2i.unadj[i] <- fit.unadj$coefficients[1]
    SE.SGLT2i.unadj[i] <- sqrt(fit.unadj$var[1,1])
    
    #adjusted analyses
    fit.adj <- coxph(f_adjusted2, temp[temp$.imp == i,])
    
    COEFS.SGLT2i.adj[i] <- fit.adj$coefficients[1]
    SE.SGLT2i.adj[i] <- sqrt(fit.adj$var[1,1])
    
    #overlap weighted analyses
    fit.ow <- coxph(f_adjusted2, temp[temp$.imp ==i,], weights = overlap2)
    
    COEFS.SGLT2i.ow[i] <- fit.ow$coefficients[1]
    SE.SGLT2i.ow[i] <- sqrt(fit.ow$var[1,1])
    
    #inverse probability of treatment weighted analyses
    fit.iptw <- coxph(f_adjusted2, temp[temp$.imp ==i,], weights = IPTW2)
    
    COEFS.SGLT2i.iptw[i] <- fit.iptw$coefficients[1]
    SE.SGLT2i.iptw[i] <- sqrt(fit.iptw$var[1,1])
    
    rm(fit.unadj)
    rm(fit.adj)
    rm(fit.ow)
    rm(fit.iptw)
  }
  
  # pool hazard ratios
  unadjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.unadj, SE.SGLT2i.unadj, n.imp)
  adjusted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.adj, SE.SGLT2i.adj, n.imp)
  overlapweighted_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.ow, SE.SGLT2i.ow, n.imp)
  iptw_SGLT2i <- pool.rubin.HR(COEFS.SGLT2i.iptw, SE.SGLT2i.iptw, n.imp)
  
  # save pooled HR and 95% confidence interval
  unadjusted_SGLT2i_string <- paste0(sprintf("%.2f", round(unadjusted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(unadjusted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(unadjusted_SGLT2i[3], 2)), ")")
  adjusted_SGLT2i_string <- paste0(sprintf("%.2f", round(adjusted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(adjusted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(adjusted_SGLT2i[3], 2)), ")")
  overlapweighted_SGLT2i_string <- paste0(sprintf("%.2f", round(overlapweighted_SGLT2i[1], 2)), " (", sprintf("%.2f", round(overlapweighted_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(overlapweighted_SGLT2i[3], 2)), ")")
  iptw_SGLT2i_string <- paste0(sprintf("%.2f", round(iptw_SGLT2i[1], 2)), " (", sprintf("%.2f", round(iptw_SGLT2i[2], 2)), ", ", sprintf("%.2f", round(iptw_SGLT2i[3], 2)), ")")
  
  # combine in dataframe that we can tabulate
  outcome_SGLT2ivsDPP4iSU_hr <- cbind(outcome=k, count[c(1:2)], followup[c(1:2)], events[c(1:2)],
                                      unadjusted_SGLT2i_string, adjusted_SGLT2i_string, overlapweighted_SGLT2i_string, iptw_SGLT2i_string)
  
  SGLT2ivsDPP4iSU_hrs2 <- rbind(SGLT2ivsDPP4iSU_hrs2, outcome_SGLT2ivsDPP4iSU_hr)
  
  SGLT2ivsDPP4iSU_hr <- rbind(
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Unadjusted", 
          HR = unadjusted_SGLT2i[1], LB = unadjusted_SGLT2i[2], UB = unadjusted_SGLT2i[3], string = unadjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Adjusted", 
          HR = adjusted_SGLT2i[1], LB = adjusted_SGLT2i[2], UB = adjusted_SGLT2i[3], string = adjusted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "Overlap-weighted", 
          HR = overlapweighted_SGLT2i[1], LB = overlapweighted_SGLT2i[2], UB = overlapweighted_SGLT2i[3], string = overlapweighted_SGLT2i_string),
    cbind(outcome = k, contrast = "SGLT2i vs DPP4i/SU", analysis = "IPTW", 
          HR = iptw_SGLT2i[1], LB = iptw_SGLT2i[2], UB = iptw_SGLT2i[3], string = iptw_SGLT2i_string))
  
  all_SGLT2ivsDPP4iSU_hrs2 <- rbind(all_SGLT2ivsDPP4iSU_hrs2, SGLT2ivsDPP4iSU_hr)
  
  rm(temp)
}
############################3 SUBGROUP ANALYSES################################################################


#create empty data frame to which we can append the hazard ratios once calculated
subgroup_SGLT2ivsDPP4iSU_hrs <- subgroup_hrs <- data.frame()

#ensure albuminuria recognised as logical variable:
cohort <- cohort %>% mutate(albuminuria = as.logical(albuminuria))

## analyses stratified by subgroup
for (k in outcomes) {
  
  print(paste0("Calculating hazard ratios for outcome ", k))
  
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  
    
    # calculate number of subjects in each group
    count <- cohort %>%
      group_by(studydrug2,albuminuria) %>%
      summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_count",
                  values_from=count)
    
    # calculate median follow up time (years) per group
    followup <- cohort %>%
      group_by(studydrug2,albuminuria) %>%
      summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_followup",
                  values_from=time)
    
    # summarise number of events per group
    events <- cohort %>%
      group_by(studydrug2,albuminuria) %>%
      summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
                drug_count=round(n()/n.imp, 0)) %>%
      mutate(events_perc=round(event_count*100/drug_count, 1),
             events=paste0(event_count, " (", events_perc, "%)")) %>%
      select(studydrug2, albuminuria, events) %>%
      pivot_wider(names_from=studydrug2,
                  names_glue="{studydrug2}_events",
                  values_from=events)
    
    
    # write formulas for adjusted and unadjusted analyses
    f2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*albuminuria"))
    
    f_adjusted2 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2*albuminuria + ", paste(covariates, collapse=" + ")))
    
    # create empty vectors to store the hazard ratios from every imputed dataset
    # for the unadjusted survival models
    COEFS.noalb.unadj <- SE.noalb.unadj <-
      COEFS.microalb.unadj <- SE.microalb.unadj <-
      # for the adjusted survival models
      COEFS.noalb.adj <- SE.noalb.adj <-
      COEFS.microalb.adj <- SE.microalb.adj <-
      # for the overlap-weighted survival models
      COEFS.noalb.ow <- SE.noalb.ow <-
      COEFS.microalb.ow <- SE.microalb.ow <-
      p_value_interaction <-
      rep(NA,n.imp)
    
    for (i in 1:n.imp) {
      print(paste0("Analyses in imputed dataset number ", i))
      
      #unadjusted analyses first
      fit.unadj <- coxph(f2, cohort[cohort$.imp == i,])
      
      #store coefficients and standard errors from this model
      COEFS.noalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"]
      SE.noalb.unadj[i] <- sqrt(fit.unadj$var[1,1])
      
      COEFS.microalb.unadj[i] <- fit.unadj$coefficients["studydrug2SGLT2i"] + fit.unadj$coefficients["studydrug2SGLT2i:albuminuriaTRUE"]
      SE.microalb.unadj[i] <- sqrt(abs(fit.unadj$var[1]) + abs(fit.unadj$var[nrow(fit.unadj$var),nrow(fit.unadj$var)]) + 2 * vcov(fit.unadj)[1,nrow(fit.unadj$var)])
      
      #adjusted analyses
      fit.adj <- coxph(f_adjusted2, cohort[cohort$.imp == i,])
      
      COEFS.noalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"]
      SE.noalb.adj[i] <- sqrt(fit.adj$var[1,1])
      
      COEFS.microalb.adj[i] <- fit.adj$coefficients["studydrug2SGLT2i"] + fit.adj$coefficients["studydrug2SGLT2i:albuminuriaTRUE"]
      SE.microalb.adj[i] <- sqrt(abs(fit.adj$var[1]) + abs(fit.adj$var[nrow(fit.adj$var),nrow(fit.adj$var)]) + 2 * vcov(fit.adj)[1,nrow(fit.adj$var)])
      
      #overlap-weighted analyses
      fit.ow <- coxph(f_adjusted2, cohort[cohort$.imp == i,], weights = overlap2)
      
      COEFS.noalb.ow[i] <- fit.ow$coefficients["studydrug2SGLT2i"]
      SE.noalb.ow[i] <- sqrt(fit.ow$var[1,1])
      
      COEFS.microalb.ow[i] <- fit.ow$coefficients["studydrug2SGLT2i"] + fit.ow$coefficients["studydrug2SGLT2i:albuminuriaTRUE"]
      SE.microalb.ow[i] <- sqrt(abs(fit.ow$var[1]) + abs(fit.ow$var[nrow(fit.ow$var),nrow(fit.ow$var)]) + 2 * vcov(fit.ow)[1,nrow(fit.ow$var)])
      
      if (k == "ckd_egfr50") {
          f3 <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug2 + ", paste(covariates, collapse = " + ")))
          fit.no_interaction <- coxph(f3, cohort[cohort$.imp == i,], weights = overlap2)
          
          chi <- 2 * fit.ow$loglik[2] - 2 * fit.no_interaction$loglik[2]
          p_value_interaction[i] <- 1 - pchisq(chi, df = 2)
        }
      
    }
    
    # pool hazard ratios
    unadjusted_noalb <- pool.rubin.HR(COEFS.noalb.unadj, SE.noalb.unadj, n.imp)
    adjusted_noalb <- pool.rubin.HR(COEFS.noalb.adj, SE.noalb.adj, n.imp)
    ow_noalb <- pool.rubin.HR(COEFS.noalb.ow, SE.noalb.ow, n.imp)
    
    unadjusted_microalb <- pool.rubin.HR(COEFS.microalb.unadj, SE.microalb.unadj, n.imp)
    adjusted_microalb <- pool.rubin.HR(COEFS.microalb.adj, SE.microalb.adj, n.imp)
    ow_microalb <- pool.rubin.HR(COEFS.microalb.ow, SE.microalb.ow, n.imp)
    
    # save pooled HR and 95% confidence interval
    unadjusted_noalb_string <- paste0(sprintf("%.2f", round(unadjusted_noalb[1], 2)), " (", sprintf("%.2f", round(unadjusted_noalb[2], 2)), ", ", sprintf("%.2f", round(unadjusted_noalb[3], 2)), ")")
    adjusted_noalb_string <- paste0(sprintf("%.2f", round(adjusted_noalb[1], 2)), " (", sprintf("%.2f", round(adjusted_noalb[2], 2)), ", ", sprintf("%.2f", round(adjusted_noalb[3], 2)), ")")
    ow_noalb_string <- paste0(sprintf("%.2f", round(ow_noalb[1], 2)), " (", sprintf("%.2f", round(ow_noalb[2], 2)), ", ", sprintf("%.2f", round(ow_noalb[3], 2)), ")")
    
    unadjusted_microalb_string <- paste0(sprintf("%.2f", round(unadjusted_microalb[1], 2)), " (", sprintf("%.2f", round(unadjusted_microalb[2], 2)), ", ", sprintf("%.2f", round(unadjusted_microalb[3], 2)), ")")
    adjusted_microalb_string <- paste0(sprintf("%.2f", round(adjusted_microalb[1], 2)), " (", sprintf("%.2f", round(adjusted_microalb[2], 2)), ", ", sprintf("%.2f", round(adjusted_microalb[3], 2)), ")")
    ow_microalb_string <- paste0(sprintf("%.2f", round(ow_microalb[1], 2)), " (", sprintf("%.2f", round(ow_microalb[2], 2)), ", ", sprintf("%.2f", round(ow_microalb[3], 2)), ")")
    
    # combine in dataframe that we can tabulate
    presegfr_noalb_hr <- cbind(outcome=k, count[1,c(2:3)], followup[1,c(2:3)], events[1,c(2:3)], contrast = "uACR <3mg/mmol",
                               unadjusted=unadjusted_noalb_string, adjusted=adjusted_noalb_string, overlapweighted=ow_noalb_string
    )
    presegfr_microalb_hr <- cbind(outcome=k, count[2,c(2:3)], followup[2,c(2:3)], events[2,c(2:3)], contrast = "uACR 3-30mg/mmol",
                                  unadjusted=unadjusted_microalb_string, adjusted=adjusted_microalb_string, overlapweighted=ow_microalb_string
    )
    
    outcome_subgroup_SGLT2ivsDPP4iSU_hrs <- rbind(presegfr_noalb_hr, presegfr_microalb_hr)
    
    temp <- rbind(
      cbind(outcome = k, contrast = "uACR <3mg/mmol", analysis = "Adjusted",
            HR = adjusted_noalb[1], LB = adjusted_noalb[2], UB = adjusted_noalb[3], string = adjusted_noalb_string),
      cbind(outcome = k, contrast = "uACR 3-30mg/mmol", analysis = "Adjusted",
            HR = adjusted_microalb[1], LB = adjusted_microalb[2], UB = adjusted_microalb[3], string = adjusted_microalb_string),
      cbind(outcome = k, contrast = "uACR <3mg/mmol", analysis = "Overlap-weighted",
            HR = ow_noalb[1], LB = ow_noalb[2], UB = ow_noalb[3], string = ow_noalb_string),
      cbind(outcome = k, contrast = "uACR 3-30mg/mmol", analysis = "Overlap-weighted",
            HR = ow_microalb[1], LB = ow_microalb[2], UB = ow_microalb[3], string = ow_microalb_string)
        )
    
    subgroup_hrs <- rbind(subgroup_hrs, temp)
    subgroup_SGLT2ivsDPP4iSU_hrs <- rbind(subgroup_SGLT2ivsDPP4iSU_hrs, outcome_subgroup_SGLT2ivsDPP4iSU_hrs)
    
}


############################4 STORE HAZARD RATIOS################################################################

# add numbers treated and events in suitable format for graphs

all_hrs$model <- paste0(all_hrs$string, " [", all_hrs$analysis, "]")
all_hrs$model <- factor(all_hrs$model, levels = unique(all_hrs$model))

all_SGLT2ivsDPP4iSU_hrs$model <- paste0(all_SGLT2ivsDPP4iSU_hrs$string, " [", all_SGLT2ivsDPP4iSU_hrs$analysis, "]")
all_SGLT2ivsDPP4iSU_hrs$model <- factor(all_SGLT2ivsDPP4iSU_hrs$model, levels = unique(all_SGLT2ivsDPP4iSU_hrs$model))

all_SGLT2ivsDPP4iSU_hrs2$model <- paste0(all_SGLT2ivsDPP4iSU_hrs2$string, " [", all_SGLT2ivsDPP4iSU_hrs2$analysis, "]")
all_SGLT2ivsDPP4iSU_hrs2$model <- factor(all_SGLT2ivsDPP4iSU_hrs2$model, levels = unique(all_SGLT2ivsDPP4iSU_hrs2$model))

# have to coerce HR and CI to class numeric as they sometimes default to character
class(all_hrs$HR) <- class(all_hrs$LB) <- class(all_hrs$UB) <-
  class(all_SGLT2ivsDPP4iSU_hrs$HR) <- class(all_SGLT2ivsDPP4iSU_hrs$LB) <- class(all_SGLT2ivsDPP4iSU_hrs$UB) <- 
  class(all_SGLT2ivsDPP4iSU_hrs2$HR) <- class(all_SGLT2ivsDPP4iSU_hrs2$LB) <- class(all_SGLT2ivsDPP4iSU_hrs2$UB) <- 
  class(subgroup_hrs$HR) <- class(subgroup_hrs$LB) <- class(subgroup_hrs$UB) <- "numeric"

# remove unadjusted HR as we do not want to plot these
all_hrs <- all_hrs[!all_hrs$analysis == "Unadjusted",]
all_SGLT2ivsDPP4iSU_hrs <- all_SGLT2ivsDPP4iSU_hrs[!all_SGLT2ivsDPP4iSU_hrs$analysis == "Unadjusted",]
all_SGLT2ivsDPP4iSU_hrs2 <- all_SGLT2ivsDPP4iSU_hrs2[!all_SGLT2ivsDPP4iSU_hrs2$analysis == "Unadjusted",]

#add event count and total count to hr dataframe
all_SGLT2ivsDPP4iSU_hrs <- all_SGLT2ivsDPP4iSU_hrs %>% left_join(SGLT2ivsDPP4iSU_hrs %>% select(1:7))
all_SGLT2ivsDPP4iSU_hrs2 <- all_SGLT2ivsDPP4iSU_hrs2 %>% left_join(SGLT2ivsDPP4iSU_hrs2 %>% select(1:7))
subgroup_hrs <- subgroup_hrs %>% left_join(subgroup_SGLT2ivsDPP4iSU_hrs %>% select(-c(adjusted, unadjusted)), by = c("outcome", "contrast"))


all_hrs <- all_hrs %>%
  separate(`DPP4i_events`, into = c("DPP4i_events_number", "DPP4i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(`SU_events`, into = c("SU_events_number", "SU_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(SGLT2i_events, into = c("SGLT2i_events_number", "SGLT2i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `DPP4i_events_percentage` = str_replace(`DPP4i_events_percentage`, "\\)", ""),
    SU_events_percentage = str_replace(SU_events_percentage, "\\)", ""),
    SGLT2i_events_percentage = str_replace(SGLT2i_events_percentage, "\\)", ""),
    `DPP4i_nN` = paste0(`DPP4i_events_number`, "/", `DPP4i_count`),
    `SU_nN` = paste0(`SU_events_number`, "/", `SU_count`),
    SGLT2i_nN = paste0(SGLT2i_events_number, "/", SGLT2i_count)
  )


all_SGLT2ivsDPP4iSU_hrs <- all_SGLT2ivsDPP4iSU_hrs %>%
  separate(`DPP4i/SU_events`, into = c("DPP4i/SU_events_number", "DPP4i/SU_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(SGLT2i_events, into = c("SGLT2i_events_number", "SGLT2i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `DPP4i/SU_events_percentage` = str_replace(`DPP4i/SU_events_percentage`, "\\)", ""),
    SGLT2i_events_percentage = str_replace(SGLT2i_events_percentage, "\\)", ""),
    `DPP4i/SU_nN` = paste0(`DPP4i/SU_events_number`, "/", `DPP4i/SU_count`),
    SGLT2i_nN = paste0(SGLT2i_events_number, "/", SGLT2i_count)
  )

all_SGLT2ivsDPP4iSU_hrs2 <- all_SGLT2ivsDPP4iSU_hrs2 %>%
  separate(`DPP4i/SU_events`, into = c("DPP4i/SU_events_number", "DPP4i/SU_events_percentage"), sep = " \\(", remove = FALSE) %>%
  separate(SGLT2i_events, into = c("SGLT2i_events_number", "SGLT2i_events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `DPP4i/SU_events_percentage` = str_replace(`DPP4i/SU_events_percentage`, "\\)", ""),
    SGLT2i_events_percentage = str_replace(SGLT2i_events_percentage, "\\)", ""),
    `DPP4i/SU_nN` = paste0(`DPP4i/SU_events_number`, "/", `DPP4i/SU_count`),
    SGLT2i_nN = paste0(SGLT2i_events_number, "/", SGLT2i_count)
  )

subgroup_hrs <- subgroup_hrs %>%
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

save(all_hrs, file=paste0(today, "_all_hrs.Rda"))
save(all_SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_all_SGLT2ivsDPP4iSU_hrs.Rda"))
save(SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_SGLT2ivsDPP4iSU_hrs.Rda"))
save(all_SGLT2ivsDPP4iSU_hrs2, file=paste0(today, "_all_SGLT2ivsDPP4iSU_hrs2.Rda"))
save(SGLT2ivsDPP4iSU_hrs2, file=paste0(today, "_SGLT2ivsDPP4iSU_hrs2.Rda"))
save(subgroup_hrs, file=paste0(today, "_subgroup_hrs.Rda"))
save(subgroup_SGLT2ivsDPP4iSU_hrs, file=paste0(today, "_subgroup_SGLT2ivsDPP4iSU_hrs.Rda"))
