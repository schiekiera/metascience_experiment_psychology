# ==============================================================================
# remove everything
rm(list = ls())

# load libraries
library(tidyverse)
library(lme4)
library(mediation)
library(utils)
library(report)
library(MuMIn)
library(writexl)
library(ggthemes)
library(jtools)
library(gridExtra)

#########################################################
#####                                               #####
##### Script: Model Results Stage 2                #####
#####                                               #####
#########################################################

# set seed
set.seed(2025)

##############
## Load data #
##############
url1 <- "https://raw.githubusercontent.com/schiekiera/metascience_experiment_psychology/main/data/ap1_experiment1_anonymized.csv"
url2 <- "https://raw.githubusercontent.com/schiekiera/metascience_experiment_psychology/main/data/ap1_experiment2_anonymized.csv"
url3 <- "https://raw.githubusercontent.com/schiekiera/metascience_experiment_psychology/main/data/ap1_experiment3_anonymized.csv"
url4 <- "https://raw.githubusercontent.com/schiekiera/metascience_experiment_psychology/main/data/ap1_experiment4_anonymized.csv"
df1<- read.csv(url1)
df2<- read.csv(url2)
df3<- read.csv(url3)
df4<- read.csv(url4)

## Add experiment column
df1$experiment <- "experiment_1"
df2$experiment <- "experiment_2"
df3$experiment <- "experiment_3"
df4$experiment <- "experiment_4"


######### exclude trials with reading times below 2 seconds
# experiment 1
nrow_before1 <- nrow(df1)
df1 <- df1 %>% filter(read1_rt > 2000)
nrow_after1 <- nrow(df1)
# experiment 2
nrow_before2 <- nrow(df2)
df2 <- df2 %>% filter(read1_rt > 2000)
nrow_after2 <- nrow(df2)
# experiment 3
nrow_before3 <- nrow(df3)
df3 <- df3 %>% filter(read1_rt > 2000)
nrow_after3 <- nrow(df3)
# experiment 4
nrow_before4 <- nrow(df4)
df4 <- df4 %>% filter(read1_rt > 2000)
nrow_after4 <- nrow(df4)

## recode binary treatment from [0,1] to [1,0]
## 0 == control, 1 == treatment
### HC: 0 == hypothesis-consistent, 1 == hypothesis-inconsistent
### Sig: 0 == significant, 1 == non-significant
df1$treatment <- ifelse(df1$treatment == 0, 1, 0)
df2$treatment <- ifelse(df2$treatment == 0, 1, 0)
df3$treatment <- ifelse(df3$treatment == 0, 1, 0)
df4$treatment <- ifelse(df4$treatment == 0, 1, 0)


## combine df1 and df2, and df 3 and df4
df1_2<-rbind(df1,df2)
df3_4<-rbind(df3,df4)


# create a new variable for the decision change
df1_2$decision_change <- df1_2$decision2_resp - df1_2$decision1_resp
df3_4$decision_change_read <- df3_4$decision2_read_resp - df3_4$decision1_read_resp
df3_4$decision_change_cite <- df3_4$decision2_cite_resp - df3_4$decision1_cite_resp
df1$decision_change <-df1$decision2_resp - df1$decision1_resp
df2$decision_change <-df2$decision2_resp - df2$decision1_resp
df3$decision_change_read <-df3$decision2_read_resp - df3$decision1_read_resp
df3$decision_change_cite <-df3$decision2_cite_resp - df3$decision1_cite_resp
df4$decision_change_read <-df4$decision2_read_resp - df4$decision1_read_resp
df4$decision_change_cite <-df4$decision2_cite_resp - df4$decision1_cite_resp

# report in manuscript
# paste XX % (SD = XX, Range = XX-XX%))
# decision_change
cat(round(mean(df1_2$decision_change, na.rm = T),2)," (SD = ",round(sd(df1_2$decision_change, na.rm = T),2),", Range = ",round(min(df1_2$decision_change, na.rm = T),2),"-",round(max(df1_2$decision_change, na.rm = T),2),")")
cat(round(mean(df3_4$decision_change_read, na.rm = T),2)," (SD = ",round(sd(df3_4$decision_change_read, na.rm = T),2),", Range = ",round(min(df3_4$decision_change_read, na.rm = T),2),"-",round(max(df3_4$decision_change_read, na.rm = T),2),")")
cat(round(mean(df3_4$decision_change_cite, na.rm = T),2)," (SD = ",round(sd(df3_4$decision_change_cite, na.rm = T),2),", Range = ",round(min(df3_4$decision_change_cite, na.rm = T),2),"-",round(max(df3_4$decision_change_cite, na.rm = T),2),")")

# for_resp
# add 1 to for_resp to match the manuscript: (1-7, “When I gave the answer, I felt" - "Very uncertain " - " Very certain "),
df1_2$for_resp<-df1_2$for_resp+1
df3_4$for_read_resp<-df3_4$for_read_resp+1
df3_4$for_cite_resp<-df3_4$for_cite_resp+1
cat(round(mean(df1_2$for_resp, na.rm = T),2)," (SD = ",round(sd(df1_2$for_resp, na.rm = T),2),", Range = ",round(min(df1_2$for_resp, na.rm = T),2),"-",round(max(df1_2$for_resp, na.rm = T),2),")")
cat(round(mean(df3_4$for_read_resp, na.rm = T),2)," (SD = ",round(sd(df3_4$for_read_resp, na.rm = T),2),", Range = ",round(min(df3_4$for_read_resp, na.rm = T),2),"-",round(max(df3_4$for_read_resp, na.rm = T),2),")")
cat(round(mean(df3_4$for_cite_resp, na.rm = T),2)," (SD = ",round(sd(df3_4$for_cite_resp, na.rm = T),2),", Range = ",round(min(df3_4$for_cite_resp, na.rm = T),2),"-",round(max(df3_4$for_cite_resp, na.rm = T),2),")")


# t-tests for differences between conditions
CLoS_Significance<-t.test(data = df1, decision_change ~ treatment)
CLoS_Hypothesis<-t.test(data = df2, decision_change ~ treatment)
CLoC_Significance<-t.test(data = df3, decision_change_cite ~ treatment)
CLoC_Hypothesis<-t.test(data = df4, decision_change_cite ~ treatment)
CLoR_Significance<-t.test(data = df3, decision_change_read ~ treatment)
CLoR_Hypothesis<-t.test(data = df4, decision_change_read ~ treatment)

# print t-tests
## significance
report(CLoS_Significance)
report(CLoR_Significance)
report(CLoC_Significance)

# hypothesis-consistency
report(CLoS_Hypothesis)
report(CLoR_Hypothesis)
report(CLoC_Hypothesis)





##################
# Stage 2 Models #
##################

# Experiment 1

# Model 1.2
# Path a: treatment → FOR
model_1.2_med <- lmer(for_resp ~ treatment + (treatment | id), data = df1)
report(model_1.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_1.2_out <- lmer(decision_change ~ for_resp + treatment + (treatment | id), data = df1)
# Model with singular fit, so drop random slope
model_1.2_out <- lmer(decision_change ~ for_resp + treatment + (1 | id), data = df1)
report(model_1.2_out)

# Fit the mediation model
model_1.2 <- mediate(
  model_1.2_med, # random slope
  model_1.2_out, # random intercept
  treat = "treatment",
  mediator = "for_resp",
  sims = 1000
)
summary(model_1.2)
sub_models_1.2<- c("random_slope", "random_intercept")

# Experiment 2

# Model 2.2
# Path a: treatment → FOR
model_2.2_med <- lmer(for_resp ~ treatment + (treatment | id), data = df2)
## Drop random slope due to singularity
model_2.2_med <- lmer(for_resp ~ treatment + (1 | id), data = df2)
report(model_2.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_2.2_out <- lmer(decision_change ~ for_resp + treatment + (treatment | id), data = df2)
report(model_2.2_out)

# Fit the mediation model
model_2.2 <- mediate(
  model_2.2_med, # random intercept
  model_2.2_out, # random slope
  treat = "treatment",
  mediator = "for_resp",
  sims = 1000
)
summary(model_2.2)

sub_models_2.2<- c("random_intercept", "random_slope")


# Experiment 3

# Model 3A.2
# Path a: treatment → FOR
model_3A.2_med <- lmer(for_read_resp ~ treatment + (treatment | id), data = df3)
report(model_3A.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_3A.2_out <- lmer(decision_change_read ~ for_read_resp + treatment + (treatment | id), data = df3)
## Model with singular fit, so drop random slope
model_3A.2_out <- lmer(decision_change_read ~ for_read_resp + treatment + (1 | id), data = df3)
report(model_3A.2_out)

# Fit the mediation model
model_3A.2 <- mediate(
  model_3A.2_med, # random slope
  model_3A.2_out, # random intercept
  treat = "treatment",
  mediator = "for_read_resp",
  sims = 1000
)
summary(model_3A.2)

sub_models_3A.2<- c("random_slope", "random_intercept")


# Model 3B.2
# Path a: treatment → FOR
model_3B.2_med <- lmer(for_cite_resp ~ treatment + (treatment | id), data = df3)
# Drop random slope due to singularity
model_3B.2_med <- lmer(for_cite_resp ~ treatment + (1 | id), data = df3)
report(model_3B.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_3B.2_out <- lmer(decision_change_cite ~ for_cite_resp + treatment + (treatment | id), data = df3)
report(model_3B.2_out)

# Fit the mediation model
model_3B.2 <- mediate(
  model_3B.2_med, # random intercept
  model_3B.2_out, # random slope
  treat = "treatment",
  mediator = "for_cite_resp",
  sims = 1000
)
summary(model_3B.2)

sub_models_3B.2<- c("random_intercept", "random_slope")

# Experiment 4

# Model 4A.2
# Path a: treatment → FOR
model_4A.2_med <- lmer(for_read_resp ~ treatment + (treatment | id), data = df4)
# Drop random slope due to singularity
model_4A.2_med <- lmer(for_read_resp ~ treatment + (1 | id), data = df4)
report(model_4A.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_4A.2_out <- lmer(decision_change_read ~ for_read_resp + treatment + (treatment | id), data = df4)
# Model with singular fit, so drop random slope
model_4A.2_out <- lmer(decision_change_read ~ for_read_resp + treatment + (1 | id), data = df4)
report(model_4A.2_out)

# Fit the mediation model
model_4A.2 <- mediate(
  model_4A.2_med, # random intercept
  model_4A.2_out, # random intercept
  treat = "treatment",
  mediator = "for_read_resp",
  sims = 1000
)
summary(model_4A.2)

sub_models_4A.2<- c("random_intercept", "random_intercept")

# Model 4B.2
# Path a: treatment → FOR
model_4B.2_med <- lmer(for_cite_resp ~ treatment + (treatment | id), data = df4)
report(model_4B.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_4B.2_out <- lmer(decision_change_cite ~ for_cite_resp + treatment + (treatment | id), data = df4)
# Drop random slope due to singularity
model_4B.2_out <- lmer(decision_change_cite ~ for_cite_resp + treatment + (1 | id), data = df4)
report(model_4B.2_out)

# Fit the mediation model
model_4B.2 <- mediate(
  model_4B.2_med, # random_slope
  model_4B.2_out, # random intercept
  treat = "treatment",
  mediator = "for_cite_resp",
  sims = 1000
)
summary(model_4B.2)

sub_models_4B.2<- c("random_slope", "random_intercept")

# combine sub_models
sub_models<-data.frame(exp1= sub_models_1.2, 
                       exp2= sub_models_2.2,
                       exp3A= sub_models_3A.2,
                       exp3B= sub_models_3B.2,
                       exp4A= sub_models_4A.2,
                       exp4B= sub_models_4B.2,
                       row.names = c("med: path a","out: path b & c"))

# Load the lmerTest package
# lmerTest extends the functionality of lme4 by providing methods for
# calculating approximate p-values for fixed effects using the Satterthwaite (or other)
# degrees-of-freedom approximations.
library(lmerTest)

# Refit models for p-value computation
print(sub_models)

## Modsub_models## Model 1.2
model_1.2_med_refit<-lmer(for_resp ~ treatment + (treatment | id), data = df1)
model_1.2_out_refit<-lmer(decision_change ~ for_resp + treatment + (1 | id), data = df1)

## Model 2.2
model_2.2_med_refit<-lmer(for_resp ~ treatment + (1 | id), data = df2)
model_2.2_out_refit<-lmer(decision_change ~ for_resp + treatment + (treatment | id), data = df2)

## Model 3A.2
model_3A.2_med_refit<-lmer(for_read_resp ~ treatment + (treatment | id), data = df3)
model_3A.2_out_refit<-lmer(decision_change_read ~ for_read_resp + treatment + (1 | id), data = df3)

## Model 3B.2
model_3B.2_med_refit<-lmer(for_cite_resp ~ treatment + (1 | id), data = df3)
model_3B.2_out_refit<-lmer(decision_change_cite ~ for_cite_resp + treatment + (treatment | id), data = df3)

## Model 4A.2
model_4A.2_med_refit<-lmer(for_read_resp ~ treatment + (1 | id), data = df4)
model_4A.2_out_refit<-lmer(decision_change_read ~ for_read_resp + treatment + (1 | id), data = df4)

## Model 4B.2
model_4B.2_med_refit<-lmer(for_cite_resp ~ treatment + (treatment | id), data = df4)
model_4B.2_out_refit<-lmer(decision_change_cite ~ for_cite_resp + treatment + (1 | id), data = df4)



# functions for processing the table and texts
# process beta
process_beta<-function(beta, p_value){
  if (p_value < 0.001) {
    beta <- paste0(round(beta, 2), "***")
  } else if (p_value < 0.01) {
    beta <- paste0(round(beta, 2), "**")
  } else if (p_value < 0.05) {
    beta <- paste0(round(beta, 2), "*")
  } else {
    beta <- round(beta, 2)
  }
  return(beta)
}

# define rounding functions
round_normal <- function(x, digits = 2) {
  round(x, digits = digits)
}
round_p_value <- function(x, digits = 3) {
  p<-round(x, digits = digits)
  if(p<0.001){
    p<- "< .001"
  }
  else {
    # turn it into a string 
    p<-as.character(p)
    # delete first character
    p<-substr(p, 2, nchar(p))
    # add "= "
    p<-paste0("= ", p)
  }
  
  return(p)
}


# define function to obtain model parameters
obtain_model_parameters <- function(model_name,
                                    mediation_model,
                                    path_a_model,
                                    path_b_c_model,
                                    for_type,
                                    output_type = "text") {

  
  # Path a
  ## Define fixed_effects_a and confidence_intervals_a
  fixed_effects_a <- summary(path_a_model)$coefficients
  confidence_intervals_a<-confint(path_a_model, method = "Wald")
  ## Extract beta, p-value and confidence intervals fot treatment
  beta_treatment_a <-  round_normal(fixed_effects_a["treatment", "Estimate"])
  p_value_treatment_a <- fixed_effects_a["treatment", "Pr(>|t|)"]
  beta_p_treatment_a<-process_beta(beta_treatment_a, p_value_treatment_a)
  ci_treatment_a<- round_normal(confidence_intervals_a["treatment",])
  ci_treatment_a<-paste0("[",paste0(ci_treatment_a, collapse = ", "),"]")

  # Path b and c
  ## Define fixed_effects_b_c
  fixed_effects_b_c <- summary(path_b_c_model)$coefficients
  confidence_intervals_b_c<-confint(path_b_c_model, method = "Wald")
  ## Extract beta, p-value and confidence intervals
  ### path b: for_resp
  beta_for_resp_b <-  round_normal(fixed_effects_b_c[for_type, "Estimate"])
  p_value_for_resp_b <- fixed_effects_b_c[for_type, "Pr(>|t|)"]
  beta_p_for_resp_b<-process_beta(beta_for_resp_b, p_value_for_resp_b)
  ci_for_resp_b<- round_normal(confidence_intervals_b_c[for_type,])
  ci_for_resp_b<-paste0("[",paste0(ci_for_resp_b, collapse = ", "),"]")
  ### path c: treatment
  beta_treatment_c <-  round_normal(fixed_effects_b_c["treatment", "Estimate"])
  p_value_treatment_c <- fixed_effects_b_c["treatment", "Pr(>|t|)"]
  beta_p_treatment_c<-process_beta(beta_treatment_c, p_value_treatment_c)
  ci_treatment_c<- round_normal(confidence_intervals_b_c["treatment",])
  ci_treatment_c<-paste0("[",paste0(ci_treatment_c, collapse = ", "),"]")
  
  # ACME: d0
  # ACME = Average Causal Mediation Effect
  acme<-round_normal(mediation_model$d0)
  acme_p<-mediation_model$d0.p
  acme_beta_p<-process_beta(acme, acme_p)
  acme_ci<-round_normal(mediation_model$d0.ci)
  acme_ci<-paste0("[",paste0(acme_ci, collapse = ", "),"]")
  # ADE: z0
  ## ADE = Average Direct Effect
  ade<-round_normal(mediation_model$z0)
  ade_p<-mediation_model$z0.p
  ade_beta_p<-process_beta(ade, ade_p)
  ade_ci<-round_normal(mediation_model$z0.ci)
  ade_ci<-paste0("[",paste0(ade_ci, collapse = ", "),"]")
  # Total effect: tau
  total_effect<-round_normal(mediation_model$tau.coef)
  total_effect_p<-mediation_model$tau.p
  total_effect_beta_p<-process_beta(total_effect, total_effect_p)
  total_effect_ci<-round_normal(mediation_model$tau.ci)
  total_effect_ci<-paste0("[",paste0(total_effect_ci, collapse = ", "),"]")
  # proportion mediated: n0
  prop_mediated<-round_normal(mediation_model$n0)
  prop_mediated_p<-mediation_model$n0.p
  prop_mediated_beta_p<-process_beta(prop_mediated, prop_mediated_p)
  prop_mediated_ci<-round_normal(mediation_model$n0.ci)
  prop_mediated_ci<-paste0("[",paste0(prop_mediated_ci, collapse = ", "),"]")

  # Create fixed effects table
  results_table <- data.frame(
    beta_p = c(acme_beta_p, beta_p_treatment_a, beta_p_for_resp_b, ade_beta_p, total_effect_beta_p, prop_mediated_beta_p),
    ci = c(acme_ci, ci_treatment_a, ci_for_resp_b, ade_ci, total_effect_ci, prop_mediated_ci),
    row.names = c("ACME", "Path a", "Path b", "ADE", "Total Effect", "Proportion Mediated")
  )
  
  if (output_type == "text") {
    # Print results
    cat("Model: ", model_name, "\n")
    cat(paste0("ADE. (Beta = ", ade, ", p ", round_p_value(ade_p), ")\n"))
    cat(paste0("ACME. (Beta = ", acme, ", p ", round_p_value(acme_p), ")\n"))
    cat(paste0("Path a: treatment → FOR. (Beta = ", beta_treatment_a, ", p ", round_p_value(p_value_treatment_a), ")\n"))
    cat(paste0("Path b: FOR → C-LoSRC. (Beta = ", beta_for_resp_b,", p ", round_p_value(p_value_for_resp_b), ")\n"))
  } else if (output_type == "table"){
    # Return results
    return(results_table)
  }
}



# Get model parameters for text
obtain_model_parameters("model_1.2", model_1.2, model_1.2_med_refit, model_1.2_out_refit, for_type = "for_resp", output_type = "text")
obtain_model_parameters("model_2.2", model_2.2, model_2.2_med_refit, model_2.2_out_refit, for_type = "for_resp", output_type = "text")
obtain_model_parameters("model_3A.2", model_3A.2, model_3A.2_med_refit, model_3A.2_out_refit, for_type = "for_read_resp", output_type = "text")
obtain_model_parameters("model_3B.2", model_3B.2, model_3B.2_med_refit, model_3B.2_out_refit, for_type = "for_cite_resp", output_type = "text")
obtain_model_parameters("model_4A.2", model_4A.2, model_4A.2_med_refit, model_4A.2_out_refit, for_type = "for_read_resp", output_type = "text")
obtain_model_parameters("model_4B.2", model_4B.2, model_4B.2_med_refit, model_4B.2_out_refit, for_type = "for_cite_resp", output_type = "text")


# Get model parameters for table
results_table1<-obtain_model_parameters("model_1.2", model_1.2, model_1.2_med_refit, model_1.2_out_refit, for_type = "for_resp", output_type = "table")
results_table2<-obtain_model_parameters("model_2.2", model_2.2, model_2.2_med_refit, model_2.2_out_refit, for_type = "for_resp", output_type = "table")
results_table3A<-obtain_model_parameters("model_3A.2", model_3A.2, model_3A.2_med_refit, model_3A.2_out_refit, for_type = "for_read_resp", output_type = "table")
results_table3B<-obtain_model_parameters("model_3B.2", model_3B.2, model_3B.2_med_refit, model_3B.2_out_refit, for_type = "for_cite_resp", output_type = "table")
results_table4A<-obtain_model_parameters("model_4A.2", model_4A.2, model_4A.2_med_refit, model_4A.2_out_refit, for_type = "for_read_resp", output_type = "table")
results_table4B<-obtain_model_parameters("model_4B.2", model_4B.2, model_4B.2_med_refit, model_4B.2_out_refit, for_type = "for_cite_resp", output_type = "table")


# combine all results
df_part1<-data.frame(cbind(results_table1,results_table3A,results_table3B))
df_part2<-data.frame(cbind(results_table2,results_table4A,results_table4B),row.names = c("ACME1", "Path a1", "Path b1", "ADE1", "Total Effect1", "Proportion Mediated1"))

# Move row names into a column
df_part1$Effect <- rownames(df_part1)
df_part2$Effect <- rownames(df_part2)

# Ensure column names match
colnames(df_part2) <- colnames(df_part1)

# Stack dataframes
df_combined <- rbind(df_part1, df_part2)

# Restore row names and remove the extra column
rownames(df_combined) <- df_combined$Effect
# make Effect the first column
df_combined <- df_combined[, c(ncol(df_combined), 1:(ncol(df_combined)-1))]

# Print result
print(df_combined)

# write to xlsx
#write_xlsx(df_combined, paste0(output_path,"summary_stage_2_model_overview.xlsx"))




