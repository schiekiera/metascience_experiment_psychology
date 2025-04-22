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



####################
# IMPORTANT NOTE:  #
####################
# Since the analyses in this script rely on personal data,
# The sensitivity analyses II are NOT executable.



#############################################################################
#####                                                                   #####
##### Script: Sensitivity Analysis II - Model Results Stage 1 & 2          #####
#####                                                                   #####
#############################################################################

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

## recode binary treatment from [0,1] to [1,0]
## 0 == control, 1 == treatment
### HC: 0 == hypothesis-consistent, 1 == hypothesis-inconsistent
### Sig: 0 == significant, 1 == non-significant
df1$treatment <- ifelse(df1$treatment == 0, 1, 0)
df2$treatment <- ifelse(df2$treatment == 0, 1, 0)
df3$treatment <- ifelse(df3$treatment == 0, 1, 0)
df4$treatment <- ifelse(df4$treatment == 0, 1, 0)


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



###########################
## Sensitivity Analysis  ##
###########################
# Exclude participants with less than 10 trials

# Exclude all participants who not come from US, NED, GER, AUS, UK OR CAN
# This means --> drop all other == 1 FOR ALL DFS
df1 <- df1 %>% filter(other!=1)
df2 <- df2 %>% filter(other!=1)
df3 <- df3 %>% filter(other!=1)
df4 <- df4 %>% filter(other!=1)

# count number of experiments
print(paste0("Number of Experiments for Experiment 1: ",length(unique(df1$id))))
print(paste0("Number of Experiments for Experiment 2: ",length(unique(df2$id))))
print(paste0("Number of Experiments for Experiment 3: ",length(unique(df3$id))))
print(paste0("Number of Experiments for Experiment 4: ",length(unique(df4$id))))

# count number of participants from each country
df1_subset<-df1[duplicated(df1$id)==FALSE,]
df2_subset<-df2[duplicated(df2$id)==FALSE,]
df3_subset<-df3[duplicated(df3$id)==FALSE,]
df4_subset<-df4[duplicated(df4$id)==FALSE,]
countries1<-df1_subset[,c("uk","us","can","aus","nd")]
countries2<-df2_subset[,c("uk","us","can","aus","nd")]
countries3<-df3_subset[,c("uk","us","can","aus","nd")]
countries4<-df4_subset[,c("uk","us","can","aus","nd")]
all_countries<-rbind(countries1,countries2,countries3,countries4)
# one hot encode germany if all other countries are 0
all_countries$ger<-ifelse(rowSums(all_countries)==0,1,0)
# sum up all countries 
sum_countries<-colSums(all_countries)
print(sum_countries)



# define function to report model results
extract_model_stats_random_slope_model <- function(model, data, predictor, outcome) {
  # Extract coefficients
  fixed_effects <- summary(model)$coefficients
  
  # Extract beta and p-value
  beta <- fixed_effects[predictor, "Estimate"]
  p_value <- fixed_effects[predictor, "Pr(>|t|)"]
  
  # Calculate standardized beta
  sd_x <- sd(data[[predictor]])
  sd_y <- sd(data[[outcome]])
  beta_standardized <- beta * (sd_x / sd_y)
  
  # Extract random effects
  random_effects <- as.data.frame(VarCorr(model))
  slope_variance <- random_effects[random_effects$grp == "id" &
                                     grepl(predictor, random_effects$var1), "vcov"]
  slope_variance <- as.numeric(slope_variance)
  
  # Extract the random slopes for treatment for each id
  random_slopes <- ranef(model)$id$treatment
  
  # Compute the subject-specific slopes for each id
  subject_specific_slopes <- beta + random_slopes
  
  # Calculate the 95% central (2.5th and 97.5th percentile) interval
  central_95 <- quantile(subject_specific_slopes, probs = c(0.025, 0.975))
  
  # Extract R2 values
  r2 <- performance::r2(model)
  marginal_r2 <- r2$R2_marginal
  conditional_r2 <- r2$R2_conditional
  
  # Print results
  cat("Beta:", round(beta, 2), "\n")
  cat("P-Value:", round(p_value, 3), "\n")
  cat("Standardized Beta:", round(beta_standardized, 2), "\n")
  cat("Slope CI Lower:", round(central_95[1], 2), "\n")
  cat("Slope CI Upper:", round(central_95[2], 2), "\n")
  cat("Slope Variance:", round(slope_variance, 2), "\n")
  cat("Marginal R2:", round(marginal_r2, 2), "\n")
  cat("Conditional R2:", round(conditional_r2, 2), "\n")
  
  # Return results as a list
  return(
    list(
      Beta = beta,
      P_Value = p_value,
      Standardized_Beta = beta_standardized,
      Slope_Variance = slope_variance,
      Slope_Lower = central_95[1],
      Slope_Upper = central_95[2],
      Marginal_R2 = marginal_r2,
      Conditional_R2 = conditional_r2
    )
  )
}


# define function to report model results
extract_model_stats_random_intercept_model <- function(model, data, predictor, outcome) {
  # Extract coefficients
  fixed_effects <- summary(model)$coefficients
  
  # Extract beta and p-value
  beta <- fixed_effects[predictor, "Estimate"]
  p_value <- fixed_effects[predictor, "Pr(>|t|)"]
  
  # Calculate standardized beta
  sd_x <- sd(data[[predictor]])
  sd_y <- sd(data[[outcome]])
  beta_standardized <- beta * (sd_x / sd_y)
  
  # Extract R2 values
  r2 <- performance::r2(model)
  marginal_r2 <- r2$R2_marginal
  conditional_r2 <- r2$R2_conditional
  
  # Print results
  cat("Beta:", round(beta, 2), "\n")
  cat("P-Value:", round(p_value, 3), "\n")
  cat("Standardized Beta:", round(beta_standardized, 2), "\n")
  cat("Marginal R2:", round(marginal_r2, 2), "\n")
  cat("Conditional R2:", round(conditional_r2, 2), "\n")
  
  # Return results as a list
  return(
    list(
      Beta = beta,
      P_Value = p_value,
      Standardized_Beta = beta_standardized,
      Marginal_R2 = marginal_r2,
      Conditional_R2 = conditional_r2
    )
  )
}




##################
# Stage 1 Models #
##################

## Model 1.1
model_1.1 <-
  lmer(decision1_resp ~ treatment + (treatment | id) + (1 |
                                                          abstract), data = df1)
results_1.1 <- extract_model_stats_random_slope_model(model_1.1,
                                                      data = df1,
                                                      predictor = "treatment",
                                                      outcome = "decision1_resp")

## Model 2.1
model_2.1 <-
  lmer(decision1_resp ~ treatment + (treatment | id) + (1 |
                                                          abstract), data = df2)
# singularity --> drop random slopes
model_2.1 <-
  lmer(decision1_resp ~ treatment + (1 | id) + (1 |
                                                   abstract), data = df2)
results_2.1 <- extract_model_stats_random_slope_model(model_2.1,
                                                      data = df2,
                                                      predictor = "treatment",
                                                      outcome = "decision1_resp")

## Model 3.1A: reading
model_3A.1 <-
  lmer(decision1_read_resp ~ treatment + (treatment | id) + (1 |
                                                               abstract),
       data = df3)

results_3A.1 <- extract_model_stats_random_slope_model(model_3A.1,
                                                       data = df3,
                                                       predictor = "treatment",
                                                       outcome = "decision1_read_resp")

## Model 3.1B: Citing
model_3B.1 <-
  lmer(decision1_cite_resp ~ treatment + (treatment | id) + (1 |
                                                               abstract),
       data = df3)

results_3B.1 <- extract_model_stats_random_slope_model(model_3B.1,
                                                       data = df3,
                                                       predictor = "treatment",
                                                       outcome = "decision1_cite_resp")

## Model 4A.1: reading
model_4A.1 <-
  lmer(decision1_read_resp ~ treatment + (treatment | id) + (1 |
                                                               abstract),
       data = df4)
results_4A.1 <- extract_model_stats_random_slope_model(model_4A.1,
                                                       data = df4,
                                                       predictor = "treatment",
                                                       outcome = "decision1_read_resp")

## Model 4B.1B: Citing
model_4B.1 <-
  lmer(decision1_cite_resp ~ treatment + (treatment |
                                            id) + (1 | abstract),
       data = df4)
# singularity --> drop random slopes
model_4B.1 <-
  lmer(decision1_cite_resp ~ treatment + (1 | id) + (1 |
                                                       abstract), data = df4)

results_4B.1 <- extract_model_stats_random_slope_model(model_4B.1,
                                                       data = df4,
                                                       predictor = "treatment",
                                                       outcome = "decision1_cite_resp")




## Effect size plots
## Reproducible beta effect sizes plot

# define function to report model results
extract_parameters_for_plotting <- function(model, data, predictor, outcome) {
  # Extract coefficients
  fixed_effects <- summary(model)$coefficients
  
  # Extract beta and p-value
  beta <- fixed_effects[predictor, "Estimate"]
  p_value <- fixed_effects[predictor, "Pr(>|t|)"]
  
  # Calculate standardized beta
  sd_x <- sd(data[[predictor]])
  sd_y <- sd(data[[outcome]])
  beta_standardized <- beta * (sd_x / sd_y)
  se<- fixed_effects[predictor, "Std. Error"]
  se_standardized<- se * (sd_x / sd_y)
  result<-data.frame(beta_std=beta_standardized, se=se_standardized)
  return(result)
}


# betas
b1 <- extract_parameters_for_plotting(model_1.1, df1, "treatment", "decision1_resp")$beta_std
b2<- extract_parameters_for_plotting(model_2.1, df2, "treatment", "decision1_resp")$beta_std
b3a <- extract_parameters_for_plotting(model_3A.1, df3, "treatment", "decision1_read_resp")$beta_std
b3b <- extract_parameters_for_plotting(model_3B.1, df3, "treatment", "decision1_cite_resp")$beta_std
b4a <- extract_parameters_for_plotting(model_4A.1, df4, "treatment", "decision1_read_resp")$beta_std
b4b <- extract_parameters_for_plotting(model_4B_random_intercept, df4, "treatment", "decision1_read_resp")$beta_std ## select random intercept model because of singular fit of the random slope model

se1 <- extract_parameters_for_plotting(model_1.1, df1, "treatment", "decision1_resp")$se
se2<- extract_parameters_for_plotting(model_2.1, df2, "treatment", "decision1_resp")$se
se3a <- extract_parameters_for_plotting(model_3A.1, df3, "treatment", "decision1_read_resp")$se
se3b <- extract_parameters_for_plotting(model_3B.1, df3, "treatment", "decision1_cite_resp")$se
se4a <- extract_parameters_for_plotting(model_4A.1, df4, "treatment", "decision1_read_resp")$se
se4b <- extract_parameters_for_plotting(model_4B_random_intercept, df4, "treatment", "decision1_read_resp")$se ## select random intercept model because of singular fit of the random slope model



d <- data.frame(
  effect = c(b1, b2, b3a, b3b, b4a, b4b),
  se = c(se1, se2, se3a, se3b, se4a, se4b),
  expNr = factor(
    c(
      "Exp.1: Publish",
      "Exp.2: Publish",
      "Exp.3a: Read",
      "Exp.3b: Cite",
      "Exp.4a: Read",
      "Exp.4b: Cite"
    )
  ),
  expFac = factor(
    c(
      "Significance",
      "Hyp.-Consistency",
      "Significance",
      "Significance",
      "Hyp.-Consistency",
      "Hyp.-Consistency"
    )
  ),
  decision = factor(c(
    "Publish", "Publish", "Cite", "Read", "Cite", "Read"
  ))
)

figure4 <- ggplot(data = d,
            aes(effect, expNr),
            group = expFac,
            color = expFac) +
  geom_point(data = d, aes(effect, expNr, shape = expFac)) +
  geom_errorbar(
    aes(xmin = effect - se, xmax = effect + se),
    width = .2,
    position = position_dodge(width = .4)
  ) +
  labs(x = "Standardized beta", y = "Experiment") +
  scale_y_discrete(limits = sort(unique(d$expNr), decreasing = TRUE))+
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(
    axis.text = element_text(size = rel(0.8)),
    axis.title = element_text(size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme_apa(legend.pos = "bottom")


# save plot 
ggsave(paste0(plot_path,"sensitivity_figure4.png"), plot = figure4, width = 15, height = 10, units = "cm")
# saves as tiff file
ggsave(paste0(plot_path,"sensitivity_figure4.tiff"), plot = figure4, width = 15, height = 10, units = "cm", dpi = 300)



extract_model_results_table <- function(model,model_type) {
  # Extract fixed effects
  fixed_effects <- summary(model)$coefficients
  
  # Treatment
  beta_treatment <- fixed_effects["treatment", "Estimate"]
  p_value_treatment <- fixed_effects["treatment", "Pr(>|t|)"]
  if (p_value_treatment < 0.001) {
    beta_treatment <- paste0(round(beta_treatment, 2), "***")
  } else if (p_value_treatment < 0.01) {
    beta_treatment <- paste0(round(beta_treatment, 2), "**")
  } else if (p_value_treatment < 0.05) {
    beta_treatment <- paste0(round(beta_treatment, 2), "*")
  } else {
    beta_treatment <- round(beta_treatment, 2)
  }
  se_treatment <- round(fixed_effects["treatment", "Std. Error"], 2)
  
  # Intercept
  beta_intercept <- fixed_effects["(Intercept)", "Estimate"]
  p_value_intercept <- fixed_effects["(Intercept)", "Pr(>|t|)"]
  if (p_value_intercept < 0.001) {
    beta_intercept <- paste0(round(beta_intercept, 2), "***")
  } else if (p_value_intercept < 0.01) {
    beta_intercept <- paste0(round(beta_intercept, 2), "**")
  } else if (p_value_intercept < 0.05) {
    beta_intercept <- paste0(round(beta_intercept, 2), "*")
  } else {
    beta_intercept <- round(beta_intercept, 2)
  }
  se_intercept <- round(fixed_effects["(Intercept)", "Std. Error"], 2)
  
  # Extract random effects
  random_effects <- as.data.frame(VarCorr(model))
  random_intercept_subject <- round(random_effects[random_effects$grp == "id" &
                                                     grepl("(Intercept)", random_effects$var1) & is.na(random_effects$var2), "sdcor"], 2)
  random_intercept_abstract <- round(random_effects[random_effects$grp == "abstract" &
                                                      grepl("(Intercept)", random_effects$var1), "sdcor"], 2)
  residual <- round(random_effects[random_effects$grp == "Residual", "sdcor"], 2)
  if(model_type == "random_slope"){
    random_slope_treatment <- round(random_effects[random_effects$grp == "id" &
                                                     grepl("treatment", random_effects$var1), "sdcor"], 2)
  }else if (model_type == "random_intercept"){
    random_slope_treatment <- NA
  }
  
  # Extract R2 values
  r2 <- performance::r2(model)
  marginal_r2 <- round(r2$R2_marginal, 2)
  conditional_r2 <- round(r2$R2_conditional, 2)
  
  # Create fixed effects table
  fixed_table <- data.frame(
    column1 = c(beta_intercept, beta_treatment),
    column2 = c(se_intercept, se_treatment),
    row.names = c("Intercept", "Treatment")
  )
  
  # Placeholder for random effects header
  place_holder1 <- data.frame(column1 = "SD", column2 = NA, row.names = "NA")
  
  # Create random effects table
  randomized_table <- data.frame(
    column1 = c(
      random_intercept_subject,
      random_intercept_abstract,
      residual,
      random_slope_treatment
    ),
    column2 = c(NA, NA, NA, NA),
    row.names = c(
      "Random Intercept (Subject)",
      "Random Intercept (Abstract)",
      "Residual",
      "Random Slope (Treatment)"
    )
  )
  
  # Placeholder for model quality header
  place_holder2 <- data.frame(column1 = "Est", column2 = NA, row.names = "NA")
  
  # Create R2 table
  r2_table <- data.frame(
    column1 = c(marginal_r2, conditional_r2),
    column2 = c(NA, NA),
    row.names = c("Marginal R2", "Conditional R2")
  )
  
  # Combine tables
  results_table <- rbind(fixed_table, place_holder1, randomized_table, place_holder2, r2_table)
  
  # Add column with row names
  results_table$column0 <- rownames(results_table)
  
  # Order columns alphabetically
  results_table <- results_table[, c("column0", "column1", "column2")]
  
  return(results_table)
}

# apply function to all models
results_table1 <- extract_model_results_table(model_1.1,"random_slope")
results_table2 <- extract_model_results_table(model_2.1,"random_slope")
results_table3A <- extract_model_results_table(model_3A.1,"random_slope")
results_table3B <- extract_model_results_table(model_3B.1,"random_slope")
results_table4A <- extract_model_results_table(model_4A.1,"random_slope")
results_table4B <- extract_model_results_table(model_4B_random_intercept,"random_intercept")

# combine all results
df_part1<-data.frame(cbind(results_table1,results_table3A,results_table3B))
df_part2<-data.frame(cbind(results_table2,results_table4A,results_table4B))
model_results_overview_stage1<-data.frame(rbind(df_part1,df_part2))
str(model_results_overview_stage1)
# write to xlsx
#write_xlsx(model_results_overview_stage1, paste0(output_path,"sensitivity_summary_stage_1_model_overview.xlsx"))



##################
# Stage 2 Models #
##################


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


# Experiment 1

# Model 1.2
# Path a: treatment → FOR
model_1.2_med <- lmer(for_resp ~ treatment + (treatment | id), data = df1)
# drop random slope due to singularity
model_1.2_med <- lmer(for_resp ~ treatment + (1 | id), data = df1)
report(model_1.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_1.2_out <- lmer(change_decision ~ for_resp + treatment + (treatment | id), data = df1)
# Model with singular fit, so drop random slope
model_1.2_out <- lmer(change_decision ~ for_resp + treatment + (1 | id), data = df1)
report(model_1.2_out)

# Fit the mediation model
model_1.2 <- mediate(
  model_1.2_med, # random intercept
  model_1.2_out, # random intercept
  treat = "treatment",
  mediator = "for_resp",
  sims = 1000
)
summary(model_1.2)
sub_models_1.2<- c("random_intercept", "random_intercept")

# Experiment 2

# Model 2.2
# Path a: treatment → FOR
model_2.2_med <- lmer(for_resp ~ treatment + (treatment | id), data = df2)
## Drop random slope due to singularity
model_2.2_med <- lmer(for_resp ~ treatment + (1 | id), data = df2)
report(model_2.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_2.2_out <- lmer(change_decision ~ for_resp + treatment + (treatment | id), data = df2)
# Drop random slope due to singularity
model_2.2_out <- lmer(change_decision ~ for_resp + treatment + (1 | id), data = df2)
report(model_2.2_out)

# Fit the mediation model
model_2.2 <- mediate(
  model_2.2_med, # random intercept
  model_2.2_out, # random intercept
  treat = "treatment",
  mediator = "for_resp",
  sims = 1000
)
summary(model_2.2)

sub_models_2.2<- c("random_intercept", "random_intercept")


# Experiment 3

# Model 3A.2
# Path a: treatment → FOR
model_3A.2_med <- lmer(for_read_resp ~ treatment + (treatment | id), data = df3)
# drop random slope due to singularity
model_3A.2_med <- lmer(for_read_resp ~ treatment + (1 | id), data = df3)
report(model_3A.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_3A.2_out <- lmer(decision_change_read ~ for_read_resp + treatment + (treatment | id), data = df3)
## Model with singular fit, so drop random slope
model_3A.2_out <- lmer(decision_change_read ~ for_read_resp + treatment + (1 | id), data = df3)
# linear instead of linear mixed model: drop also random intercept model due to singularity
model_3A.2_out <- lm(decision_change_read ~ for_read_resp + treatment, data = df3)
report(model_3A.2_out)


# Fit the mediation model
model_3A.2 <- mediate(
  model_3A.2_med, # random slope
  model_3A.2_out, # no random effect
  treat = "treatment",
  mediator = "for_read_resp",
  sims = 1000
)
summary(model_3A.2)

sub_models_3A.2<- c("random_intercept", "no_random_effect")

# Model 3B.2
# Path a: treatment → FOR
model_3B.2_med <- lmer(for_cite_resp ~ treatment + (treatment | id), data = df3)
# Drop random slope due to singularity
model_3B.2_med <- lmer(for_cite_resp ~ treatment + (1 | id), data = df3)
report(model_3B.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_3B.2_out <- lmer(decision_change_cite ~ for_cite_resp + treatment + (treatment | id), data = df3)
# drop random slope due to singularity
model_3B.2_out <- lmer(decision_change_cite ~ for_cite_resp + treatment + (1 | id), data = df3)
# linear instead of linear mixed model: drop also random intercept model due to singularity
model_3B.2_out <- lm(decision_change_cite ~ for_cite_resp + treatment, data = df3)
report(model_3B.2_out)

# Fit the mediation model
model_3B.2 <- mediate(
  model_3B.2_med, # random intercept
  model_3B.2_out, # no random effect
  treat = "treatment",
  mediator = "for_cite_resp",
  sims = 1000
)
summary(model_3B.2)

sub_models_3B.2<- c("random_intercept", "no_random_effect")

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
model_1.2_med_refit<-lmer(for_resp ~ treatment + (1 | id), data = df1)
model_1.2_out_refit<-lmer(change_decision ~ for_resp + treatment + (1 | id), data = df1)

## Model 2.2
model_2.2_med_refit<-lmer(for_resp ~ treatment + (1 | id), data = df2)
model_2.2_out_refit<-lmer(change_decision ~ for_resp + treatment + (1 | id), data = df2)

## Model 3A.2
model_3A.2_med_refit<-lmer(for_read_resp ~ treatment + (1 | id), data = df3)
model_3A.2_out_refit<-lm(decision_change_read ~ for_read_resp + treatment, data = df3)

## Model 3B.2
model_3B.2_med_refit<-lmer(for_cite_resp ~ treatment + (1 | id), data = df3)
model_3B.2_out_refit<-lm(decision_change_cite ~ for_cite_resp + treatment, data = df3)

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