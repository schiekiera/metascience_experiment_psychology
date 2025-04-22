# ==============================================================================
# remove everything
rm(list = ls())

# load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
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
##### Script: Model Results Stage 1                #####
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

results_4B.1 <- extract_model_stats_random_slope_model(model_4B.1,
                                                       data = df4,
                                                       predictor = "treatment",
                                                       outcome = "decision1_cite_resp")
# Singularity --> drop random slopes
model_4B_random_intercept <-
  lmer(decision1_cite_resp ~ treatment + (1 | id) + (1 |
                                                       abstract), data = df4)
results_4B_random_intercept <- extract_model_stats_random_intercept_model(
  model_4B_random_intercept,
  data = df4,
  predictor = "treatment",
  outcome = "decision1_cite_resp"
)


## qq plots
# QQ plot for model_1.1
#png(paste0(plot_path,"qqplot_model_1.1.png"), width = 1000, height = 1000, res = 300)
qqnorm(resid(model_1.1),main="Model 1.1")
#dev.off()

# QQ plot for model_2.1
#png(paste0(plot_path,"qqplot_model_2.1.png"), width = 1000, height = 1000, res = 300)
qqnorm(resid(model_2.1),main="Model 2.1")
#dev.off()

# QQ plot for model_3A.1
#png(paste0(plot_path,"qqplot_model_3A.1.png"), width = 1000, height = 1000, res = 300)
qqnorm(resid(model_3A.1),main="Model 3A.1")
#dev.off()

# QQ plot for model_3B.1
#png(paste0(plot_path,"qqplot_model_3B.1.png"), width = 1000, height = 1000, res = 300)
qqnorm(resid(model_3B.1),main="Model 3B.1")
#dev.off()

# QQ plot for model_4A.1
#png(paste0(plot_path,"qqplot_model_4A.1.png"), width = 1000, height = 1000, res = 300)
qqnorm(resid(model_4A.1),main="Model 4A.1")
#dev.off()

# QQ plot for model_4B_random_intercept
#png(paste0(plot_path,"qqplot_model_4B.1.png"), width = 1000, height = 1000, res = 300)
qqnorm(resid(model_4B_random_intercept),main="Model 4B.1")
#dev.off()




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


# create data frame for effect size plot
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
#ggsave(paste0(plot_path,"figure4.png"), plot = figure4, width = 15, height = 10, units = "cm")
# saves as tiff file
#ggsave(paste0(plot_path,"figure4.tiff"), plot = figure4, width = 15, height = 10, units = "cm", dpi = 300)


# define function to extract model results table
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

# print model results
print(model_results_overview_stage1)

# write to xlsx
#write_xlsx(model_results_overview_stage1, paste0(output_path,"summary_stage_1_model_overview.xlsx"))