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



############################################################
#####                                                  #####
##### Script: Exploratory Results Stage 1 & 2          #####
#####                                                  #####
############################################################

####################
# IMPORTANT NOTE:  #
####################
# Since the analyses in this script rely on personal data,
# The exploratory analyses are NOT executable.

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

## postdoc and professor variable
# write function to apply this transformation to dataframes
transform_position <- function(df) {
  df$postdoc<-ifelse(df$position == "Post-Doctoral Researcher" | df$position == "Other: PhD/Doctorate Degree", 1, 0)
  df$professor<-ifelse(df$position == "Junior Professor" | df$position == "Professor", 1, 0)
  return(df)
}
df1_2 <- transform_position(df1_2)
df3_4 <- transform_position(df3_4)
df1 <- transform_position(df1)
df2 <- transform_position(df2)
df3 <- transform_position(df3)
df4 <- transform_position(df4)

# or operator | is not vectorized

######################
## Center Predictors #
######################

# Recode pressure_to_publish
recode_pressure_to_publish <- function(df) {
  # Check if the column exists in the dataset
  if (!"pressure_to_publish" %in% colnames(df)) {
    stop("Error: The dataset does not contain a column named 'pressure_to_publish'.")
  }
  
  # Recode the variable
  df$pressure_to_publish_numeric <- as.numeric(factor(df$pressure_to_publish, 
                                                      levels = c("Not at all pressured", "Neither pressured nor unpressured", 
                                                                 "Somewhat pressured", "Pressured", "Very pressured")))
  
  return(df)
}

# apply function to df1, df2, df3, and df4
df1 <- recode_pressure_to_publish(df1)
df2 <- recode_pressure_to_publish(df2)
df3 <- recode_pressure_to_publish(df3)
df4 <- recode_pressure_to_publish(df4)

# Grand Mean Centering
grand_mean_center <- function(df, variables) {
  for (var in variables) {
    if (!var %in% colnames(df)) {
      stop(paste("Error: The dataset does not contain a column named", var))
    }
    df[[paste0(var, "_centered")]] <- df[[var]] - mean(df[[var]], na.rm = TRUE)
  }
  return(df)
}

# List of variables to grand mean center
vars_to_center <- c("pressure_to_publish_numeric", "bfi_2_S_score", "familiarity_PB", "trial")

# Apply the function to df1, df2, df3, and df4
df1 <- grand_mean_center(df1, vars_to_center)
df2 <- grand_mean_center(df2, vars_to_center)
df3 <- grand_mean_center(df3, vars_to_center)
df4 <- grand_mean_center(df4, vars_to_center)

# summary function
round_coef <- function(model) {
  coef <- summary(model)$coefficients
  coef[, 1] <- round(coef[, 1], 2)
  coef[, 2] <- round(coef[, 2], 2)
  coef[, 4] <- round(coef[, 4], 2)
  coef[, 5] <- round(coef[, 5], 3)
  
  return(coef)
}



##################################
# Stage 1 Models #
##################

## Experiment 1
model_1.1.4 <- lmer(
  decision1_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df1
)
round_coef(model_1.1.4)

# Experiment 2
model_2.1.4 <- lmer(
  decision1_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df2
)
round_coef(model_2.1.4)

# Experiment 3
# Reading
model_3A.1.4 <- lmer(
  decision1_read_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df3
)
round_coef(model_3A.1.4)

# Citing
model_3B.1.4 <- lmer(
  decision1_cite_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df3
)
round_coef(model_3B.1.4)


# Experiment 4
# Reading
model_4A.1.4 <- lmer(
  decision1_read_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df4
)
round_coef(model_4A.1.4)

# Citing
model_4B.1.4 <- lmer(
  decision1_cite_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df4
)
# drop random slope due to singularity
model_4B.1.4 <- lmer(
  decision1_cite_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (1 | id) + (1 | abstract),
  data = df4
)
round_coef(model_4B.1.4)





##################
# Stage 2 Models #
##################
# Fit models using df_clean to assure that both data.frames have the same size for mediation

# Experiment1
df_clean1 <- df1 %>% 
  filter(!is.na(for_resp), !is.na(change_decision), !is.na(treatment),
         !is.na(familiarity_PB_centered), !is.na(pressure_to_publish_numeric_centered),
         !is.na(bfi_2_S_score_centered), !is.na(nd), !is.na(uk), !is.na(us), !is.na(can),
         !is.na(aus), !is.na(other), !is.na(trial_centered))

# Experiment 2
df_clean2 <- df2 %>% 
  filter(!is.na(for_resp), !is.na(change_decision), !is.na(treatment),
         !is.na(familiarity_PB_centered), !is.na(pressure_to_publish_numeric_centered),
         !is.na(bfi_2_S_score_centered), !is.na(nd), !is.na(uk), !is.na(us), !is.na(can),
         !is.na(aus), !is.na(other), !is.na(trial_centered))

# Experiment 3 - Reading
df_clean3_read <- df3 %>% 
  filter(!is.na(for_read_resp), !is.na(decision_change_read), !is.na(treatment),
         !is.na(familiarity_PB_centered), !is.na(pressure_to_publish_numeric_centered),
         !is.na(bfi_2_S_score_centered), !is.na(nd), !is.na(uk), !is.na(us), !is.na(can),
         !is.na(aus), !is.na(other), !is.na(trial_centered))

# Experiment 3 - Citing
df_clean3_cite <- df3 %>% 
  filter(!is.na(for_cite_resp), !is.na(decision_change_cite), !is.na(treatment),
         !is.na(familiarity_PB_centered), !is.na(pressure_to_publish_numeric_centered),
         !is.na(bfi_2_S_score_centered), !is.na(nd), !is.na(uk), !is.na(us), !is.na(can),
         !is.na(aus), !is.na(other), !is.na(trial_centered))

# Experiment 4 - Reading
df_clean4_read <- df4 %>% 
  filter(!is.na(for_read_resp), !is.na(decision_change_read), !is.na(treatment),
         !is.na(familiarity_PB_centered), !is.na(pressure_to_publish_numeric_centered),
         !is.na(bfi_2_S_score_centered), !is.na(nd), !is.na(uk), !is.na(us), !is.na(can),
         !is.na(aus), !is.na(other), !is.na(trial_centered))

# Experiment 4 - Citing
df_clean4_cite <- df4 %>% 
  filter(!is.na(for_cite_resp), !is.na(decision_change_cite), !is.na(treatment),
         !is.na(familiarity_PB_centered), !is.na(pressure_to_publish_numeric_centered),
         !is.na(bfi_2_S_score_centered), !is.na(nd), !is.na(uk), !is.na(us), !is.na(can),
         !is.na(aus), !is.na(other), !is.na(trial_centered))


################
## Fit Models ##
################

# Experiment 1
## med
model_1.2.2_med <- lmer(for_resp ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + (treatment | id),
                        data = df_clean1)
report(model_1.2.2_med)

## out
model_1.2.2_out <- lmer(change_decision ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (treatment | id),
                        data = df_clean1)
### Model with singular fit, so drop random slope
model_1.2.2_out <- lmer(change_decision ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (1 | id),
                        data = df_clean1)
report(model_1.2.2_out)


## Fit mediation model
model_1.2.2 <- mediation::mediate(
  model_1.2.2_med,
  model_1.2.2_out,
  treat = "treatment",
  mediator = "for_resp",
  sims = 1000
)

# Experiment 2
## med
model_2.2.2_med <- lmer(for_resp ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + (treatment | id),
                        data = df_clean2)
### model with singular fit, so drop random slope
model_2.2.2_med <- lmer(for_resp ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + (1 | id),
                        data = df_clean2)
report(model_2.2.2_med)

## out
model_2.2.2_out <- lmer(change_decision ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (treatment | id),
                        data = df_clean2)
report(model_2.2.2_out)

## Fit mediation model
model_2.2.2 <- mediation::mediate(
  model_2.2.2_med,
  model_2.2.2_out,
  treat = "treatment",
  mediator = "for_resp",
  sims = 1000
)


# Experiment 3 - Reading
## med
model_3A.2.2_med <- lmer(for_read_resp ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered +  (treatment | id) ,
                        data = df_clean3_read)
report(model_3A.2.2_med)

## out
model_3A.2.2_out <- lmer(decision_change_read ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_read_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (treatment | id),
                        data = df_clean3_read)
### model with singular fit, so drop random slope
model_3A.2.2_out <- lmer(decision_change_read ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_read_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (1 | id),
                        data = df_clean3_read)
report(model_3A.2.2_out)

## Fit mediation model
model_3A.2.2 <- mediation::mediate(
  model_3A.2.2_med,
  model_3A.2.2_out,
  treat = "treatment",
  mediator = "for_read_resp",
  sims = 1000
)

# Experiment 3 - Citing
## med
model_3B.2.2_med <- lmer(for_cite_resp ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + (treatment | id),
                        data = df_clean3_cite)
report(model_3B.2.2_med)

## out
model_3B.2.2_out <- lmer(decision_change_cite ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_cite_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (treatment | id),
                        data = df_clean3_cite)
### Model with singular fit, so drop random slope
model_3B.2.2_out <- lmer(decision_change_cite ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_cite_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (1 | id),
                        data = df_clean3_cite)
### Still singular fit, so drop random intercept as well --> full linear fixed regression model
model_3B.2.2_out <- lm(decision_change_cite ~ treatment * familiarity_PB_centered +
                         treatment * pressure_to_publish_numeric_centered + for_cite_resp + bfi_2_S_score_centered +
                         nd + uk + us + can + aus + other +
                         trial_centered,
                       data = df_clean3_cite)
report(model_3B.2.2_out)

## Fit mediation model
model_3B.2.2 <- mediation::mediate(
  model_3B.2.2_med,
  model_3B.2.2_out,
  treat = "treatment",
  mediator = "for_cite_resp",
  sims = 1000
)


# Experiment 4 - Reading
## med
model_4A.2.2_med <- lmer(for_read_resp ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + (treatment | id),
                        data = df_clean4_read)
### Model with singular fit, so drop random slope
model_4A.2.2_med <- lmer(for_read_resp ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + (1 | id),
                        data = df_clean4_read)
report(model_4A.2.2_med)

## out
model_4A.2.2_out <- lmer(decision_change_read ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_read_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (treatment | id),
                        data = df_clean4_read)
### Model with singular fit, so drop random slope
model_4A.2.2_out <- lmer(decision_change_read ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_read_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (1 | id),
                        data = df_clean4_read)
report(model_4A.2.2_out)

## Fit mediation model
model_4A.2.2 <- mediation::mediate(
  model_4A.2.2_med,
  model_4A.2.2_out,
  treat = "treatment",
  mediator = "for_read_resp",
  sims = 1000
)

# Experiment 4 - Citing
## med
model_4B.2.2_med <- lmer(for_cite_resp ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + (treatment | id),
                        data = df_clean4_cite)
report(model_4B.2.2_med)

## out
model_4B.2.2_out <- lmer(decision_change_cite ~ treatment * familiarity_PB_centered +
                          treatment * pressure_to_publish_numeric_centered + for_cite_resp + bfi_2_S_score_centered +
                          postdoc + professor +nd + uk + us + can + aus + other +
                          trial_centered + (treatment | id),
                        data = df_clean4_cite)
report(model_4B.2.2_out)

## Fit mediation model
model_4B.2.2 <- mediation::mediate(
  model_4B.2.2_med,
  model_4B.2.2_out,
  treat = "treatment",
  mediator = "for_cite_resp",
  sims = 1000
)

# Print verbose summaries for all models
summary(model_1.2.2, verbose = TRUE)
summary(model_2.2.2, verbose = TRUE)
summary(model_3A.2.2, verbose = TRUE)
summary(model_3B.2.2, verbose = TRUE)
summary(model_4A.2.2, verbose = TRUE)
summary(model_4B.2.2, verbose = TRUE)



# print significant results
print_significant_rows <- function(model) {
  # Extract significant rows
  significant_rows <- model[model[, "Pr(>|t|)"] < 0.05 & row.names(model)!="(Intercept)", , drop = FALSE]
  
  # If no significant rows, return NA
  if (nrow(significant_rows) == 0) {
    return(NA)
  }
  else {
    # round everything to 2 decimal places but p-values 3 decimals
    significant_rows <- round(significant_rows, 2)
    significant_rows[, "Pr(>|t|)"] <- round(significant_rows[, "Pr(>|t|)"], 3)
    return(significant_rows)
  }
}

library(lmerTest)
s_1.2.2_med<-summary(model_1.2.2_med)$coefficients
s_1.2.2_out<-summary(model_1.2.2_out)$coefficients
s_2.2.2_med<-summary(model_2.2.2_med)$coefficients
s_2.2.2_out<-summary(model_2.2.2_out)$coefficients
s_3A.2.2_med<-summary(model_3A.2.2_med)$coefficients
s_3A.2.2_out<-summary(model_3A.2.2_out)$coefficients
s_3B.2.2_med<-summary(model_3B.2.2_med)$coefficients
s_3B.2.2_out<-summary(model_3B.2.2_out)$coefficients
s_4A.2.2_med<-summary(model_4A.2.2_med)$coefficients
s_4A.2.2_out<-summary(model_4A.2.2_out)$coefficients
s_4B.2.2_med<-summary(model_4B.2.2_med)$coefficients
s_4B.2.2_out<-summary(model_4B.2.2_out)$coefficients


# print significant results
print_significant_rows(s_1.2.2_med)
print_significant_rows(s_1.2.2_out)
print_significant_rows(s_2.2.2_med)
print_significant_rows(s_2.2.2_out)
print_significant_rows(s_3A.2.2_med)
print_significant_rows(s_3A.2.2_out)
print_significant_rows(s_3B.2.2_med)
print_significant_rows(s_3B.2.2_out)
print_significant_rows(s_4A.2.2_med)
print_significant_rows(s_4A.2.2_out)
print_significant_rows(s_4B.2.2_med)


