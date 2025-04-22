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
##### Script: Appendix Results                      #####
#####                                               #####
#########################################################

####################
# IMPORTANT NOTE:  #
####################
# Since the analyses in this script rely on personal data (covariate models),
# The analyses are NOT fully executable.

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


#################################################
##                                             ##
## R E S U L T S  T A B L E  F U N C T I O N S ##
##                                             ##
#################################################


# 1.  extract model results table
extract_model_results_table <- function(model, model_type) {
  # Extract fixed effects
  fixed_effects <- summary(model)$coefficients
  
  # Format beta and SE values with significance stars
  format_beta <- function(estimate, p_value) {
    stars <- ifelse(p_value < 0.001, "***",
                    ifelse(p_value < 0.01, "**",
                           ifelse(p_value < 0.05, "*", "")))
    paste0(round(estimate, 2), stars)
  }
  
  # Initialize fixed effects table
  fixed_table <- data.frame(
    column0 = character(),
    column1 = character(),
    column2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through all fixed effects
  for (effect_name in rownames(fixed_effects)) {
    beta <- format_beta(fixed_effects[effect_name, "Estimate"],
                        fixed_effects[effect_name, "Pr(>|t|)"])
    se <- round(fixed_effects[effect_name, "Std. Error"], 2)
    fixed_table <- rbind(fixed_table, data.frame(
      column0 = effect_name,
      column1 = beta,
      column2 = se,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add placeholder row to split fixed and random effects
  place_holder1 <- data.frame(column0 = "Random Effects", column1 = "SD", column2 = NA)
  
  # Extract random effects
  random_effects <- as.data.frame(VarCorr(model))
  
  # Create random effects table
  random_rows <- list()
  if ("id" %in% random_effects$grp) {
    intercept_subject <- random_effects[random_effects$grp == "id" & random_effects$var1 == "(Intercept)", "sdcor"]
    if (length(intercept_subject) > 0) {
      random_rows[["Random Intercept (Subject)"]] <- round(intercept_subject, 2)
    }
    if (model_type == "random_slope") {
      slope <- random_effects[random_effects$grp == "id" & grepl("treatment", random_effects$var1), "sdcor"]
      if (length(slope) > 0) {
        random_rows[["Random Slope (Treatment)"]] <- round(slope, 2)
      }
    }
  }
  if ("abstract" %in% random_effects$grp) {
    intercept_abstract <- random_effects[random_effects$grp == "abstract" & random_effects$var1 == "(Intercept)", "sdcor"]
    if (length(intercept_abstract) > 0) {
      random_rows[["Random Intercept (Abstract)"]] <- round(intercept_abstract, 2)
    }
  }
  if ("Residual" %in% random_effects$grp) {
    residual <- random_effects[random_effects$grp == "Residual", "sdcor"]
    if (length(residual) > 0) {
      random_rows[["Residual"]] <- round(residual, 2)
    }
  }
  
  # Format random effects into a dataframe
  random_table <- do.call(rbind, lapply(names(random_rows), function(name) {
    data.frame(column0 = name, column1 = random_rows[[name]], column2 = NA)
  }))
  
  # Placeholder for model quality section
  place_holder2 <- data.frame(column0 = "Model Quality", column1 = "Est", column2 = NA)
  
  # Model fit metrics
  r2_vals <- performance::r2(model)
  r2_table <- data.frame(
    column0 = c("Marginal R2", "Conditional R2"),
    column1 = round(c(r2_vals$R2_marginal, r2_vals$R2_conditional), 2),
    column2 = NA
  )
  
  # Combine all sections
  results_table <- rbind(fixed_table, place_holder1, random_table, place_holder2, r2_table)
  
  return(results_table)
}

# 2.  align model tables for export
align_model_tables_for_export <- function(model_list, model_types, model_labels) {
  # Extract model results
  model_results <- mapply(extract_model_results_table, model_list, model_types, SIMPLIFY = FALSE)
  
  # Get union of all unique rownames across all models
  all_rownames <- unique(unlist(lapply(model_results, function(df) df$column0)))
  
  # Align each model table to full rowname set
  aligned_tables <- lapply(seq_along(model_results), function(i) {
    df <- model_results[[i]]
    label <- model_labels[i]
    
    beta_vec <- vapply(all_rownames, function(r) {
      matches <- df[df$column0 == r, "column1"]
      if (length(matches) > 0) as.character(matches[1]) else ""
    }, character(1))
    
    se_vec <- vapply(all_rownames, function(r) {
      matches <- df[df$column0 == r, "column2"]
      if (length(matches) > 0) as.character(matches[1]) else ""
    }, character(1))
    
    aligned <- data.frame(
      rownames = all_rownames,
      beta = beta_vec,
      se = se_vec,
      stringsAsFactors = FALSE
    )
    
    # Rename columns
    colnames(aligned) <- c(
      paste0("Parameter_", label),
      paste0("Beta_", label),
      paste0("SE_", label)
    )
    
    return(aligned)
  })
  
  # Combine all model tables column-wise
  final_df <- do.call(cbind, aligned_tables)
  return(final_df)
}

# function to extract AIC, BIC, and logLik
get_aic_bic_loglik <- function(model) {
  aic <- round(AIC(model),1)
  bic <- round(BIC(model),1)
  loglik <- round(as.numeric(logLik(model)),1)
  return(c(aic, bic, loglik))
}

# extract model comparisons
extract_model_comparisons <- function(model_list, model_labels) {
  comparison_table <- data.frame(
    Comparison = character(),
    Chisq = numeric(),
    Df = numeric(),
    p = character(),
    AIC = numeric(),
    BIC = numeric(),
    logLik = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 2:length(model_list)) {
    m1 <- model_list[[i - 1]]
    m2 <- model_list[[i]]
    
    cmp <- anova(m1, m2)
    
    # Safely extract comparison info
    chisq <- round(cmp$Chisq[2], 2)
    df <- cmp$Df[2]
    pval <- cmp$`Pr(>Chisq)`[2]
    
    # Format p-value nicely
    p <- ifelse(pval < 0.001, "<.001", format(round(pval, 3), nsmall = 3))
    
    # AIC/BIC/logLik from larger model
    aic <- round(AIC(m2), 1)
    bic <- round(BIC(m2), 1)
    loglik <- round(logLik(m2), 2)
    
    comparison_table <- rbind(comparison_table, data.frame(
      Comparison = paste(model_labels[i - 1], "vs.", model_labels[i]),
      Chisq = chisq,
      Df = df,
      p = p,
      AIC = aic,
      BIC = bic,
      logLik = loglik
    ))
  }
  
  return(comparison_table)
}


# append model comparisons
append_model_comparison<-function(df, model_list, model_labels){
  comparison_df <- extract_model_comparisons(model_list[1:3], model_labels)
  matrix<-matrix(NA, nrow = 7, ncol = 9)
  matrix[1:3,2]<-get_aic_bic_loglik(model_list[[1]])
  matrix[1:3,4]<-get_aic_bic_loglik(model_list[[2]])
  matrix[1:3,6]<-get_aic_bic_loglik(model_list[[3]])
  matrix[1:3,8]<-get_aic_bic_loglik(model_list[[4]])
  matrix[5:7,4]<-data.frame(t(comparison_df))$X1[2:4]
  matrix[5:7,6]<-data.frame(t(comparison_df))$X2[2:4]
  df_comp<-data.frame(matrix)
  colnames(df_comp)<-colnames(df)
  df_comp[,1]<-c("AIC", "BIC", "logLik", "Comparison", "chi", "df", "p")
  df<-rbind(df,df_comp)
  return(df)
}


##################
# Stage 1 Models #
##################

############
# Model 1 ##
############

# Model 1.1.1
# Null model with Fixed Effects only (without predictor)
model_1.1.1 <-
  lmer(decision1_resp ~ 1 + (1 | id) + (1 |
                                                  abstract), data = df1)
## Model 1.1.2
# Adding treatment as a predictor
model_1.1.2 <-
  lmer(decision1_resp ~ treatment + (1 | id) + (1 |
                                                          abstract), data = df1)
## Model 1.1.3
## Adding a random slope for treatment within participants (id)
model_1.1.3 <-
  lmer(decision1_resp ~ treatment + (treatment | id) + (1 |
                                                          abstract), data = df1)
## Model 1.1.4
## Exploratory covariate model
model_1.1.4 <- lmer(
  decision1_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df1
)

############
# Model 2 ##
############

# Model 2.1.1
# Null model with Fixed Effects only (without predictor)
model_2.1.1 <-
  lmer(decision1_resp ~ 1 + (1 | id) + (1 |
                                                  abstract), data = df2)

## Model 2.1.2
# Adding treatment as a predictor
model_2.1.2 <-
  lmer(decision1_resp ~ treatment + (1 | id) + (1 |
                                                  abstract), data = df2)

## Model 2.1.3
## Adding a random slope for treatment within participants (id)
model_2.1.3 <-
  lmer(decision1_resp ~ treatment + (treatment | id) + (1 |
                                                          abstract), data = df2)

## Model 2.1.4
## Exploratory covariate model
model_2.1.4 <- lmer(
  decision1_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df2
)

summary(model_2.1.4)
#############
# Model 3A ##
#############

# Model 3A.1.1
# Null model with Fixed Effects only (without predictor)
model_3A.1.1 <-
  lmer(decision1_read_resp ~ 1 + (1 | id) + (1 |
                                          abstract), data = df3)

## Model 3A.1.2
# Adding treatment as a predictor
model_3A.1.2 <-
  lmer(decision1_read_resp ~ treatment + (1 | id) + (1 |
                                                  abstract), data = df3)

## Model 3.1A: reading
## Adding a random slope for treatment within participants (id)
model_3A.1.3 <-
  lmer(decision1_read_resp ~ treatment + (treatment | id) + (1 |
                                                               abstract),
       data = df3)

## Model 3A.1.4
## Exploratory covariate model
model_3A.1.4 <- lmer(
  decision1_read_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus  + other +
    (treatment | id) + (1 | abstract),
  data = df3
)


#############
# Model 3B ##
#############

# 3B.1.1 Null
model_3B.1.1 <- lmer(decision1_cite_resp ~ 1 + (1 | id) + (1 | abstract), data = df3)

# 3B.1.2 Fixed only
model_3B.1.2 <- lmer(decision1_cite_resp ~ treatment + (1 | id) + (1 | abstract), data = df3)

# 3B.1.3 Random slope — already exists
model_3B.1.3 <-
  lmer(decision1_cite_resp ~ treatment + (treatment | id) + (1 |
                                                               abstract),
       data = df3)

# 3B.1.4 Covariates
model_3B.1.4 <- lmer(
  decision1_cite_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus + other +
    (treatment | id) + (1 | abstract),
  data = df3
)


##############
# Model 4A1 ##
##############

# 4A.1.1 Null model
model_4A.1.1 <- lmer(decision1_read_resp ~ 1 + (1 | id) + (1 | abstract), data = df4)

# 4A.1.2 Random intercept
model_4A.1.2 <- lmer(decision1_read_resp ~ treatment + (1 | id) + (1 | abstract), data = df4)

# 4A.1.3 Random slope
model_4A.1.3 <-
  lmer(decision1_read_resp ~ treatment + (treatment | id) + (1 |
                                                               abstract),
       data = df4)

# 4A.1.4 Covariates
model_4A.1.4 <- lmer(
  decision1_read_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus + other +
    (treatment | id) + (1 | abstract),
  data = df4
)

summary(model_4A.1.4)
##############
# Model 4B1 ##
##############

# 4B.1.1 Null
model_4B.1.1 <- lmer(decision1_cite_resp ~ 1 + (1 | id) + (1 | abstract), data = df4)

# 4B.1.3 Random intercept
model_4B.1.2 <- lmer(decision1_cite_resp ~ treatment + (1 | id) + (1 | abstract), data = df4)

# 4B.1.3 Random slope
model_4B.1.3 <-
  lmer(decision1_cite_resp ~ treatment + (treatment |
                                            id) + (1 | abstract),
       data = df4)

# 4B.1.3 Covariates --> drop random slope
model_4B.1.4 <- lmer(
  decision1_cite_resp ~ treatment * familiarity_PB_centered +
    treatment * pressure_to_publish_numeric_centered +
    postdoc + professor + trial_centered +
    nd + uk + us + can + aus + other +
    (1 | id) + (1 | abstract),
  data = df4
)

# sort vector
sort_vector <- c(0, 15, 16, 17, 18, 20, 21, 22, 1, 19, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
length(sort_vector)

# Model 1
models1 <- list(model_1.1.1, model_1.1.2, model_1.1.3, model_1.1.4)
types1 <- c("random_intercept", "random_intercept", "random_slope", "random_slope")
labels1 <- c("1.1.1", "1.1.2", "1.1.3", "1.1.4")
export_df1 <- align_model_tables_for_export(models1, types1, labels1)
export_df1$sort_vector <- sort_vector
export_df1 <- export_df1[order(export_df1$sort_vector), ]
export_df1<-export_df1[,c(-4,-7,-10,-13),]
export_df1<-append_model_comparison(export_df1, models1, labels1)
#write_xlsx(export_df1, file.path(output_path, "Appendix/Stage1/appendix_model1_summary_stage_1_model_overview.xlsx"))


# Model 2
models2 <- list(model_2.1.1, model_2.1.2, model_2.1.3, model_2.1.4)
types2 <- c("random_intercept", "random_intercept", "random_slope", "random_slope")
labels2 <- c("2.1.1", "2.1.2", "2.1.3", "2.1.4")
export_df2 <- align_model_tables_for_export(models2, types2, labels2)
export_df2$sort_vector <- sort_vector
export_df2 <- export_df2[order(export_df2$sort_vector), ]
export_df2<-export_df2[,c(-4,-7,-10,-13),]
export_df2<-append_model_comparison(export_df2, models2, labels2)
#write_xlsx(export_df2, file.path(output_path, "Appendix/Stage1/appendix_model2_summary_stage_1_model_overview.xlsx"))

# Model 3A
models3A <- list(model_3A.1.1, model_3A.1.2, model_3A.1.3, model_3A.1.4)
types3A <- c("random_intercept", "random_intercept", "random_slope", "random_slope")
labels3A <- c("3A.1.1", "3A.1.2", "3A.1.3", "3A.1.4")
export_df3A <- align_model_tables_for_export(models3A, types3A, labels3A)
export_df3A$sort_vector <- sort_vector
export_df3A <- export_df3A[order(export_df3A$sort_vector), ]
export_df3A<-export_df3A[,c(-4,-7,-10,-13),]
export_df3A<-append_model_comparison(export_df3A, models3A, labels3A)
#write_xlsx(export_df3A, file.path(output_path, "Appendix/Stage1/appendix_model3A_summary_stage_1_model_overview.xlsx"))

# Model 3B
models3B <- list(model_3B.1.1, model_3B.1.2, model_3B.1.3, model_3B.1.4)
types3B <- c("random_intercept", "random_intercept", "random_slope", "random_slope")
labels3B <- c("3B.1.1", "3B.1.2", "3B.1.3", "3B.1.4")
export_df3B <- align_model_tables_for_export(models3B, types3B, labels3B)
export_df3B$sort_vector <- sort_vector
export_df3B <- export_df3B[order(export_df3B$sort_vector), ]
export_df3B<-export_df3B[,c(-4,-7,-10,-13),]
export_df3B<-append_model_comparison(export_df3B, models3B, labels3B)
#write_xlsx(export_df3B, file.path(output_path, "Appendix/Stage1/appendix_model3B_summary_stage_1_model_overview.xlsx"))

# Model 4A
models4A <- list(model_4A.1.1, model_4A.1.2, model_4A.1.3, model_4A.1.4)
types4A <- c("random_intercept", "random_intercept", "random_slope", "random_slope")
labels4A <- c("4A.1.1", "4A.1.2", "4A.1.3", "4A.1.4")
export_df4A <- align_model_tables_for_export(models4A, types4A, labels4A)
export_df4A$sort_vector <- sort_vector
export_df4A <- export_df4A[order(export_df4A$sort_vector), ]
export_df4A<-export_df4A[,c(-4,-7,-10,-13),]
export_df4A<-append_model_comparison(export_df4A, models4A, labels4A)
#write_xlsx(export_df4A, file.path(output_path, "Appendix/Stage1/appendix_model4A_summary_stage_1_model_overview.xlsx"))

# Model 4B
models4B <- list(model_4B.1.1, model_4B.1.2, model_4B.1.3, model_4B.1.4)
types4B <- c("random_intercept", "random_intercept", "random_slope", "random_slope")
labels4B <- c("4B.1.1", "4B.1.2", "4B.1.3", "4B.1.4")
export_df4B <- align_model_tables_for_export(models4B, types4B, labels4B)
export_df4B$sort_vector <- sort_vector
export_df4B <- export_df4B[order(export_df4B$sort_vector), ]
export_df4B<-export_df4B[,c(-4,-7,-10,-13),]
export_df4B<-append_model_comparison(export_df4B, models4B, labels4B)
#write_xlsx(export_df4B, file.path(output_path, "Appendix/Stage1/appendix_model4B_summary_stage_1_model_overview.xlsx"))



################
## Functions  ##
################

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
obtain_model_parameters_main <- function(model_name,
                                    mediation_model,
                                    path_a_model,
                                    path_b_c_model,
                                    for_type) {
  
  
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
  # Check if it's a linear model or mixed model and handle accordingly
  if(inherits(path_b_c_model, "lm")) {
    beta_for_resp_b <-  round_normal(fixed_effects_b_c[for_type, 1])
    p_value_for_resp_b <- fixed_effects_b_c[for_type, 4]
  } else {
    beta_for_resp_b <-  round_normal(fixed_effects_b_c[for_type, "Estimate"])
    p_value_for_resp_b <- fixed_effects_b_c[for_type, "Pr(>|t|)"]
  }
  beta_p_for_resp_b<-process_beta(beta_for_resp_b, p_value_for_resp_b)
  ci_for_resp_b<- round_normal(confidence_intervals_b_c[for_type,])
  ci_for_resp_b<-paste0("[",paste0(ci_for_resp_b, collapse = ", "),"]")
  ### path c: treatment
  if(inherits(path_b_c_model, "lm")) {
    beta_treatment_c <-  round_normal(fixed_effects_b_c["treatment", 1])
    p_value_treatment_c <- fixed_effects_b_c["treatment", 4]
  } else {
    beta_treatment_c <-  round_normal(fixed_effects_b_c["treatment", "Estimate"])
    p_value_treatment_c <- fixed_effects_b_c["treatment", "Pr(>|t|)"]
  }
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
    beta_p = c(beta_p_treatment_a, NA,ade_beta_p, acme_beta_p, beta_p_for_resp_b, total_effect_beta_p, prop_mediated_beta_p),
    ci = c(ci_treatment_a, NA, ade_ci, acme_ci, ci_for_resp_b, total_effect_ci, prop_mediated_ci),
    row.names = c("Path a", "empty1","ADE", "ACME", "Path b", "Total Effect", "Proportion Mediated")
  )
  
  return(results_table)
}

# define function to obtain exploratory model parameters
obtain_model_parameters_exploratory <- function(model_name,
                                           mediation_model,
                                           path_a_model,
                                           path_b_c_model,
                                           for_type) {
  
  # Path a
  ## Define fixed_effects_a and confidence_intervals_a
  fixed_effects_a <- summary(path_a_model)$coefficients
  confidence_intervals_a<-confint(path_a_model, method = "Wald")
  ## Extract beta, p-value and confidence intervals for treatment
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
  # Check if it's a linear model or mixed model and handle accordingly
  if(inherits(path_b_c_model, "lm")) {
    beta_for_resp_b <-  round_normal(fixed_effects_b_c[for_type, 1])
    p_value_for_resp_b <- fixed_effects_b_c[for_type, 4]
  } else {
    beta_for_resp_b <-  round_normal(fixed_effects_b_c[for_type, "Estimate"])
    p_value_for_resp_b <- fixed_effects_b_c[for_type, "Pr(>|t|)"]
  }
  beta_p_for_resp_b<-process_beta(beta_for_resp_b, p_value_for_resp_b)
  ci_for_resp_b<- round_normal(confidence_intervals_b_c[for_type,])
  ci_for_resp_b<-paste0("[",paste0(ci_for_resp_b, collapse = ", "),"]")
  ### path c: treatment
  if(inherits(path_b_c_model, "lm")) {
    beta_treatment_c <-  round_normal(fixed_effects_b_c["treatment", 1])
    p_value_treatment_c <- fixed_effects_b_c["treatment", 4]
  } else {
    beta_treatment_c <-  round_normal(fixed_effects_b_c["treatment", "Estimate"])
    p_value_treatment_c <- fixed_effects_b_c["treatment", "Pr(>|t|)"]
  }
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
  
  # Extract additional parameters from mediator model (path_a_model)
  ## Pressure to Publish
  beta_pressure_med <- round_normal(fixed_effects_a["pressure_to_publish_numeric_centered", "Estimate"])
  p_value_pressure_med <- fixed_effects_a["pressure_to_publish_numeric_centered", "Pr(>|t|)"]
  beta_p_pressure_med <- process_beta(beta_pressure_med, p_value_pressure_med)
  ci_pressure_med <- round_normal(confidence_intervals_a["pressure_to_publish_numeric_centered",])
  ci_pressure_med <- paste0("[",paste0(ci_pressure_med, collapse = ", "),"]")
  
  ## Familiarity with PB
  beta_familiarity_med <- round_normal(fixed_effects_a["familiarity_PB_centered", "Estimate"])
  p_value_familiarity_med <- fixed_effects_a["familiarity_PB_centered", "Pr(>|t|)"]
  beta_p_familiarity_med <- process_beta(beta_familiarity_med, p_value_familiarity_med)
  ci_familiarity_med <- round_normal(confidence_intervals_a["familiarity_PB_centered",])
  ci_familiarity_med <- paste0("[",paste0(ci_familiarity_med, collapse = ", "),"]")
  
  ## Pressure to Publish * treatment
  beta_pressure_treat_med <- round_normal(fixed_effects_a["treatment:pressure_to_publish_numeric_centered", "Estimate"])
  p_value_pressure_treat_med <- fixed_effects_a["treatment:pressure_to_publish_numeric_centered", "Pr(>|t|)"]
  beta_p_pressure_treat_med <- process_beta(beta_pressure_treat_med, p_value_pressure_treat_med)
  ci_pressure_treat_med <- round_normal(confidence_intervals_a["treatment:pressure_to_publish_numeric_centered",])
  ci_pressure_treat_med <- paste0("[",paste0(ci_pressure_treat_med, collapse = ", "),"]")
  
  ## Familiarity with PB * treatment
  beta_familiarity_treat_med <- round_normal(fixed_effects_a["treatment:familiarity_PB_centered", "Estimate"])
  p_value_familiarity_treat_med <- fixed_effects_a["treatment:familiarity_PB_centered", "Pr(>|t|)"]
  beta_p_familiarity_treat_med <- process_beta(beta_familiarity_treat_med, p_value_familiarity_treat_med)
  ci_familiarity_treat_med <- round_normal(confidence_intervals_a["treatment:familiarity_PB_centered",])
  ci_familiarity_treat_med <- paste0("[",paste0(ci_familiarity_treat_med, collapse = ", "),"]")
  
  # Extract additional parameters from outcome model (path_b_c_model)
  ## Get coefficients based on model type
  if(inherits(path_b_c_model, "lm")) {
    ## Pressure to Publish
    beta_pressure_out <- round_normal(fixed_effects_b_c["pressure_to_publish_numeric_centered", 1])
    p_value_pressure_out <- fixed_effects_b_c["pressure_to_publish_numeric_centered", 4]
    
    ## Familiarity with PB
    beta_familiarity_out <- round_normal(fixed_effects_b_c["familiarity_PB_centered", 1])
    p_value_familiarity_out <- fixed_effects_b_c["familiarity_PB_centered", 4]
    
    ## Pressure to Publish * treatment
    beta_pressure_treat_out <- round_normal(fixed_effects_b_c["treatment:pressure_to_publish_numeric_centered", 1])
    p_value_pressure_treat_out <- fixed_effects_b_c["treatment:pressure_to_publish_numeric_centered", 4]
    
    ## Familiarity with PB * treatment
    beta_familiarity_treat_out <- round_normal(fixed_effects_b_c["treatment:familiarity_PB_centered", 1])
    p_value_familiarity_treat_out <- fixed_effects_b_c["treatment:familiarity_PB_centered", 4]
    
    ## Postdoc Level
    beta_postdoc_out <- round_normal(fixed_effects_b_c["postdoc", 1])
    p_value_postdoc_out <- fixed_effects_b_c["postdoc", 4]
    
    ## Professor Level
    beta_professor_out <- round_normal(fixed_effects_b_c["professor", 1])
    p_value_professor_out <- fixed_effects_b_c["professor", 4]
    
    ## Abstract Position
    beta_abstract_out <- round_normal(fixed_effects_b_c["trial_centered", 1])
    p_value_abstract_out <- fixed_effects_b_c["trial_centered", 4]
    
    ## Conscientiousness
    beta_conscientiousness_out <- round_normal(fixed_effects_b_c["bfi_2_S_score_centered", 1])
    p_value_conscientiousness_out <- fixed_effects_b_c["bfi_2_S_score_centered", 4]
  } else {
    ## Pressure to Publish
    beta_pressure_out <- round_normal(fixed_effects_b_c["pressure_to_publish_numeric_centered", "Estimate"])
    p_value_pressure_out <- fixed_effects_b_c["pressure_to_publish_numeric_centered", "Pr(>|t|)"]
    
    ## Familiarity with PB
    beta_familiarity_out <- round_normal(fixed_effects_b_c["familiarity_PB_centered", "Estimate"])
    p_value_familiarity_out <- fixed_effects_b_c["familiarity_PB_centered", "Pr(>|t|)"]
    
    ## Pressure to Publish * treatment
    beta_pressure_treat_out <- round_normal(fixed_effects_b_c["treatment:pressure_to_publish_numeric_centered", "Estimate"])
    p_value_pressure_treat_out <- fixed_effects_b_c["treatment:pressure_to_publish_numeric_centered", "Pr(>|t|)"]
    
    ## Familiarity with PB * treatment
    beta_familiarity_treat_out <- round_normal(fixed_effects_b_c["treatment:familiarity_PB_centered", "Estimate"])
    p_value_familiarity_treat_out <- fixed_effects_b_c["treatment:familiarity_PB_centered", "Pr(>|t|)"]
    
    ## Postdoc Level
    beta_postdoc_out <- round_normal(fixed_effects_b_c["postdoc", "Estimate"])
    p_value_postdoc_out <- fixed_effects_b_c["postdoc", "Pr(>|t|)"]
    
    ## Professor Level
    beta_professor_out <- round_normal(fixed_effects_b_c["professor", "Estimate"])
    p_value_professor_out <- fixed_effects_b_c["professor", "Pr(>|t|)"]
    
    ## Abstract Position
    beta_abstract_out <- round_normal(fixed_effects_b_c["trial_centered", "Estimate"])
    p_value_abstract_out <- fixed_effects_b_c["trial_centered", "Pr(>|t|)"]
    
    ## Conscientiousness
    beta_conscientiousness_out <- round_normal(fixed_effects_b_c["bfi_2_S_score_centered", "Estimate"])
    p_value_conscientiousness_out <- fixed_effects_b_c["bfi_2_S_score_centered", "Pr(>|t|)"]
  }
  
  # Process outcome model parameters
  beta_p_pressure_out <- process_beta(beta_pressure_out, p_value_pressure_out)
  ci_pressure_out <- round_normal(confidence_intervals_b_c["pressure_to_publish_numeric_centered",])
  ci_pressure_out <- paste0("[",paste0(ci_pressure_out, collapse = ", "),"]")
  
  beta_p_familiarity_out <- process_beta(beta_familiarity_out, p_value_familiarity_out)
  ci_familiarity_out <- round_normal(confidence_intervals_b_c["familiarity_PB_centered",])
  ci_familiarity_out <- paste0("[",paste0(ci_familiarity_out, collapse = ", "),"]")
  
  beta_p_pressure_treat_out <- process_beta(beta_pressure_treat_out, p_value_pressure_treat_out)
  ci_pressure_treat_out <- round_normal(confidence_intervals_b_c["treatment:pressure_to_publish_numeric_centered",])
  ci_pressure_treat_out <- paste0("[",paste0(ci_pressure_treat_out, collapse = ", "),"]")
  
  beta_p_familiarity_treat_out <- process_beta(beta_familiarity_treat_out, p_value_familiarity_treat_out)
  ci_familiarity_treat_out <- round_normal(confidence_intervals_b_c["treatment:familiarity_PB_centered",])
  ci_familiarity_treat_out <- paste0("[",paste0(ci_familiarity_treat_out, collapse = ", "),"]")
  
  beta_p_postdoc_out <- process_beta(beta_postdoc_out, p_value_postdoc_out)
  ci_postdoc_out <- round_normal(confidence_intervals_b_c["postdoc",])
  ci_postdoc_out <- paste0("[",paste0(ci_postdoc_out, collapse = ", "),"]")
  
  beta_p_professor_out <- process_beta(beta_professor_out, p_value_professor_out)
  ci_professor_out <- round_normal(confidence_intervals_b_c["professor",])
  ci_professor_out <- paste0("[",paste0(ci_professor_out, collapse = ", "),"]")
  
  beta_p_abstract_out <- process_beta(beta_abstract_out, p_value_abstract_out)
  ci_abstract_out <- round_normal(confidence_intervals_b_c["trial_centered",])
  ci_abstract_out <- paste0("[",paste0(ci_abstract_out, collapse = ", "),"]")
  
  beta_p_conscientiousness_out <- process_beta(beta_conscientiousness_out, p_value_conscientiousness_out)
  ci_conscientiousness_out <- round_normal(confidence_intervals_b_c["bfi_2_S_score_centered",])
  ci_conscientiousness_out <- paste0("[",paste0(ci_conscientiousness_out, collapse = ", "),"]")
  
  # Create fixed effects table with all parameters
  results_table <- data.frame(
    beta_p = c(beta_p_treatment_a, NA, ade_beta_p, acme_beta_p, beta_p_for_resp_b, total_effect_beta_p, prop_mediated_beta_p, 
               NA, beta_p_pressure_med, beta_p_familiarity_med, beta_p_pressure_treat_med, beta_p_familiarity_treat_med,
               NA, beta_p_pressure_out, beta_p_familiarity_out, beta_p_pressure_treat_out, beta_p_familiarity_treat_out,
               beta_p_postdoc_out, beta_p_professor_out, beta_p_abstract_out, beta_p_conscientiousness_out),
    ci = c(ci_treatment_a, NA, ade_ci, acme_ci, ci_for_resp_b, total_effect_ci, prop_mediated_ci,
           NA, ci_pressure_med, ci_familiarity_med, ci_pressure_treat_med, ci_familiarity_treat_med,
           NA, ci_pressure_out, ci_familiarity_out, ci_pressure_treat_out, ci_familiarity_treat_out,
           ci_postdoc_out, ci_professor_out, ci_abstract_out, ci_conscientiousness_out),
    row.names = c("Path a", "empty1","ADE", "ACME", "Path b", "Total Effect", "Proportion Mediated", 
                 "empty2", "Pressure to Publish_med", "Familiarity w. PB_med", "Pressure to Publish * treatment_med", "Familiarity w. PB * treatment_med",
                 "empty3", "Pressure to Publish_out", "Familiarity w. PB_out", "Pressure to Publish * treatment_out", "Familiarity w. PB * treatment_out",
                 "Postdoc Level_out", "Professorial Level_out", "Abstract Position_out", "Conscientiousness_out")
  )
  
  return(results_table)
}





########
## DF ##
########

df_list <- list()

for (i in 1:4) {
  path_ = paste0(path, i, "/")
  
  
  # List all files in the current folder
  files <-
    file.info(list.files(path_, full.names = T))
  
  # re-order data by date
  files <- files[order(files$mtime, decreasing = TRUE), ]
  
  
  # read all files and append to one data.frame "data"
  df <- data.frame()
  
  for (identifier in 1:nrow(files)) {
    df <- rbind(df, read.csv(rownames(files)[identifier]))
  }
  
  df$experiment <- paste0("experiment_", i)
  
  df_list[[i]] <- df
  
}

df1 <- df_list[[1]]
df2 <- df_list[[2]]
df3 <- df_list[[3]]
df4 <- df_list[[4]]


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

# create a new variable for the decision change
df1$decision_change <-df1$decision2_resp - df1$decision1_resp
df2$decision_change <-df2$decision2_resp - df2$decision1_resp
df3$decision_change_read <-df3$decision2_read_resp - df3$decision1_read_resp
df3$decision_change_cite <-df3$decision2_cite_resp - df3$decision1_cite_resp
df4$decision_change_read <-df4$decision2_read_resp - df4$decision1_read_resp
df4$decision_change_cite <-df4$decision2_cite_resp - df4$decision1_cite_resp

######################
# Exploratory models #
######################

## postdoc and professor variable
# write function to apply this transformation to dataframes
transform_position <- function(df) {
  df$postdoc<-ifelse(df$position == "Post-Doctoral Researcher" | df$position == "Other: PhD/Doctorate Degree", 1, 0)
  df$professor<-ifelse(df$position == "Junior Professor" | df$position == "Professor", 1, 0)
  return(df)
}

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




##################
# Stage 2 Models #
##################


##########
## Main ##
##########


# Experiment 1

# Model 1.2
# Path a: treatment → FOR
model_1.2_med <- lmer(for_resp ~ treatment + (treatment | id), data = df1)
report(model_1.2_med)

# Path c and b: treatment + FOR → C-LoSRC
model_1.2_out <- lmer(change_decision ~ for_resp + treatment + (treatment | id), data = df1)
# Model with singular fit, so drop random slope
model_1.2_out <- lmer(change_decision ~ for_resp + treatment + (1 | id), data = df1)
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
model_2.2_out <- lmer(change_decision ~ for_resp + treatment + (treatment | id), data = df2)
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
print(sub_models)


#################
## Exploratory ##
#################

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
### Usually fit lm, but for consistency with other models, fit lmer
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
                           trial_centered + (1 | id),
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





# Load the lmerTest package
library(lmerTest)

########################################
# Refit models for p-value computation #
########################################



## Model 1.2
model_1.2_med_refit<-lmer(for_resp ~ treatment + (treatment | id), data = df1)
model_1.2_out_refit<-lmer(change_decision ~ for_resp + treatment + (1 | id), data = df1)

## Model 2.2
model_2.2_med_refit<-lmer(for_resp ~ treatment + (1 | id), data = df2)
model_2.2_out_refit<-lmer(change_decision ~ for_resp + treatment + (treatment | id), data = df2)

## Model 3A.2
model_3A.2_med_refit<-lmer(for_read_resp ~ treatment + (treatment | id), data = df3)
model_3A.2_out_refit<-lm(decision_change_read ~ for_read_resp + treatment, data = df3)

## Model 3B.2
model_3B.2_med_refit<-lmer(for_cite_resp ~ treatment + (1 | id), data = df3)
model_3B.2_out_refit<-lmer(decision_change_cite ~ for_cite_resp + treatment + (treatment | id), data = df3)

## Model 4A.2
model_4A.2_med_refit<-lmer(for_read_resp ~ treatment + (1 | id), data = df4)
model_4A.2_out_refit<-lmer(decision_change_read ~ for_read_resp + treatment + (1 | id), data = df4)

## Model 4B.2
model_4B.2_med_refit<-lmer(for_cite_resp ~ treatment + (treatment | id), data = df4)
model_4B.2_out_refit<-lmer(decision_change_cite ~ for_cite_resp + treatment + (1 | id), data = df4)




#################
## Exploratory ##
#################

# Model 1.2
model_1.2_med_refit<-lmer(for_resp ~ treatment * familiarity_PB_centered +
                            treatment * pressure_to_publish_numeric_centered + (treatment | id),
                          data = df_clean1)

model_1.2_out_refit<-lmer(change_decision ~ treatment * familiarity_PB_centered +
                            treatment * pressure_to_publish_numeric_centered + for_resp + bfi_2_S_score_centered +
                            postdoc + professor +nd + uk + us + can + aus + other +
                            trial_centered + (1 | id),
                          data = df_clean1)

## Model 2.2
model_2.2_med_refit<-lmer(for_resp ~ treatment * familiarity_PB_centered +
                            treatment * pressure_to_publish_numeric_centered + (1 | id),
                          data = df_clean2)

model_2.2_out_refit<-lmer(change_decision ~ treatment * familiarity_PB_centered +
                            treatment * pressure_to_publish_numeric_centered + for_resp + bfi_2_S_score_centered +
                            postdoc + professor +nd + uk + us + can + aus + other +
                            trial_centered + (treatment | id),
                          data = df_clean2)

# Experiment 3 - Reading
## med
model_3A.2_med_refit<-lmer(for_read_resp ~ treatment * familiarity_PB_centered +
                             treatment * pressure_to_publish_numeric_centered + (treatment | id),
                           data = df_clean3_read)

## out
model_3A.2_out_refit<-lmer(decision_change_read ~ treatment * familiarity_PB_centered +
                             treatment * pressure_to_publish_numeric_centered + for_read_resp + bfi_2_S_score_centered +
                             postdoc + professor +nd + uk + us + can + aus + other +
                             trial_centered + (1 | id),
                           data = df_clean3_read)

# Experiment 3 - Citing
## med
model_3B.2_med_refit<-lmer(for_cite_resp ~ treatment * familiarity_PB_centered +
                             treatment * pressure_to_publish_numeric_centered + (treatment | id),
                           data = df_clean3_cite)

## out
model_3B.2_out_refit<-lmer(decision_change_cite ~ treatment * familiarity_PB_centered +
                             treatment * pressure_to_publish_numeric_centered + for_cite_resp + bfi_2_S_score_centered +
                             postdoc + professor +nd + uk + us + can + aus + other +
                             trial_centered + (1 | id),
                           data = df_clean3_read)



# Experiment 4 - Reading
## med
model_4A.2_med_refit<-lmer(for_read_resp ~ treatment * familiarity_PB_centered +
                             treatment * pressure_to_publish_numeric_centered + (1 | id),
                           data = df_clean4_read)

## out
model_4A.2_out_refit<-lmer(decision_change_read ~ treatment * familiarity_PB_centered +
                             treatment * pressure_to_publish_numeric_centered + for_read_resp + bfi_2_S_score_centered +
                             postdoc + professor +nd + uk + us + can + aus + other +
                             trial_centered + (1 | id),
                           data = df_clean4_read)

# Experiment 4 - Citing
## med
model_4B.2_med_refit<-lmer(for_cite_resp ~ treatment * familiarity_PB_centered +
                             treatment * pressure_to_publish_numeric_centered + (treatment | id),
                           data = df_clean4_cite)

## out
model_4B.2_out_refit<-lmer(decision_change_cite ~ treatment * familiarity_PB_centered +
                             treatment * pressure_to_publish_numeric_centered + for_cite_resp + bfi_2_S_score_centered +
                             postdoc + professor +nd + uk + us + can + aus + other +
                             trial_centered + (treatment | id),
                           data = df_clean4_cite)


###############################
## S T O R E   R E S U L T S ##
###############################

# Define main and exploratory rownames
main_rownames = c("Path a", "empty1", "ADE", "ACME", "Path b", "Total Effect", "Proportion Mediated")
exploratory_rownames = c("Path a", "empty1", "ADE", "ACME", "Path b", "Total Effect", "Proportion Mediated", 
                         "empty2", "Pressure to Publish_med", "Familiarity w. PB_med", 
                         "Pressure to Publish * treatment_med", "Familiarity w. PB * treatment_med",
                         "empty3", "Pressure to Publish_out", "Familiarity w. PB_out", 
                         "Pressure to Publish * treatment_out", "Familiarity w. PB * treatment_out",
                         "Postdoc Level_out", "Professorial Level_out", "Abstract Position_out", "Conscientiousness_out")

# Get model parameters for table
results_table1<-obtain_model_parameters_main("model_1.2", model_1.2, model_1.2_med_refit, model_1.2_out_refit, for_type = "for_resp")
results_table2<-obtain_model_parameters_main("model_2.2", model_2.2, model_2.2_med_refit, model_2.2_out_refit, for_type = "for_resp")
results_table3A<-obtain_model_parameters_main("model_3A.2", model_3A.2, model_3A.2_med_refit, model_3A.2_out_refit, for_type = "for_read_resp")
results_table3B<-obtain_model_parameters_main("model_3B.2", model_3B.2, model_3B.2_med_refit, model_3B.2_out_refit, for_type = "for_cite_resp")
results_table4A<-obtain_model_parameters_main("model_4A.2", model_4A.2, model_4A.2_med_refit, model_4A.2_out_refit, for_type = "for_read_resp")
results_table4B<-obtain_model_parameters_main("model_4B.2", model_4B.2, model_4B.2_med_refit, model_4B.2_out_refit, for_type = "for_cite_resp")

# add empty rows to main results
add_empty_rows<-function(results_table){
  rownames_exp<-exploratory_rownames[8:21]
  empty_vector<-rep(NA, 14)
  exploratory_empty<-data.frame(beta_p=empty_vector, ci=empty_vector)
  rownames(exploratory_empty)<-rownames_exp
  results_table<-rbind(results_table, exploratory_empty)
  return(results_table)
}
results_table1<-add_empty_rows(results_table1)
results_table2<-add_empty_rows(results_table2)
results_table3A<-add_empty_rows(results_table3A)
results_table3B<-add_empty_rows(results_table3B)
results_table4A<-add_empty_rows(results_table4A)
results_table4B<-add_empty_rows(results_table4B)

# Get exploratory model parameters for table
results_table1_exp<-obtain_model_parameters_exploratory("model_1.2.2", model_1.2.2, model_1.2_med_refit, model_1.2_out_refit, for_type = "for_resp")
results_table2_exp<-obtain_model_parameters_exploratory("model_2.2.2", model_2.2.2, model_2.2_med_refit, model_2.2_out_refit, for_type = "for_resp")
results_table3A_exp<-obtain_model_parameters_exploratory("model_3A.2.2", model_3A.2.2, model_3A.2_med_refit, model_3A.2_out_refit, for_type = "for_read_resp")
results_table3B_exp<-obtain_model_parameters_exploratory("model_3B.2.2", model_3B.2.2, model_3B.2_med_refit, model_3B.2_out_refit, for_type = "for_cite_resp")
results_table4A_exp<-obtain_model_parameters_exploratory("model_4A.2.2", model_4A.2.2, model_4A.2_med_refit, model_4A.2_out_refit, for_type = "for_read_resp")
results_table4B_exp<-obtain_model_parameters_exploratory("model_4B.2.2", model_4B.2.2, model_4B.2_med_refit, model_4B.2_out_refit, for_type = "for_cite_resp")

# cbind main and exploratory results
results_table1<-cbind(results_table1, results_table1_exp)
results_table2<-cbind(results_table2, results_table2_exp)
results_table3A<-cbind(results_table3A, results_table3A_exp)
results_table3B<-cbind(results_table3B, results_table3B_exp)
results_table4A<-cbind(results_table4A, results_table4A_exp)
results_table4B<-cbind(results_table4B, results_table4B_exp)

# add rownames
results_table1$rownames<-main_rownames
results_table1<-results_table1[,c(5,1:4)]
results_table2$rownames<-main_rownames
results_table2<-results_table2[,c(5,1:4)]
results_table3A$rownames<-main_rownames
results_table3A<-results_table3A[,c(5,1:4)]
results_table3B$rownames<-main_rownames
results_table3B<-results_table3B[,c(5,1:4)]
results_table4A$rownames<-main_rownames
results_table4A<-results_table4A[,c(5,1:4)]
results_table4B$rownames<-main_rownames
results_table4B<-results_table4B[,c(5,1:4)]

# write to excel
#write_xlsx(results_table1, file.path(output_path, "Appendix/Stage2/appendix_model1_summary_stage_2_model_overview.xlsx"))
#write_xlsx(results_table2, file.path(output_path, "Appendix/Stage2/appendix_model2_summary_stage_2_model_overview.xlsx"))
#write_xlsx(results_table3A, file.path(output_path, "Appendix/Stage2/appendix_model3A_summary_stage_2_model_overview.xlsx"))
#write_xlsx(results_table3B, file.path(output_path, "Appendix/Stage2/appendix_model3B_summary_stage_2_model_overview.xlsx"))
#write_xlsx(results_table4A, file.path(output_path, "Appendix/Stage2/appendix_model4A_summary_stage_2_model_overview.xlsx"))
#write_xlsx(results_table4B, file.path(output_path, "Appendix/Stage2/appendix_model4B_summary_stage_2_model_overview.xlsx"))