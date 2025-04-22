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
##### Script: Descriptive Statistics                #####
#####                                               #####
#########################################################

####################
# IMPORTANT NOTE:  #
####################
# Since the analyses in this script rely on personal data,
# The demographic analyses are NOT executable, all other analyses are executable.

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
nrow_before1<-nrow(df1)
df1<-df1 %>% filter(read1_rt>2000)
nrow_after1<-nrow(df1)
# experiment 2
nrow_before2<-nrow(df2)
df2<-df2 %>% filter(read1_rt>2000)
nrow_after2<-nrow(df2)
# experiment 3
nrow_before3<-nrow(df3)
df3<-df3 %>% filter(read1_rt>2000)
nrow_after3<-nrow(df3)
# experiment 4
nrow_before4<-nrow(df4)
df4<-df4 %>% filter(read1_rt>2000)
nrow_after4<-nrow(df4)
# print exclusion
print(paste0("Experiment 1: ",nrow_before1-nrow_after1," trials were excluded due to reading times below 2 seconds."))
print(paste0("Experiment 2: ",nrow_before2-nrow_after2," trials were excluded due to reading times below 2 seconds."))
print(paste0("Experiment 3: ",nrow_before3-nrow_after3," trials were excluded due to reading times below 2 seconds."))
print(paste0("Experiment 4: ",nrow_before4-nrow_after4," trials were excluded due to reading times below 2 seconds."))

# sum up all exclusions
sum_exclusions<-sum(nrow_before1-nrow_after1,nrow_before2-nrow_after2,nrow_before3-nrow_after3,nrow_before4-nrow_after4)
print(paste0("In total, ",sum_exclusions," trials were excluded due to reading times below 2 seconds."))


# ==============================================================================
#                           SAMPLE CHARACTERISTICS
# ==============================================================================

# Calculate sample sizes for each experiment
n_participants1 <- length(unique(df1$id))
n_participants2 <- length(unique(df2$id))
n_participants3 <- length(unique(df3$id))
n_participants4 <- length(unique(df4$id))

# Calculate number of trials
n_trials1 <- nrow(df1)
n_trials2 <- nrow(df2)
n_trials3 <- nrow(df3)
n_trials4 <- nrow(df4)

# Print sample characteristics
cat(
  "For experiment 1, we have ",
  n_participants1,
  " participants and ",
  n_trials1,
  " trials.\nFor experiment 2, we have ",
  n_participants2,
  " participants and ",
  n_trials2,
  " trials.\nFor experiment 3, we have ",
  n_participants3,
  " participants and ",
  n_trials3,
  " trials.\nFor experiment 4, we have ",
  n_participants4,
  " participants and ",
  n_trials4,
  " trials."
)



# ==============================================================================
#                           DATA MERGING
# ==============================================================================

# Find common variables across all experiments
c1 <- colnames(df1)
c2 <- colnames(df2)
c3 <- colnames(df3)
c4 <- colnames(df4)
cnames <- unique(c(c1, c2, c3, c4))

# Keep only variables present in all experiments
cnames <- cnames[cnames %in% c1]
cnames <- cnames[cnames %in% c2]
cnames <- cnames[cnames %in% c3]
cnames <- cnames[cnames %in% c4]

# Combine all experiments into one dataset
df <- rbind(df1[,cnames],
           df2[,cnames],
           df3[,cnames],
           df4[,cnames])

# Create unique participant IDs across experiments
df$id_experiments <- paste0(df$experiment,"_",df$id)
table(df$id_experiments)
df_all_trials <- df
df <- df %>% distinct(id_experiments, .keep_all = TRUE)

# ==============================================================================
#                           HELPER FUNCTIONS
# ==============================================================================

# Function to calculate and print absolute and relative frequencies
abs_and_rel <- function(x) {
    abs <- table(x)
    rel <- round(table(x)/length(x)*100, 0)
    # Sort frequencies in descending order
    abs <- abs[order(abs, decreasing = TRUE)]
    rel <- rel[order(rel, decreasing = TRUE)]
    for (i in 1:length(abs)) {
        cat(paste0(names(abs)[i],": ",abs[i]," (",rel[i],"%)\n"))
    }
}

# Function to create formatted frequency tables
abs_and_rel_table <- function(x) {
    abs <- table(x)
    rel <- round(table(x) / length(x) * 100, 0)
    
    # Sort frequencies by category names
    abs <- abs[order(names(abs), decreasing = FALSE)]
    rel <- rel[order(names(rel), decreasing = FALSE)]
    
    # Ensure all categories are included
    all_levels <- unique(x)
    all_levels <- sort(all_levels)
    
    # Create formatted strings (e.g., "47 (62.7)")
    store <- vector("character", length(all_levels))
    names_store <- all_levels
    
    for (i in seq_along(all_levels)) {
        level <- all_levels[i]
        count <- ifelse(level %in% names(abs), abs[level], 0)
        percent <- ifelse(level %in% names(rel), rel[level], 0)
        store[i] <- paste0(count, " (", percent, ")")
    }
    
    # Convert to dataframe
    local_df <- as.data.frame(t(store), stringsAsFactors = FALSE)
    colnames(local_df) <- names_store
    
    return(local_df)
}





# ==============================================================================
#                           VISUALIZATION
# ==============================================================================

# Create boxplots for each experiment
# Experiment 1: Effect of significance on submission likelihood
plot_exp1 <- ggplot(aes(y = decision1_resp, x = factor(treatment_string)), data = df1, fill = factor(treatment_string)) +
  geom_boxplot() + theme_apa() + 
  xlab("Significance of results") + 
  ylab("Likelihood of submitting an abstract for publication") + 
  ggtitle("Experiment 1")

# Experiment 2: Effect of hypothesis-consistency on submission likelihood
plot_exp2 <- ggplot(aes(y = decision1_resp, x = factor(treatment_string, levels = rev(levels(factor(treatment_string))))), data = df2) +
  geom_boxplot() + theme_apa() +
  xlab("Hypothesis-consistency of results") + 
  ylab("Likelihood of submitting an abstract for publication") + 
  ggtitle("Experiment 2")

# Experiment 3a: Effect of significance on citation likelihood
plot_exp3a <- ggplot(aes(y = decision1_cite_resp, x = factor(treatment_string)), data = df3) +
  geom_boxplot() + theme_apa() +
  xlab("Significance of results") + 
  ylab("Likelihood of citing an abstract") + 
  ggtitle("Experiment 3")

# Experiment 3b: Effect of significance on reading likelihood
plot_exp3b <- ggplot(aes(y = decision1_read_resp, x = factor(treatment_string)), data = df3) +
  geom_boxplot() + theme_apa() +
  xlab("Significance of results") + 
  ylab("Likelihood of reading an abstract") + 
  ggtitle("Experiment 3")

# Experiment 4a: Effect of hypothesis-consistency on citation likelihood
plot_exp4a <- ggplot(aes(y = decision1_cite_resp, x = factor(treatment_string, levels = rev(levels(factor(treatment_string))))), data = df4) +
  geom_boxplot() + theme_apa() +
  xlab("Hypothesis-consistency of results") + 
  ylab("Likelihood of citing an Abstract") + 
  ggtitle("Experiment 4")

# Experiment 4b: Effect of hypothesis-consistency on reading likelihood
plot_exp4b <- ggplot(aes(y = decision1_read_resp, x = factor(treatment_string, levels = rev(levels(factor(treatment_string))))), data = df4) +
  geom_boxplot() + theme_apa() +
  xlab("Hypothesis-consistency of results") + 
  ylab("Likelihood of reading an abstract") + 
  ggtitle("Experiment 4") 

# Arrange plots in a 3x2 grid
figure3 <- grid.arrange(plot_exp1, plot_exp2, plot_exp3a, plot_exp4a, plot_exp3b, plot_exp4b, ncol = 2)

# Save plots in different formats
#ggsave(paste0(plot_path,"figure3.png"), plot = figure3, width = 25, height = 30, units = "cm")
#ggsave(paste0(plot_path,"figure3.tiff"), plot = figure3, width = 25, height = 30, units = "cm", dpi = 300)



# ==============================================================================
#                           STATISTICAL TESTS
# ==============================================================================

# Conduct t-tests for differences between conditions
# LoS = Likelihood of Submission
# LoC = Likelihood of Citation
# LoR = Likelihood of Reading

# Effect of significance
LoS_Significance <- t.test(data = df1, decision1_resp ~ treatment)
LoC_Significance <- t.test(data = df3, decision1_cite_resp ~ treatment)
LoR_Significance <- t.test(data = df3, decision1_read_resp ~ treatment)

# Effect of hypothesis-consistency
LoS_Hypothesis <- t.test(data = df2, decision1_resp ~ treatment)
LoC_Hypothesis <- t.test(data = df4, decision1_cite_resp ~ treatment)
LoR_Hypothesis <- t.test(data = df4, decision1_read_resp ~ treatment)

# Print statistical test results
print("Effects of significance:")
report(LoS_Significance)
report(LoR_Significance)
report(LoC_Significance)

print("Effects of hypothesis-consistency:")
report(LoS_Hypothesis)
report(LoR_Hypothesis)
report(LoC_Hypothesis)


# ==============================================================================
#                           SUMMARY STATISTICS
# ==============================================================================

# Calculate means and standard deviations by condition
# Format: Mean (SD)

# Experiment 1: Significance effect on submission
res1 <- df1 %>% 
  group_by(treatment_string) %>% 
  summarise(mean_sd_pub = paste0(round(mean(decision1_resp),2)," (", round(sd(decision1_resp),2), ")"))

# Experiment 2: Hypothesis-consistency effect on submission
res2 <- df2 %>% 
  group_by(treatment_string) %>% 
  summarise(mean_sd_pub = paste0(round(mean(decision1_resp),2)," (", round(sd(decision1_resp),2), ")"))

# Experiment 3: Significance effect on citation and reading
res3a <- df3 %>% 
  group_by(treatment_string) %>% 
  summarise(mean_sd_cite = paste0(round(mean(decision1_cite_resp),2)," (", round(sd(decision1_cite_resp),2), ")"))
res3b <- df3 %>% 
  group_by(treatment_string) %>% 
  summarise(mean_sd_read = paste0(round(mean(decision1_read_resp),2)," (", round(sd(decision1_read_resp),2), ")"))

# Experiment 4: Hypothesis-consistency effect on citation and reading
res4a <- df4 %>% 
  group_by(treatment_string) %>% 
  summarise(mean_sd_cite = paste0(round(mean(decision1_cite_resp),2)," (", round(sd(decision1_cite_resp),2), ")"))
res4b <- df4 %>% 
  group_by(treatment_string) %>% 
  summarise(mean_sd_read = paste0(round(mean(decision1_read_resp),2)," (", round(sd(decision1_read_resp),2), ")"))

# Combine results into a single table
merge_results <- merge(rbind(res1,res2), rbind(res3b,res4b), by="treatment_string")
merge_results <- merge(merge_results, rbind(res3a,res4a), by="treatment_string")

# Sort results
merge_results$sort <- c(3,4,2,1)
merge_results <- merge_results[order(merge_results$sort),]
merge_results <- merge_results[,-ncol(merge_results)]

# Transpose results table for better readability
merge_results <- data.frame(t(merge_results))
merge_results$rowname <- rownames(merge_results)
merge_results <- merge_results[,c(ncol(merge_results),1:(ncol(merge_results)-1))]

# Export results to Excel
#write_xlsx(merge_results, paste0(output_path,"summary_table_mean_differences.xlsx"))



## NOT REPRODUCIBLE
# ==============================================================================
#                           DEMOGRAPHIC ANALYSIS
# ==============================================================================

# Analyze gender distribution
print("Gender distribution across all experiments:")
abs_and_rel(df$gender)

# Analyze age distribution
print("Age distribution across all experiments:")
abs_and_rel(df$age)

# Recode and analyze academic positions
df$position_table <- df$position
df$position_table[df$position %in% c("Junior Professor","Professor")] <- "3_professor"
df$position_table[df$position %in% c("Other: PhD/Doctorate Degree","Post-Doctoral Researcher")] <- "2_postdoctoral"
df$position_table[df$position %in% c("PhD Candidate","Other: No PhD/Doctorate Degree")] <- "1_predoctoral"
print("Academic position distribution:")
abs_and_rel(df$position_table)

# gender
df$gender_table<-df$gender
df$gender_table[df$gender %in% "Female"]<-"1_Female"
df$gender_table[df$gender %in% "Diverse"]<-"2_Diverse"
df$gender_table[df$gender %in% "Male"]<-"3_Male"


# Define the expected column order
expected_columns <- c("X1_Female", "X2_Diverse", "X3_Male", 
                      "X18.25", "X26.35", "X36.45", "X46.55", "X56.or.older", 
                      "X1_predoctoral", "X2_postdoctoral", "X3_professor")

# Function to ensure each overview table has all expected columns
ensure_columns <- function(df, expected_cols) {
  # Find missing columns
  missing_cols <- setdiff(expected_cols, colnames(df))
  
  # Add missing columns with "0 (0)"
  for (col in missing_cols) {
    df[[col]] <- "0 (0)"
  }
  
  # Reorder columns
  df <- df[, expected_cols, drop=FALSE]
  
  return(df)
}

# Generate each experiment overview
g1 <- abs_and_rel_table(df[df$experiment == "experiment_1",]$gender_table)
a1 <- abs_and_rel_table(df[df$experiment == "experiment_1",]$age)
p1 <- abs_and_rel_table(df[df$experiment == "experiment_1",]$position_table)
overview1 <- data.frame(g1, a1, p1)
overview1 <- ensure_columns(overview1, expected_columns)

g2 <- abs_and_rel_table(df[df$experiment == "experiment_2",]$gender_table)
a2 <- abs_and_rel_table(df[df$experiment == "experiment_2",]$age)
p2 <- abs_and_rel_table(df[df$experiment == "experiment_2",]$position_table)
overview2 <- data.frame(g2, a2, p2)
overview2 <- ensure_columns(overview2, expected_columns)

g3 <- abs_and_rel_table(df[df$experiment == "experiment_3",]$gender_table)
a3 <- abs_and_rel_table(df[df$experiment == "experiment_3",]$age)
p3 <- abs_and_rel_table(df[df$experiment == "experiment_3",]$position_table)
overview3 <- data.frame(g3, a3, p3)
overview3 <- ensure_columns(overview3, expected_columns)

g4 <- abs_and_rel_table(df[df$experiment == "experiment_4",]$gender_table)
a4 <- abs_and_rel_table(df[df$experiment == "experiment_4",]$age)
p4 <- abs_and_rel_table(df[df$experiment == "experiment_4",]$position_table)
overview4 <- data.frame(g4, a4, p4)
overview4 <- ensure_columns(overview4, expected_columns)

# Combine all experiments into a single final table
final_overview <- rbind(overview1, overview2, overview3, overview4)
rownames(final_overview) <- c("Experiment 1", "Experiment 2", "Experiment 3", "Experiment 4")

# flipping the data frame
final_overview <- t(final_overview)
final_overview<-data.frame(final_overview)
final_overview$rowname<-rownames(final_overview)
final_overview<-final_overview[,c(ncol(final_overview),1:(ncol(final_overview)-1))]

# Print final combined table
print(final_overview)

# write to xlsx
#write_xlsx(data.frame(final_overview), paste0(output_path,"descriptive_statistics_across_experiments.xlsx"))

