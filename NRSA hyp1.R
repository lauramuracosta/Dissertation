library(lme4)    # For multilevel modeling
library(tidyr)   # For data manipulation
library(dplyr)   # For data manipulation
library(psych)   # For descriptive statistics
library(ggplot2) # For plotting
library(lmerTest) # This package adds p-values to lmer output

# Set working directory to where the CSV file is located
setwd("/Users/lacosta3/Desktop/NRSA Data Cleaning Done")

# Read the data
data <- read.csv('df_hypo_1_MLM.csv', na.strings = c("", "NA"))

# Examine the structure of the data
str(data)
head(data)
summary(data[, c("ID", "T1NA", "T1C", "T1DSC", "T2NA", "T2C", "T2SU", "T3NA", "T3C", "T3DSC")])

# Check for missing values in key variables
missing_count <- colSums(is.na(data[, c("T1NA", "T1C", "T1DSC", "T2NA", "T2C", "T2SU", "T3NA", "T3C", "T3DSC")]))
print(missing_count)

# For hypothesis 1 we will examine:
# 1. Path 1: Discrimination at Time 1 -> Negative Affect at Time 2
# 2. Path 2: Negative Affect at Time 2 -> Craving at Time 3
# 3. Path 3: Craving at Time 2 -> Substance Use at Time 2
# 4. Path 4: Discrimination at Time 1 -> Substance Use at Time 2

# Multilevel approach using lme4

# Path 1: Discrimination (T1D) -> Negative Affect (T2NA)
path1_lmer <- lmer(T2NA ~ T1DSC + (1|ID), data = data, REML = FALSE)
summary(path1_lmer)

# let's create a null model to see if discrimination as a predictor fit the data significantly better than the model without it
path1_null <- lmer(T2NA ~ (1 | ID), data = data, REML = FALSE)

# compare the models
anova(path1_null, path1_lmer)

# Path 2: Negative Affect (T2NA) -> Craving (T3C)
path2_lmer <- lmer(T3C ~ T2NA + (1|ID), data = data, REML = FALSE)
summary(path2_lmer)

# Path 3: Craving (T2C) -> Substance Use (T2SU) # For substance use count outcome
path3_count <- glmer(T2SU ~ T2C + (1|ID), data = data, family = poisson)
summary(path3_count)

# Path 4: Discrimination (T1D) -> Substance Use (T2SU) # For substance use count outcome
path4_count <- glmer(T2SU ~ T1DSC + (1|ID), data = data, family = poisson)
summary(path4_count)

# Calculate person means
data <- data %>%
  group_by(ID) %>%
  mutate(
    person_mean_DSC = mean(c(T1DSC, T3DSC), na.rm = TRUE),
    T1DSC_centered = T1DSC - person_mean_DSC,
    T3DSC_centered = T3DSC - person_mean_DSC
  ) %>%
  ungroup()

# Use these in the models
path1MC_lmer <- lmer(T2NA ~ T1DSC_centered + person_mean_DSC + (1|ID), data = data, REML = FALSE)
summary(path1MC_lmer)

# Function to calculate design effect for linear mixed models
calc_deff_lmer <- function(model) {
  # Extract variance components
  vc <- VarCorr(model)
  var_intercept <- as.numeric(attr(vc[[1]], "stddev")^2)
  var_residual <- attr(vc, "sc")^2
  
  # Calculate ICC
  icc <- var_intercept / (var_intercept + var_residual)
  
  # Get sample size information
  n_obs <- nobs(model)
  n_groups <- length(unique(model@flist[[1]]))
  avg_cluster_size <- n_obs / n_groups
  
  # Calculate design effect
  deff <- 1 + (avg_cluster_size - 1) * icc
  
  return(list(
    icc = icc,
    avg_cluster_size = avg_cluster_size,
    design_effect = deff,
    se_inflation = sqrt(deff),
    effective_sample_size = n_obs / deff
  ))
}

# Function to calculate design effect for GLMMs (Poisson)
calc_deff_glmer <- function(model) {
  # Extract random intercept variance
  vc <- VarCorr(model)
  var_intercept <- as.numeric(attr(vc[[1]], "stddev")^2)
  
  # For Poisson models, use π²/3 as theoretical level-1 variance
  var_residual <- (pi^2) / 3
  
  # Calculate ICC
  icc <- var_intercept / (var_intercept + var_residual)
  
  # Get sample size information
  n_obs <- nobs(model)
  n_groups <- length(unique(model@flist[[1]]))
  avg_cluster_size <- n_obs / n_groups
  
  # Calculate design effect
  deff <- 1 + (avg_cluster_size - 1) * icc
  
  return(list(
    icc = icc,
    avg_cluster_size = avg_cluster_size,
    design_effect = deff,
    se_inflation = sqrt(deff),
    effective_sample_size = n_obs / deff
  ))
}

# Calculate design effects for each model
deff_path1 <- calc_deff_lmer(path1_lmer)
deff_path2 <- calc_deff_lmer(path2_lmer)
deff_path3 <- calc_deff_glmer(path3_count)
deff_path4 <- calc_deff_glmer(path4_count)

# Print results
print("Design Effect for Path 1 (Discrimination → Negative Affect):")
print(deff_path1)

print("Design Effect for Path 2 (Negative Affect → Craving):")
print(deff_path2)

print("Design Effect for Path 3 (Craving → Substance Use):")
print(deff_path3)

print("Design Effect for Path 4 (Discrimination → Substance Use):")
print(deff_path4)


# For linear mixed models (Path 1 and Path 2)
# Path 1
var_int_p1 <- as.numeric(VarCorr(path1_lmer)$ID[1])
var_res_p1 <- sigma(path1_lmer)^2
icc_p1 <- var_int_p1 / (var_int_p1 + var_res_p1)
avg_size_p1 <- nobs(path1_lmer) / length(unique(path1_lmer@flist[[1]]))
deff_p1 <- 1 + (avg_size_p1 - 1) * icc_p1

# Path 2
var_int_p2 <- as.numeric(VarCorr(path2_lmer)$ID[1])
var_res_p2 <- sigma(path2_lmer)^2
icc_p2 <- var_int_p2 / (var_int_p2 + var_res_p2)
avg_size_p2 <- nobs(path2_lmer) / length(unique(path2_lmer@flist[[1]]))
deff_p2 <- 1 + (avg_size_p2 - 1) * icc_p2

# For Poisson models (Path 3 and Path 4)
# Path 3
var_int_p3 <- as.numeric(VarCorr(path3_count)$ID[1])
var_res_p3 <- (pi^2) / 3  # Theoretical level-1 variance for Poisson
icc_p3 <- var_int_p3 / (var_int_p3 + var_res_p3)
avg_size_p3 <- nobs(path3_count) / length(unique(path3_count@flist[[1]]))
deff_p3 <- 1 + (avg_size_p3 - 1) * icc_p3

# Path 4
var_int_p4 <- as.numeric(VarCorr(path4_count)$ID[1])
var_res_p4 <- (pi^2) / 3  # Theoretical level-1 variance for Poisson
icc_p4 <- var_int_p4 / (var_int_p4 + var_res_p4)
avg_size_p4 <- nobs(path4_count) / length(unique(path4_count@flist[[1]]))
deff_p4 <- 1 + (avg_size_p4 - 1) * icc_p4

# Print results
cat("Design Effect for Path 1:", round(deff_p1, 2), "\n")
cat("Design Effect for Path 2:", round(deff_p2, 2), "\n")
cat("Design Effect for Path 3:", round(deff_p3, 2), "\n")
cat("Design Effect for Path 4:", round(deff_p4, 2), "\n")

#Adding in demographics as moderators

# Read the demographic data
demo_data <- read.csv("Level 2 demos", na.strings = c("", "NA"))

# Ensure ID columns are of the same type in both datasets for merging
demo_data$ID <- as.character(demo_data$ID)
data$ID <- as.character(data$ID)

# Merge the datasets
merged_data <- data %>%
  left_join(demo_data %>% select(ID, GEN, ETH, EDU, EMP), by = "ID")

# Visually inspect to see if anything is off
write_csv(merged_data, "testing.csv")

# Check for successful merge and examine distributions
summary(merged_data[c("GEN", "ETH", "EDU", "EMP")])

# Recode demographic variables as factors with meaningful labels
merged_data <- merged_data %>%
  mutate(
    GEN_factor = factor(GEN, levels = c(0, 1), labels = c("Male", "Female")),
    ETH_factor = factor(ETH, levels = c(0, 1, 2), 
                        labels = c("Black only", "Hispanic only", "Both Black and Hispanic")),
    EDU_factor = factor(EDU, levels = c(0, 1, 2, 3, 4), 
                        labels = c("Less than HS", "HS/GED", "Some college", 
                                   "2-year degree", "4-year degree")),
    EMP_factor = factor(EMP, levels = c(0, 1), labels = c("Not employed", "Employed"))
  )

#===============================================================================
# PATH 1: DISCRIMINATION → NEGATIVE AFFECT
#===============================================================================

# Gender moderation
path1_gender_mod <- lmer(T2NA ~ T1DSC * GEN_factor + (1 + T1DSC|ID), 
                         data = merged_data, REML = FALSE)

# Ethnicity moderation 
path1_eth_mod <- lmer(T2NA ~ T1DSC * ETH_factor + (1 + T1DSC|ID), 
                      data = merged_data, REML = FALSE)

# Education moderation
path1_edu_mod <- lmer(T2NA ~ T1DSC * EDU_factor + (1 + T1DSC|ID), 
                      data = merged_data, REML = FALSE)

# Employment moderation
path1_emp_mod <- lmer(T2NA ~ T1DSC * EMP_factor + (1 + T1DSC|ID), 
                      data = merged_data, REML = FALSE)

# Summarize results for Path 1
summary(path1_gender_mod)
summary(path1_eth_mod)
summary(path1_edu_mod)
summary(path1_emp_mod)

#===============================================================================
# PATH 2: NEGATIVE AFFECT → CRAVING
#===============================================================================

# Gender moderation
path2_gender_mod <- lmer(T3C ~ T2NA * GEN_factor + (1 + T2NA|ID), 
                         data = merged_data, REML = FALSE)

# Ethnicity moderation
path2_eth_mod <- lmer(T3C ~ T2NA * ETH_factor + (1 + T2NA|ID), 
                      data = merged_data, REML = FALSE)

# Education moderation
path2_edu_mod <- lmer(T3C ~ T2NA * EDU_factor + (1 + T2NA|ID), 
                      data = merged_data, REML = FALSE)

# Employment moderation
path2_emp_mod <- lmer(T3C ~ T2NA * EMP_factor + (1 + T2NA|ID), 
                      data = merged_data, REML = FALSE)

# Summarize results for Path 2
summary(path2_gender_mod)
summary(path2_eth_mod)
summary(path2_edu_mod)
summary(path2_emp_mod)

# Let's simplify due to the convergence issues
# For gender moderation that had convergence issues
path2_gender_mod_simple <- lmer(T3C ~ T2NA * GEN_factor + (1|ID), 
                                data = merged_data, REML = FALSE)
summary(path2_gender_mod_simple)

# For education moderation that had singular fit
path2_edu_mod_simple <- lmer(T3C ~ T2NA * EDU_factor + (1|ID), 
                             data = merged_data, REML = FALSE)
summary(path2_edu_mod_simple)

# For employment moderation that had convergence issues
path2_emp_mod_simple <- lmer(T3C ~ T2NA * EMP_factor + (1|ID), 
                             data = merged_data, REML = FALSE)
summary(path2_emp_mod_simple)
#===============================================================================
# PATH 3: CRAVING → SUBSTANCE USE (POISSON)
#===============================================================================

# Try with bobyqa optimizer and increased iterations
control_params <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

# Gender moderation
path3_gender_mod <- try(glmer(T2SU ~ T2C * GEN_factor + (1 + T2C|ID), 
                              family = poisson, data = merged_data, 
                              control = control_params))

# Ethnicity moderation
path3_eth_mod <- try(glmer(T2SU ~ T2C * ETH_factor + (1 + T2C|ID), 
                           family = poisson, data = merged_data, 
                           control = control_params))

# Education moderation
path3_edu_mod <- try(glmer(T2SU ~ T2C * EDU_factor + (1 + T2C|ID), 
                           family = poisson, data = merged_data, 
                           control = control_params))

# Employment moderation
path3_emp_mod <- try(glmer(T2SU ~ T2C * EMP_factor + (1 + T2C|ID), 
                           family = poisson, data = merged_data, 
                           control = control_params))

# Convergence issues again, simplify the random effects structure
# Using only random intercepts:
if(inherits(path3_gender_mod, "try-error")) {
  path3_gender_mod <- glmer(T2SU ~ T2C * GEN_factor + (1|ID), 
                            family = poisson, data = merged_data, 
                            control = control_params)
}

# Repeat for other demographic variables if needed
if(inherits(path3_eth_mod, "try-error")) {
  path3_eth_mod <- glmer(T2SU ~ T2C * ETH_factor + (1|ID), 
                         family = poisson, data = merged_data, 
                         control = control_params)
}

if(inherits(path3_edu_mod, "try-error")) {
  path3_edu_mod <- glmer(T2SU ~ T2C * EDU_factor + (1|ID), 
                         family = poisson, data = merged_data, 
                         control = control_params)
}

if(inherits(path3_emp_mod, "try-error")) {
  path3_emp_mod <- glmer(T2SU ~ T2C * EMP_factor + (1|ID), 
                         family = poisson, data = merged_data, 
                         control = control_params)
}

# Summarize results for Path 3
if(!inherits(path3_gender_mod, "try-error")) summary(path3_gender_mod)
if(!inherits(path3_eth_mod, "try-error")) summary(path3_eth_mod)
if(!inherits(path3_edu_mod, "try-error")) summary(path3_edu_mod)
if(!inherits(path3_emp_mod, "try-error")) summary(path3_emp_mod)


#===============================================================================
# PATH 4: DISCRIMINATION → SUBSTANCE USE (POISSON)
#===============================================================================

# Gender moderation
path4_gender_mod <- try(glmer(T2SU ~ T1DSC * GEN_factor + (1 + T1DSC|ID), 
                              family = poisson, data = merged_data, 
                              control = control_params))

# Ethnicity moderation
path4_eth_mod <- try(glmer(T2SU ~ T1DSC * ETH_factor + (1 + T1DSC|ID), 
                           family = poisson, data = merged_data, 
                           control = control_params))

# Education moderation
path4_edu_mod <- try(glmer(T2SU ~ T1DSC * EDU_factor + (1 + T1DSC|ID), 
                           family = poisson, data = merged_data, 
                           control = control_params))

# Employment moderation
path4_emp_mod <- try(glmer(T2SU ~ T1DSC * EMP_factor + (1 + T1DSC|ID), 
                           family = poisson, data = merged_data, 
                           control = control_params))

# Simplify if convergence issues persist
if(inherits(path4_gender_mod, "try-error")) {
  path4_gender_mod <- glmer(T2SU ~ T1DSC * GEN_factor + (1|ID), 
                            family = poisson, data = merged_data, 
                            control = control_params)
}

if(inherits(path4_eth_mod, "try-error")) {
  path4_eth_mod <- glmer(T2SU ~ T1DSC * ETH_factor + (1|ID), 
                         family = poisson, data = merged_data, 
                         control = control_params)
}

if(inherits(path4_edu_mod, "try-error")) {
  path4_edu_mod <- glmer(T2SU ~ T1DSC * EDU_factor + (1|ID), 
                         family = poisson, data = merged_data, 
                         control = control_params)
}

if(inherits(path4_emp_mod, "try-error")) {
  path4_emp_mod <- glmer(T2SU ~ T1DSC * EMP_factor + (1|ID), 
                         family = poisson, data = merged_data, 
                         control = control_params)
}

# Summarize results for Path 4
if(!inherits(path4_gender_mod, "try-error")) summary(path4_gender_mod)
if(!inherits(path4_eth_mod, "try-error")) summary(path4_eth_mod)
if(!inherits(path4_edu_mod, "try-error")) summary(path4_edu_mod)
if(!inherits(path4_emp_mod, "try-error")) summary(path4_emp_mod)

#===============================================================================
# TRYING OUT A MORE SIMPLE MODERATION CODE
# THESE ALL YEILD THE EXACT SAME RESULTS AS ABOVE just in a more organized format that makes patterns easier to see
#===============================================================================

extract_moderation_summary <- function(model, moderator_name) {
  if(inherits(model, "try-error")) {
    return(paste("Model for", moderator_name, "failed to converge"))
  }
  
  sum_model <- summary(model)
  
  if(inherits(model, "lmerMod")) {
    # Extract coefficients table
    coef_table <- sum_model$coefficients
    
    # Find interaction terms
    interaction_rows <- grep(":", rownames(coef_table))
    
    if(length(interaction_rows) > 0) {
      result_df <- data.frame(
        Moderator = moderator_name,
        Interaction = rownames(coef_table)[interaction_rows],
        Estimate = coef_table[interaction_rows, "Estimate"],
        SE = coef_table[interaction_rows, "Std. Error"],
        t_value = coef_table[interaction_rows, "t value"],
        p_value = coef_table[interaction_rows, "Pr(>|t|)"]
      )
      return(result_df)
    } else {
      return(paste("No interaction terms found in", moderator_name, "model"))
    }
  } else if(inherits(model, "glmerMod")) {
    # Extract coefficients for GLMMs
    coef_table <- sum_model$coefficients
    
    # Find interaction terms
    interaction_rows <- grep(":", rownames(coef_table))
    
    if(length(interaction_rows) > 0) {
      result_df <- data.frame(
        Moderator = moderator_name,
        Interaction = rownames(coef_table)[interaction_rows],
        Estimate = coef_table[interaction_rows, "Estimate"],
        SE = coef_table[interaction_rows, "Std. Error"],
        z_value = coef_table[interaction_rows, "z value"],
        p_value = coef_table[interaction_rows, "Pr(>|z|)"]
      )
      return(result_df)
    } else {
      return(paste("No interaction terms found in", moderator_name, "model"))
    }
  } else {
    return(paste("Unknown model type for", moderator_name))
  }
}

# Extract and compile all moderation results
path1_moderation_results <- rbind(
  extract_moderation_summary(path1_gender_mod, "Gender"),
  extract_moderation_summary(path1_eth_mod, "Ethnicity"),
  extract_moderation_summary(path1_edu_mod, "Education"),
  extract_moderation_summary(path1_emp_mod, "Employment")
)

path2_moderation_results <- rbind(
  extract_moderation_summary(path2_gender_mod, "Gender"),
  extract_moderation_summary(path2_eth_mod, "Ethnicity"),
  extract_moderation_summary(path2_edu_mod, "Education"),
  extract_moderation_summary(path2_emp_mod, "Employment")
)

path3_moderation_results <- rbind(
  if(!inherits(path3_gender_mod, "try-error")) extract_moderation_summary(path3_gender_mod, "Gender"),
  if(!inherits(path3_eth_mod, "try-error")) extract_moderation_summary(path3_eth_mod, "Ethnicity"),
  if(!inherits(path3_edu_mod, "try-error")) extract_moderation_summary(path3_edu_mod, "Education"),
  if(!inherits(path3_emp_mod, "try-error")) extract_moderation_summary(path3_emp_mod, "Employment")
)

path4_moderation_results <- rbind(
  if(!inherits(path4_gender_mod, "try-error")) extract_moderation_summary(path4_gender_mod, "Gender"),
  if(!inherits(path4_eth_mod, "try-error")) extract_moderation_summary(path4_eth_mod, "Ethnicity"),
  if(!inherits(path4_edu_mod, "try-error")) extract_moderation_summary(path4_edu_mod, "Education"),
  if(!inherits(path4_emp_mod, "try-error")) extract_moderation_summary(path4_emp_mod, "Employment")
)

# Print summarized results for all paths
cat("PATH 1: DISCRIMINATION → NEGATIVE AFFECT - MODERATION RESULTS\n")
print(path1_moderation_results)

cat("\nPATH 2: NEGATIVE AFFECT → CRAVING - MODERATION RESULTS\n")
print(path2_moderation_results)

cat("\nPATH 3: CRAVING → SUBSTANCE USE - MODERATION RESULTS\n")
print(path3_moderation_results)

cat("\nPATH 4: DISCRIMINATION → SUBSTANCE USE - MODERATION RESULTS\n")
print(path4_moderation_results)
