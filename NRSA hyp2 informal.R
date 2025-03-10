library(lme4)    # For multilevel modeling
library(tidyr)   # For data manipulation
library(dplyr)   # For data manipulation
library(psych)   # For descriptive statistics
library(ggplot2) # For plotting
library(lmerTest) # This package adds p-values to lmer output
library(mediation) # For mediation analysis

# Set working directory to where the CSV files are located
setwd("/Users/lacosta3/Desktop/NRSA Data Cleaning Done")

# Read the data from two separate files
pte_hsb_data <- read.csv("PTE and HSB_level 1.csv")
ldsc_data <- read.csv("LDSC_lv2.csv")

# Replace ALL -99 values with NA throughout both datasets
pte_hsb_data[pte_hsb_data == -99] <- NA
ldsc_data[ldsc_data == -99] <- NA

# Verify the conversion
cat("Checking for any remaining -99 values in the datasets:\n")
any_missing_pte <- any(pte_hsb_data == -99, na.rm = TRUE)
any_missing_ldsc <- any(ldsc_data == -99, na.rm = TRUE)
cat("Any -99 values remaining in PTE/HSB data:", any_missing_pte, "\n")
cat("Any -99 values remaining in LDSC data:", any_missing_ldsc, "\n")

# Merge level-1 and level-2 data by ID
data <- merge(pte_hsb_data, ldsc_data, by = "ID", all.x = TRUE)

# Create person-mean centered variables for PTEI
# Get person means for treatment effectiveness
person_means <- aggregate(PTEI ~ ID, data = data, mean, na.rm = TRUE)
colnames(person_means) <- c("ID", "PTEI_mean")

# Merge person means back to the dataset
data <- merge(data, person_means, by = "ID", all.x = TRUE)

# Create within-person centered variables for PTEI
data$PTEI_within <- data$PTEI - data$PTEI_mean

# Create grand-mean centered level-2 predictor (LDSC)
data$LDSC_centered <- scale(data$LDSC, center = TRUE, scale = FALSE)

# Analyze descriptive statistics
# Descriptives for key variables
cat("\nDescriptive Statistics:\n")
for(var in c("LDSC", "PTEI", "HSBI")) {
  cat(var, "- Mean:", mean(data[[var]], na.rm = TRUE), 
      "SD:", sd(data[[var]], na.rm = TRUE), "\n")
}

# Count observations per participant
obs_per_id <- table(data$ID)
cat("\nNumber of unique participants:", length(obs_per_id), "\n")
cat("Average observations per participant:", mean(obs_per_id), "\n")
cat("Range of observations per participant:", range(obs_per_id), "\n")

# Hypothesis 2 - Test path 5: Lifetime Discrimination → Treatment Effectiveness
# Model with random intercept only
model_h2_path5a <- lmer(PTEI ~ LDSC_centered + 
                          (1 | ID), 
                        data = data, 
                        REML = FALSE)  # Using ML for FIML estimation

cat("\nModel for Path 5 (Discrimination → Treatment Effectiveness) - Random Intercept:\n")
print(summary(model_h2_path5a))

# Sort data by ID and date
data <- data[order(data$ID, data$SDT),]  # Sort by ID and date
data$time <- unlist(tapply(data$ID, data$ID, function(x) 1:length(x)))

# Model with random slope for time
# Changed PTEF to PTEI since we're using PTEI in this analysis
model_h2_path5b <- lmer(PTEI ~ LDSC_centered + time + 
                          (time | ID), 
                        data = data, 
                        REML = FALSE)  # Using ML for FIML estimation

# Print results
cat("\nModel for Path 5 (Discrimination → Treatment Effectiveness) - Random Intercept and Slope:\n")
print(summary(model_h2_path5b))

# Compare models to see if random slope improves fit
anova_result <- anova(model_h2_path5a, model_h2_path5b)
cat("\nModel Comparison for Random Intercept vs. Random Intercept and Slope:\n")
print(anova_result)

# Hypothesis 2 - Test path 6: Treatment Effectiveness → Treatment Use
# We'll use a between-person approach since treatment utilization is at level 2
# Create person-level dataset with aggregated variables
person_level <- aggregate(cbind(LDSC, PTEI, HSBI) ~ ID, 
                          data = data, 
                          mean,
                          na.rm = TRUE)

# After creating person_level dataset, convert HSBI back to binary
person_level$HSBI_binary <- ifelse(person_level$HSBI >= 0.5, 1, 0)

# Then use the binary version in your model
model_h2_path6 <- glm(HSBI_binary ~ PTEI, 
                      data = person_level,
                      family = binomial(link = "logit"))

# Print results
cat("\nModel for Path 6 (Treatment Effectiveness → Treatment Utilization):\n")
print(summary(model_h2_path6))

# Calculate odds ratios for easier interpretation
cat("\nOdds Ratios for Treatment Utilization Model:\n")
OR <- exp(coef(model_h2_path6))
print(OR)

# Complete mediation analysis for Hypothesis 2 (LDSC → PTEI → HSBI)
med_model <- mediate(
  model.m = lm(PTEI ~ LDSC, data = person_level, na.action = na.exclude),
  model.y = glm(HSBI_binary ~ PTEI + LDSC, 
                data = person_level, 
                family = binomial(link = "logit"),
                na.action = na.exclude),
  treat = "LDSC",
  mediator = "PTEI",
  boot = TRUE,
  sims = 1000)

# Print mediation results
cat("\nMediation Analysis Results (LDSC → PTEI → HSBI):\n")
print(summary(med_model))