# Set  working directory to where the CSV file is located
setwd("/Users/lacosta3/Desktop/NRSA Data Cleaning Done")

# Hypothesis 3 Analysis
# This analysis examines whether daily and lifetime discrimination decrease formal help seeking,
# increase perceived barriers to treatment, and affect substance use patterns

# Load necessary libraries
library(lme4)
library(lmerTest)
library(mediation)
library(dplyr)
library(tidyr)
library(ggplot2)
library(psych)

# ======== DATA PREPARATION ========

# Load datasets
df_level1 <- read.csv("DDSC and HSA_level 1.csv")
df_level2_demo <- read.csv("Level 2 demos.csv")
df_weekly <- read.csv("PTE and HSB_level 1.csv")

# Convert dates to proper format
# Convert dates with explicit format specifications

# Try to convert dates with error handling
tryCatch({
  # First attempt with standard format
  df_level1$SDa <- as.Date(df_level1$SDa)
}, error = function(e) {
  # If that fails, try with explicit format 
  df_level1$SDa <- as.Date(df_level1$SDa, format = "%m/%d/%Y")
})

tryCatch({
  # First attempt with standard format
  df_level1$SDT <- as.POSIXct(df_level1$SDT)
}, error = function(e) {
  # If that fails, try with explicit format 
  df_level1$SDT <- as.POSIXct(df_level1$SDT, format = "%m/%d/%Y %H:%M:%S")
})

tryCatch({
  # First attempt with standard format
  df_weekly$SDT <- as.POSIXct(df_weekly$SDT)
}, error = function(e) {
  # If that fails, try with explicit format 
  df_weekly$SDT <- as.POSIXct(df_weekly$SDT, format = "%m/%d/%Y %H:%M:%S")
})

# Merge level 1 and level 2 data
df_merged <- merge(df_level1, df_level2_demo, by = "ID")

# Merge with weekly data for help-seeking outcomes

# Load the additional level 2 file containing LDSC variable
df_ldsc <- read.csv("LDSC_lv2.csv")
print("LDSC data file loaded successfully")


# Now merge in the LDSC data
df_merged <- merge(df_merged, df_ldsc, by = "ID", all.x = TRUE)
print("LDSC variable merged into the dataset")

# Function to match daily with weekly data using sequence numbers instead of dates
match_by_sequence <- function(daily_data, weekly_data) {
  # This will hold our matched results
  matched_data <- data.frame()
  
  # Get unique IDs
  unique_ids <- unique(daily_data$ID)
  
  for (id in unique_ids) {
    # Subset data for this participant
    daily_subset <- daily_data[daily_data$ID == id, ]
    weekly_subset <- weekly_data[weekly_data$ID == id, ]
    
    if (nrow(weekly_subset) == 0) {
      next  # Skip if no weekly data for this ID
    }
    
    # For each daily observation, find the closest following weekly observation
    for (i in 1:nrow(daily_subset)) {
      daily_seq <- daily_subset$seq_num[i]
      
      # Find weekly observations that follow this daily observation
      # Assuming that 7 sequence steps roughly corresponds to a week
      # Adjust this parameter based on your data collection frequency
      seq_diff <- weekly_subset$seq_num - daily_seq
      valid_weeks <- which(seq_diff > 0 & seq_diff <= 7)
      
      if (length(valid_weeks) > 0) {
        # Find the closest following weekly observation
        closest_week <- valid_weeks[which.min(seq_diff[valid_weeks])]
        
        # Create row with merged data
        merged_row <- data.frame(
          # Daily data
          ID = daily_subset$ID[i],
          SDT_daily = daily_subset$SDT[i],
          SDa = daily_subset$SDa[i],
          DDSC = daily_subset$DDSC[i],
          DESHELP = daily_subset$DESHELP[i],
          PROBREC = daily_subset$PROBREC[i],
          TXREADY = daily_subset$TXREADY[i],
          PRESSTX = daily_subset$PRESSTX[i],
          TXNEEDS = daily_subset$TXNEEDS[i],
          # Include LDSC explicitly
          LDSC = daily_subset$LDSC[i],
          # Weekly data
          SDT_weekly = weekly_subset$SDT[closest_week],
          PTEF = weekly_subset$PTEF[closest_week],
          PTEI = weekly_subset$PTEI[closest_week],
          HSBANY = weekly_subset$HSBANY[closest_week],
          HSBF = weekly_subset$HSBF[closest_week],
          HSBI = weekly_subset$HSBI[closest_week],
          HSBO = weekly_subset$HSBO[closest_week],
          # Get help-seeking variables for the following week if available
          HSBF2 = ifelse(closest_week + 1 <= nrow(weekly_subset), 
                         weekly_subset$HSBF[closest_week + 1], NA),
          HSBI2 = ifelse(closest_week + 1 <= nrow(weekly_subset), 
                         weekly_subset$HSBI[closest_week + 1], NA)
        )
        
        # Add any demographic variables that were in daily_subset but not explicitly listed above
        for (col in names(daily_subset)) {
          if (!(col %in% names(merged_row)) && col != "seq_num") {
            merged_row[[col]] <- daily_subset[[col]][i]
          }
        }
        
        # Add this row to our result dataset
        if (nrow(matched_data) == 0) {
          matched_data <- merged_row
        } else {
          matched_data <- rbind(matched_data, merged_row)
        }
      }
    }
  }
  
  return(matched_data)
}

# Apply the matching function and handle empty results
df_analysis <- match_by_sequence(df_merged, df_weekly)

# Check if any rows were matched
if (nrow(df_analysis) == 0) {
  warning("No matching observations were found between daily and weekly data. This could be because:
          1. There are no overlapping IDs between the datasets
          2. The sequence matching criteria might be too strict
          3. The ordering of observations might not reflect the temporal sequence")
  
  # Create a simple merged dataset instead, without trying to match by time
  message("Creating a simple merged dataset instead...")
  df_analysis <- merge(df_merged, 
                       df_weekly[, c("ID", "PTEF", "PTEI", "HSBANY", "HSBF", "HSBI", "HSBO", "HSBF2", "HSBI2")], 
                       by = "ID", all.x = TRUE)
}

# Verify that LDSC is in the final analysis dataset
if ("LDSC" %in% names(df_analysis)) {
  message("LDSC variable is present in the final analysis dataset")
} else {
  warning("LDSC variable is missing from the final analysis dataset")
  # If missing, add it directly from the original LDSC file
  df_analysis <- merge(df_analysis, df_ldsc, by = "ID", all.x = TRUE)
}

# Create person-mean centered variables for level 1 predictors
# First check if we have DDSC values to center
if (!("DDSC" %in% names(df_analysis))) {
  warning("DDSC variable not found in the analysis dataset")
} else if (all(is.na(df_analysis$DDSC))) {
  warning("All DDSC values are NA, cannot calculate centered values")
} else {
  # Using base R approach to avoid potential issues
  # Create a temp data frame with just ID and DDSC (removing NAs)
  temp_data <- df_analysis[!is.na(df_analysis$DDSC), c("ID", "DDSC")]
  
  if (nrow(temp_data) > 0) {
    # Calculate person means
    person_means <- aggregate(DDSC ~ ID, data = temp_data, FUN = mean, na.rm = TRUE)
    names(person_means)[2] <- "DDSC_mean"
    
    # Merge the means back to the original dataframe
    df_analysis <- merge(df_analysis, person_means, by = "ID", all.x = TRUE)
    
    # Calculate the centered values
    df_analysis$DDSC_centered <- df_analysis$DDSC - df_analysis$DDSC_mean
    
    message("Successfully created person-mean centered DDSC variable")
  } else {
    warning("No valid DDSC values found for centering")
    # Add empty columns to avoid errors in later code
    df_analysis$DDSC_mean <- NA
    df_analysis$DDSC_centered <- NA
  }
}

# ======== DESCRIPTIVE STATISTICS ========

# Summary of key variables
descriptives <- describe(df_analysis[, c("DDSC", "LDSC", "PTEF", "PTEI", "HSBF", "HSBI", "DESHELP", "PROBREC", "TXREADY")])
print(descriptives)

# Correlation matrix
cor_matrix <- cor(df_analysis[, c("DDSC", "LDSC", "PTEF", "PTEI", "DESHELP", "PROBREC", "TXREADY")], 
                  use = "pairwise.complete.obs")
print(cor_matrix)

# Check frequency of help-seeking behaviors
help_seeking_freq <- df_analysis %>%
  summarise(
    formal_help = sum(HSBF, na.rm = TRUE) / sum(!is.na(HSBF)) * 100,
    informal_help = sum(HSBI, na.rm = TRUE) / sum(!is.na(HSBI)) * 100,
    any_help = sum(HSBANY, na.rm = TRUE) / sum(!is.na(HSBANY)) * 100
  )
print(help_seeking_freq)

# ======== HYPOTHESIS TESTING ========

# Part 1: Examine if daily discrimination predicts decreased formal help seeking

# Create binary version of HSBF
df_analysis$HSBF_binary <- ifelse(df_analysis$HSBF == -99, NA, df_analysis$HSBF)

# Check values
table(df_analysis$HSBF_binary, useNA = "always")

# Model 1: Daily discrimination predicting formal help seeking
model_formal <- glmer(HSBF ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                      data = df_analysis, 
                      family = binomial)
summary(model_formal)



# Model 1: Daily discrimination predicting formal help seeking
# Create binary version of HSBI
df_analysis$HSBI_binary <- ifelse(df_analysis$HSBI == -99, NA, df_analysis$HSBI)

# Check values
table(df_analysis$HSBI_binary, useNA = "always")

# Now run the model with the fixed variable
model_informal <- glmer(HSBI_binary ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                        data = df_analysis, 
                        family = binomial)

# Model 2: Daily discrimination predicting informal help seeking
model_informal <- glmer(HSBI ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                        data = df_analysis, 
                        family = binomial)
summary(model_informal)

# Part 2: Test if discrimination affects perceived treatment effectiveness (potential mediator)

# Model 3: Daily discrimination predicting perceived treatment effectiveness for formal treatment
model_ptef <- lmer(PTEF ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                   data = df_analysis)
summary(model_ptef)

# Model 4: Daily discrimination predicting perceived treatment effectiveness for informal treatment
model_ptei <- lmer(PTEI ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                   data = df_analysis)
summary(model_ptei)

# Part 3: Test if perceived treatment effectiveness predicts help seeking

# Model 5: Perceived treatment effectiveness predicting formal help seeking
model_ptef_help <- glmer(HSBF_binary ~ PTEF + DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                         data = df_analysis, 
                         family = binomial)

summary(model_ptef_help)

# Model 6: Perceived treatment effectiveness predicting informal help seeking
model_ptei_help <- glmer(HSBI_binary ~ PTEI + DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                         data = df_analysis, 
                         family = binomial)
summary(model_ptei_help)

# Part 4: Test if discrimination predicts treatment readiness and other barriers

# Model 7: Daily discrimination predicting desire for help
model_deshelp <- lmer(DESHELP ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                      data = df_analysis)
summary(model_deshelp)

# Model 8: Daily discrimination predicting problem recognition
model_probrec <- lmer(PROBREC ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                      data = df_analysis)
summary(model_probrec)

# Model 9: Daily discrimination predicting treatment readiness
model_txready <- lmer(TXREADY ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                      data = df_analysis)
summary(model_txready)

# Part 5: Test if treatment readiness and other barriers predict help seeking

# Model 10: Treatment readiness predicting formal help seeking
model_ready_formal <- glmer(HSBF_binary ~ TXREADY + DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                            data = df_analysis, 
                            family = binomial)
summary(model_ready_formal)

# Model 11: Problem recognition predicting formal help seeking
model_prob_formal <- glmer(HSBF_binary ~ PROBREC + DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                           data = df_analysis, 
                           family = binomial)
summary(model_prob_formal)

# Part 6: Look at subsequent help seeking (lagged dependent variables)

# Create binary versions for lagged help-seeking variables
df_analysis$HSBF2_binary <- ifelse(df_analysis$HSBF2 == -99, NA, df_analysis$HSBF2)
df_analysis$HSBI2_binary <- ifelse(df_analysis$HSBI2 == -99, NA, df_analysis$HSBI2)

# Check the values
table(df_analysis$HSBF2_binary, useNA = "always")
table(df_analysis$HSBI2_binary, useNA = "always")

# Update the models to use the binary versions
model_formal_lag <- glmer(HSBF2_binary ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                          data = df_analysis, 
                          family = binomial)

model_informal_lag <- glmer(HSBI2_binary ~ DDSC_centered + DDSC_mean + LDSC + GEN + ETH + (1|ID), 
                            data = df_analysis, 
                            family = binomial)

# ======== MODERATION ANALYSES ========

# Test if gender moderates the relationship between discrimination and help seeking

# Model 14: Gender moderation for formal help seeking
model_gender_formal <- glmer(HSBF_binary ~ DDSC_centered*GEN + DDSC_mean + LDSC + ETH + (1|ID), 
                             data = df_analysis, 
                             family = binomial)
summary(model_gender_formal)

# Model 15: Gender moderation for informal help seeking
model_gender_informal <- glmer(HSBI_binary ~ DDSC_centered*GEN + DDSC_mean + LDSC + ETH + (1|ID), 
                               data = df_analysis, 
                               family = binomial)
summary(model_gender_informal)

# Test if ethnicity moderates the relationship between discrimination and help seeking

# Model 16: Ethnicity moderation for formal help seeking
model_eth_formal <- glmer(HSBF_binary ~ DDSC_centered*ETH + DDSC_mean + LDSC + GEN + (1|ID), 
                          data = df_analysis, 
                          family = binomial)
summary(model_eth_formal)

# Model 17: Ethnicity moderation for informal help seeking
model_eth_informal <- glmer(HSBI_binary ~ DDSC_centered*ETH + DDSC_mean + LDSC + GEN + (1|ID), 
                            data = df_analysis, 
                            family = binomial)
summary(model_eth_informal)

# ======== VISUALIZATION ========

# Create plots to visualize key findings

# Plot 1: Daily discrimination and formal help seeking
ggplot(df_analysis, aes(x = DDSC, y = HSBF_binary)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Relationship between Daily Discrimination and\nFormal Help Seeking",
       x = "Daily Discrimination Score",
       y = "Probability of Formal Help Seeking") +
  theme_minimal()

# Plot 2: Daily discrimination and perceived treatment effectiveness
ggplot(df_analysis, aes(x = DDSC, y = PTEF)) +
  geom_smooth(method = "lm") +
  labs(title = "Relationship between Daily Discrimination and\nPerceived Treatment Effectiveness",
       x = "Daily Discrimination Score",
       y = "Perceived Treatment Effectiveness (Formal)") +
  theme_minimal()

# Plot 3: Perceived treatment effectiveness and formal help seeking
ggplot(df_analysis, aes(x = PTEF, y = HSBF_binary)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Relationship between Perceived Treatment Effectiveness and\nFormal Help Seeking",
       x = "Perceived Treatment Effectiveness (Formal)",
       y = "Probability of Formal Help Seeking") +
  theme_minimal()

# Plot 4: Compare formal vs informal help seeking by discrimination level
df_analysis %>%
  mutate(DDSC_cat = cut(DDSC, breaks = 3, labels = c("Low", "Medium", "High"))) %>%
  group_by(DDSC_cat) %>%
  summarise(
    formal = mean(HSBF_binary, na.rm = TRUE),
    informal = mean(HSBI_binary, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(formal, informal), names_to = "Help_Type", values_to = "Proportion") %>%
  ggplot(aes(x = DDSC_cat, y = Proportion, fill = Help_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of Help Seeking by \nDiscrimination Level",
       x = "Daily Discrimination Category",
       y = "Proportion Seeking Help",
       fill = "Type of Help") +
  theme_minimal()