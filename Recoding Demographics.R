# Set working directory to where the CSV file is located
setwd("/Users/lacosta3/Desktop/NRSA Data Cleaning Done")

# Read the data
data <- read.csv('orig_Level 2_demo.csv', na.strings = c("", "NA"))

# Load necessary libraries
library(dplyr)

# First, convert all "MISSING_DATA" and "UNASKED_QUESTION" to NA throughout the dataset
data_clean <- data %>%
  mutate(across(everything(), 
                ~ifelse(. %in% c("MISSING_DATA", "UNASKED_QUESTION"), NA, .)))

# Now create all new variables in one transformation
data_transformed <- data_clean %>%
  mutate(
    # Gender variable: 0 = Male, 1 = Female, NA = Other categories
    GEN = case_when(
      GEN_1 == 1 ~ 0,           # Male
      GEN_5 == 1 ~ 1,           # Female
      GEN_6 == 1 ~ NA_real_,    # Transgender Male (missing)
      GEN_7 == 1 ~ NA_real_,    # Transgender Female (missing)
      GEN_8 == 1 ~ NA_real_,    # Non-binary (missing)
      GEN_9 == 1 ~ NA_real_,    # Prefer not to respond (missing)
      TRUE ~ NA_real_           # If none of the above conditions are met
    ),
    
    # Ethnicity variable: 0 = Black only, 1 = Hispanic only, 2 = Both, NA = Neither
    ETH = case_when(
      ETH_6 == 1 & ETH_7 == 1 ~ 2,      # Both Black and Hispanic
      ETH_6 == 1 & ETH_7 == 0 ~ 0,      # Black only
      ETH_6 == 0 & ETH_7 == 1 ~ 1,      # Hispanic only
      TRUE ~ NA_real_                   # Neither Black nor Hispanic (missing)
    ),
    
    # Education variable recode
    EDU = case_when(
      EDU == 1 ~ 0,    # Less than high school → 0
      EDU == 6 ~ 1,    # Completed high school or GED → 1
      EDU == 7 ~ 2,    # Some college → 2
      EDU == 8 ~ 3,    # Completed 2-year degree → 3
      EDU == 9 ~ 4,    # Completed 4-year degree → 4
      TRUE ~ NA_real_  # If none of the above conditions are met
    ),
    
    # Employment variable: 1 = employed, 0 = not employed
    EMP = case_when(
      EMP_1 == 1 ~ 1,    # Employed = 1
      EMP_6 == 0 ~ 0,    # Not employed = 0
      TRUE ~ NA_real_    # If neither condition is met (e.g., missing data)
    )
  )

# View the first few rows to verify the new variables
head(data_transformed)

# Save the transformed data to a new file
write.csv(data_transformed, "Level 2 demos.csv", row.names = FALSE)