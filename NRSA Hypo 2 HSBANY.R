# Set  working directory to where the CSV file is located
setwd("/Users/lacosta3/Desktop/NRSA Data Cleaning Done")

# Using HSBANY (general help-seeking) instead of specific HSBF/HSBI

# Load required packages
library(lme4)      # For multilevel modeling
library(lmerTest)  # For p-values in mixed models
library(dplyr)     # For data manipulation
library(tidyr)     # For data reshaping
library(ggplot2)   # For visualization
library(mediation) # For mediation analysis

#--------------------------------------------------------------
# STEP 1: Import and prepare datasets
#--------------------------------------------------------------

# Import datasets
df_lv1 <- read.csv("DDSC and HSA_level 1.csv")    # Level-1 data (daily discrimination, treatment attitudes)
df_ldsc <- read.csv("LDSC_lv2.csv")      # Level-2 data (lifetime discrimination)
df_demo <- read.csv("Level 2 demos.csv") # Level-2 data (demographics)
df_hlm <- read.csv("PTE and HSB_level 1.csv")    # Level-1 data (treatment effectiveness and utilization)

# Replace -99 with NA in all dataframes
df_lv1[df_lv1 == -99] <- NA
df_ldsc[df_ldsc == -99] <- NA
df_demo[df_demo == -99] <- NA
df_hlm[df_hlm == -99] <- NA

#--------------------------------------------------------------
# STEP 2: Create Level-2 dataset (person-level variables)
#--------------------------------------------------------------

# Create base level-2 dataset from LDSC
df_level2 <- df_ldsc %>%
  select(ID, LDSC) %>%
  # Center LDSC for analysis
  mutate(LDSC_centered = scale(LDSC, center = TRUE, scale = FALSE)[,1])

# Add gender (start with gender to filter to participants with complete data)
if ("GEN" %in% names(df_demo)) {
  df_level2 <- df_level2 %>%
    inner_join(df_demo %>% 
                 select(ID, GEN) %>%
                 # GEN: 0=Male, 1=Female
                 mutate(
                   GEN_factor = factor(GEN, 
                                       levels = c(0, 1),
                                       labels = c("Male", "Female")),
                   GEN_num = GEN  # Numeric version for modeling
                 ),
               by = "ID")
} else {
  df_level2 <- df_level2 %>%
    mutate(GEN_num = NA)
}

# Print summary of level-2 dataset
cat("\nLevel-2 dataset summary:\n")
cat("Number of participants:", nrow(df_level2), "\n")

#--------------------------------------------------------------
# STEP 3: Prepare Level-1 datasets
#--------------------------------------------------------------

# Prepare treatment effectiveness dataset
df_effectiveness <- df_hlm %>%
  select(ID, SDT, PTEF, PTEI, HSBANY) %>%  # Now including HSBANY
  arrange(ID, SDT) %>%
  group_by(ID) %>%
  mutate(
    time = row_number(),
    # Person-mean centering
    PTEF_mean = mean(PTEF, na.rm = TRUE),
    PTEF_within = PTEF - PTEF_mean,
    PTEI_mean = mean(PTEI, na.rm = TRUE),
    PTEI_within = PTEI - PTEI_mean
  ) %>%
  ungroup()

# Prepare daily discrimination dataset
df_ddsc <- df_lv1 %>%
  select(ID, SDT, DDSC) %>%
  arrange(ID, SDT) %>%
  group_by(ID) %>%
  mutate(
    time = row_number(),
    # Person-mean centering
    DDSC_mean = mean(DDSC, na.rm = TRUE),
    DDSC_within = DDSC - DDSC_mean
  ) %>%
  ungroup()

# Create dataset for Path 5 analysis (Discrimination → Treatment Effectiveness)
# Try joining on ID and SDT first
df_path5 <- df_effectiveness %>%
  inner_join(df_ddsc, by = c("ID", "SDT"))

# If no matches, try alternative join strategy
if(nrow(df_path5) == 0) {
  cat("\nNo exact ID-SDT matches found. Trying alternative matching approach.\n")
  
  # Try matching by sequence within participant
  common_ids <- intersect(unique(df_effectiveness$ID), unique(df_ddsc$ID))
  
  # Build dataset by matching observations by order within participant
  df_path5 <- data.frame()
  for(id in common_ids) {
    # Get all observations for this ID from both datasets
    eff_data <- df_effectiveness %>% filter(ID == id)
    disc_data <- df_ddsc %>% filter(ID == id)
    
    # Match by sequence position if possible
    n_obs <- min(nrow(eff_data), nrow(disc_data))
    
    if(n_obs > 0) {
      # Use the first n_obs from each dataset
      matched_data <- cbind(
        eff_data[1:n_obs, ],
        disc_data[1:n_obs, setdiff(names(disc_data), c("ID", "SDT", "time"))]
      )
      df_path5 <- rbind(df_path5, matched_data)
    }
  }
}

# Add level-2 variables to the Path 5 dataset
df_path5 <- df_path5 %>%
  left_join(df_level2, by = "ID")

# Create dataset for Path 6 analysis (Treatment Effectiveness → Help-Seeking (HSBANY))
df_path6 <- df_hlm %>%
  group_by(ID) %>%
  summarize(
    PTEF = mean(PTEF, na.rm = TRUE),
    PTEI = mean(PTEI, na.rm = TRUE),
    # Convert to binary for logistic regression
    HSBANY = as.numeric(mean(HSBANY, na.rm = TRUE) > 0.5),
    # Also calculate proportion of weeks with help-seeking as an alternative
    prop_HSBANY = mean(HSBANY, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Add level-2 variables
  left_join(df_level2, by = "ID")

# Create dataset for mediation analysis
df_mediation <- df_path6

# Create dataset for Ad Hoc analyses (attitudes)
df_attitudes <- df_lv1 %>%
  select(ID, SDT, DESHELP, PROBREC, TXREADY, PRESSTX, TXNEEDS) %>%
  arrange(ID, SDT) %>%
  group_by(ID) %>%
  mutate(
    time = row_number(),
    # Person-mean centering
    DESHELP_mean = mean(DESHELP, na.rm = TRUE),
    DESHELP_within = DESHELP - DESHELP_mean,
    PROBREC_mean = mean(PROBREC, na.rm = TRUE),
    PROBREC_within = PROBREC - PROBREC_mean,
    TXREADY_mean = mean(TXREADY, na.rm = TRUE),
    TXREADY_within = TXREADY - TXREADY_mean,
    PRESSTX_mean = mean(PRESSTX, na.rm = TRUE),
    PRESSTX_within = PRESSTX - PRESSTX_mean,
    TXNEEDS_mean = mean(TXNEEDS, na.rm = TRUE),
    TXNEEDS_within = TXNEEDS - TXNEEDS_mean
  ) %>%
  ungroup()

# Create matched dataset for Ad Hoc analyses
# Try to match on ID and SDT first
df_adhoc <- df_attitudes %>%
  inner_join(df_ddsc, by = c("ID", "SDT"))

# If no matches, try alternative join strategy
if(nrow(df_adhoc) == 0) {
  # Try matching by sequence within participant
  common_ids <- intersect(unique(df_attitudes$ID), unique(df_ddsc$ID))
  
  # Build dataset by matching observations by order within participant
  df_adhoc <- data.frame()
  for(id in common_ids) {
    # Get all observations for this ID from both datasets
    att_data <- df_attitudes %>% filter(ID == id)
    disc_data <- df_ddsc %>% filter(ID == id)
    
    # Match by sequence position if possible
    n_obs <- min(nrow(att_data), nrow(disc_data))
    
    if(n_obs > 0) {
      # Use the first n_obs from each dataset
      matched_data <- cbind(
        att_data[1:n_obs, ],
        disc_data[1:n_obs, setdiff(names(disc_data), c("ID", "SDT", "time"))]
      )
      df_adhoc <- rbind(df_adhoc, matched_data)
    }
  }
}

# Add level-2 variables
df_adhoc <- df_adhoc %>%
  left_join(df_level2, by = "ID")

# Print dataset summaries
cat("\nPath 5 dataset (Discrimination → Effectiveness):", nrow(df_path5), "observations,", length(unique(df_path5$ID)), "participants\n")
cat("Path 6 dataset (Effectiveness → Help-Seeking):", nrow(df_path6), "participants\n")
cat("Ad Hoc dataset (Discrimination → Attitudes):", nrow(df_adhoc), "observations,", length(unique(df_adhoc$ID)), "participants\n")

#--------------------------------------------------------------
# STEP 4: Test Path 5 (Discrimination → Treatment Effectiveness)
#--------------------------------------------------------------

# Check if we have enough data for multilevel modeling
if(nrow(df_path5) > 20 && length(unique(df_path5$ID)) > 5) {
  # Model 1a: LDSC and DDSC predicting formal treatment effectiveness (PTEF)
  model_1a <- try(lmer(PTEF ~ DDSC_within + DDSC_mean + LDSC_centered + 
                         (1 | ID), 
                       data = df_path5, 
                       REML = FALSE))
  
  if(!inherits(model_1a, "try-error")) {
    cat("\nPath 5a: Effects on Formal Treatment Effectiveness (PTEF)\n")
    print(summary(model_1a))
  } else {
    cat("\nUnable to fit multilevel model for PTEF. Using simpler model.\n")
    
    # Alternative: person-level model
    df_path5_person <- df_path5 %>%
      group_by(ID) %>%
      summarize(
        PTEF_mean = mean(PTEF, na.rm = TRUE),
        DDSC_mean = mean(DDSC, na.rm = TRUE),
        LDSC_centered = first(LDSC_centered),
        GEN_num = first(GEN_num),
        .groups = "drop"
      )
    
    model_1a_person <- lm(PTEF_mean ~ DDSC_mean + LDSC_centered, data = df_path5_person)
    cat("\nPath 5a (Person-Level): Effects on Formal Treatment Effectiveness (PTEF)\n")
    print(summary(model_1a_person))
  }
  
  # Model 1b: LDSC and DDSC predicting informal treatment effectiveness (PTEI)
  model_1b <- try(lmer(PTEI ~ DDSC_within + DDSC_mean + LDSC_centered + 
                         (1 | ID), 
                       data = df_path5, 
                       REML = FALSE))
  
  if(!inherits(model_1b, "try-error")) {
    cat("\nPath 5b: Effects on Informal Treatment Effectiveness (PTEI)\n")
    print(summary(model_1b))
  } else {
    cat("\nUnable to fit multilevel model for PTEI. Using simpler model.\n")
    
    # Alternative: person-level model
    model_1b_person <- lm(PTEI_mean ~ DDSC_mean + LDSC_centered, data = df_path5_person)
    cat("\nPath 5b (Person-Level): Effects on Informal Treatment Effectiveness (PTEI)\n")
    print(summary(model_1b_person))
  }
  
} else {
  cat("\nInsufficient data for multilevel modeling. Using person-level analysis.\n")
  
  # Create person-level dataset
  df_path5_person <- df_path5 %>%
    group_by(ID) %>%
    summarize(
      PTEF_mean = mean(PTEF, na.rm = TRUE),
      PTEI_mean = mean(PTEI, na.rm = TRUE),
      DDSC_mean = mean(DDSC, na.rm = TRUE),
      LDSC_centered = first(LDSC_centered),
      GEN_num = first(GEN_num),
      .groups = "drop"
    )
  
  # Person-level models
  model_1a_person <- lm(PTEF_mean ~ DDSC_mean + LDSC_centered, data = df_path5_person)
  cat("\nPath 5a (Person-Level): Effects on Formal Treatment Effectiveness (PTEF)\n")
  print(summary(model_1a_person))
  
  model_1b_person <- lm(PTEI_mean ~ DDSC_mean + LDSC_centered, data = df_path5_person)
  cat("\nPath 5b (Person-Level): Effects on Informal Treatment Effectiveness (PTEI)\n")
  print(summary(model_1b_person))
}

#--------------------------------------------------------------
# STEP A: Directly test the relationship between discrimination and help-seeking
#--------------------------------------------------------------

# Test if discrimination directly predicts help-seeking (not mediated by treatment effectiveness)
discrimination_helpseek <- glm(HSBANY ~ LDSC_centered, 
                               data = df_path6,
                               family = binomial(link = "logit"))

cat("\nDirect effect of Lifetime Discrimination on Any Help-Seeking Behavior:\n")
print(summary(discrimination_helpseek))
cat("Odds Ratio:", exp(coef(discrimination_helpseek)["LDSC_centered"]), "\n\n")

# Alternative model using the proportion of weeks with help-seeking (continuous outcome)
discrimination_helpseek_prop <- lm(prop_HSBANY ~ LDSC_centered, data = df_path6)
cat("Direct effect of Lifetime Discrimination on Proportion of Help-Seeking:\n")
print(summary(discrimination_helpseek_prop))

#--------------------------------------------------------------
# STEP 5: Test Path 6 (Treatment Effectiveness → Any Help-Seeking)
#--------------------------------------------------------------

# Model 4a: PTEF predicting any help-seeking behavior (HSBANY)
model_4a <- try(glm(HSBANY ~ PTEF, 
                    data = df_path6,
                    family = binomial(link = "logit")))

if(!inherits(model_4a, "try-error")) {
  cat("\nPath 6a: Effect of Formal Treatment Effectiveness on Any Help-Seeking (HSBANY)\n")
  print(summary(model_4a))
  cat("Odds Ratio:", exp(coef(model_4a)["PTEF"]), "\n")
}

# Model 4b: PTEI predicting any help-seeking behavior (HSBANY)
model_4b <- try(glm(HSBANY ~ PTEI, 
                    data = df_path6,
                    family = binomial(link = "logit")))

if(!inherits(model_4b, "try-error")) {
  cat("\nPath 6b: Effect of Informal Treatment Effectiveness on Any Help-Seeking (HSBANY)\n")
  print(summary(model_4b))
  cat("Odds Ratio:", exp(coef(model_4b)["PTEI"]), "\n")
}

# Model 4c: Both PTEF and PTEI predicting any help-seeking behavior (HSBANY)
model_4c <- try(glm(HSBANY ~ PTEF + PTEI, 
                    data = df_path6,
                    family = binomial(link = "logit")))

if(!inherits(model_4c, "try-error")) {
  cat("\nPath 6c: Effect of Both Treatment Effectiveness Types on Any Help-Seeking (HSBANY)\n")
  print(summary(model_4c))
  cat("Odds Ratio PTEF:", exp(coef(model_4c)["PTEF"]), "\n")
  cat("Odds Ratio PTEI:", exp(coef(model_4c)["PTEI"]), "\n")
}

#--------------------------------------------------------------
# STEP 6: Test Mediation (LDSC → PTEF/PTEI → HSBANY)
#--------------------------------------------------------------

# Model 5a: LDSC → PTEF → HSBANY
med_model_formal <- try(mediate(
  model.m = lm(PTEF ~ LDSC_centered, data = df_mediation),
  model.y = glm(HSBANY ~ PTEF + LDSC_centered, 
                data = df_mediation, 
                family = binomial(link = "logit")),
  treat = "LDSC_centered",
  mediator = "PTEF",
  boot = TRUE,
  sims = 1000))

if(!inherits(med_model_formal, "try-error")) {
  cat("\nMediation Analysis: LDSC → PTEF → HSBANY\n")
  print(summary(med_model_formal))
}

# Model 5b: LDSC → PTEI → HSBANY
med_model_informal <- try(mediate(
  model.m = lm(PTEI ~ LDSC_centered, data = df_mediation),
  model.y = glm(HSBANY ~ PTEI + LDSC_centered, 
                data = df_mediation, 
                family = binomial(link = "logit")),
  treat = "LDSC_centered",
  mediator = "PTEI",
  boot = TRUE,
  sims = 1000))

if(!inherits(med_model_informal, "try-error")) {
  cat("\nMediation Analysis: LDSC → PTEI → HSBANY\n")
  print(summary(med_model_informal))
}

#--------------------------------------------------------------
# STEP 7: Test Gender Moderation
#--------------------------------------------------------------

# Only test gender moderation if we have gender data
if("GEN_num" %in% names(df_path5)) {
  # For multilevel models (if available)
  if(exists("model_1a") && !inherits(model_1a, "try-error")) {
    # Gender moderation for PTEF
    model_2a <- try(lmer(PTEF ~ DDSC_within * GEN_num + DDSC_mean * GEN_num + 
                           LDSC_centered * GEN_num + 
                           (1 | ID), 
                         data = df_path5, 
                         REML = FALSE))
    
    if(!inherits(model_2a, "try-error")) {
      cat("\nGender Moderation: Effects on Formal Treatment Effectiveness (PTEF)\n")
      print(summary(model_2a))
    }
  }
}

# Gender moderation of the direct effect of discrimination on help-seeking
discrimination_helpseek_gender <- glm(HSBANY ~ LDSC_centered * GEN_num, 
                                      data = df_path6,
                                      family = binomial(link = "logit"))

cat("\nGender Moderation of Direct Effect of Discrimination on Help-Seeking:\n")
print(summary(discrimination_helpseek_gender))

#--------------------------------------------------------------
# STEP 8: Ad Hoc Analyses - Treatment Attitudes → Any Help-Seeking
#--------------------------------------------------------------

# Create person-level dataset for treatment attitudes
df_attitudes_person <- df_adhoc %>%
  group_by(ID) %>%
  summarize(
    DESHELP = mean(DESHELP, na.rm = TRUE),
    PROBREC = mean(PROBREC, na.rm = TRUE),
    TXREADY = mean(TXREADY, na.rm = TRUE),
    PRESSTX = mean(PRESSTX, na.rm = TRUE),
    TXNEEDS = mean(TXNEEDS, na.rm = TRUE),
    .groups = "drop"
  )

# Merge with treatment utilization data
df_attitudes_util <- df_attitudes_person %>%
  inner_join(df_path6 %>% select(ID, HSBANY), by = "ID")

# Model for any help-seeking behavior
model_7a <- try(glm(HSBANY ~ DESHELP + PROBREC + TXREADY + PRESSTX + TXNEEDS, 
                    data = df_attitudes_util,
                    family = binomial(link = "logit")))

if(!inherits(model_7a, "try-error")) {
  cat("\nTreatment Attitudes Predicting Any Help-Seeking Behavior (HSBANY)\n")
  print(summary(model_7a))
}

# Test each attitude predictor individually
for(attitude in c("DESHELP", "PROBREC", "TXREADY", "PRESSTX", "TXNEEDS")) {
  formula_str <- paste("HSBANY ~", attitude)
  model <- try(glm(as.formula(formula_str), 
                   data = df_attitudes_util,
                   family = binomial(link = "logit")))
  
  if(!inherits(model, "try-error")) {
    cat("\nEffect of", attitude, "on Any Help-Seeking Behavior:\n")
    print(summary(model))
    cat("Odds Ratio:", exp(coef(model)[2]), "\n")
  }
}

#--------------------------------------------------------------
# STEP 9: Summarize Key Findings
#--------------------------------------------------------------

cat("\n\n==========================================\n")
cat("SUMMARY OF KEY FINDINGS (USING HSBANY)\n")
cat("==========================================\n\n")

# Summarize direct effect of discrimination on help-seeking
cat("Direct Effect of Discrimination on Help-Seeking:\n")
if(!inherits(discrimination_helpseek, "try-error")) {
  beta <- coef(discrimination_helpseek)["LDSC_centered"]
  se <- summary(discrimination_helpseek)$coefficients["LDSC_centered", "Std. Error"]
  p <- summary(discrimination_helpseek)$coefficients["LDSC_centered", "Pr(>|z|)"]
  or <- exp(beta)
  
  cat("- Lifetime discrimination → Help-seeking: β =", round(beta, 2), 
      ", p =", round(p, 3), ", OR =", round(or, 2), "\n")
  
  if(p < 0.05) {
    cat("  * SIGNIFICANT: Lifetime discrimination directly affects help-seeking behavior\n")
  } else {
    cat("  * NOT SIGNIFICANT: No direct effect of lifetime discrimination on help-seeking\n")
  }
}

# Summarize treatment effectiveness effects on help-seeking
cat("\nEffects of Treatment Effectiveness on Help-Seeking:\n")
if(!inherits(model_4a, "try-error")) {
  beta <- coef(model_4a)["PTEF"]
  se <- summary(model_4a)$coefficients["PTEF", "Std. Error"]
  p <- summary(model_4a)$coefficients["PTEF", "Pr(>|z|)"]
  or <- exp(beta)
  
  cat("- Formal effectiveness → Help-seeking: β =", round(beta, 2), 
      ", p =", round(p, 3), ", OR =", round(or, 2), "\n")
  
  if(p < 0.05) {
    cat("  * SIGNIFICANT: Perceived formal treatment effectiveness affects help-seeking\n")
  } else {
    cat("  * NOT SIGNIFICANT: No effect of formal treatment effectiveness on help-seeking\n")
  }
}

if(!inherits(model_4b, "try-error")) {
  beta <- coef(model_4b)["PTEI"]
  se <- summary(model_4b)$coefficients["PTEI", "Std. Error"]
  p <- summary(model_4b)$coefficients["PTEI", "Pr(>|z|)"]
  or <- exp(beta)
  
  cat("- Informal effectiveness → Help-seeking: β =", round(beta, 2), 
      ", p =", round(p, 3), ", OR =", round(or, 2), "\n")
  
  if(p < 0.05) {
    cat("  * SIGNIFICANT: Perceived informal treatment effectiveness affects help-seeking\n")
  } else {
    cat("  * NOT SIGNIFICANT: No effect of informal treatment effectiveness on help-seeking\n")
  }
}

# Summarize gender moderation
if(!inherits(discrimination_helpseek_gender, "try-error")) {
  beta <- coef(discrimination_helpseek_gender)["LDSC_centered:GEN_num"]
  se <- summary(discrimination_helpseek_gender)$coefficients["LDSC_centered:GEN_num", "Std. Error"]
  p <- summary(discrimination_helpseek_gender)$coefficients["LDSC_centered:GEN_num", "Pr(>|z|)"]
  
  cat("\nGender Moderation:\n")
  cat("- LDSC × Gender interaction: β =", round(beta, 2), ", p =", round(p, 3), "\n")
  
  if(p < 0.05) {
    cat("  * SIGNIFICANT: Gender moderates the effect of discrimination on help-seeking\n")
  } else {
    cat("  * NOT SIGNIFICANT: No gender moderation detected\n")
  }
}

# Overall conclusion about hypothesis
cat("\nOverall Conclusion:\n")
cat("Based on using HSBANY (any help-seeking behavior) as the outcome:\n")

path5_significant <- FALSE
path6_significant <- FALSE

# Check Path 5 significance
if(exists("model_1a") && !inherits(model_1a, "try-error")) {
  p_ldsc <- coef(summary(model_1a))["LDSC_centered", "Pr(>|t|)"]
  if(p_ldsc < 0.05) path5_significant <- TRUE
}

# Check Path 6 significance
if(exists("model_4a") && !inherits(model_4a, "try-error")) {
  p_ptef <- summary(model_4a)$coefficients["PTEF", "Pr(>|z|)"]
  if(p_ptef < 0.05) path6_significant <- TRUE
}

if(exists("model_4b") && !inherits(model_4b, "try-error")) {
  p_ptei <- summary(model_4b)$coefficients["PTEI", "Pr(>|z|)"]
  if(p_ptei < 0.05) path6_significant <- TRUE
}

# Direct effect significance
direct_effect_significant <- FALSE
if(!inherits(discrimination_helpseek, "try-error")) {
  p_direct <- summary(discrimination_helpseek)$coefficients["LDSC_centered", "Pr(>|z|)"]
  if(p_direct < 0.05) direct_effect_significant <- TRUE
}

# Overall conclusion
if(path5_significant && path6_significant) {
  cat("- FULL SUPPORT for hypothesis: Discrimination affects treatment effectiveness, which affects help-seeking\n")
} else if(path5_significant && !path6_significant) {
  cat("- PARTIAL SUPPORT: Discrimination affects treatment effectiveness, but this doesn't affect help-seeking\n")
} else if(!path5_significant && path6_significant) {
  cat("- PARTIAL SUPPORT: Treatment effectiveness affects help-seeking, but isn't affected by discrimination\n")
} else if(direct_effect_significant) {
  cat("- ALTERNATIVE PATHWAY: Discrimination directly affects help-seeking without mediation through effectiveness\n")
} else {
  cat("- NO SUPPORT for hypothesis: Neither path is significant\n")
}

cat("\nNote on Statistical Power: These analyses may still be limited by sample size and missing data.\n")