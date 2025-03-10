# Set  working directory to where the CSV file is located
setwd("/Users/lacosta3/Desktop/NRSA Data Cleaning Done")

# Testing the pathway: Discrimination → Treatment Effectiveness → Treatment Utilization

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

# Add demographics (start with gender to filter to participants with complete data)
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

# Add other demographic variables if available
if ("ETH" %in% names(df_demo)) {
  df_level2 <- df_level2 %>%
    left_join(df_demo %>% 
                select(ID, ETH) %>%
                mutate(ETH_factor = factor(ETH,
                                           levels = c(0, 1, 2),
                                           labels = c("Black only", "Hispanic only", "Both Black and Hispanic"))),
              by = "ID")
}

if ("EDU" %in% names(df_demo)) {
  df_level2 <- df_level2 %>%
    left_join(df_demo %>% 
                select(ID, EDU) %>%
                mutate(EDU_factor = factor(EDU,
                                           levels = c(0, 1, 2, 3, 4),
                                           labels = c("Less than high school", 
                                                      "Completed high school/GED",
                                                      "Some college",
                                                      "Completed 2-year degree",
                                                      "Completed 4-year degree"))),
              by = "ID")
}

if ("EMP" %in% names(df_demo)) {
  df_level2 <- df_level2 %>%
    left_join(df_demo %>% 
                select(ID, EMP) %>%
                mutate(EMP_factor = factor(EMP,
                                           levels = c(0, 1),
                                           labels = c("Not employed", "Employed"))),
              by = "ID")
}

# Print summary of level-2 dataset
cat("\nLevel-2 dataset summary:\n")
cat("Number of participants:", nrow(df_level2), "\n")

#--------------------------------------------------------------
# STEP 3: Prepare Level-1 datasets
#--------------------------------------------------------------

# Prepare treatment effectiveness dataset
df_effectiveness <- df_hlm %>%
  select(ID, SDT, PTEF, PTEI) %>%
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

# Create dataset for Path 6 analysis (Treatment Effectiveness → Treatment Utilization)
df_path6 <- df_hlm %>%
  group_by(ID) %>%
  summarize(
    PTEF = mean(PTEF, na.rm = TRUE),
    PTEI = mean(PTEI, na.rm = TRUE),
    # Convert to binary for logistic regression
    HSBF = as.numeric(mean(HSBF, na.rm = TRUE) > 0.5),
    HSBI = as.numeric(mean(HSBI, na.rm = TRUE) > 0.5),
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
cat("Path 6 dataset (Effectiveness → Utilization):", nrow(df_path6), "participants\n")
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
# STEP 5: Test Path 6 (Treatment Effectiveness → Treatment Utilization)
#--------------------------------------------------------------

# Model 4a: PTEF predicting formal treatment utilization (HSBF)
model_4a <- try(glm(HSBF ~ PTEF, 
                    data = df_path6,
                    family = binomial(link = "logit")))

if(!inherits(model_4a, "try-error")) {
  cat("\nPath 6a: Effect of Formal Treatment Effectiveness on Utilization (HSBF)\n")
  print(summary(model_4a))
  cat("Odds Ratio:", exp(coef(model_4a)["PTEF"]), "\n")
}

# Model 4b: PTEI predicting informal treatment utilization (HSBI)
model_4b <- try(glm(HSBI ~ PTEI, 
                    data = df_path6,
                    family = binomial(link = "logit")))

if(!inherits(model_4b, "try-error")) {
  cat("\nPath 6b: Effect of Informal Treatment Effectiveness on Utilization (HSBI)\n")
  print(summary(model_4b))
  cat("Odds Ratio:", exp(coef(model_4b)["PTEI"]), "\n")
}

#--------------------------------------------------------------
# STEP 6: Test Mediation (LDSC → PTEF/PTEI → HSBF/HSBI)
#--------------------------------------------------------------

# Model 5a: LDSC → PTEF → HSBF
med_model_formal <- try(mediate(
  model.m = lm(PTEF ~ LDSC_centered, data = df_mediation),
  model.y = glm(HSBF ~ PTEF + LDSC_centered, 
                data = df_mediation, 
                family = binomial(link = "logit")),
  treat = "LDSC_centered",
  mediator = "PTEF",
  boot = TRUE,
  sims = 1000))

if(!inherits(med_model_formal, "try-error")) {
  cat("\nMediation Analysis: LDSC → PTEF → HSBF\n")
  print(summary(med_model_formal))
}

# Model 5b: LDSC → PTEI → HSBI
med_model_informal <- try(mediate(
  model.m = lm(PTEI ~ LDSC_centered, data = df_mediation),
  model.y = glm(HSBI ~ PTEI + LDSC_centered, 
                data = df_mediation, 
                family = binomial(link = "logit")),
  treat = "LDSC_centered",
  mediator = "PTEI",
  boot = TRUE,
  sims = 1000))

if(!inherits(med_model_informal, "try-error")) {
  cat("\nMediation Analysis: LDSC → PTEI → HSBI\n")
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
    
    # Gender moderation for PTEI
    model_2b <- try(lmer(PTEI ~ DDSC_within * GEN_num + DDSC_mean * GEN_num + 
                           LDSC_centered * GEN_num + 
                           (1 | ID), 
                         data = df_path5, 
                         REML = FALSE))
    
    if(!inherits(model_2b, "try-error")) {
      cat("\nGender Moderation: Effects on Informal Treatment Effectiveness (PTEI)\n")
      print(summary(model_2b))
    }
  } else if(exists("model_1a_person")) {
    # Person-level gender moderation
    model_2a_person <- try(lm(PTEF_mean ~ DDSC_mean * GEN_num + LDSC_centered * GEN_num, 
                              data = df_path5_person))
    
    if(!inherits(model_2a_person, "try-error")) {
      cat("\nGender Moderation (Person-Level): Effects on PTEF\n")
      print(summary(model_2a_person))
    }
    
    model_2b_person <- try(lm(PTEI_mean ~ DDSC_mean * GEN_num + LDSC_centered * GEN_num, 
                              data = df_path5_person))
    
    if(!inherits(model_2b_person, "try-error")) {
      cat("\nGender Moderation (Person-Level): Effects on PTEI\n")
      print(summary(model_2b_person))
    }
  }
}

#--------------------------------------------------------------
# STEP 8: Ad Hoc Analyses - Discrimination and Treatment Attitudes
#--------------------------------------------------------------

# Function to test discrimination effects on treatment attitudes
test_attitude_model <- function(attitude_var, attitude_label) {
  # Formula for within-person centered variable
  formula_str <- paste0(attitude_var, " ~ DDSC_within + DDSC_mean + LDSC_centered + (1 | ID)")
  
  # Try to fit the model
  model <- try(lmer(as.formula(formula_str), data = df_adhoc, REML = FALSE))
  
  if(!inherits(model, "try-error")) {
    cat("\nDiscrimination Effects on", attitude_label, "\n")
    print(summary(model))
    return(model)
  } else {
    cat("\nUnable to fit multilevel model for", attitude_label, "\n")
    return(NULL)
  }
}

# Test discrimination effects on each treatment attitude
attitude_vars <- c("DESHELP", "PROBREC", "TXREADY", "PRESSTX", "TXNEEDS")
attitude_labels <- c("Desire for Help", "Problem Recognition", "Treatment Readiness", 
                     "Pressure to be in Treatment", "Treatment Needs")

attitude_models <- list()
for(i in 1:length(attitude_vars)) {
  attitude_models[[i]] <- test_attitude_model(attitude_vars[i], attitude_labels[i])
}

#--------------------------------------------------------------
# STEP 9: Ad Hoc Analyses - Treatment Attitudes → Treatment Utilization
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
  inner_join(df_path6 %>% select(ID, HSBF, HSBI), by = "ID")

# Model for formal treatment utilization
model_7a <- try(glm(HSBF ~ DESHELP + PROBREC + TXREADY + PRESSTX + TXNEEDS, 
                    data = df_attitudes_util,
                    family = binomial(link = "logit")))

if(!inherits(model_7a, "try-error")) {
  cat("\nTreatment Attitudes Predicting Formal Treatment Utilization (HSBF)\n")
  print(summary(model_7a))
}

# Model for informal treatment utilization
model_7b <- try(glm(HSBI ~ DESHELP + PROBREC + TXREADY + PRESSTX + TXNEEDS, 
                    data = df_attitudes_util,
                    family = binomial(link = "logit")))

if(!inherits(model_7b, "try-error")) {
  cat("\nTreatment Attitudes Predicting Informal Treatment Utilization (HSBI)\n")
  print(summary(model_7b))
}

#--------------------------------------------------------------
# STEP 10: Export Results
#--------------------------------------------------------------

# Create results directory
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# Function to safely export model results
export_model <- function(model, filename, type = "mixed") {
  if(exists(model) && !inherits(get(model), "try-error")) {
    if(type == "mixed") {
      write.csv(coef(summary(get(model))), file.path(results_dir, filename))
    } else if(type == "lm") {
      write.csv(summary(get(model))$coefficients, file.path(results_dir, filename))
    } else if(type == "mediation" && !is.null(get(model)$d0)) {
      med_results <- data.frame(
        Effect = c("Indirect", "Direct", "Total", "Proportion Mediated"),
        Estimate = c(get(model)$d0, get(model)$z0, get(model)$tau.coef, get(model)$n0),
        p_value = c(get(model)$d0.p, get(model)$z0.p, get(model)$tau.p, NA)
      )
      write.csv(med_results, file.path(results_dir, filename), row.names = FALSE)
    }
  }
}

# Export key model results
if(exists("model_1a") && !inherits(model_1a, "try-error")) {
  export_model("model_1a", "results_path5_formal.csv")
} else if(exists("model_1a_person")) {
  export_model("model_1a_person", "results_path5_formal_person.csv", "lm")
}

if(exists("model_1b") && !inherits(model_1b, "try-error")) {
  export_model("model_1b", "results_path5_informal.csv")
} else if(exists("model_1b_person")) {
  export_model("model_1b_person", "results_path5_informal_person.csv", "lm")
}

if(exists("model_4a")) export_model("model_4a", "results_path6_formal.csv", "lm")
if(exists("model_4b")) export_model("model_4b", "results_path6_informal.csv", "lm")
if(exists("med_model_formal")) export_model("med_model_formal", "results_mediation_formal.csv", "mediation")
if(exists("med_model_informal")) export_model("med_model_informal", "results_mediation_informal.csv", "mediation")

cat("\nAnalysis complete. Results exported to 'results' directory.\n")