# =============================================================================
# LONGITUDINAL ANALYSIS OF METABOLIC SYNDROME IN PLHIV ON DTG-BASED ART
# =============================================================================
# VISEND Trial Data (2026) - 144-Week Analysis from Zambia
# =============================================================================
#
# Author:      Mpanji Siwingwa
# Date:        28 September 2026
# Version:     3.0
# Repository:  https://github.com/MpanjiSiwingwa/VISEND-MetSyn-Longitudinal-Analysis
# License:     MIT
#
# Description:
#   This script analyzes longitudinal metabolic syndrome outcomes among
#   people living with HIV receiving dolutegravir-based and PI-based ART
#   regimens over 144 weeks in the VISEND trial in Zambia.
#
# Correspondence: mpanjisiwingwa@gmail.com
# =============================================================================


# =============================================================================
# 0. Setup
# =============================================================================

set.seed(123)


# -----------------------------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------------------------

# Data manipulation and cleaning
library(dplyr)
library(tidyr)
library(forcats)
library(janitor)
library(purrr)
library(stringr)

# Statistical modelling
library(geepack)
library(emmeans)
library(ggeffects)
library(margins)
library(epitools)
library(splines)

# Model tidying and results
library(broom)
library(knitr)
library(gtsummary)
library(gt)

# Data visualization
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)
library(patchwork)

# Publication and document outputs
library(flextable)
library(officer)
library(gtExtras)
library(sysfonts)
library(webshot)


# =============================================================================
# 1. Load and Clean Dataset
# =============================================================================

# -----------------------------------------------------------------------------
# Data configuration
# -----------------------------------------------------------------------------

data_file <- "MetSyn_dataset_2025_cleaned_080626.csv"


# -----------------------------------------------------------------------------
# Check that the input dataset exists
# -----------------------------------------------------------------------------

if (!file.exists(data_file)) {
  stop(
    paste0(
      "Input dataset not found: ", data_file, "\n\n",
      "Please place the authorized VISEND dataset in the expected ",
      "data directory or provide the correct file path."
    )
  )
}


# -----------------------------------------------------------------------------
# Read dataset
# -----------------------------------------------------------------------------

data <- read.csv(
  data_file,
  na.strings = c("", "NA"),
  strip.white = TRUE,
  check.names = FALSE
)


# -----------------------------------------------------------------------------
# Define variables required for the analysis
# -----------------------------------------------------------------------------

required_variables <- c(
  "Trial_number",
  "Event_Name",
  "Drug_Code.",
  "Sex",
  "Age..yrs.",
  "CD4_count.cells.µl.",
  "Viral_Load.cp.ml.",
  "Blood_Sugar..mmol.L.",
  "Triglycerides..mmol.L.",
  "Cholesterol_HDL_.mmol.L.",
  "Bp_Systolic..mmHg.",
  "Bp_Diastolic..mmHg.",
  "Waist_Circumference.cm.",
  "Alcohol_Consuption",
  "Tobbacco_Use",
  "Diabetes_Mellitus_status",
  "Educational_Level",
  "Previous_Regimen"
)


# -----------------------------------------------------------------------------
# Check that all required variables are present
# -----------------------------------------------------------------------------

missing_variables <- setdiff(
  required_variables,
  names(data)
)

if (length(missing_variables) > 0) {
  stop(
    paste0(
      "The following required variables are missing from the dataset:\n",
      paste(missing_variables, collapse = ", ")
    )
  )
}


# -----------------------------------------------------------------------------
# Keep only variables required for the analysis
# -----------------------------------------------------------------------------

data <- data %>%
  select(all_of(required_variables))


# -----------------------------------------------------------------------------
# Inspect imported data
# -----------------------------------------------------------------------------

cat("\nVariables successfully loaded:\n")
print(names(data))

cat("\nNumber of rows:", nrow(data), "\n")
cat("Number of participants:", n_distinct(data$Trial_number), "\n")

cat("\nOriginal visit labels:\n")
print(sort(unique(data$Event_Name)))


# =============================================================================
# 1.1 Standardize Visit Time Points
# =============================================================================

# -----------------------------------------------------------------------------
# Standardize raw visit labels
# -----------------------------------------------------------------------------

data <- data %>%
  mutate(
    Event_Name = str_squish(Event_Name),
    
    Event_Name = case_when(
      
      str_to_lower(Event_Name) == "baseline" ~ "Baseline",
      
      str_to_lower(Event_Name) %in% c(
        "week_24",
        "week 24"
      ) ~ "Week_24",
      
      str_to_lower(Event_Name) %in% c(
        "week_48",
        "week 48"
      ) ~ "Week_48",
      
      # VISEND analysis time point:
      # raw Week 64 is analysed as Week 72
      str_to_lower(Event_Name) %in% c(
        "week_64",
        "week 64",
        "week_72",
        "week 72"
      ) ~ "Week_72",
      
      str_to_lower(Event_Name) %in% c(
        "week_96",
        "week 96"
      ) ~ "Week_96",
      
      # VISEND analysis time point:
      # raw Week 112 is analysed as Week 120
      str_to_lower(Event_Name) %in% c(
        "week_112",
        "week _112",
        "week 112",
        "week_120",
        "week 120"
      ) ~ "Week_120",
      
      str_to_lower(Event_Name) %in% c(
        "week_144",
        "week_ 144",
        "week 144"
      ) ~ "Week_144",
      
      TRUE ~ Event_Name
    )
  )


# -----------------------------------------------------------------------------
# Create ordered analysis visit variable
# -----------------------------------------------------------------------------

data <- data %>%
  mutate(
    Event_Name = factor(
      Event_Name,
      levels = c(
        "Baseline",
        "Week_24",
        "Week_48",
        "Week_72",
        "Week_96",
        "Week_120",
        "Week_144"
      )
    )
  )


# -----------------------------------------------------------------------------
# Create human-readable visit labels
# -----------------------------------------------------------------------------

data <- data %>%
  mutate(
    clean_event = factor(
      Event_Name,
      levels = c(
        "Baseline",
        "Week_24",
        "Week_48",
        "Week_72",
        "Week_96",
        "Week_120",
        "Week_144"
      ),
      labels = c(
        "Baseline",
        "Week 24",
        "Week 48",
        "Week 72",
        "Week 96",
        "Week 120",
        "Week 144"
      )
    )
  )


# -----------------------------------------------------------------------------
# Verify standardized visit labels
# -----------------------------------------------------------------------------

cat("\nStandardized visit labels:\n")
print(table(data$Event_Name, useNA = "ifany"))

cat("\nOrdered analysis visits:\n")
print(levels(data$Event_Name))


# =============================================================================
# 2. Define Metabolic Syndrome
# =============================================================================
# NCEP ATP III criteria
# =============================================================================

flag_metabolic_syndrome <- function(data) {
  
  data %>%
    mutate(
      
      # -----------------------------------------------------------------------
      # Individual MetS components
      # -----------------------------------------------------------------------
      
      glucose_abnormal =
        Blood_Sugar..mmol.L. >= 5.6,
      
      triglycerides_abnormal =
        Triglycerides..mmol.L. >= 1.69,
      
      hdl_abnormal =
        if_else(
          Sex == "Male",
          Cholesterol_HDL_.mmol.L. < 1.0,
          Cholesterol_HDL_.mmol.L. < 1.3
        ),
      
      blood_pressure_abnormal =
        Bp_Systolic..mmHg. >= 130 |
        Bp_Diastolic..mmHg. >= 85,
      
      waist_abnormal =
        if_else(
          Sex == "Male",
          Waist_Circumference.cm. >= 94,
          Waist_Circumference.cm. >= 80
        ),
      
      
      # -----------------------------------------------------------------------
      # Number of available MetS components
      # -----------------------------------------------------------------------
      
      MetSyn_components_available =
        rowSums(
          !is.na(
            cbind(
              glucose_abnormal,
              triglycerides_abnormal,
              hdl_abnormal,
              blood_pressure_abnormal,
              waist_abnormal
            )
          )
        ),
      
      
      # -----------------------------------------------------------------------
      # Number of abnormal MetS components
      # -----------------------------------------------------------------------
      
      MetSyn_components_abnormal =
        rowSums(
          cbind(
            glucose_abnormal,
            triglycerides_abnormal,
            hdl_abnormal,
            blood_pressure_abnormal,
            waist_abnormal
          ),
          na.rm = TRUE
        ),
      
      
      # -----------------------------------------------------------------------
      # Final MetS classification
      # -----------------------------------------------------------------------
      
      Metabolic_Syndrome = case_when(
        MetSyn_components_available < 5 ~ NA_character_,
        MetSyn_components_abnormal >= 3 ~ "MetSyn",
        TRUE ~ "No MetSyn"
      )
    )
}


# -----------------------------------------------------------------------------
# Apply MetS definition
# -----------------------------------------------------------------------------

data_full <- flag_metabolic_syndrome(data)


# =============================================================================
# 3. Exclude Participants with Metabolic Syndrome at Baseline
# =============================================================================

# -----------------------------------------------------------------------------
# Identify baseline MetS participants
# -----------------------------------------------------------------------------

baseline_metsyn_ids <- data_full %>%
  filter(
    Event_Name == "Baseline",
    Metabolic_Syndrome == "MetSyn"
  ) %>%
  pull(Trial_number)


# -----------------------------------------------------------------------------
# Exclude baseline MetS participants and all their follow-up observations
# -----------------------------------------------------------------------------

data_clean <- data_full %>%
  filter(
    !Trial_number %in% baseline_metsyn_ids
  )


# -----------------------------------------------------------------------------
# Verify missing baseline information after exclusion
# -----------------------------------------------------------------------------

baseline_missing_after <- data_clean %>%
  filter(Event_Name == "Baseline") %>%
  summarise(
    missing_alcohol =
      sum(is.na(Alcohol_Consuption)),
    
    missing_tobacco =
      sum(is.na(Tobbacco_Use)),
    
    missing_diabetes =
      sum(is.na(Diabetes_Mellitus_status)),
    
    total_missing =
      sum(
        is.na(Alcohol_Consuption) |
          is.na(Tobbacco_Use) |
          is.na(Diabetes_Mellitus_status)
      ),
    
    total_participants =
      n_distinct(Trial_number)
  )

print(baseline_missing_after)


# -----------------------------------------------------------------------------
# Summary after baseline MetS exclusion
# -----------------------------------------------------------------------------

rows_before <- nrow(data_full)
rows_after <- nrow(data_clean)
rows_removed <- rows_before - rows_after

n_participants_before <-
  n_distinct(data_full$Trial_number)

n_participants_after <-
  n_distinct(data_clean$Trial_number)

n_participants_removed <-
  n_participants_before - n_participants_after


cat("\n--- Summary After Baseline MetS Exclusion ---\n")
cat("Original total rows:     ", rows_before, "\n")
cat("Rows removed:            ", rows_removed, "\n")
cat("Rows remaining:          ", rows_after, "\n")
cat("Participants removed:    ", n_participants_removed, "\n")
cat("Participants remaining:  ", n_participants_after, "\n")


# =============================================================================
# 4. Descriptive Follow-up Summaries
# =============================================================================

# -----------------------------------------------------------------------------
# Sample size by time point and regimen
# -----------------------------------------------------------------------------

sample_size_by_week <- data_clean %>%
  group_by(Event_Name, Drug_Code.) %>%
  summarise(
    n = n_distinct(Trial_number),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Drug_Code.,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    Total = rowSums(
      across(where(is.numeric))
    )
  ) %>%
  arrange(Event_Name)

print(sample_size_by_week)


# -----------------------------------------------------------------------------
# MetS by event and regimen
# -----------------------------------------------------------------------------

metSyn_by_event_regimen <- data_clean %>%
  group_by(
    Event_Name,
    Drug_Code.,
    Metabolic_Syndrome
  ) %>%
  summarise(
    n = n_distinct(Trial_number),
    .groups = "drop"
  ) %>%
  filter(
    Metabolic_Syndrome == "MetSyn"
  ) %>%
  pivot_wider(
    names_from = Drug_Code.,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    Total = rowSums(
      across(where(is.numeric))
    )
  )

print(metSyn_by_event_regimen)


# -----------------------------------------------------------------------------
# MetS prevalence at baseline after exclusion
# -----------------------------------------------------------------------------

baseline_prevalence <- data_clean %>%
  filter(Event_Name == "Baseline") %>%
  count(Metabolic_Syndrome) %>%
  mutate(
    prevalence = round(
      100 * n / sum(n),
      1
    )
  )

print(baseline_prevalence)


# -----------------------------------------------------------------------------
# Number of observations by time point
# -----------------------------------------------------------------------------

weekly_summary <- data_clean %>%
  group_by(Event_Name) %>%
  summarise(
    samples_remaining = n(),
    .groups = "drop"
  ) %>%
  arrange(Event_Name)

print(weekly_summary)

# =============================================================================
# 5. Follow-up Summary
# =============================================================================

expected_n <- 899

followup_summary <- data_clean %>%
  group_by(Event_Name) %>%
  summarise(
    participants_with_data = n(),
    
    new_MetSyn_cases =
      sum(New_MetSyn_Case, na.rm = TRUE),
    
    incidence =
      round(
        100 * new_MetSyn_cases / participants_with_data,
        1
      ),
    
    missing_data =
      sum(is.na(New_MetSyn_Case)),
    
    missed_appointments =
      expected_n - participants_with_data,
    
    .groups = "drop"
  ) %>%
  arrange(Event_Name)

print(followup_summary)

# =============================================================================
# 6. Variable Re-coding and Data Preparation
# =============================================================================

# -----------------------------------------------------------------------------
# Complete data set (includes baseline)
# -----------------------------------------------------------------------------

metabolic_data_clean <- data_clean %>%
  rename(
    Regimen_Type = Drug_Code.,
    ID = Trial_number,
    Alcohol_Consumption = Alcohol_Consuption
  ) %>%
  mutate(
    
    # -------------------------------------------------------------------------
    # Recode ART regimen
    # -------------------------------------------------------------------------
    
    Regimen_Type = case_when(
      Regimen_Type %in% c(
        "AZT+3TC+LPVr",
        "AZT+3TC+ATVr"
      ) ~ "PI_based",
      
      Regimen_Type == "TAFED" ~ "TAFED",
      
      Regimen_Type == "TLD" ~ "TLD",
      
      TRUE ~ NA_character_
    ),
    
    Regimen_Type = factor(
      Regimen_Type,
      levels = c(
        "PI_based",
        "TAFED",
        "TLD"
      )
    ),
    
    
    # -------------------------------------------------------------------------
    # MetS binary outcome
    # -------------------------------------------------------------------------
    
    metabolic_num = case_when(
      Metabolic_Syndrome == "MetSyn" ~ 1,
      Metabolic_Syndrome == "No MetSyn" ~ 0,
      TRUE ~ NA_real_
    ),
    
    
    # -------------------------------------------------------------------------
    # Standardized continuous variables
    # -------------------------------------------------------------------------
    
    Age_scaled =
      as.numeric(scale(Age..yrs.)),
    
    CD4_scaled =
      as.numeric(scale(CD4_count.cells.µl.)),
    
    ViralLoad_scaled =
      as.numeric(scale(Viral_Load.cp.ml.)),
    
    
    # -------------------------------------------------------------------------
    # Participant ID
    # -------------------------------------------------------------------------
    
    ID = as.factor(ID)
  )


# -----------------------------------------------------------------------------
# Verify recoded variables
# -----------------------------------------------------------------------------

cat("\nART regimen distribution:\n")
print(table(
  metabolic_data_clean$Regimen_Type,
  useNA = "ifany"
))

cat("\nMetS outcome distribution:\n")
print(table(
  metabolic_data_clean$metabolic_num,
  useNA = "ifany"
))

cat("\nVisit distribution:\n")
print(table(
  metabolic_data_clean$Event_Name,
  useNA = "ifany"
))


# =============================================================================
# 6.1 Baseline Descriptive Statistics
# =============================================================================

# -----------------------------------------------------------------------------
# Calculate median and IQR for age
# -----------------------------------------------------------------------------

age_summary <- metabolic_data_clean %>%
  filter(Event_Name == "Baseline") %>%
  summarise(
    median_age = median(
      Age..yrs.,
      na.rm = TRUE
    ),
    
    IQR_lower = quantile(
      Age..yrs.,
      0.25,
      na.rm = TRUE
    ),
    
    IQR_upper = quantile(
      Age..yrs.,
      0.75,
      na.rm = TRUE
    )
  )


# -----------------------------------------------------------------------------
# Print publication-ready age summary
# -----------------------------------------------------------------------------

cat(
  "Median age:",
  age_summary$median_age,
  "years (IQR:",
  age_summary$IQR_lower,
  "–",
  age_summary$IQR_upper,
  ")\n"
)


# =============================================================================
# 6.2 Test Distribution and Summarize Baseline Variables
# =============================================================================

# -----------------------------------------------------------------------------
# Filter to baseline data
# -----------------------------------------------------------------------------

baseline_data <- metabolic_data_clean %>%
  filter(Event_Name == "Baseline")


# -----------------------------------------------------------------------------
# Test distribution using Shapiro-Wilk
# -----------------------------------------------------------------------------

shapiro_age <- shapiro.test(
  baseline_data$Age..yrs.
)

shapiro_cd4 <- shapiro.test(
  baseline_data$CD4_count.cells.µl.
)

shapiro_vl <- shapiro.test(
  baseline_data$Viral_Load.cp.ml.
)


# -----------------------------------------------------------------------------
# Summaries for Age
# -----------------------------------------------------------------------------

age_summary <- baseline_data %>%
  summarise(
    median_age =
      median(
        Age..yrs.,
        na.rm = TRUE
      ),
    
    IQR_lower =
      quantile(
        Age..yrs.,
        0.25,
        na.rm = TRUE
      ),
    
    IQR_upper =
      quantile(
        Age..yrs.,
        0.75,
        na.rm = TRUE
      ),
    
    mean_age =
      mean(
        Age..yrs.,
        na.rm = TRUE
      ),
    
    sd_age =
      sd(
        Age..yrs.,
        na.rm = TRUE
      )
  )


# -----------------------------------------------------------------------------
# Summaries for CD4 count
# -----------------------------------------------------------------------------

cd4_summary <- baseline_data %>%
  summarise(
    median_cd4 =
      median(
        CD4_count.cells.µl.,
        na.rm = TRUE
      ),
    
    IQR_lower =
      quantile(
        CD4_count.cells.µl.,
        0.25,
        na.rm = TRUE
      ),
    
    IQR_upper =
      quantile(
        CD4_count.cells.µl.,
        0.75,
        na.rm = TRUE
      ),
    
    mean_cd4 =
      mean(
        CD4_count.cells.µl.,
        na.rm = TRUE
      ),
    
    sd_cd4 =
      sd(
        CD4_count.cells.µl.,
        na.rm = TRUE
      )
  )


# -----------------------------------------------------------------------------
# Summaries for Viral Load
# -----------------------------------------------------------------------------

vl_summary <- baseline_data %>%
  summarise(
    median_vl =
      median(
        Viral_Load.cp.ml.,
        na.rm = TRUE
      ),
    
    IQR_lower =
      quantile(
        Viral_Load.cp.ml.,
        0.25,
        na.rm = TRUE
      ),
    
    IQR_upper =
      quantile(
        Viral_Load.cp.ml.,
        0.75,
        na.rm = TRUE
      ),
    
    mean_vl =
      mean(
        Viral_Load.cp.ml.,
        na.rm = TRUE
      ),
    
    sd_vl =
      sd(
        Viral_Load.cp.ml.,
        na.rm = TRUE
      )
  )


# =============================================================================
# 6.3 Print Baseline Results
# =============================================================================

cat(
  "Age: median",
  age_summary$median_age,
  "years (IQR",
  age_summary$IQR_lower,
  "–",
  age_summary$IQR_upper,
  "); mean",
  round(age_summary$mean_age, 1),
  "±",
  round(age_summary$sd_age, 1),
  "(Shapiro-Wilk p =",
  signif(shapiro_age$p.value, 3),
  ")\n"
)


cat(
  "CD4 count: median",
  cd4_summary$median_cd4,
  "cells/µL (IQR",
  cd4_summary$IQR_lower,
  "–",
  cd4_summary$IQR_upper,
  "); mean",
  round(cd4_summary$mean_cd4, 1),
  "±",
  round(cd4_summary$sd_cd4, 1),
  "(Shapiro-Wilk p =",
  signif(shapiro_cd4$p.value, 3),
  ")\n"
)


cat(
  "Viral load: median",
  vl_summary$median_vl,
  "cp/mL (IQR",
  vl_summary$IQR_lower,
  "–",
  vl_summary$IQR_upper,
  "); mean",
  round(vl_summary$mean_vl, 1),
  "±",
  round(vl_summary$sd_vl, 1),
  "(Shapiro-Wilk p =",
  signif(shapiro_vl$p.value, 3),
  ")\n"
)

# ============================================================================
# Figure 1: Baseline Distributions and Normality Assessment
# Figure Completeness and Quality
# ============================================================================


# Plot A: Age Distribution
p1 <- ggplot(baseline_data, aes(x = Age..yrs.)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  labs(title = "A: Age Distribution", x = "Age (years)", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# Plot B: Q-Q Plot for Age
p2 <- ggqqplot(baseline_data$Age..yrs., title = "B: Q-Q Plot: Age") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# Plot C: CD4 Count Distribution
p3 <- ggplot(baseline_data, aes(x = CD4_count.cells.µl.)) +
  geom_histogram(binwidth = 50, fill = "darkgreen", color = "white") +
  labs(title = "C: CD4 Count Distribution", x = "CD4 (cells/µL)", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# Plot D: Q-Q Plot for CD4 Count
p4 <- ggqqplot(baseline_data$CD4_count.cells.µl., title = "D: Q-Q Plot: CD4 Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# Create log-transformed viral load
baseline_data$log_vl <- log10(baseline_data$Viral_Load.cp.ml. + 1)

# Plot E: Log10 Viral Load Distribution
p5 <- ggplot(baseline_data, aes(x = log_vl)) +
  geom_histogram(binwidth = 0.2, fill = "firebrick", color = "white") +
  labs(title = "E: Log10 Viral Load Distribution", 
       x = "Log10 Viral Load (cp/mL)", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# Plot F: Q-Q Plot for Log10 Viral Load
p6 <- ggqqplot(baseline_data$log_vl, title = "F: Q-Q Plot: Log10 Viral Load") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# Combine all plots (no auto-labels since titles have A-F)
combined_plot <- ggarrange(p1, p2, p3, p4, p5, p6,
                           ncol = 2, nrow = 3,
                           labels = NULL)

# Add main title and margins
combined_plot_titled <- annotate_figure(
  combined_plot,
  top = text_grob(
    "Figure S1. Distribution and Normality Assessment of Baseline Clinical Variables",
    face = "bold",
    size = 14
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))

# Display the plot
print(combined_plot_titled)

# =============================================================================
# Table 1: Baseline Characteristics of the Study Population by ART Regimen
# =============================================================================

# -----------------------------------------------------------------------------
# Step 1: Create baseline dataset
# -----------------------------------------------------------------------------

baseline_data <- data_clean %>%
  filter(Event_Name == "Baseline") %>%
  mutate(
    
    # -------------------------------------------------------------------------
    # Age categories
    # -------------------------------------------------------------------------
    
    Age_group = case_when(
      Age..yrs. >= 40 ~ "≥40 years",
      Age..yrs. < 40 ~ "<40 years",
      TRUE ~ NA_character_
    ),
    
    
    # -------------------------------------------------------------------------
    # CD4 categories
    # -------------------------------------------------------------------------
    
    CD4_group = case_when(
      CD4_count.cells.µl. <= 100 ~ "≤100",
      CD4_count.cells.µl. <= 350 ~ "101–350",
      CD4_count.cells.µl. > 350 ~ ">350",
      TRUE ~ NA_character_
    ),
    
    
    # -------------------------------------------------------------------------
    # Viral load categories
    # -------------------------------------------------------------------------
    
    ViralLoad_group = case_when(
      Viral_Load.cp.ml. >= 1000 ~ "≥1000",
      Viral_Load.cp.ml. < 1000 ~ "<1000",
      TRUE ~ NA_character_
    )
  )


# -----------------------------------------------------------------------------
# Step 2: Rename variables for presentation
# -----------------------------------------------------------------------------

baseline_data <- baseline_data %>%
  rename(
    Alcohol_Consumption = Alcohol_Consuption,
    Tobacco_Use = Tobbacco_Use
  )


# -----------------------------------------------------------------------------
# Step 3: Create MetS component indicators
# -----------------------------------------------------------------------------

baseline_data <- baseline_data %>%
  mutate(
    
    BloodSugar_high = factor(
      Blood_Sugar..mmol.L. >= 5.6,
      levels = c(FALSE, TRUE),
      labels = c("Normal", "Elevated")
    ),
    
    Triglycerides_high = factor(
      Triglycerides..mmol.L. >= 1.69,
      levels = c(FALSE, TRUE),
      labels = c("Normal", "Elevated")
    ),
    
    HDL_low = factor(
      case_when(
        Sex == "Male" ~ Cholesterol_HDL_.mmol.L. < 1.0,
        Sex == "Female" ~ Cholesterol_HDL_.mmol.L. < 1.3,
        TRUE ~ NA
      ),
      levels = c(FALSE, TRUE),
      labels = c("Normal", "Low")
    ),
    
    Hypertension = factor(
      Bp_Systolic..mmHg. >= 130 |
        Bp_Diastolic..mmHg. >= 85,
      levels = c(FALSE, TRUE),
      labels = c(
        "No Hypertension",
        "Hypertension"
      )
    ),
    
    Abdominal_obesity = factor(
      case_when(
        Sex == "Male" ~ Waist_Circumference.cm. >= 94,
        Sex == "Female" ~ Waist_Circumference.cm. >= 80,
        TRUE ~ NA
      ),
      levels = c(FALSE, TRUE),
      labels = c(
        "No Abdominal Obesity",
        "Abdominal Obesity"
      )
    )
  )


# -----------------------------------------------------------------------------
# Step 4: Custom Fisher's exact test with larger workspace
# -----------------------------------------------------------------------------

fisher_large <- function(data, variable, by, ...) {
  
  test_result <- stats::fisher.test(
    table(
      data[[by]],
      data[[variable]]
    ),
    workspace = 2e7
  )
  
  tibble::tibble(
    p.value = test_result$p.value
  )
}


# -----------------------------------------------------------------------------
# Step 5: Generate Table 1
# -----------------------------------------------------------------------------

table1 <- baseline_data %>%
  select(
    Regimen_Type,
    Age_group,
    Sex,
    CD4_group,
    ViralLoad_group,
    Diabetes_Mellitus_status,
    Tobacco_Use,
    Alcohol_Consumption,
    Educational_Level,
    BloodSugar_high,
    Triglycerides_high,
    HDL_low,
    Hypertension,
    Abdominal_obesity
  ) %>%
  
  tbl_summary(
    by = Regimen_Type,
    
    label = list(
      Age_group ~ "Age",
      Sex ~ "Sex",
      CD4_group ~ "CD4 Count",
      ViralLoad_group ~ "Viral Load",
      Diabetes_Mellitus_status ~ "Diabetes Mellitus Status",
      Tobacco_Use ~ "Tobacco Use",
      Alcohol_Consumption ~ "Alcohol Use",
      Educational_Level ~ "Level of Education",
      BloodSugar_high ~ "Blood Sugar ≥5.6 mmol/L",
      Triglycerides_high ~ "Triglycerides ≥1.69 mmol/L",
      HDL_low ~ "Low HDL Cholesterol",
      Hypertension ~ "Hypertension (SBP ≥130 or DBP ≥85 mmHg)",
      Abdominal_obesity ~ "Abdominal Obesity"
    ),
    
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    missing = "no"
  ) %>%
  
  add_p(
    test = list(
      Tobacco_Use ~ fisher.test,
      Alcohol_Consumption ~ fisher_large,
      Diabetes_Mellitus_status ~ fisher.test,
      Educational_Level ~ chisq.test,
      BloodSugar_high ~ fisher.test,
      Triglycerides_high ~ fisher.test,
      HDL_low ~ fisher.test,
      Hypertension ~ fisher.test,
      Abdominal_obesity ~ fisher.test
    )
  ) %>%
  
  modify_caption(
    "**Table 1: Baseline Characteristics of the Study Population by ART Regimen Group**"
  ) %>%
  
  bold_labels()


# -----------------------------------------------------------------------------
# Step 6: Print Table 1
# -----------------------------------------------------------------------------

print(table1)

# =============================================================================
# Table 1A. Sociodemographic and Clinical Characteristics
# =============================================================================

table1A <- baseline_data %>%
  select(
    Age_group,
    Sex,
    CD4_group,
    ViralLoad_group,
    Diabetes_Mellitus_status,
    Tobacco_Use,
    Alcohol_Consumption,
    Educational_Level,
    Previous_Regimen,
    Regimen_Type
  ) %>%
  tbl_summary(
    by = Regimen_Type,
    label = list(
      Age_group ~ "Age",
      Sex ~ "Sex",
      CD4_group ~ "CD4 Count (cells/µL)",
      ViralLoad_group ~ "Viral Load",
      Diabetes_Mellitus_status ~ "Diabetes Mellitus",
      Tobacco_Use ~ "Tobacco Use",
      Alcohol_Consumption ~ "Alcohol Use",
      Educational_Level ~ "Level of Education",
      Previous_Regimen ~ "Previous ART Regimen"
    ),
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no"
  ) %>%
  add_p() %>%
  modify_caption(
    "**Table 1A. Baseline Sociodemographic and Clinical Characteristics by ART Regimen**"
  ) %>%
  bold_labels()

print(table1A)

# 🔥 STEP 2: TABLE 1B (Metabolic & Cardiovascular Risk)

table1B <- baseline_data %>%
  select(
    Regimen_Type,
    BloodSugar_high,
    Triglycerides_high,
    HDL_low,
    Hypertension,
    Abdominal_obesity
  ) %>%
  tbl_summary(
    by = Regimen_Type,
    label = list(
      BloodSugar_high ~ "Elevated fasting glucose (≥5.6 mmol/L)",
      Triglycerides_high ~ "Elevated Triglycerides (≥1.69 mmol/L)",
      HDL_low ~ "Low HDL Cholesterol",
      Hypertension ~ "Hypertension (SBP ≥130 or DBP ≥85 mmHg)",
      Abdominal_obesity ~ "Abdominal Obesity"
    ),
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no"
  ) %>%
  add_p(test = list(
    BloodSugar_high ~ "fisher.test",
    Triglycerides_high ~ "fisher.test",
    HDL_low ~ "fisher.test",
    Hypertension ~ "fisher.test",
    Abdominal_obesity ~ "fisher.test"
  )) %>%
  modify_caption("**Table 1B. Baseline Metabolic and Cardiovascular Risk Factors by ART Regimen**") %>%
  bold_labels()

# print Table1B
print(table1B)


# =============================================================================
# 7. Non-Linear Time Trend Analysis
# =============================================================================
# Natural cubic spline and piecewise GEE analysis of MetS trajectories
# =============================================================================

# -----------------------------------------------------------------------------
# Prepare spline analysis dataset
# -----------------------------------------------------------------------------

dat_spline <- metabolic_data_clean %>%
  mutate(
    time_weeks = case_when(
      clean_event == "Baseline" ~ 0,
      clean_event == "Week 24" ~ 24,
      clean_event == "Week 48" ~ 48,
      clean_event == "Week 72" ~ 72,
      clean_event == "Week 96" ~ 96,
      clean_event == "Week 120" ~ 120,
      clean_event == "Week 144" ~ 144,
      TRUE ~ NA_real_
    ),
    
    Regimen_Type = fct_drop(
      factor(Regimen_Type)
    )
  ) %>%
  filter(
    clean_event != "Baseline",
    !is.na(metabolic_num),
    !is.na(time_weeks)
  )


# -----------------------------------------------------------------------------
# Set TAFED as reference group
# -----------------------------------------------------------------------------

dat_spline$Regimen_Type <-
  relevel(
    dat_spline$Regimen_Type,
    ref = "TAFED"
  )


# =============================================================================
# 7.1 Linear Time Trend GEE
# =============================================================================

gee_linear <- geeglm(
  metabolic_num ~
    Regimen_Type * time_weeks +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = dat_spline,
  family = binomial,
  corstr = "exchangeable"
)


# =============================================================================
# 7.2 Natural Cubic Spline GEE
# =============================================================================

gee_spline <- geeglm(
  metabolic_num ~
    Regimen_Type * ns(time_weeks, df = 3) +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = dat_spline,
  family = binomial,
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# Test for non-linearity
# -----------------------------------------------------------------------------

nonlinear_p <- anova(
  gee_linear,
  gee_spline,
  test = "Wald"
)

cat("\n========== TEST FOR NON-LINEAR TIME TREND ==========\n")
print(nonlinear_p)


# =============================================================================
# 7.3 Piecewise GEE Analysis
# =============================================================================
# Compare the time trend during:
#   Early period: 0-48 weeks
#   Late period: 48-144 weeks
# =============================================================================

dat_spline <- dat_spline %>%
  mutate(
    period = if_else(
      time_weeks <= 48,
      "Early (0-48 weeks)",
      "Late (48-144 weeks)"
    ),
    
    time_period = if_else(
      time_weeks <= 48,
      time_weeks,
      time_weeks - 48
    )
  )


# -----------------------------------------------------------------------------
# Early period: 0-48 weeks
# -----------------------------------------------------------------------------

slope_early <- geeglm(
  metabolic_num ~
    Regimen_Type * time_period +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = dat_spline %>%
    filter(period == "Early (0-48 weeks)"),
  family = binomial,
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# Late period: 48-144 weeks
# -----------------------------------------------------------------------------

slope_late <- geeglm(
  metabolic_num ~
    Regimen_Type * time_period +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = dat_spline %>%
    filter(period == "Late (48-144 weeks)"),
  family = binomial,
  corstr = "exchangeable"
)


# =============================================================================
# 7.4 Extract TAFED-Specific Piecewise Slopes
# =============================================================================

early_slope <- coef(slope_early)["time_period"]
late_slope <- coef(slope_late)["time_period"]

early_p <- summary(slope_early)$coefficients[
  "time_period",
  "Pr(>|W|)"
]

late_p <- summary(slope_late)$coefficients[
  "time_period",
  "Pr(>|W|)"
]


cat("\n========== PIECEWISE SLOPES FOR TAFED GROUP ==========\n")

cat(
  "Early period (0-48 weeks): ",
  round(exp(early_slope), 3),
  " OR per week, p = ",
  round(early_p, 4),
  "\n",
  sep = ""
)

cat(
  "Late period (48-144 weeks): ",
  round(exp(late_slope), 3),
  " OR per week, p = ",
  round(late_p, 4),
  "\n",
  sep = ""
)


# =============================================================================
# 7.5 Save Spline Analysis Results
# =============================================================================

spline_results <- list(
  linear_model = gee_linear,
  spline_model = gee_spline,
  non_linearity_test = nonlinear_p,
  early_slope = early_slope,
  late_slope = late_slope
)

saveRDS(
  spline_results,
  "spline_analysis_results.rds"
)

# =============================================================================
# Supplementary Table S1. Piecewise GEE Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# Create piecewise time variables
# -----------------------------------------------------------------------------

dat_pw <- dat_spline %>%
  mutate(
    early_time = pmin(time_weeks, 48),
    late_time = pmax(time_weeks - 48, 0)
  )


# -----------------------------------------------------------------------------
# Set PI-based ART as the reference group
# -----------------------------------------------------------------------------

dat_pw$Regimen_Type <- relevel(
  factor(dat_pw$Regimen_Type),
  ref = "PI_based"
)


# -----------------------------------------------------------------------------
# Fit piecewise GEE model
# -----------------------------------------------------------------------------

pw_model <- geeglm(
  metabolic_num ~
    Regimen_Type * (early_time + late_time) +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = dat_pw,
  family = binomial,
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# Estimate time slopes by ART regimen
# -----------------------------------------------------------------------------

early_tab <- emtrends(
  pw_model,
  specs = "Regimen_Type",
  var = "early_time"
) %>%
  summary(infer = TRUE) %>%
  mutate(
    Period = "Weeks 0–48",
    OR = exp(early_time.trend),
    LCL = exp(lower.CL),
    UCL = exp(upper.CL)
  )


late_tab <- emtrends(
  pw_model,
  specs = "Regimen_Type",
  var = "late_time"
) %>%
  summary(infer = TRUE) %>%
  mutate(
    Period = "Weeks 48–144",
    OR = exp(late_time.trend),
    LCL = exp(lower.CL),
    UCL = exp(upper.CL)
  )


# -----------------------------------------------------------------------------
# Combine results
# -----------------------------------------------------------------------------

supp_table <- bind_rows(
  early_tab %>%
    select(Regimen_Type, Period, OR, LCL, UCL, p.value),
  
  late_tab %>%
    select(Regimen_Type, Period, OR, LCL, UCL, p.value)
)


# -----------------------------------------------------------------------------
# Format publication table
# -----------------------------------------------------------------------------

supp_table2 <- supp_table %>%
  mutate(
    `Adjusted OR (95% CI)` =
      sprintf(
        "%.2f (%.2f–%.2f)",
        OR,
        LCL,
        UCL
      ),
    
    `p-value` = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  select(
    Regimen = Regimen_Type,
    `Follow-up period` = Period,
    `Adjusted OR (95% CI)`,
    `p-value`
  )


# -----------------------------------------------------------------------------
# Export to Word
# -----------------------------------------------------------------------------

ft <- flextable(supp_table2) %>%
  autofit()

doc <- read_docx() %>%
  body_add_par(
    "Supplementary Table S1. Piecewise GEE analysis of temporal changes in metabolic syndrome risk by ART regimen",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft) %>%
  body_add_par(
    "Abbreviations: OR, odds ratio; CI, confidence interval. Models were adjusted for sex, age, and baseline CD4 count. ORs represent the relative change in odds of metabolic syndrome per week within each follow-up period.",
    style = "Normal"
  )

print(
  doc,
  target = "Supplementary_Table_S1.docx"
)

#-----------------------------------------------------------
# Figure 1. Non-linear trajectories of metabolic syndrome
#-----------------------------------------------------------

# Generate model-based predictions for all ART regimens
pred_all <- ggpredict(
  gee_spline,
  terms = c("time_weeks [0:144 by=1]", "Regimen_Type")
)

# Create publication-quality plot
figure1 <- ggplot(
  pred_all,
  aes(x = x, y = predicted, color = group, fill = group)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.15,
    linetype = 0
  ) +
  geom_line(linewidth = 1.2) +
  geom_vline(
    xintercept = 48,
    linetype = "dashed",
    color = "gray40",
    alpha = 0.7
  ) +
  annotate(
    "text",
    x = 24,
    y = 0.25,
    label = "Early phase\n(0–48 weeks)",
    size = 3
  ) +
  annotate(
    "text",
    x = 100,
    y = 0.25,
    label = "Late phase\n(48–144 weeks)",
    size = 3
  ) +
  labs(
    title = "Non-linear Trajectories of Metabolic Syndrome by ART Regimen",
    x = "Weeks on Treatment",
    y = "Predicted Probability of Metabolic Syndrome",
    color = "ART Regimen",
    fill = "ART Regimen"
  ) +
  scale_color_manual(
    values = c(
      "TLD" = "#2E86C1",
      "TAFED" = "#E67E22",
      "PI_based" = "#95A5A6"
    )
  ) +
  scale_fill_manual(
    values = c(
      "TLD" = "#2E86C1",
      "TAFED" = "#E67E22",
      "PI_based" = "#95A5A6"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

# Save Figure 1
ggsave(
  "Figure1_NonLinear_Metabolic_Syndrome.png",
  figure1,
  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# Figure 2. Metabolic Syndrome Prevalence Over Time
# =============================================================================

# ---- 1) Configure participant ID column ----
id_var <- "ID"

if (!id_var %in% names(metabolic_data_clean)) {
  stop(
    "The participant ID column '", id_var,
    "' was not found. Set id_var to your ID column name."
  )
}


# ---- 2) Define time ordering and numeric mapping ----
level_order <- c(
  "Baseline",
  "Week 24",
  "Week 48",
  "Week 72",
  "Week 96",
  "Week 120",
  "Week 144"
)

time_map <- c(
  Baseline = 0,
  "Week 24" = 24,
  "Week 48" = 48,
  "Week 72" = 72,
  "Week 96" = 96,
  "Week 120" = 120,
  "Week 144" = 144
)


# ---- 3) Prepare analysis dataset ----
dat <- metabolic_data_clean %>%
  mutate(
    clean_event = factor(clean_event, levels = level_order),
    time_weeks = unname(time_map[as.character(clean_event)]),
    Regimen_Type = fct_drop(as.factor(Regimen_Type))
  ) %>%
  filter(
    !is.na(metabolic_num),
    !is.na(time_weeks),
    !is.na(Regimen_Type)
  )


# ---- 4) Fit GEE model ----
gee_fit <- geeglm(
  metabolic_num ~ time_weeks * Regimen_Type,
  family = binomial(link = "logit"),
  id = dat[[id_var]],
  corstr = "exchangeable",
  data = dat
)


# ---- 5) Estimate temporal trends by ART regimen ----
trend_tests <- emtrends(
  gee_fit,
  ~ Regimen_Type,
  var = "time_weeks"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  transmute(
    Regimen_Type = as.character(Regimen_Type),
    slope_logit = time_weeks.trend,
    se_slope = SE,
    OR_per_week = exp(slope_logit),
    CI_lower = exp(slope_logit - 1.96 * se_slope),
    CI_upper = exp(slope_logit + 1.96 * se_slope),
    p_value = p.value,
    p_value_fmt = ifelse(
      p.value < 0.001,
      "<0.001",
      sprintf("%.3f", p.value)
    )
  )


# ---- 6) Overall temporal trend ----
overall_time_test <- emtrends(
  gee_fit,
  ~ 1,
  var = "time_weeks"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  transmute(
    label = paste0(
      "Overall temporal trend p = ",
      ifelse(
        p.value < 0.001,
        "<0.001",
        sprintf("%.3f", p.value)
      )
    )
  ) %>%
  pull(label)


# ---- 7) Calculate observed prevalence at each timepoint ----
trend_data <- dat %>%
  group_by(Regimen_Type, clean_event) %>%
  summarise(
    proportion = mean(metabolic_num),
    .groups = "drop"
  ) %>%
  filter(!is.na(proportion))


# ---- 8) Determine positions for p-value annotations ----
endpoints <- trend_data %>%
  group_by(Regimen_Type) %>%
  slice_max(
    order_by = as.numeric(clean_event),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

annot_df <- endpoints %>%
  left_join(trend_tests, by = "Regimen_Type") %>%
  arrange(Regimen_Type) %>%
  mutate(
    p_label = paste0("p = ", p_value_fmt),
    y_pos = proportion + seq(
      0.12,
      0.04,
      length.out = n()
    ),
    x_pos = as.numeric(clean_event) + 0.3
  )


# ---- 9) Create Figure 2 ----
figure2 <- ggplot(
  trend_data,
  aes(
    x = clean_event,
    y = proportion,
    color = Regimen_Type,
    group = Regimen_Type
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(
    data = annot_df,
    aes(
      x = x_pos,
      y = y_pos,
      label = p_label,
      color = Regimen_Type
    ),
    fontface = "bold",
    size = 3.2,
    show.legend = FALSE
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  scale_color_manual(
    values = c(
      "TLD" = "#2E86C1",
      "TAFED" = "#E67E22",
      "PI_based" = "#95A5A6"
    )
  ) +
  labs(
    title = "Figure 2. Metabolic Syndrome Prevalence Over Time",
    subtitle = paste0(
      "Stratified by ART Regimen. ",
      overall_time_test
    ),
    x = "Timepoint",
    y = "Proportion with Metabolic Syndrome",
    color = "ART Regimen",
    caption = paste(
      "GEE with exchangeable correlation. ",
      "P-values represent tests of the temporal trend within each ART regimen."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(
      color = "gray80",
      linetype = "dashed"
    ),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(
      hjust = 0,
      size = 10,
      face = "italic"
    ),
    legend.position = "bottom",
    plot.margin = margin(
      t = 20,
      r = 50,
      b = 30,
      l = 20
    )
  )


# ---- 10) Save Figure 2 ----
ggsave(
  "Figure2_Metabolic_Syndrome_Prevalence.png",
  figure2,
  width = 8,
  height = 6,
  dpi = 300
)



# ============================================================================
# COMBINED FIGURE 2: Metabolic Syndrome Outcomes Over Time
# ============================================================================


pred_all <- ggpredict(gee_spline, terms = c("time_weeks [0:144 by=1]", "Regimen_Type"))

plot_spline <- ggplot(pred_all, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, linetype = 0) +
  geom_vline(xintercept = 48, linetype = "dashed", color = "gray40", alpha = 0.6) +
  annotate("text", x = 24, y = 0.25, label = "Early phase\n(0-48 weeks)", size = 3) +
  annotate("text", x = 100, y = 0.25, label = "Late phase\n(48-144 weeks)", size = 3) +
  labs(x = "Weeks on Treatment",
       y = "Predicted Probability of Metabolic Syndrome",
       color = "Regimen", fill = "Regimen",
       title = "Non-linear trajectories of metabolic syndrome risk by regimen") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("TLD" = "#2E86C1", "TAFED" = "#E67E22", "PI_based" = "#95A5A6")) +
  scale_fill_manual(values = c("TLD" = "#2E86C1", "TAFED" = "#E67E22", "PI_based" = "#95A5A6")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))


# ---- 3) Panel B: Observed incidence with p-values ----
plot_incidence <- ggplot(trend_data, aes(x = clean_event, y = proportion, color = Regimen_Type, group = Regimen_Type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(
    data = annot_df,
    aes(x = x_pos, y = y_pos, label = p_label, color = Regimen_Type),
    fontface = "bold", size = 3.2, show.legend = FALSE
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Timepoint",
    y = "Proportion with Metabolic Syndrome",
    title = "B. Observed incidence with p-values",
    color = "Regimen"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 4) Combine panels vertically ----
combined_plot <- plot_spline / plot_incidence +
  plot_annotation(
    title = "Figure 2: Predicted and observed incidence of metabolic syndrome over 144 weeks by ART regimen",
    caption = "Panel A shows predicted probabilities from GEE spline models (Wald test p = 0.0012).\nPanel B shows observed incidence with regimen-specific p-values. Shaded areas = 95% CI."
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)  # <-- makes title bold
  ) &
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))


# ---- 5) Save publication-quality figure ----
ggsave("Figure2_Combined.png", combined_plot, width = 7, height = 9, dpi = 300)
ggsave("Figure2_Combined.tiff", combined_plot, width = 7, height = 9, dpi = 300, compression = "lzw")



# =============================================================================
# Component-Specific Analysis
# Section 4: Drivers of MetS Increase
# =============================================================================

# -----------------------------------------------------------------------------
# Step 1: Create standardized covariates
# -----------------------------------------------------------------------------

data_full <- data_full %>%
  mutate(
    Age_scaled = as.numeric(scale(Age..yrs.)),
    CD4_scaled = as.numeric(scale(`CD4_count.cells.µl.`))
  )


# -----------------------------------------------------------------------------
# Step 2: Prepare component-specific dataset
# -----------------------------------------------------------------------------

components_long <- data_full %>%
  select(
    Trial_number, Event_Name, Sex, Age_scaled, CD4_scaled,
    Waist_Circumference.cm., Bp_Systolic..mmHg.,
    Bp_Diastolic..mmHg., Blood_Sugar..mmol.L.,
    Triglycerides..mmol.L., Cholesterol_HDL_.mmol.L.,
    Drug_Code.
  ) %>%
  mutate(
    Event_Name_clean = case_when(
      Event_Name == "Week 24"  ~ "week_24",
      Event_Name == "Week 48"  ~ "week_48",
      Event_Name == "Week 72"  ~ "week_72",
      Event_Name == "Week 96"  ~ "week_96",
      Event_Name == "Week 120" ~ "week_120",
      Event_Name == "Week 144" ~ "week_144",
      TRUE ~ tolower(Event_Name)
    ),
    
    time_weeks = case_when(
      Event_Name_clean == "week_24"  ~ 24,
      Event_Name_clean == "week_48"  ~ 48,
      Event_Name_clean == "week_72"  ~ 72,
      Event_Name_clean == "week_96"  ~ 96,
      Event_Name_clean == "week_120" ~ 120,
      Event_Name_clean == "week_144" ~ 144,
      TRUE ~ NA_real_
    ),
    
    Regimen_Type = case_when(
      Drug_Code. %in% c("AZT+3TC+LPVr", "AZT+3TC+ATVr") ~ "PI_based",
      Drug_Code. == "TAFED" ~ "TAFED",
      Drug_Code. == "TLD"   ~ "TLD",
      TRUE ~ Drug_Code.
    ),
    
    WC_abnormal = case_when(
      Sex == "Male"   ~ Waist_Circumference.cm. >= 94,
      Sex == "Female" ~ Waist_Circumference.cm. >= 80,
      TRUE ~ NA
    ),
    
    BP_abnormal =
      Bp_Systolic..mmHg. >= 130 |
      Bp_Diastolic..mmHg. >= 85,
    
    Glucose_abnormal =
      Blood_Sugar..mmol.L. >= 5.6,
    
    TG_abnormal =
      Triglycerides..mmol.L. >= 1.69,
    
    HDL_abnormal = case_when(
      Sex == "Male"   ~ Cholesterol_HDL_.mmol.L. < 1.0,
      Sex == "Female" ~ Cholesterol_HDL_.mmol.L. < 1.3,
      TRUE ~ NA
    )
  ) %>%
  filter(
    !is.na(time_weeks),
    !is.na(Regimen_Type)
  )


# -----------------------------------------------------------------------------
# Step 3: Convert components to long format
# -----------------------------------------------------------------------------

components_long_form <- components_long %>%
  pivot_longer(
    cols = c(
      WC_abnormal,
      BP_abnormal,
      Glucose_abnormal,
      TG_abnormal,
      HDL_abnormal
    ),
    names_to = "Component",
    values_to = "Abnormal"
  ) %>%
  mutate(
    Component = recode(
      Component,
      WC_abnormal = "Waist Circumference",
      BP_abnormal = "Blood Pressure",
      Glucose_abnormal = "Glucose",
      TG_abnormal = "Triglycerides",
      HDL_abnormal = "HDL Cholesterol"
    )
  )


# -----------------------------------------------------------------------------
# Step 4: Run adjusted GEE models
# -----------------------------------------------------------------------------

components_to_analyze <- c(
  "Waist Circumference",
  "Blood Pressure",
  "Glucose",
  "Triglycerides",
  "HDL Cholesterol"
)

component_results_list <- list()

for (comp in components_to_analyze) {
  
  comp_data <- components_long_form %>%
    filter(
      Component == comp,
      !is.na(Abnormal)
    )
  
  model <- geeglm(
    Abnormal ~ time_weeks + Sex + Age_scaled + CD4_scaled,
    id = Trial_number,
    data = comp_data,
    family = binomial(link = "logit"),
    corstr = "exchangeable"
  )
  
  time_row <- broom::tidy(model) %>%
    filter(term == "time_weeks")
  
  component_results_list[[comp]] <- data.frame(
    Component = comp,
    OR_per_week = exp(time_row$estimate),
    CI_lower = exp(
      time_row$estimate - 1.96 * time_row$std.error
    ),
    CI_upper = exp(
      time_row$estimate + 1.96 * time_row$std.error
    ),
    p_value = time_row$p.value,
    Significant = ifelse(
      time_row$p.value < 0.05,
      "Yes",
      "No"
    )
  )
}


# -----------------------------------------------------------------------------
# Step 5: Format component-specific results
# -----------------------------------------------------------------------------

component_results <- bind_rows(component_results_list) %>%
  mutate(
    OR_CI = sprintf(
      "%.2f (%.2f–%.2f)",
      OR_per_week,
      CI_lower,
      CI_upper
    ),
    p_value_fmt = ifelse(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    )
  )

supplementary_table2 <- component_results %>%
  select(
    `Metabolic syndrome component` = Component,
    `Adjusted OR per week (95% CI)` = OR_CI,
    `p-value` = p_value_fmt,
    `Significant temporal trend` = Significant
  )


# -----------------------------------------------------------------------------
# Step 6: Export Supplementary Table S2
# -----------------------------------------------------------------------------

ft <- flextable(supplementary_table2) %>%
  autofit()

ft <- set_caption(
  ft,
  "Supplementary Table S2. Component-specific GEE analyses of metabolic syndrome components over follow-up."
)

doc <- read_docx() %>%
  body_add_flextable(value = ft) %>%
  body_add_par(
    paste(
      "Models were adjusted for age, sex, and baseline CD4 count.",
      "Odds ratios represent the change in odds of each metabolic",
      "abnormality per week of follow-up."
    ),
    style = "Normal"
  )

print(
  doc,
  target = "Supplementary_Table_S2.docx"
)


# -----------------------------------------------------------------------------
# Step 7: Calculate component prevalence over time
# -----------------------------------------------------------------------------

component_trends <- components_long_form %>%
  group_by(
    Component,
    Regimen_Type,
    time_weeks
  ) %>%
  summarise(
    proportion = mean(Abnormal, na.rm = TRUE),
    n = sum(!is.na(Abnormal)),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# Step 8: Plot component trajectories
# -----------------------------------------------------------------------------

component_plot <- ggplot(
  component_trends,
  aes(
    x = time_weeks,
    y = proportion,
    color = Regimen_Type,
    group = Regimen_Type
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  facet_wrap(
    ~ Component,
    scales = "free_y",
    ncol = 3
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_continuous(
    breaks = c(24, 48, 72, 96, 120, 144)
  ) +
  scale_color_manual(
    values = c(
      "TLD" = "#2E86C1",
      "TAFED" = "#E67E22",
      "PI_based" = "#95A5A6"
    )
  ) +
  labs(
    title = "Figure S2. Trajectories of Individual Metabolic Syndrome Components",
    subtitle = "By ART Regimen Over 144 Weeks",
    x = "Week",
    y = "Proportion with Abnormal Component",
    color = "ART Regimen",
    caption = paste(
      "Waist circumference cut-offs: 94 cm for males and 80 cm for females;",
      "glucose ≥5.6 mmol/L; triglycerides ≥1.69 mmol/L;",
      "HDL <1.0 mmol/L for males and <1.3 mmol/L for females;",
      "blood pressure ≥130/85 mmHg."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(component_plot)


# -----------------------------------------------------------------------------
# Step 9: Save Figure S2
# -----------------------------------------------------------------------------

ggsave(
  "Figure_S2_Component_Trajectories.png",
  component_plot,
  width = 12,
  height = 8,
  dpi = 600
)


# -----------------------------------------------------------------------------
# Step 10: Save component-specific results
# -----------------------------------------------------------------------------

write.csv(
  component_results,
  "Component_Analysis_Results.csv",
  row.names = FALSE
)


# =============================================================================
# 8. Viral Load Follow-Up Analysis
# =============================================================================
# Purpose:
#   To describe viral suppression over follow-up and assess the extent of
#   missing viral load data by ART regimen and timepoint.
# =============================================================================

-----------------------------------------------------------------------------
  # 8.1 Create viral suppression dataset
  # -----------------------------------------------------------------------------

vl_followup <- data_full %>%
  select(
    Trial_number,
    Event_Name,
    Viral_Load.cp.ml.,
    Drug_Code.
  ) %>%
  mutate(
    # Standardize follow-up timepoints
    Event_Name_clean = case_when(
      grepl("baseline", Event_Name, ignore.case = TRUE) ~ "Baseline",
      grepl("24", Event_Name) ~ "Week 24",
      grepl("48", Event_Name) ~ "Week 48",
      grepl("64", Event_Name) ~ "Week 72",
      grepl("72", Event_Name) ~ "Week 72",
      grepl("96", Event_Name) ~ "Week 96",
      grepl("112", Event_Name) ~ "Week 120",
      grepl("120", Event_Name) ~ "Week 120",
      grepl("144", Event_Name) ~ "Week 144",
      TRUE ~ NA_character_
    ),
    
    # Define viral suppression
    viral_suppressed = case_when(
      is.na(Viral_Load.cp.ml.) ~ NA_character_,
      Viral_Load.cp.ml. < 1000 ~ "Suppressed (<1000)",
      Viral_Load.cp.ml. >= 1000 ~ "Not Suppressed (≥1000)"
    ),
    
    # Standardize ART regimen
    Regimen_Type = case_when(
      Drug_Code. %in% c("AZT+3TC+LPVr", "AZT+3TC+ATVr") ~ "PI_based",
      Drug_Code. == "TAFED" ~ "TAFED",
      Drug_Code. == "TLD" ~ "TLD",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Event_Name_clean),
    !is.na(Regimen_Type)
  ) %>%
  mutate(
    Event_Name_clean = factor(
      Event_Name_clean,
      levels = c(
        "Baseline",
        "Week 24",
        "Week 48",
        "Week 72",
        "Week 96",
        "Week 120",
        "Week 144"
      )
    ),
    Regimen_Type = factor(
      Regimen_Type,
      levels = c("PI_based", "TAFED", "TLD")
    )
  )


# -----------------------------------------------------------------------------
# 8.2 Calculate viral suppression rates
# -----------------------------------------------------------------------------

vl_summary_table <- vl_followup %>%
  filter(!is.na(viral_suppressed)) %>%
  group_by(
    Regimen_Type,
    Event_Name_clean,
    viral_suppressed
  ) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(
    Regimen_Type,
    Event_Name_clean
  ) %>%
  mutate(
    total = sum(n),
    pct = 100 * n / total
  ) %>%
  filter(
    viral_suppressed == "Suppressed (<1000)"
  ) %>%
  transmute(
    Regimen_Type,
    Timepoint = Event_Name_clean,
    Suppression_Rate = round(pct, 1),
    N = total
  ) %>%
  arrange(
    Regimen_Type,
    Timepoint
  )


# Display results
cat("\n========== VIRAL SUPPRESSION RATES OVER TIME ==========\n")
print(vl_summary_table, n = 30)


# -----------------------------------------------------------------------------
# 8.3 Plot viral suppression over time
# -----------------------------------------------------------------------------

vl_plot <- ggplot(
  vl_followup %>%
    filter(!is.na(viral_suppressed)),
  aes(
    x = Event_Name_clean,
    fill = viral_suppressed
  )
) +
  geom_bar(
    position = "fill",
    width = 0.7
  ) +
  facet_wrap(
    ~ Regimen_Type,
    ncol = 3
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_fill_manual(
    values = c(
      "Suppressed (<1000)" = "#2E8B57",
      "Not Suppressed (≥1000)" = "#CD5C5C"
    )
  ) +
  labs(
    title = "Supplementary Figure S3. Viral Suppression Over Time",
    subtitle = "By ART Regimen (VL <1000 copies/mL = suppressed)",
    x = "Timepoint",
    y = "Proportion of Participants",
    fill = "Viral Load Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold"
    )
  )


# Display plot
print(vl_plot)


# -----------------------------------------------------------------------------
# 8.4 Save viral suppression plot
# -----------------------------------------------------------------------------

ggsave(
  "Figure_S3_VL_Suppression_OverTime.png",
  vl_plot,
  width = 10,
  height = 6,
  dpi = 600
)


# -----------------------------------------------------------------------------
# 8.5 Save viral suppression summary
# -----------------------------------------------------------------------------

write.csv(
  vl_summary_table,
  "Viral_Suppression_Summary.csv",
  row.names = FALSE
)

# =============================================================================
# 9. GEE Model: Main Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# 9.1 Prepare GEE analysis dataset
# -----------------------------------------------------------------------------

gee_data <- metabolic_data_clean %>%
  filter(
    clean_event != "Baseline"
  ) %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Alcohol_Consumption,
    Age_scaled,
    CD4_scaled
  ) %>%
  na.omit() %>%
  mutate(
    Regimen_Type = fct_relevel(
      Regimen_Type,
      "PI_based"
    ),
    Sex = fct_relevel(
      Sex,
      "male"
    ),
    clean_event = fct_relevel(
      droplevels(clean_event),
      "Week 24"
    ),
    across(
      where(is.factor),
      fct_drop
    )
  )


# -----------------------------------------------------------------------------
# 9.2 Fit main GEE model
# -----------------------------------------------------------------------------

gee_model <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = gee_data,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)

print(summary(gee_model))


# -----------------------------------------------------------------------------
# 9.3 Create publication-ready results table
# -----------------------------------------------------------------------------

coef_table <- broom::tidy(
  gee_model,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  filter(
    term != "(Intercept)"
  ) %>%
  mutate(
    OR = round(estimate, 2),
    `95% CI` = paste0(
      "(",
      round(conf.low, 2),
      ", ",
      round(conf.high, 2),
      ")"
    ),
    `p-value` = ifelse(
      p.value < 0.001,
      "<0.001",
      sprintf("%.3f", p.value)
    )
  ) %>%
  select(
    term,
    OR,
    `95% CI`,
    `p-value`
  )


# Display results
print(
  knitr::kable(
    coef_table,
    caption = "Adjusted GEE Odds Ratios for Metabolic Syndrome"
  )
)

# =============================================================================
# 9.2 Sensitivity Analyses: Previous ART Regimen and Baseline Viral Load
# =============================================================================
#
# Purpose:
#   To assess whether previous ART regimen and baseline viral load
#   influence the association between randomized ART regimen
#   and metabolic syndrome.
#
# Analyses:
#   Model 1: Primary model + previous ART regimen
#   Model 2: Primary model + baseline viral load
#   Model 3: Primary model + previous ART regimen + baseline viral load
#   Model 4: Participants with baseline VL <1000 copies/mL
#   Model 5: Participants with baseline VL >=1000 copies/mL
# =============================================================================


# -----------------------------------------------------------------------------
# 9.2.1 Prepare baseline viral load information
# -----------------------------------------------------------------------------

analysis_baseline_vl <- metabolic_data_clean %>%
  filter(Event_Name == "baseline") %>%
  select(
    ID,
    Baseline_VL = Viral_Load.cp.ml.
  ) %>%
  mutate(
    ID = as.character(ID),
    Baseline_VL_log10 = log10(Baseline_VL + 1),
    Baseline_VL_group = case_when(
      is.na(Baseline_VL) ~ NA_character_,
      Baseline_VL < 1000 ~ "<1000 copies/mL",
      Baseline_VL >= 1000 ~ ">=1000 copies/mL"
    )
  )


# -----------------------------------------------------------------------------
# 9.2.2 Create longitudinal dataset for sensitivity analyses
# -----------------------------------------------------------------------------

analysis_data <- metabolic_data_clean %>%
  filter(Event_Name != "baseline") %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Age_scaled,
    CD4_scaled,
    Previous_Regimen
  ) %>%
  mutate(
    ID = as.character(ID),
    Regimen_Type = fct_relevel(
      factor(Regimen_Type),
      "PI_based"
    ),
    Sex = fct_relevel(
      factor(Sex),
      "male"
    ),
    clean_event = fct_relevel(
      droplevels(factor(clean_event)),
      "Week 24"
    ),
    Previous_Regimen = factor(Previous_Regimen)
  ) %>%
  left_join(
    analysis_baseline_vl,
    by = "ID"
  )


# -----------------------------------------------------------------------------
# 9.2.3 Check previous ART regimen and baseline viral load
# -----------------------------------------------------------------------------

cat("\n========== PREVIOUS ART REGIMEN ==========\n")
print(
  table(
    analysis_data$Previous_Regimen,
    useNA = "ifany"
  )
)

cat("\n========== BASELINE VIRAL LOAD ==========\n")
print(
  table(
    analysis_data$Baseline_VL_group,
    useNA = "ifany"
  )
)


# =============================================================================
# 9.2.4 Model 1: Primary model + Previous ART regimen
# =============================================================================

analysis_prevART <- analysis_data %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Age_scaled,
    CD4_scaled,
    Previous_Regimen
  ) %>%
  na.omit() %>%
  droplevels()

analysis_prevART$Previous_Regimen <- relevel(
  analysis_prevART$Previous_Regimen,
  ref = "TLE"
)

gee_prevART <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled +
    Previous_Regimen,
  id = ID,
  data = analysis_prevART,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)

cat("\n============================================================\n")
cat("MODEL 1: PRIMARY MODEL + PREVIOUS ART REGIMEN\n")
cat("============================================================\n")

print(summary(gee_prevART))


# =============================================================================
# 9.2.5 Model 2: Primary model + baseline viral load
# =============================================================================

analysis_vl <- analysis_data %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Age_scaled,
    CD4_scaled,
    Baseline_VL_log10
  ) %>%
  na.omit() %>%
  droplevels()

gee_vl <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled +
    Baseline_VL_log10,
  id = ID,
  data = analysis_vl,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)

cat("\n============================================================\n")
cat("MODEL 2: PRIMARY MODEL + BASELINE VIRAL LOAD\n")
cat("============================================================\n")

print(summary(gee_vl))


# =============================================================================
# 9.2.6 Model 3: Primary model + Previous ART + baseline viral load
# =============================================================================

analysis_full <- analysis_data %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Age_scaled,
    CD4_scaled,
    Previous_Regimen,
    Baseline_VL_log10
  ) %>%
  na.omit() %>%
  droplevels()

analysis_full$Previous_Regimen <- relevel(
  analysis_full$Previous_Regimen,
  ref = "TLE"
)

gee_full <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled +
    Previous_Regimen +
    Baseline_VL_log10,
  id = ID,
  data = analysis_full,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)

cat("\n============================================================\n")
cat("MODEL 3: PRIMARY MODEL + PREVIOUS ART + BASELINE VL\n")
cat("============================================================\n")

print(summary(gee_full))


# =============================================================================
# 9.2.7 Model 4: Participants with baseline VL <1000 copies/mL
# =============================================================================

analysis_vl_low <- analysis_data %>%
  filter(
    !is.na(Baseline_VL),
    Baseline_VL < 1000
  ) %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Age_scaled,
    CD4_scaled
  ) %>%
  na.omit() %>%
  droplevels()

gee_vl_low <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = analysis_vl_low,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)

cat("\n============================================================\n")
cat("MODEL 4: BASELINE VL <1000 copies/mL\n")
cat("============================================================\n")

print(summary(gee_vl_low))


# =============================================================================
# 9.2.8 Model 5: Participants with baseline VL >=1000 copies/mL
# =============================================================================

analysis_vl_high <- analysis_data %>%
  filter(
    !is.na(Baseline_VL),
    Baseline_VL >= 1000
  ) %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Age_scaled,
    CD4_scaled
  ) %>%
  na.omit() %>%
  droplevels()

gee_vl_high <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = analysis_vl_high,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)

cat("\n============================================================\n")
cat("MODEL 5: BASELINE VL >=1000 copies/mL\n")
cat("============================================================\n")

print(summary(gee_vl_high))


# -----------------------------------------------------------------------------
# 9.2.9 Sample sizes for sensitivity analyses
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SAMPLE SIZES FOR SENSITIVITY ANALYSES\n")
cat("============================================================\n")

cat(
  "\nModel 1 - Previous ART:\n",
  "Observations =", nrow(analysis_prevART),
  "\nParticipants =", dplyr::n_distinct(analysis_prevART$ID),
  "\n"
)

cat(
  "\nModel 2 - Baseline VL:\n",
  "Observations =", nrow(analysis_vl),
  "\nParticipants =", dplyr::n_distinct(analysis_vl$ID),
  "\n"
)

cat(
  "\nModel 3 - Previous ART + Baseline VL:\n",
  "Observations =", nrow(analysis_full),
  "\nParticipants =", dplyr::n_distinct(analysis_full$ID),
  "\n"
)

cat(
  "\nModel 4 - Baseline VL <1000:\n",
  "Observations =", nrow(analysis_vl_low),
  "\nParticipants =", dplyr::n_distinct(analysis_vl_low$ID),
  "\n"
)

cat(
  "\nModel 5 - Baseline VL >=1000:\n",
  "Observations =", nrow(analysis_vl_high),
  "\nParticipants =", dplyr::n_distinct(analysis_vl_high$ID),
  "\n"
)

# =============================================================================
# 9.3 Extract TAFED and TLD Effects from Sensitivity Models
# =============================================================================

extract_regimen_effects <- function(model, model_name, n_participants) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(
      term %in% c(
        "Regimen_TypeTAFED",
        "Regimen_TypeTLD"
      )
    ) %>%
    mutate(
      Model = model_name,
      N = n_participants,
      
      Regimen = case_when(
        term == "Regimen_TypeTAFED" ~ "TAFED vs PI-based",
        term == "Regimen_TypeTLD" ~ "TLD vs PI-based"
      ),
      
      `OR (95% CI)` = paste0(
        sprintf("%.2f", estimate),
        " (",
        sprintf("%.2f", conf.low),
        "–",
        sprintf("%.2f", conf.high),
        ")"
      ),
      
      `p-value` = ifelse(
        p.value < 0.001,
        "<0.001",
        sprintf("%.3f", p.value)
      )
    ) %>%
    select(
      Model,
      N,
      Regimen,
      `OR (95% CI)`,
      `p-value`
    )
}


# -----------------------------------------------------------------------------
# 9.3.1 Extract regimen effects from all sensitivity models
# -----------------------------------------------------------------------------

regimen_results <- bind_rows(
  
  extract_regimen_effects(
    gee_prevART,
    "Primary + previous ART",
    dplyr::n_distinct(analysis_prevART$ID)
  ),
  
  extract_regimen_effects(
    gee_vl,
    "Primary + baseline VL",
    dplyr::n_distinct(analysis_vl$ID)
  ),
  
  extract_regimen_effects(
    gee_full,
    "Primary + previous ART + baseline VL",
    dplyr::n_distinct(analysis_full$ID)
  ),
  
  extract_regimen_effects(
    gee_vl_low,
    "Baseline VL <1000",
    dplyr::n_distinct(analysis_vl_low$ID)
  ),
  
  extract_regimen_effects(
    gee_vl_high,
    "Baseline VL >=1000",
    dplyr::n_distinct(analysis_vl_high$ID)
  )
)


# -----------------------------------------------------------------------------
# 9.3.2 Display sensitivity analysis results
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("TAFED AND TLD ESTIMATES FROM SENSITIVITY MODELS\n")
cat("============================================================\n")

print(regimen_results)


# =============================================================================
# 9.4 Export Supplementary Table S4
# =============================================================================

ft <- flextable(regimen_results) %>%
  autofit()

doc <- read_docx() %>%
  body_add_par(
    "Supplementary Table S4. Sensitivity analyses examining previous ART regimen and baseline viral load",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft) %>%
  body_add_par(
    paste(
      "All models used a binomial distribution with logit link and",
      "exchangeable working correlation structure. The primary model",
      "was adjusted for follow-up time, sex, age, and CD4 count.",
      "Sensitivity models additionally accounted for previous ART",
      "regimen and/or baseline viral load."
    ),
    style = "Normal"
  )

print(
  doc,
  target = "Supplementary_Table_S4_Sensitivity_Analyses.docx"
)

message("✔ Sensitivity analysis table completed.")


# =============================================================================
# 9.5 Stratified GEE Analysis by ART Regimen
# =============================================================================


# -----------------------------------------------------------------------------
# 9.5.1 Prepare stratified GEE dataset
# -----------------------------------------------------------------------------

gee_data <- metabolic_data_clean %>%
  filter(
    clean_event != "Baseline"
  ) %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Age_scaled
  ) %>%
  na.omit() %>%
  mutate(
    Regimen_Type = fct_relevel(
      Regimen_Type,
      "PI_based"
    ),
    Sex = fct_relevel(
      Sex,
      "male"
    ),
    clean_event = fct_relevel(
      fct_drop(clean_event),
      "Week 24"
    )
  )


# Split data by ART regimen
gee_split <- split(
  gee_data,
  gee_data$Regimen_Type
)


# -----------------------------------------------------------------------------
# 9.5.2 Function to run stratified GEE models
# -----------------------------------------------------------------------------

run_stratified_gee <- function(data, group_name) {
  
  model <- geeglm(
    metabolic_num ~
      clean_event +
      Sex +
      Age_scaled,
    id = ID,
    data = data,
    family = binomial(link = "logit"),
    corstr = "exchangeable"
  )
  
  broom::tidy(
    model,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      term != "(Intercept)"
    ) %>%
    mutate(
      OR = round(estimate, 2),
      CI = paste0(
        "(",
        round(conf.low, 2),
        ", ",
        round(conf.high, 2),
        ")"
      ),
      p = ifelse(
        p.value < 0.001,
        "<0.001",
        sprintf("%.3f", p.value)
      ),
      regimen = group_name,
      term = case_when(
        term == "clean_eventWeek 48" ~ "Week 48 vs Week 24¹",
        term == "clean_eventWeek 72" ~ "Week 72 vs Week 24¹",
        term == "clean_eventWeek 96" ~ "Week 96 vs Week 24¹",
        term == "clean_eventWeek 120" ~ "Week 120 vs Week 24¹",
        term == "clean_eventWeek 144" ~ "Week 144 vs Week 24¹",
        term == "Sexfemale" ~ "Female vs Male",
        term == "Age_scaled" ~ "Age (per 1 SD increase)",
        TRUE ~ term
      )
    ) %>%
    select(
      regimen,
      term,
      OR,
      CI,
      p
    )
}


# -----------------------------------------------------------------------------
# 9.5.3 Run stratified GEE models
# -----------------------------------------------------------------------------

results <- purrr::imap_dfr(
  gee_split,
  run_stratified_gee
)


# -----------------------------------------------------------------------------
# 9.5.4 Format results for publication
# -----------------------------------------------------------------------------

TableX <- results %>%
  pivot_wider(
    names_from = regimen,
    values_from = c(OR, CI, p),
    names_glue = "{regimen}_{.value}"
  ) %>%
  arrange(
    match(
      term,
      c(
        "Week 48 vs Week 24¹",
        "Week 72 vs Week 24¹",
        "Week 96 vs Week 24¹",
        "Week 120 vs Week 24¹",
        "Week 144 vs Week 24¹",
        "Female vs Male",
        "Age (per 1 SD increase)"
      )
    )
  )


# -----------------------------------------------------------------------------
# 9.5.5 Export stratified GEE results
# -----------------------------------------------------------------------------

doc <- read_docx() %>%
  body_add_par(
    "Table X. Stratified GEE analysis of metabolic syndrome by ART regimen",
    style = "heading 1"
  ) %>%
  body_add_flextable(
    flextable(TableX) %>%
      autofit()
  )

print(
  doc,
  target = "TableX_MetSyn_GEE_Stratified_Final.docx"
)

message(
  "✔ Stratified GEE results saved to TableX_MetSyn_GEE_Stratified_Final.docx"
)


# =============================================================================
# 9.6 Stratified GEE Analysis by Baseline Viral Load
# =============================================================================

# -----------------------------------------------------------------------------
# 9.6.1 Create viral suppression indicator for follow-up observations
# -----------------------------------------------------------------------------

gee_data_3 <- metabolic_data %>%
  mutate(
    ID = as.character(ID),
    
    ViralLoad_suppressed = case_when(
      is.na(Viral_Load.cp.ml.) ~ NA_integer_,
      Viral_Load.cp.ml. < 1000 ~ 1L,
      Viral_Load.cp.ml. >= 1000 ~ 0L
    )
  ) %>%
  select(
    ID,
    metabolic_num,
    Regimen_Type,
    clean_event,
    Sex,
    Alcohol_Consumption,
    Age_scaled,
    CD4_scaled,
    ViralLoad_suppressed
  )


# -----------------------------------------------------------------------------
# 9.6.2 Create complete-case dataset
# -----------------------------------------------------------------------------

gee_data_3_model <- gee_data_3 %>%
  filter(
    !is.na(metabolic_num),
    !is.na(Regimen_Type),
    !is.na(clean_event),
    !is.na(Sex),
    !is.na(Age_scaled),
    !is.na(CD4_scaled),
    !is.na(ViralLoad_suppressed),
    !is.na(Alcohol_Consumption)
  ) %>%
  mutate(
    Regimen_Type = droplevels(factor(Regimen_Type)),
    clean_event = droplevels(factor(clean_event)),
    Sex = droplevels(factor(Sex)),
    Alcohol_Consumption = droplevels(factor(Alcohol_Consumption))
  )


# -----------------------------------------------------------------------------
# 9.6.3 Determine baseline viral load status
# -----------------------------------------------------------------------------

baseline_vl_status <- metabolic_data_clean %>%
  filter(
    Event_Name == "baseline"
  ) %>%
  select(
    ID,
    Baseline_VL = Viral_Load.cp.ml.
  ) %>%
  mutate(
    ID = as.character(ID),
    
    Baseline_VL_group = case_when(
      is.na(Baseline_VL) ~ NA_character_,
      Baseline_VL < 1000 ~ "Suppressed at baseline",
      Baseline_VL >= 1000 ~ "Unsuppressed at baseline"
    )
  )


# -----------------------------------------------------------------------------
# 9.6.4 Merge baseline VL status with follow-up data
# -----------------------------------------------------------------------------

gee_data_stratified <- gee_data_3_model %>%
  left_join(
    baseline_vl_status,
    by = "ID"
  ) %>%
  filter(
    !is.na(Baseline_VL_group)
  ) %>%
  mutate(
    Regimen_Type = droplevels(Regimen_Type),
    clean_event = droplevels(clean_event),
    Sex = droplevels(Sex),
    Alcohol_Consumption = droplevels(Alcohol_Consumption)
  )


# -----------------------------------------------------------------------------
# 9.6.5 Fit GEE model among participants suppressed at baseline
# -----------------------------------------------------------------------------

data_vl_suppressed <- gee_data_stratified %>%
  filter(
    Baseline_VL_group == "Suppressed at baseline"
  ) %>%
  mutate(
    across(where(is.factor), droplevels)
  )

gee_suppressed <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled +
    Alcohol_Consumption +
    ViralLoad_suppressed,
  id = ID,
  data = data_vl_suppressed,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# 9.6.6 Fit GEE model among participants unsuppressed at baseline
# -----------------------------------------------------------------------------

data_vl_unsuppressed <- gee_data_stratified %>%
  filter(
    Baseline_VL_group == "Unsuppressed at baseline"
  ) %>%
  mutate(
    across(where(is.factor), droplevels)
  )

gee_unsuppressed <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled +
    Alcohol_Consumption +
    ViralLoad_suppressed,
  id = ID,
  data = data_vl_unsuppressed,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# 9.6.7 Display model summaries
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("STRATIFIED BY BASELINE VL: SUPPRESSED GROUP\n")
cat("============================================================\n")

print(summary(gee_suppressed))


cat("\n============================================================\n")
cat("STRATIFIED BY BASELINE VL: UNSUPPRESSED GROUP\n")
cat("============================================================\n")

print(summary(gee_unsuppressed))


# =============================================================================
# 9.7 Direct Comparison of TAFED vs TLD
# Addresses Reviewer 1 Comment #9: Clarify TAFED vs TLD differences
# =============================================================================

# -----------------------------------------------------------------------------
# 9.7.1 Restrict analysis to DTG-based regimens
# -----------------------------------------------------------------------------

gee_data_dtg_only <- gee_data %>%
  filter(
    Regimen_Type %in% c("TAFED", "TLD")
  ) %>%
  mutate(
    Regimen_Type = droplevels(Regimen_Type),
    Regimen_Type = forcats::fct_relevel(
      Regimen_Type,
      "TAFED"
    )
  )


# -----------------------------------------------------------------------------
# 9.7.2 Fit GEE model: TLD vs TAFED
# -----------------------------------------------------------------------------

gee_tafed_vs_tld <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = gee_data_dtg_only,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# 9.7.3 Extract direct TLD vs TAFED comparison
# -----------------------------------------------------------------------------

tafed_vs_tld_result <- broom::tidy(
  gee_tafed_vs_tld,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  filter(
    term == "Regimen_TypeTLD"
  ) %>%
  mutate(
    Comparison = "TLD vs TAFED",
    `OR (95% CI)` = paste0(
      sprintf("%.2f", estimate),
      " (",
      sprintf("%.2f", conf.low),
      "–",
      sprintf("%.2f", conf.high),
      ")"
    ),
    `p-value` = ifelse(
      p.value < 0.001,
      "<0.001",
      sprintf("%.3f", p.value)
    )
  ) %>%
  select(
    Comparison,
    `OR (95% CI)`,
    `p-value`
  )

print(tafed_vs_tld_result)


# -----------------------------------------------------------------------------
# 9.7.4 Display model summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("DIRECT COMPARISON: TLD vs TAFED\n")
cat("============================================================\n")

print(tafed_vs_tld_result)

# =============================================================================
# 9.8 Stratified GEE Analysis by ART Regimen
# =============================================================================

# -----------------------------------------------------------------------------
# 9.8.1 Split data by regimen
# -----------------------------------------------------------------------------

gee_split <- split(
  gee_data,
  gee_data$Regimen_Type
)


# -----------------------------------------------------------------------------
# 9.8.2 Function to fit regimen-specific GEE models
# -----------------------------------------------------------------------------

run_stratified_gee <- function(data, group_name) {
  
  model <- geeglm(
    metabolic_num ~
      clean_event +
      Sex +
      Age_scaled,
    id = ID,
    data = data,
    family = binomial(link = "logit"),
    corstr = "exchangeable"
  )
  
  broom::tidy(
    model,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      term != "(Intercept)"
    ) %>%
    mutate(
      regimen = group_name,
      
      term_clean = recode(
        term,
        "clean_eventWeek 48"  = "Week 48 vs Week 24¹",
        "clean_eventWeek 72"  = "Week 72 vs Week 24¹",
        "clean_eventWeek 96"  = "Week 96 vs Week 24¹",
        "clean_eventWeek 120" = "Week 120 vs Week 24¹",
        "clean_eventWeek 144" = "Week 144 vs Week 24¹",
        "Sexfemale"           = "Female vs Male",
        "Age_scaled"          = "Age (per 1 SD)",
        .default = term
      ),
      
      label = factor(
        term_clean,
        levels = c(
          "Week 48 vs Week 24¹",
          "Week 72 vs Week 24¹",
          "Week 96 vs Week 24¹",
          "Week 120 vs Week 24¹",
          "Week 144 vs Week 24¹",
          "Female vs Male",
          "Age (per 1 SD)"
        )
      ),
      
      group = case_when(
        str_detect(term_clean, "^Week") ~ "Time Point",
        term_clean == "Female vs Male" ~ "Demographics",
        term_clean == "Age (per 1 SD)" ~ "Clinical",
        TRUE ~ "Other"
      )
    )
}


# -----------------------------------------------------------------------------
# 9.8.3 Run regimen-specific GEE models
# -----------------------------------------------------------------------------

figure3 <- purrr::imap_dfr(
  gee_split,
  run_stratified_gee
)


# -----------------------------------------------------------------------------
# 9.8.4 Create stratified forest plot
# -----------------------------------------------------------------------------

forest_plot <- ggplot(
  figure3,
  aes(
    x = estimate,
    y = label
  )
) +
  geom_point(
    size = 3
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    width = 0.2,
    orientation = "y"
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  facet_grid(
    group ~ regimen,
    scales = "free_y"
  ) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = NULL,
    title = "Figure 3. Adjusted Odds Ratios for Metabolic Syndrome",
    subtitle = "Stratified by ART regimen",
    caption = paste0(
      "Odds ratios and 95% confidence intervals from regimen-specific ",
      "GEE logistic regression models. ",
      "Week comparisons are versus Week 24; ",
      "age was standardised (mean = 44.1 years, SD = 10.2 years)."
    )
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 10
    )
  )


# -----------------------------------------------------------------------------
# 9.8.5 Save Figure 3
# -----------------------------------------------------------------------------

ggsave(
  filename = "Figure3_ForestPlot_Regimen.png",
  plot = forest_plot,
  width = 11,
  height = 6,
  dpi = 600,
  units = "in"
)

# =============================================================================
# 9.9 GEE Results Table (Table 3): Crude and Adjusted Odds Ratios
# =============================================================================

# -----------------------------------------------------------------------------
# 9.9.1 Prepare analysis dataset
# -----------------------------------------------------------------------------

table3_data <- gee_data_3 %>%
  mutate(
    ID = as.character(ID),
    Regimen_Type = droplevels(factor(Regimen_Type)),
    clean_event = droplevels(factor(clean_event)),
    Sex = droplevels(factor(Sex)),
    Alcohol_Consumption = droplevels(factor(Alcohol_Consumption))
  )


# -----------------------------------------------------------------------------
# 9.9.2 Crude GEE model
# Reference regimen: PI-based
# -----------------------------------------------------------------------------

gee_table3_crude <- geeglm(
  metabolic_num ~ Regimen_Type,
  id = ID,
  data = table3_data,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# 9.9.3 Adjusted GEE model
# -----------------------------------------------------------------------------

table3_adj_data <- table3_data %>%
  filter(
    !is.na(metabolic_num),
    !is.na(Regimen_Type),
    !is.na(clean_event),
    !is.na(Sex),
    !is.na(Age_scaled),
    !is.na(CD4_scaled),
    !is.na(Alcohol_Consumption),
    !is.na(ViralLoad_suppressed)
  )

gee_table3_adj <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled +
    Alcohol_Consumption +
    ViralLoad_suppressed,
  id = ID,
  data = table3_adj_data,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# 9.9.4 Extract regimen effects
# -----------------------------------------------------------------------------

extract_regimen_effects <- function(model, model_name) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(
      term %in% c(
        "Regimen_TypeTAFED",
        "Regimen_TypeTLD"
      )
    ) %>%
    mutate(
      Model = model_name,
      
      Regimen = case_when(
        term == "Regimen_TypeTAFED" ~ "TAFED vs PI-based",
        term == "Regimen_TypeTLD"   ~ "TLD vs PI-based"
      ),
      
      `OR (95% CI)` = paste0(
        sprintf("%.2f", estimate),
        " (",
        sprintf("%.2f", conf.low),
        "–",
        sprintf("%.2f", conf.high),
        ")"
      ),
      
      `p-value` = ifelse(
        p.value < 0.001,
        "<0.001",
        sprintf("%.3f", p.value)
      )
    ) %>%
    select(
      Model,
      Regimen,
      `OR (95% CI)`,
      `p-value`
    )
}


# -----------------------------------------------------------------------------
# 9.9.5 Create Table 3
# -----------------------------------------------------------------------------

table3_results <- bind_rows(
  extract_regimen_effects(
    gee_table3_crude,
    "Crude"
  ),
  
  extract_regimen_effects(
    gee_table3_adj,
    "Adjusted"
  )
)

print(table3_results)

# -----------------------------------------------------------------------------
# 9.9.6 Add direct TAFED vs TLD comparison
# -----------------------------------------------------------------------------

table3_tafed_vs_tld <- broom::tidy(
  gee_tafed_vs_tld,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  filter(
    term == "Regimen_TypeTLD"
  ) %>%
  mutate(
    Model = "Direct comparison",
    Regimen = "TLD vs TAFED",
    
    `OR (95% CI)` = paste0(
      sprintf("%.2f", estimate),
      " (",
      sprintf("%.2f", conf.low),
      "–",
      sprintf("%.2f", conf.high),
      ")"
    ),
    
    `p-value` = ifelse(
      p.value < 0.001,
      "<0.001",
      sprintf("%.3f", p.value)
    )
  ) %>%
  select(
    Model,
    Regimen,
    `OR (95% CI)`,
    `p-value`
  )

print(table3_tafed_vs_tld)


# =============================================================================
# 9.10 Figure 4: Predicted Probability of Metabolic Syndrome
# =============================================================================

# -----------------------------------------------------------------------------
# 9.10.1 Overall predicted probabilities by regimen
# -----------------------------------------------------------------------------

emm_overall <- emmeans(
  gee_model,
  ~ Regimen_Type,
  type = "response"
)

emm_overall_df <- as.data.frame(emm_overall)

panelA <- ggplot(
  emm_overall_df,
  aes(x = Regimen_Type, y = prob)
) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.8
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "A. Overall",
    x = "ART Regimen",
    y = "Predicted Probability (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )


# -----------------------------------------------------------------------------
# 9.10.2 Predicted probabilities by sex and regimen
# -----------------------------------------------------------------------------

emm_sex <- emmeans(
  gee_model,
  ~ Sex * Regimen_Type,
  type = "response"
)

emm_sex_df <- as.data.frame(emm_sex)

panelB <- ggplot(
  emm_sex_df,
  aes(
    x = Regimen_Type,
    y = prob,
    color = Sex,
    group = Sex
  )
) +
  geom_point(size = 4) +
  geom_line(linewidth = 1.2) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.8
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "B. By Sex",
    x = "ART Regimen",
    y = "Predicted Probability (%)",
    color = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )


# -----------------------------------------------------------------------------
# 9.10.3 Combine panels
# -----------------------------------------------------------------------------

figure4_combined <- panelA + panelB +
  patchwork::plot_annotation(
    title = "Figure 4. Predicted Probability of Metabolic Syndrome by ART Regimen",
    subtitle = "Estimated Marginal Means from GEE Model",
    caption = paste0(
      "Models adjusted for follow-up time, age, and CD4 count. ",
      "Confidence intervals represent 95% Wald estimates."
    )
  )


# -----------------------------------------------------------------------------
# 9.10.4 Save Figure 4
# -----------------------------------------------------------------------------

ggsave(
  filename = "Figure4_PredictedProbability_MetS.png",
  plot = figure4_combined,
  width = 12,
  height = 6,
  dpi = 600,
  units = "in"
)

print(figure4_combined)

# =============================================================================
# 9.11 Table 4: Sex-Stratified GEE Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# 9.11.1 Female model
# -----------------------------------------------------------------------------

gee_female <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = gee_data %>%
    filter(Sex == "female") %>%
    droplevels(),
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# 9.11.2 Male model
# -----------------------------------------------------------------------------

gee_male <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = gee_data %>%
    filter(Sex == "male") %>%
    droplevels(),
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# 9.11.3 Extract results
# -----------------------------------------------------------------------------

female_results <- broom::tidy(
  gee_female,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  mutate(Sex = "Female")

male_results <- broom::tidy(
  gee_male,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  mutate(Sex = "Male")


# -----------------------------------------------------------------------------
# 9.11.4 Format Table 4
# -----------------------------------------------------------------------------

table4_results <- bind_rows(
  female_results,
  male_results
) %>%
  filter(
    term %in% c(
      "Regimen_TypeTAFED",
      "Regimen_TypeTLD",
      "Age_scaled",
      "CD4_scaled"
    )
  ) %>%
  mutate(
    Predictor = case_when(
      term == "Regimen_TypeTAFED" ~ "TAFED vs PI-based",
      term == "Regimen_TypeTLD"   ~ "TLD vs PI-based",
      term == "Age_scaled"        ~ "Age (per 1 SD)",
      term == "CD4_scaled"        ~ "CD4 (per 1 SD)"
    ),
    
    `OR (95% CI)` = paste0(
      sprintf("%.2f", estimate),
      " (",
      sprintf("%.2f", conf.low),
      "–",
      sprintf("%.2f", conf.high),
      ")"
    ),
    
    `p-value` = ifelse(
      p.value < 0.001,
      "<0.001",
      sprintf("%.3f", p.value)
    )
  ) %>%
  select(
    Sex,
    Predictor,
    `OR (95% CI)`,
    `p-value`
  )

print(table4_results)


# -----------------------------------------------------------------------------
# 9.11.5 Export Table 4
# -----------------------------------------------------------------------------

ft_table4 <- flextable(table4_results) %>%
  autofit() %>%
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header")

read_docx() %>%
  body_add_par(
    "Table 4. Sex-Stratified GEE Estimates for Metabolic Syndrome",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft_table4) %>%
  body_add_par(
    "Adjusted models include visit week, age (per 1 SD), and CD4 count (per 1 SD).",
    style = "Normal"
  ) %>%
  print(
    target = "Table4_SexStratified_GEE.docx"
  )

# =============================================================================
# 9.12 Sensitivity Analysis: Excluding Week 144
# =============================================================================

gee_data_no144 <- gee_data %>%
  filter(
    clean_event != "Week 144"
  ) %>%
  droplevels()


gee_sensitivity_no144 <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = gee_data_no144,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# Compare with primary model
# -----------------------------------------------------------------------------

sensitivity_no144 <- bind_rows(
  broom::tidy(
    gee_model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    mutate(Analysis = "Primary analysis"),
  
  broom::tidy(
    gee_sensitivity_no144,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    mutate(Analysis = "Excluding Week 144")
) %>%
  filter(
    term %in% c(
      "Regimen_TypeTAFED",
      "Regimen_TypeTLD",
      "Sexfemale",
      "Age_scaled",
      "CD4_scaled"
    )
  )

print(sensitivity_no144)

# =============================================================================
# 9.13 Sensitivity Analysis: Baseline Viral Load <1000 copies/mL
# =============================================================================

# -----------------------------------------------------------------------------
# 9.13.1 Identify participants with baseline VL <1000 copies/mL
# -----------------------------------------------------------------------------

ids_vl_lt1000 <- metabolic_data_clean %>%
  filter(
    clean_event == "Baseline",
    !is.na(Viral_Load.cp.ml.),
    Viral_Load.cp.ml. < 1000
  ) %>%
  pull(ID) %>%
  as.character()


# -----------------------------------------------------------------------------
# 9.13.2 Create restricted longitudinal dataset
# -----------------------------------------------------------------------------

gee_data_vl_lt1000 <- metabolic_data_clean %>%
  mutate(
    ID = as.character(ID)
  ) %>%
  filter(
    ID %in% ids_vl_lt1000,
    clean_event != "Baseline"
  ) %>%
  mutate(
    Regimen_Type = droplevels(factor(Regimen_Type)),
    clean_event = droplevels(factor(clean_event)),
    Sex = droplevels(factor(Sex))
  )

# -----------------------------------------------------------------------------
# 9.13.3 Set reference categories
# -----------------------------------------------------------------------------

gee_data_vl_lt1000 <- gee_data_vl_lt1000 %>%
  mutate(
    Regimen_Type = forcats::fct_relevel(
      Regimen_Type,
      "PI_based"
    ),
    clean_event = forcats::fct_relevel(
      clean_event,
      "Week 24"
    ),
    Sex = droplevels(factor(Sex))
  )


# -----------------------------------------------------------------------------
# 9.13.4 Fit sensitivity GEE
# -----------------------------------------------------------------------------

gee_vl_lt1000 <- geeglm(
  metabolic_num ~
    Regimen_Type +
    clean_event +
    Sex +
    Age_scaled +
    CD4_scaled,
  id = ID,
  data = gee_data_vl_lt1000,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)


# -----------------------------------------------------------------------------
# 9.13.5 Extract results
# -----------------------------------------------------------------------------

table_s7_vl_lt1000 <- broom::tidy(
  gee_vl_lt1000,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  filter(
    term != "(Intercept)"
  ) %>%
  mutate(
    Term = recode(
      term,
      "Regimen_TypeTAFED"   = "TAFED vs PI-based",
      "Regimen_TypeTLD"     = "TLD vs PI-based",
      "Sexfemale"           = "Female vs Male",
      "Age_scaled"          = "Age (per 1 SD)",
      "CD4_scaled"          = "CD4 (per 1 SD)",
      "clean_eventWeek 48"  = "Week 48 vs Week 24",
      "clean_eventWeek 72"  = "Week 72 vs Week 24",
      "clean_eventWeek 96"  = "Week 96 vs Week 24",
      "clean_eventWeek 120" = "Week 120 vs Week 24",
      "clean_eventWeek 144" = "Week 144 vs Week 24",
      .default = term
    ),
    
    OR = round(estimate, 2),
    
    `95% CI` = paste0(
      sprintf("%.2f", conf.low),
      "–",
      sprintf("%.2f", conf.high)
    ),
    
    `p-value` = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  select(
    Term,
    OR,
    `95% CI`,
    `p-value`
  )

print(table_s7_vl_lt1000)

# =============================================================================
# END OF ANALYSIS PIPELINE
# =============================================================================

