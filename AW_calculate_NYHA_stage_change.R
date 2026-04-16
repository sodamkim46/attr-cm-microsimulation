

# Install and load survminer (restart R session before running this)
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
library(survival)
library(survminer)

# STEP 1: Find each individual's time to death
# apply() scans each row (individual) to find the first cycle where state == "Death"
time_to_death <- apply(res_st_male$m_States, MARGIN = 1, FUN = function(individual_row) {
  death_positions <- which(individual_row == "Death")
  if (length(death_positions) == 0) {
    return(NA)  # Never died = censored
  }
  return(death_positions[1] - 1)  # Subtract 1 to convert index back to cycle number
})

# STEP 2: Create time and event variables for survival analysis
total_cycles      <- ncol(res_st_male$m_States) - 1
event_indicator   <- ifelse(!is.na(time_to_death), 1, 0)
time_to_death_clean <- ifelse(!is.na(time_to_death), time_to_death, total_cycles)

# STEP 3: Build a clean data frame
survival_data <- data.frame(
  individual = rownames(res_st_male$m_States),
  time       = time_to_death_clean,
  event      = event_indicator
)
head(survival_data)  # Sanity check

# STEP 4: Fit the Kaplan-Meier curve
km_fit <- survfit(Surv(time, event) ~ 1, data = survival_data)
summary(km_fit)

# STEP 5: Plot
ggsurvplot(
  fit              = km_fit,
  data             = survival_data,
  conf.int         = TRUE,
  risk.table       = TRUE,
  xlab             = "Cycle (Time)",
  ylab             = "Survival Probability",
  title            = "Kaplan-Meier Curve: Time to Death",
  surv.median.line = "hv",
  ggtheme          = theme_minimal()
)


######Proportion same or better at cycle 6########
# STEP 1: Assign a numeric severity rank to each NYHA state
# Lower number = better health state
nyha_rank <- c(
  "NYHA1" = 1,
  "NYHA2" = 2,
  "NYHA3" = 3,
  "NYHA4" = 4
)

# STEP 2: Pull each individual's state at baseline (cycle 0) and cycle 6
state_at_baseline <- res_st_male$m_States[ , "cycle_0"]
state_at_cycle6   <- res_st_male$m_States[ , "cycle_6"]

# STEP 3: Filter to only individuals alive at cycle 6
alive_at_cycle6  <- state_at_cycle6 != "Death"
baseline_alive   <- state_at_baseline[alive_at_cycle6]
cycle6_alive     <- state_at_cycle6[alive_at_cycle6]
n_alive          <- sum(alive_at_cycle6)
cat("Number alive at cycle 6:", n_alive, "\n")

# STEP 4: Convert health state names to numeric severity ranks
rank_baseline <- nyha_rank[baseline_alive]
rank_cycle6   <- nyha_rank[cycle6_alive]

# STEP 5: Determine who is same or better at cycle 6
# Same or better = cycle 6 rank <= baseline rank
same_or_better <- rank_cycle6 <= rank_baseline

# STEP 6: Calculate and report the proportion
n_same_or_better <- sum(same_or_better)
proportion       <- n_same_or_better / n_alive
cat("Number same or better at cycle 6:", n_same_or_better, "\n")
cat("Proportion same or better at cycle 6:", round(proportion, 3), "\n")
cat("Percentage same or better at cycle 6:", round(proportion * 100, 1), "%\n")

# STEP 7 (Optional): Cross-tab of baseline vs. cycle 6 health states
# Useful for seeing exactly how people moved between states
state_comparison_table <- table(
  Baseline = baseline_alive,
  Cycle_6  = cycle6_alive
)
print(state_comparison_table) #50.9%


################################################################################
### Analysis for Standard of Care (SoC) Arm
################################################################################
# This section calculates the same metrics for the SoC arm, using the results
# generated from the 'Validation using Excel model.R' script (res_soc_male).
# Ensure that script has been run so 'res_soc_male' is in the environment.

cat("\n--- SoC Arm Analysis ---\n")

# --- Kaplan-Meier Curve for SoC Arm ---
cat("\n--- Kaplan-Meier Curve: Time to Death (SoC Arm) ---\n")
# STEP 1: Find each individual's time to death for SoC arm
time_to_death_soc <- apply(res_soc_male$m_States, MARGIN = 1, FUN = function(individual_row) {
  death_positions <- which(individual_row == "Death")
  if (length(death_positions) == 0) {
    return(NA)  # Never died = censored
  }
  return(death_positions[1] - 1)  # Subtract 1 to convert index back to cycle number
})

# STEP 2: Create time and event variables for survival analysis for SoC arm
total_cycles_soc      <- ncol(res_soc_male$m_States) - 1
event_indicator_soc   <- ifelse(!is.na(time_to_death_soc), 1, 0)
time_to_death_clean_soc <- ifelse(!is.na(time_to_death_soc), time_to_death_soc, total_cycles_soc)

# STEP 3: Build a clean data frame for SoC arm
survival_data_soc <- data.frame(
  individual = rownames(res_soc_male$m_States),
  time       = time_to_death_clean_soc,
  event      = event_indicator_soc
)
head(survival_data_soc)  # Sanity check

# STEP 4: Fit the Kaplan-Meier curve for SoC arm
km_fit_soc <- survfit(Surv(time, event) ~ 1, data = survival_data_soc)
summary(km_fit_soc)

# STEP 5: Plot for SoC arm
ggsurvplot(
  fit              = km_fit_soc,
  data             = survival_data_soc,
  conf.int         = TRUE,
  risk.table       = TRUE,
  xlab             = "Cycle (Time)",
  ylab             = "Survival Probability",
  title            = "Kaplan-Meier Curve: Time to Death (SoC Arm)",
  surv.median.line = "hv",
  ggtheme          = theme_minimal()
)

######Proportion same or better at cycle 6 for SoC Arm########

################################################################################
### Analysis for Standard of Care (SoC) Arm
################################################################################
# This section calculates the same metric for the SoC arm, using the results
# generated from the 'Validation using Excel model.R' script (res_soc_male).
# Ensure that script has been run so 'res_soc_male' is in the environment.

######Proportion same or better at cycle 6 for SoC Arm########

# STEP 2: Pull each individual's state at baseline (cycle 0) and cycle 6 from SoC results
state_at_baseline_soc <- res_soc_male$m_States[ , "cycle_0"]
state_at_cycle6_soc   <- res_soc_male$m_States[ , "cycle_6"]

# STEP 3: Filter to only individuals alive at cycle 6
alive_at_cycle6_soc  <- state_at_cycle6_soc != "Death"
baseline_alive_soc   <- state_at_baseline_soc[alive_at_cycle6_soc]
cycle6_alive_soc     <- state_at_cycle6_soc[alive_at_cycle6_soc]
n_alive_soc          <- sum(alive_at_cycle6_soc)
cat("\n--- SoC Arm Analysis (Cycle 6) ---\n")
cat("Number alive at cycle 6 (SoC):", n_alive_soc, "\n")

# STEP 4: Convert health state names to numeric severity ranks
rank_baseline_soc <- nyha_rank[baseline_alive_soc]
rank_cycle6_soc   <- nyha_rank[cycle6_alive_soc]

# STEP 5: Determine who is same or better at cycle 6
same_or_better_soc <- rank_cycle6_soc <= rank_baseline_soc

# STEP 6: Calculate and report the proportion
n_same_or_better_soc <- sum(same_or_better_soc)
proportion_soc       <- n_same_or_better_soc / n_alive_soc
cat("Number same or better at cycle 6 (SoC):", n_same_or_better_soc, "\n")
cat("Proportion same or better at cycle 6 (SoC):", round(proportion_soc, 3), "\n")
cat("Percentage same or better at cycle 6 (SoC):", round(proportion_soc * 100, 1), "%\n")

# STEP 7 (Optional): Cross-tab of baseline vs. cycle 6 health states for SoC
state_comparison_table_soc <- table(
  Baseline_SoC = baseline_alive_soc,
  Cycle_6_SoC  = cycle6_alive_soc
)
print(state_comparison_table_soc) #35.9%