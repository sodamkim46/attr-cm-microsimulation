# Mortality Fitting Script using 2022 Actuarial Life Table Logic
# Fits a Gompertz model using flexsurv

rm(list = ls())

if (!require("flexsurv")) {
  install.packages("flexsurv", repos = "https://cloud.r-project.org")
}

## load packages
library(survival)
library(flexsurv)
library(ggplot2)
library(dplyr)
library(kableExtra)
library(survminer)

# 1. Load 2022 Actuarial Life Table & Simulate Cohort -----------------------
# Load the raw CSV file.
# We assume it has columns 'Exact.age', 'Male.Death.probability', 'Female.Death.probability'
# Please ensure your CSV file is in the same directory as this script.
tryCatch({
  df_raw_lifetable <- read.csv("2022 Actuarial Life Table.csv")
}, error = function(e) {
  stop("Could not find '2022 Actuarial Life Table.csv'. Please ensure the file is in the correct directory and named correctly.")
})

set.seed(123)

ages_from_csv <- df_raw_lifetable$Exact.age

# Function to generate synthetic cohort for a given qx vector
generate_cohort <- function(qx_vector, ages, sex_label, n_sim = 50000) {
  # Construct the survival function S(t)
  # S(t) = product_{i=0 to t-1} (1 - qx_i)
  surv_probs <- c(1, cumprod(1 - qx_vector))
  
  # Corresponding time points (0, 1, 2, ...)
  time_points <- c(min(ages), ages + 1)
  
  # The CDF is F(t) = 1 - S(t)
  cdf_death  <- 1 - surv_probs
  
  # Sample random death ages
  random_uniform <- runif(n_sim)
  death_ages <- approx(x = cdf_death, y = time_points, xout = random_uniform, rule = 2)$y
  
  data.frame(
    ID = 1:n_sim,
    time = death_ages,
    status = 1, # All are events (death)
    Sex = sex_label
  )
}

# Generate cohorts for Male and Female
df_male   <- generate_cohort(df_raw_lifetable$Male.Death.probability, ages_from_csv, "Male")
df_female <- generate_cohort(df_raw_lifetable$Female.Death.probability, ages_from_csv, "Female")

df_life_table <- rbind(df_male, df_female)

# Filter for adult population (e.g., age > 50) relevant for ATTR-CM
# to ensure the fit is optimized for the older population.
df_adults <- df_life_table %>% filter(time > 50)

# 2. Plot Actual Mortality (Kaplan-Meier) by Sex ----------------------------

# Create survival object
surv_obj <- Surv(time = df_adults$time, event = df_adults$status)

# Fit KM
km_fit <- survfit(surv_obj ~ Sex, data = df_adults)

# Plot
ggsurvplot(
  km_fit, 
  data = df_adults,
  title = "Actual Survival Curve by Sex (2022 Actuarial Table - Adults > 50)",
  xlab = "Age (Years)",
  ylab = "Survival Probability",
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  legend.labs = c("Female", "Male")
)

# 3. Fit Parametric Model (Gompertz) Separately -----------------------------

# Fit Male
fit_male <- flexsurvreg(
  formula = Surv(time, status) ~ 1, 
  data = df_adults %>% filter(Sex == "Male"), 
  dist = "gompertz"
)

# Fit Female
fit_female <- flexsurvreg(
  formula = Surv(time, status) ~ 1, 
  data = df_adults %>% filter(Sex == "Female"), 
  dist = "gompertz"
)

# Print Results
cat("\n--- Male Fit ---\n")
print(fit_male)
cat("\n--- Female Fit ---\n")
print(fit_female)

# Plot the fits
plot(fit_male, main = "Male: Gompertz Fit vs Actual", xlab = "Age", ylab = "Survival")
plot(fit_female, main = "Female: Gompertz Fit vs Actual", xlab = "Age", ylab = "Survival")

# 4. Extract and Save Parameters --------------------------------------------
params <- list(
  male = list(
    gom_shape = fit_male$res["shape", "est"],
    gom_rate  = fit_male$res["rate", "est"]
  ),
  female = list(
    gom_shape = fit_female$res["shape", "est"],
    gom_rate  = fit_female$res["rate", "est"]
  )
)

print(params)
