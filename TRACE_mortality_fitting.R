# clear R's session memory (Global Environment)
rm(list = ls())

#install.packages("flexsurv")
#install.packages("survminer")

library(survival)
library(flexsurv)
library(ggplot2)
library(dplyr)
library(kableExtra)
library(survminer)

# to caclulate cumulative hazard every cycle with the continuous age

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

# Filter for adult population (e.g., age > 50) relevant for ATTR-CM -> THIS MAY BE WRONG! LET'S USE THE WHOLE POPULATION
# to ensure the fit is optimized for the older population.
# df_adults <- df_life_table %>% filter(time > 50)

# 2. Plot Actual Mortality (Kaplan-Meier) by Sex ----------------------------

# Create survival object
surv_obj <- Surv(time = df_life_table$time, event = df_life_table$status)

# Fit KM
km_fit <- survfit(surv_obj ~ Sex, data = df_life_table)

# Plot
ggsurvplot(
  km_fit, 
  data = df_life_table,
  title = "Actual Survival Curve by Sex (2022 Actuarial Table)",
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
  data = df_life_table %>% filter(Sex == "Male"), 
  dist = "gompertz"
)

# Fit Female
fit_female <- flexsurvreg(
  formula = Surv(time, status) ~ 1, 
  data = df_life_table %>% filter(Sex == "Female"), 
  dist = "gompertz"
)

# Print Results
cat("\n--- Male Fit ---\n")
print(fit_male)
cat("\n--- Female Fit ---\n")
print(fit_female)

# Plot the fits
plot(fit_male, main = "Male: Gompertz Fit vs Actual", xlab = "Age", ylab = "Survival")
# add legend
legend("bottomleft", legend = c("Actual (Kaplan-Meier)", "Fitted (Gompertz)"), 
       col = c("black", "red"), lty = 1, lwd = 2, bty = "n")

plot(fit_female, main = "Female: Gompertz Fit vs Actual", xlab = "Age", ylab = "Survival")
# add legend
legend("bottomleft", legend = c("Actual (Kaplan-Meier)", "Fitted (Gompertz)"), 
       col = c("black", "red"), lty = 1, lwd = 2, bty = "n")

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

# 5. Left-truncated Gompertz fit (age >= 60) --------------------------------
# The whole-life fit above overestimates mortality at ages 70-85 (the ATTR-CM range).
# Refit by maximum likelihood on the life table, conditional on survival to age 60:
# each age x contributes d_x deaths in [x, x+1) and l_{x+1} survivors, with
# h(a) = rate * exp(shape * a), H(x, x+1) = (rate / shape) * (exp(shape * (x + 1)) - exp(shape * x))
fit_age_min <- 60

fit_gompertz_lt <- function(ages, qx, age_min = fit_age_min) {
  keep <- ages >= age_min & qx < 1
  x  <- ages[keep]; q <- qx[keep]
  lx <- cumprod(c(1, 1 - q))[seq_along(q)] # survivors at start of each age, from age_min
  dx <- lx * q # ㅔ proportion of individual who die
  negll <- function(par) {
    shape <- exp(par[1]); rate <- exp(par[2]) # shape, rate always
    H <- (rate / shape) * (exp(shape * (x + 1)) - exp(shape * x))
    p_die <- pmin(pmax(1 - exp(-H), 1e-12), 1 - 1e-12) # with these parameters -> p_die?
    -sum(dx * log(p_die) + (lx - dx) * log(1 - p_die)) # LIKELIHOOD
  }
  opt <- optim(c(log(0.09), log(1e-5)), negll, control = list(maxit = 5000, reltol = 1e-12))# find the optimal params
  list(shape = exp(opt$par[1]), rate = exp(opt$par[2]))
}

params_gompertz_lt <- list(
  male   = fit_gompertz_lt(df_raw_lifetable$Exact.age, df_raw_lifetable$Male.Death.probability),
  female = fit_gompertz_lt(df_raw_lifetable$Exact.age, df_raw_lifetable$Female.Death.probability)
)
print(params_gompertz_lt)

# Check: 3.5-year survival, life table vs left-truncated Gompertz
surv_gompertz <- function(a, g, t) exp(-(g$rate / g$shape) * (exp(g$shape * (a + t)) - exp(g$shape * a)))
surv_lifetable <- function(a, qx, t) {
  l <- approx(df_raw_lifetable$Exact.age, log(cumprod(c(1, 1 - qx))[seq_along(qx)]), xout = c(a, a + t))$y
  exp(l[2] - l[1])
}
for (a in c(65, 71, 75, 77.5, 80, 85, 90)) {
  cat(sprintf("Age %.1f | Male: LT %.3f Gompertz %.3f | Female: LT %.3f Gompertz %.3f\n", a,
              surv_lifetable(a, df_raw_lifetable$Male.Death.probability, 3.5),   surv_gompertz(a, params_gompertz_lt$male, 3.5),
              surv_lifetable(a, df_raw_lifetable$Female.Death.probability, 3.5), surv_gompertz(a, params_gompertz_lt$female, 3.5)))
}
