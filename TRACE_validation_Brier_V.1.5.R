# Brier score: model survival vs observed real-world survival
# Definition follows the DOC-M model (Kim et al., Med Decis Making 2023;43(7-8):930-948):
# mean squared difference between predicted and observed aggregate survival.

source("TRACE_microsim_V.1.5.R") # base case results; it starts with rm(list = ls()), so it must come first

### Block 1. Calibration data: Masri tafamidis-treated cohort ###
# Masri et al., JACC CardioOncol 2025;7(3):282-293, Figure 1 (n = 624; 5 US centers; time zero = tafamidis start).
# The figure is a raster image: values read from its pixels (about 0.2 percentage points per pixel), at each 6-month tick.
# n_risk is approximate: the paper reports it every 5 months, so it is interpolated to 6-month ticks.
# Months 6-54 are the calibration points; months 60 and 65 have few patients at risk and are not used to fit.
df_masri <- data.frame(
  month  = c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 65),
  n_risk = c(591, 553, 513, 476, 446, 403, 336, 255, 141, 78, 7),
  km     = c(0.9491, 0.8911, 0.8312, 0.7682, 0.7253, 0.6803, 0.6503, 0.6224, 0.5794, 0.5674, 0.5425),
  lower  = c(0.9301, 0.8661, 0.8022, 0.7343, 0.6903, 0.6444, 0.6124, 0.5844, 0.5385, 0.5225, 0.4745),
  upper  = c(0.9640, 0.9121, 0.8601, 0.7982, 0.7582, 0.7143, 0.6863, 0.6603, 0.6224, 0.6084, 0.6044)
)
df_masri$use <- ifelse(df_masri$month <= 54, "calibration", "not fitted")
stopifnot(all(df_masri$km[df_masri$month <= 54] == df_calib_masri$km)) # same points as the base case calibration

# Check: the only values published as numbers are at month 65: 54.1% (95% CI 47.4-60.4)
v_masri_published <- c(km = 0.541, lower = 0.474, upper = 0.604)
v_masri_extracted <- unlist(df_masri[df_masri$month == 65, c("km", "lower", "upper")])
stopifnot(all(abs(v_masri_extracted - v_masri_published) < 0.002))

print(df_masri, row.names = FALSE)
cat("Masri: largest difference from the published values:", max(abs(v_masri_extracted - v_masri_published)), "\n")

### Block 2. Validation data: THAOS tafamidis-treated cohort (not used to fit) ###
# Garcia-Pavia et al., J Card Fail 2025;31(3):525-533, Kaplan-Meier figure (n = 587).
# Values read from the vector paths of the figure in the PDF, at each 6-month tick.
# The lower bound at month 6 is hidden by the overlapping band of the untreated group.
df_thaos <- data.frame(
  month  = c(6, 12, 18, 24, 30, 36, 42),
  n_risk = c(549, 475, 401, 317, 232, 145, 65),
  km     = c(0.9774, 0.9476, 0.9221, 0.9018, 0.8444, 0.8144, 0.7681),
  lower  = c(NA, 0.9279, 0.8955, 0.8717, 0.8047, 0.7693, 0.7094),
  upper  = c(0.9865, 0.9650, 0.9421, 0.9251, 0.8766, 0.8515, 0.8166)
)

# Check: the extraction reproduces the values published in the paper (within 0.1 percentage point)
v_published <- c(km_30 = 0.844, lower_30 = 0.805, upper_30 = 0.877,
                 km_42 = 0.768, lower_42 = 0.709, upper_42 = 0.817)
v_extracted <- c(km_30 = df_thaos$km[df_thaos$month == 30], lower_30 = df_thaos$lower[df_thaos$month == 30],
                 upper_30 = df_thaos$upper[df_thaos$month == 30], km_42 = df_thaos$km[df_thaos$month == 42],
                 lower_42 = df_thaos$lower[df_thaos$month == 42], upper_42 = df_thaos$upper[df_thaos$month == 42])
stopifnot(all(abs(v_extracted - v_published) < 0.001))

print(df_thaos, row.names = FALSE)
cat("THAOS: largest difference from the published values:", max(abs(v_extracted - v_published)), "\n")

### Block 3. Model survival next to the observed survival ###
# Stabilizer arm, sex-weighted. Both cohorts were tafamidis-treated, so they are compared with this arm.
v_surv_st <- prop_male * colMeans(res_st_male$m_States != "Death") +
  prop_female * colMeans(res_st_female$m_States != "Death")
v_month_model <- (0:num_cycles) * cycle_length * 12

# approx(): model value at each observed month; month 65 falls between two cycles, so it is interpolated
df_masri$model <- approx(v_month_model, v_surv_st, xout = df_masri$month)$y
df_thaos$model <- approx(v_month_model, v_surv_st, xout = df_thaos$month)$y
df_masri$diff  <- df_masri$model - df_masri$km
df_thaos$diff  <- df_thaos$model - df_thaos$km

cat("\nModel vs Masri (theta =", round(calibrated_theta, 4), ")\n")
print(df_masri[, c("month", "n_risk", "km", "model", "diff", "use")], row.names = FALSE, digits = 3)
cat("\nModel vs THAOS, before adjusting THAOS to the model population\n")
print(df_thaos[, c("month", "n_risk", "km", "model", "diff")], row.names = FALSE, digits = 3)

### Block 4. Move THAOS to the model population, then calculate the Brier score ###
# THAOS has a different age and sex mix, so survival is multiplied by the ratio of expected (general population) survival,
# the same relative survival transport as in the base case script. Masri is the model population: no adjustment.
calc_expected_mix <- function(v_age_male, v_age_female, p_male, t) {
  p_male * calc_expected_surv(v_age_male, params_gompertz$male, t) +
    (1 - p_male) * calc_expected_surv(v_age_female, params_gompertz$female, t)
}
v_ratio <- sapply(df_thaos$month / 12, function(t) {
  calc_expected_mix(m_indi_features_male[, "age"], m_indi_features_female[, "age"], prop_male, t) /
    calc_expected_mix(v_age_ref_male, v_age_ref_female, ref_prop_male, t)
})
df_thaos$target       <- df_thaos$km * v_ratio
df_thaos$target_lower <- df_thaos$lower * v_ratio
df_thaos$target_upper <- df_thaos$upper * v_ratio
df_thaos$diff         <- df_thaos$model - df_thaos$target

# Brier score (DOC-M definition) with the root mean squared error, which is on the survival scale
summarize_fit <- function(label, v_model, v_observed) {
  v_diff <- v_model - v_observed
  data.frame(comparison = label, n_points = length(v_diff), brier = mean(v_diff^2),
             rmse_pp = 100 * sqrt(mean(v_diff^2)), max_abs_diff_pp = 100 * max(abs(v_diff)))
}
is_fit <- df_masri$use == "calibration"
df_fit <- rbind(
  summarize_fit("Masri, months 6-54 (calibration points)", df_masri$model[is_fit], df_masri$km[is_fit]),
  summarize_fit("Masri, months 60 and 65 (not fitted)", df_masri$model[!is_fit], df_masri$km[!is_fit]),
  summarize_fit("THAOS, months 6-42 (not fitted)", df_thaos$model, df_thaos$target)
)

cat("\nModel vs THAOS, moved to the model population\n")
print(df_thaos[, c("month", "target", "target_lower", "target_upper", "model", "diff")], row.names = FALSE, digits = 3)
cat("\nFit summary (theta =", round(calibrated_theta, 4), ")\n")
print(df_fit, row.names = FALSE, digits = 3)

### Block 5. Survival by starting NYHA class: Masri 65-month values (not used to fit) ###
# Masri et al. 2025, Figure 2A. NYHA IV: 6 patients at baseline. Few patients remain at risk at 65 months.
df_nyha <- data.frame(
  nyha  = c("NYHA1", "NYHA2", "NYHA3", "NYHA4"),
  n     = c(84, 324, 210, 6),
  km    = c(0.830, 0.620, 0.332, 0.000),
  lower = c(0.721, 0.554, 0.204, 0.000),
  upper = c(0.899, 0.679, 0.465, 0.000)
)
surv_by_start <- function(res, k) colMeans(res$m_States[res$m_States[, 1] == k, ] != "Death")
df_nyha$model <- sapply(df_nyha$nyha, function(k) {
  v_surv_k <- prop_male * surv_by_start(res_st_male, k) + prop_female * surv_by_start(res_st_female, k)
  approx(v_month_model, v_surv_k, xout = 65)$y
})
df_nyha$diff <- df_nyha$model - df_nyha$km
df_nyha$in_ci <- df_nyha$model >= df_nyha$lower & df_nyha$model <= df_nyha$upper

cat("\nModel vs Masri at 65 months, by starting NYHA class\n")
print(df_nyha, row.names = FALSE, digits = 3)

### Block 6. Save the results ###
dir_val <- "Validation_output"
dir.create(dir_val, showWarnings = FALSE)
write.csv(df_masri, file.path(dir_val, "Brier_Masri_V.1.5.csv"), row.names = FALSE)
write.csv(df_thaos, file.path(dir_val, "Brier_THAOS_V.1.5.csv"), row.names = FALSE)
write.csv(df_fit,   file.path(dir_val, "Brier_summary_V.1.5.csv"), row.names = FALSE)
write.csv(df_nyha,  file.path(dir_val, "Brier_NYHA_V.1.5.csv"), row.names = FALSE)
