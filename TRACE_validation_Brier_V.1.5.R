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

### Block 6. External validation: run the model in each cohort's own population ###
# Pre-specified: each cohort's age, sex and NYHA mix go in, theta stays at the base case (no re-fit), results are reported as they are.
# Model survival at months 6 to 60, sex-weighted
si_params_calib <- list(v_states_costs_arm = v_states_costs_si, l_trans_probs_arm = l_trans_probs_si_tv,
                        p_disc = p_disc_si, v_hosp_probs = v_hosp_probs_si, arm_name = "silencer")
get_cohort_curve <- function(age, sd, male, v_nyha, arm_params) {
  set.seed(seed + 1)
  com <- common_params
  com$m_indi_features_male   <- cbind(age = rnorm(num_i, age, sd))
  com$m_indi_features_female <- cbind(age = rnorm(num_i, age, sd))
  com$v_starting_states      <- sample(v_states_names, size = num_i, replace = TRUE, prob = c(v_nyha, 0))
  com$prop_male              <- male
  v_surv <- get_survival_curve(calibrated_theta, 10, arm_params, com)
  setNames(v_surv, (1:10) * 6)
}
compare_cohort <- function(df_obs, v_curve) {
  df_obs$model <- v_curve[as.character(df_obs$month)]
  df_obs$diff  <- df_obs$model - df_obs$km
  df_obs$in_ci <- df_obs$model >= df_obs$lower & df_obs$model <= df_obs$upper
  df_obs
}

### Block 7. External validation: trial survival (not used to fit) ###
# The figures show no confidence band, so the 95% CI comes from the at-risk table: Greenwood formula, one interval per table column.
# The HELIOS-B curve is IPTW-adjusted and its at-risk count falls fast after month 36, so its late CIs are approximate.
add_greenwood_ci <- function(df_km, v_month, v_n, v_events) {
  v_d <- diff(v_events)
  v_n_start <- v_n[-length(v_n)]
  v_var <- cumsum(v_d / (v_n_start * (v_n_start - v_d)))
  v_se <- df_km$km * sqrt(v_var[match(df_km$month, v_month[-1])])
  df_km$lower <- df_km$km - 1.96 * v_se
  df_km$upper <- df_km$km + 1.96 * v_se
  df_km
}

# ATTRibute-CM open-label extension, continuous acoramidis. Judge et al., Circulation 2025;151:601-611, Figure 1 (n = 409).
# Raster figure, read from its pixels. At months 6 and 18 a step falls on the tick: resolved with the cumulative events in the table.
df_attribute <- add_greenwood_ci(
  data.frame(cohort = "ATTRibute-CM (acoramidis)", month = (1:7) * 6,
             km = c(0.9806, 0.9415, 0.8923, 0.8412, 0.8069, 0.7972, 0.7622)),
  v_month  = (0:14) * 3,
  v_n      = c(409, 407, 401, 393, 385, 369, 365, 358, 344, 336, 297, 260, 247, 243, 216),
  v_events = c(0, 2, 8, 16, 24, 40, 44, 51, 65, 73, 79, 79, 82, 85, 93)
)
# HELIOS-B overall population, vutrisiran. Fontana et al., NEJM 2025, Figure 2A (n = 326). Vector paths of the PDF.
df_helios <- add_greenwood_ci(
  data.frame(cohort = "HELIOS-B (vutrisiran)", month = (1:7) * 6,
             km = c(0.9850, 0.9501, 0.9134, 0.8931, 0.8579, 0.8359, 0.8133)),
  v_month  = (0:7) * 6,
  v_n      = c(326, 321, 308, 296, 289, 277, 198, 33),
  v_events = c(0, 5, 18, 30, 37, 49, 56, 60)
)
# ATTR-ACT pooled tafamidis. Maurer et al., NEJM 2018, Figure 2B (n = 264). Vector paths of the PDF.
# Not independent: the transition probabilities come from this trial. Reference only.
df_attract <- add_greenwood_ci(
  data.frame(cohort = "ATTR-ACT (tafamidis; not independent)", month = (1:5) * 6,
             km = c(0.9500, 0.8902, 0.8206, 0.7607, 0.7008)),
  v_month  = (0:10) * 3,
  v_n      = c(264, 259, 252, 244, 235, 222, 216, 209, 200, 193, 99),
  v_events = c(0, 5, 12, 20, 29, 42, 48, 55, 64, 71, 78)
)

# Checks against the published numbers
# ATTRibute 30 months: 80.7% (76.9-84.5). ATTR-ACT 30 months: 186 of 264 alive; the published curve moves in steps of about 1 point.
v_attribute_30 <- unlist(df_attribute[df_attribute$month == 30, c("km", "lower", "upper")])
stopifnot(all(abs(v_attribute_30 - c(0.807, 0.769, 0.845)) < 0.001))
stopifnot(abs(df_attract$km[df_attract$month == 30] - 186 / 264) < 0.005)

# Population: age, SD, male, NYHA I-IV at baseline. HELIOS-B SD not reported, ATTRibute value assumed.
df_trials <- rbind(
  compare_cohort(df_attribute, get_cohort_curve(77.4, 6.5, 0.912, c(0.121, 0.696, 0.183, 0), st_params_calib)),
  compare_cohort(df_helios,    get_cohort_curve(77.0, 6.5, 0.920, c(0.150, 0.770, 0.080, 0), si_params_calib)),
  compare_cohort(df_attract,   get_cohort_curve(74.5, 7.2, 0.913, c(0.084, 0.596, 0.320, 0), st_params_calib))
)
cat("\nModel vs trial survival (theta =", round(calibrated_theta, 4), ", not re-fitted)\n")
print(df_trials, row.names = FALSE, digits = 3)

### Block 8. External validation: other real-world tafamidis cohorts (not used to fit) ###
# Rule set before running: tafamidis-treated, n >= 100, age and sex reported, overall survival reported at 30 months or later.
# Mayo is 85% treated and is kept so that no high-survival cohort is dropped. Masri and THAOS are in Blocks 1 and 2.
# Screened out: Greece (n = 75), Takashio 2023 (to 24 months), Sarkar 2024 and HEAR (baseline inputs not available), early access alone (within Fumagalli).
# NYHA not reported: Masri mix assumed (TriNetX, Mayo). Charite reports I/II together: split by the Masri ratio. Age from a median: SD = IQR / 1.349.
v_nyha_masri <- c(0.134, 0.519, 0.337, 0.010)
v_nyha_berlin <- c(0.641 * v_nyha_masri[1:2] / sum(v_nyha_masri[1:2]), 0.359, 0)
l_realworld <- list(
  # Fumagalli et al., Cardiol Ther 2026;15(Suppl 1):61-80, Table 1 and Figure 2. Early access cohort + THAOS 2019-2023; time zero = enrollment or tafamidis start.
  list(cohort = "Fumagalli <80 y (n = 1188)", age = 72.3, sd = 6.12, male = 0.902, nyha = c(0.182, 0.560, 0.249, 0.009),
       month = c(30, 42), km = c(0.770, 0.709), lower = c(0.734, 0.663), upper = c(0.802, 0.750)),
  list(cohort = "Fumagalli 80-85 y (n = 505)", age = 82.1, sd = 1.62, male = 0.875, nyha = c(0.099, 0.523, 0.370, 0.008),
       month = c(30, 42), km = c(0.662, 0.572), lower = c(0.599, 0.495), upper = c(0.717, 0.642)),
  list(cohort = "Fumagalli >85 y (n = 231)", age = 88.0, sd = 2.22, male = 0.848, nyha = c(0.093, 0.528, 0.346, 0.033),
       month = c(30, 42), km = c(0.563, 0.429), lower = c(0.463, 0.301), upper = c(0.651, 0.551)),
  # Debonnaire et al., J Comp Eff Res 2026;15(10):e260114, Tables 1-3. TriNetX new users, reweighted to each trial (ReplicATTR); survival = 1 - 42-month mortality.
  list(cohort = "TriNetX, ATTRibute-like (n = 617)", age = 77.3, sd = 6.5, male = 0.914, nyha = v_nyha_masri,
       month = 42, km = 0.646, lower = 0.580, upper = 0.704),
  list(cohort = "TriNetX, ATTR-ACT-like (n = 624)", age = 74.5, sd = 7.2, male = 0.913, nyha = v_nyha_masri,
       month = 42, km = 0.631, lower = 0.566, upper = 0.689),
  # Porcari et al., Eur J Heart Fail 2025;27(12):2952-2960. Italy, wild-type only; time zero = diagnosis, so time before treatment start counts as survival.
  list(cohort = "Italy, treated (n = 965)", age = 78.5, sd = 6.5, male = 0.896, nyha = c(0.177, 0.761, 0.062, 0),
       month = 30, km = 0.897, lower = 0.866, upper = 0.921),
  # Nies et al., ESC Heart Fail 2026;13(5):xvag211. Germany, 8 centers; baseline = first visit in 2020-2021. No CI reported.
  list(cohort = "Germany, treated (n = 259)", age = 79, sd = 7 / 1.349, male = 0.866, nyha = c(0.111, 0.438, 0.443, 0.009),
       month = c(30, 36), km = c(0.88, 0.86), lower = c(NA, NA), upper = c(NA, NA)),
  # Muller et al., J Am Heart Assoc 2024;13(14):e033478. Berlin; time zero = tafamidis start; median follow-up 18 months, so 60 months rests on few patients. No CI reported.
  list(cohort = "Berlin, treated (n = 167)", age = 80, sd = 7 / 1.349, male = 0.85, nyha = v_nyha_berlin,
       month = 60, km = 0.702, lower = NA, upper = NA),
  # Skov et al., J Am Heart Assoc 2026;15(16):e048585. Mayo Clinic, 85% on tafamidis; start near diagnosis.
  list(cohort = "Mayo, 85% treated (n = 441)", age = 76.3, sd = 8.3, male = 0.932, nyha = v_nyha_masri,
       month = 36, km = 0.787, lower = 0.747, upper = 0.831)
)
df_realworld <- do.call(rbind, lapply(l_realworld, function(k) {
  compare_cohort(data.frame(cohort = k$cohort, month = k$month, km = k$km, lower = k$lower, upper = k$upper),
                 get_cohort_curve(k$age, k$sd, k$male, k$nyha, st_params_calib))
}))

# Fumagalli strata pooled, weighted by patients at baseline; CI from the stratum CIs, strata treated as independent
v_fuma_n <- c(1188, 505, 231)
df_fuma <- df_realworld[grepl("Fumagalli", df_realworld$cohort), ]
df_fuma_pooled <- do.call(rbind, lapply(c(30, 42), function(m) {
  df_m <- df_fuma[df_fuma$month == m, ]
  v_w <- v_fuma_n / sum(v_fuma_n)
  v_se <- (df_m$upper - df_m$lower) / (2 * 1.96)
  km <- sum(v_w * df_m$km)
  se <- sqrt(sum((v_w * v_se)^2))
  data.frame(cohort = "Fumagalli, all ages (n = 1924)", month = m, km = km, lower = km - 1.96 * se, upper = km + 1.96 * se,
             model = sum(v_w * df_m$model))
}))
df_fuma_pooled$diff  <- df_fuma_pooled$model - df_fuma_pooled$km
df_fuma_pooled$in_ci <- df_fuma_pooled$model >= df_fuma_pooled$lower & df_fuma_pooled$model <= df_fuma_pooled$upper
df_realworld <- rbind(df_realworld, df_fuma_pooled)

cat("\nModel vs other real-world cohorts (theta =", round(calibrated_theta, 4), ", not re-fitted)\n")
print(df_realworld, row.names = FALSE, digits = 3)

df_fit <- rbind(df_fit,
  summarize_fit("ATTRibute-CM OLE, months 6-42", df_trials$model[df_trials$cohort == df_attribute$cohort[1]], df_attribute$km),
  summarize_fit("HELIOS-B, months 6-42", df_trials$model[df_trials$cohort == df_helios$cohort[1]], df_helios$km),
  summarize_fit("ATTR-ACT, months 6-30 (not independent)", df_trials$model[df_trials$cohort == df_attract$cohort[1]], df_attract$km),
  summarize_fit("Fumagalli, three age strata, months 30 and 42", df_fuma$model, df_fuma$km)
)

### Block 9. Save the results ###
dir_val <- "Validation_output"
dir.create(dir_val, showWarnings = FALSE)
write.csv(df_masri, file.path(dir_val, "Brier_Masri_V.1.5.csv"), row.names = FALSE)
write.csv(df_thaos, file.path(dir_val, "Brier_THAOS_V.1.5.csv"), row.names = FALSE)
write.csv(df_fit,   file.path(dir_val, "Brier_summary_V.1.5.csv"), row.names = FALSE)
write.csv(df_nyha,  file.path(dir_val, "Brier_NYHA_V.1.5.csv"), row.names = FALSE)
write.csv(df_trials, file.path(dir_val, "External_trials_V.1.5.csv"), row.names = FALSE)
write.csv(df_realworld, file.path(dir_val, "External_realworld_V.1.5.csv"), row.names = FALSE)
