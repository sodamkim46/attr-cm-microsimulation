# Scenario analysis: mortality calibration (Silencer vs Stabilizer, base case population)
# Base case: one theta fitted to the Masri survival curve, months 6-54
# Scenarios: THAOS 42-month relative survival target (one theta, the previous base case);
#            two-part theta fitted to the Masri curve, split at 18 months (as in the ICER 2024 model)
#            or at 36 months (chosen after seeing the curve, so exploratory)
# After the last observed month, the late theta applies for life. Everything else stays at base case.

# Step 1. Run the base case model and load the SA helper functions
source("TRACE_microsim_V.1.5.R")
source("TRACE_SA_functions_V.1.5.R") # base_inputs, run_owsa_case(), calibrate_theta()

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

### Step 2. Two-part theta ###
make_theta <- function(theta_early, theta_late, split_month) {
  n_early <- split_month / 12 / cycle_length
  c(rep(theta_early, n_early), rep(theta_late, num_cycles - n_early))
}

fit_two_part <- function(split_month) {
  calc_error <- function(v_par) {
    if (any(v_par < 0)) return(1e6) # theta >= 0
    calc_calib_error(make_theta(v_par[1], v_par[2], split_month), df_calib_masri$km)
  }
  optim(c(calibrated_theta, calibrated_theta), calc_error, method = "Nelder-Mead", control = list(maxit = 36, reltol = 1e-3))
}

### Step 3. Run the scenarios ###
run_calib_case <- function(label, theta_early, theta_late, split_month) {
  v_theta <- if (is.na(split_month)) theta_early else make_theta(theta_early, theta_late, split_month)
  inp <- base_inputs
  inp$theta <- v_theta
  r <- run_owsa_case(inp)
  v_surv <- get_survival_curve(v_theta, 11, st_params_calib, common_params) # up to 66 months
  data.frame(
    scenario = label, theta_early = theta_early, theta_late = theta_late, split_month = split_month,
    weighted_error = sum(((v_surv[1:9] - df_calib_masri$km) / df_calib_masri$se)^2),
    surv_42 = v_surv[7], surv_54 = v_surv[9], surv_60 = v_surv[10],
    surv_65 = approx(c(60, 66), v_surv[10:11], xout = 65)$y,
    LYs_st = r["LYs_st"], LYs_si = r["LYs_si"], inc_Dcosts = r["inc_Dcosts"], inc_Dqalys = r["inc_Dqalys"], ICER = r["ICER"],
    list_price_100k = r["list_price_100k"], list_price_150k = r["list_price_150k"], row.names = NULL
  )
}

cat("\n--- SCENARIO: mortality calibration ---\n")
theta_thaos <- calibrate_theta(ref_surv_obs)
cat("THAOS 42-month theta:", theta_thaos, "\n")
fit_18 <- fit_two_part(18)
fit_36 <- fit_two_part(36)
cat("Two-part theta, split at 18 months:", fit_18$par, "| split at 36 months:", fit_36$par, "\n")

calib_results <- rbind(
  run_calib_case("Base case: Masri curve, one theta", calibrated_theta, NA, NA),
  run_calib_case("THAOS 42-month relative survival, one theta", theta_thaos, NA, NA),
  run_calib_case("Masri curve, two-part theta, split at 18 months", fit_18$par[1], fit_18$par[2], 18),
  run_calib_case("Masri curve, two-part theta, split at 36 months", fit_36$par[1], fit_36$par[2], 36)
)
print(calib_results, row.names = FALSE, digits = 5)
write.csv(calib_results, file.path(dir_out, "Scenario_calibration_V.1.5.csv"), row.names = FALSE)
