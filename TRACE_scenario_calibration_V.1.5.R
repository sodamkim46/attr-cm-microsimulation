# Scenario analysis: mortality calibration (Silencer vs Stabilizer, base case population)
# Base case: one theta fitted to the Masri survival curve, months 6-54
# Scenarios: THAOS 42-month relative survival target (one theta, the previous base case);
#            two-part theta fitted to the Masri curve, split at 18 months (as in the ICER 2024 model)
#            or at 36 months (chosen after seeing the curve, so exploratory);
#            one theta fitted by unweighted least squares (Vanni et al. 2011, Eq. 1) instead of chi-squared;
#            NYHA moves slowed to 75% or 50% in both transition matrices, theta re-fitted (exploratory ranges)
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

### Step 3. Unweighted least squares, and slower NYHA transitions ###
fit_least_squares <- function() {
  calc_error <- function(theta_val) {
    mean((get_survival_curve(theta_val, nrow(df_calib_masri), st_params_calib, common_params) - df_calib_masri$km)^2)
  }
  optimize(calc_error, interval = c(0, 3), tol = 0.002)$minimum
}

# Inputs with both transition matrices scaled: scale_trans_probs() multiplies the moves between NYHA classes
get_scaled_params <- function(mult) {
  arm <- st_params_calib; arm$l_trans_probs_arm <- scale_trans_probs(l_trans_probs_st_tv, mult)
  com <- common_params; com$l_trans_probs_soc <- scale_trans_probs(l_trans_probs_soc_tv, mult)
  list(arm = arm, com = com)
}

fit_slower_transitions <- function(mult) {
  l_par <- get_scaled_params(mult)
  calc_error <- function(theta_val) {
    v_model <- get_survival_curve(theta_val, nrow(df_calib_masri), l_par$arm, l_par$com)
    sum(((v_model - df_calib_masri$km) / df_calib_masri$se)^2)
  }
  optimize(calc_error, interval = c(0, 3), tol = 0.002)$minimum
}

### Step 4. Run the scenarios ###
run_calib_case <- function(label, theta_early, theta_late, split_month, transition_speed = 1) {
  v_theta <- if (is.na(split_month)) theta_early else make_theta(theta_early, theta_late, split_month)
  inp <- base_inputs
  inp$theta <- v_theta
  inp$trans_mult_tx <- transition_speed
  inp$trans_mult_soc <- transition_speed
  r <- run_owsa_case(inp)
  l_par <- get_scaled_params(transition_speed)
  v_surv <- get_survival_curve(v_theta, 11, l_par$arm, l_par$com) # up to 66 months
  data.frame(
    scenario = label, theta_early = theta_early, theta_late = theta_late, split_month = split_month,
    transition_speed = transition_speed,
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
theta_ls <- fit_least_squares()
theta_75 <- fit_slower_transitions(0.75)
theta_50 <- fit_slower_transitions(0.5)
cat("Least squares theta:", theta_ls, "| NYHA moves at 75%:", theta_75, "| at 50%:", theta_50, "\n")

calib_results <- rbind(
  run_calib_case("Base case: Masri curve, one theta", calibrated_theta, NA, NA),
  run_calib_case("THAOS 42-month relative survival, one theta", theta_thaos, NA, NA),
  run_calib_case("Masri curve, two-part theta, split at 18 months", fit_18$par[1], fit_18$par[2], 18),
  run_calib_case("Masri curve, two-part theta, split at 36 months", fit_36$par[1], fit_36$par[2], 36),
  run_calib_case("Masri curve, one theta, unweighted least squares", theta_ls, NA, NA),
  run_calib_case("Masri curve, one theta, NYHA moves at 75%", theta_75, NA, NA, transition_speed = 0.75),
  run_calib_case("Masri curve, one theta, NYHA moves at 50%", theta_50, NA, NA, transition_speed = 0.5)
)
print(calib_results, row.names = FALSE, digits = 5)
write.csv(calib_results, file.path(dir_out, "Scenario_calibration_V.1.5.csv"), row.names = FALSE)
