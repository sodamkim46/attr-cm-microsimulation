# Scenario analysis: stopping treatment raises mortality beyond the change in NYHA mix (Silencer vs Stabilizer, base case population)
# Base case: on-treatment and discontinued patients share one theta, so mortality differs only through NYHA class.
# Scenario: theta_off = r x theta_on. r is chosen so that all-on vs all-off in the model reproduces the ATTR-ACT
# tafamidis vs placebo mortality effect at 30 months, and theta_on is refitted so the stabilizer arm still fits the Masri curve.
# Targets (Maurer et al., NEJM 2018): Cox HR 0.70, and the 30-month ratio log(0.705) / log(0.571) = 0.625.
# ATTR-ACT compares treatment with placebo from baseline, not stopping later, so this is an upper bound on the effect of stopping.

# Step 1. Run the base case model and load the SA helper functions
source("TRACE_microsim_V.1.5.R")
source("TRACE_SA_functions_V.1.5.R") # base_inputs, run_owsa_case()

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

### Step 2. Fit theta_on and r ###
# theta_on that fits the Masri curve, for a given r
fit_theta_on <- function(r) {
  optimize(function(th) {
    v_model <- get_survival_curve(th, nrow(df_calib_masri), st_params_calib, common_params, theta_off_val = r * th)
    sum(((v_model - df_calib_masri$km) / df_calib_masri$se)^2)
  }, interval = c(0, 3), tol = 0.002)$minimum
}

# 30-month ratio, all on treatment vs all stopped at cycle 1
calc_on_off <- function(th_on, r) {
  st_on <- st_params_calib; st_on$p_disc <- 0
  st_off <- st_params_calib; st_off$p_disc <- 1
  s_on  <- get_survival_curve(th_on, 5, st_on,  common_params, theta_off_val = r * th_on)[5]
  s_off <- get_survival_curve(th_on, 5, st_off, common_params, theta_off_val = r * th_on)[5]
  unlist(list(s_on = unname(s_on), s_off = unname(s_off), ratio = unname(log(s_on) / log(s_off))))
}

solve_r <- function(target) {
  uniroot(function(r) calc_on_off(fit_theta_on(r), r)["ratio"] - target, interval = c(1.5, 5), tol = 0.02)$root
}

### Step 3. Run each case ###
run_case <- function(label, r) {
  th_on <- if (r == 1) calibrated_theta else fit_theta_on(r)
  inp <- base_inputs
  inp$theta <- th_on
  inp$theta_off <- r * th_on
  res <- run_owsa_case(inp)
  v_on_off <- calc_on_off(th_on, r)
  v_curve <- get_survival_curve(th_on, nrow(df_calib_masri), st_params_calib, common_params, theta_off_val = r * th_on)
  cat(sprintf("  %-34s r = %.3f theta_on = %.4f theta_off = %.4f on/off ratio = %.3f ICER = %.0f\n",
              label, r, th_on, r * th_on, v_on_off["ratio"], res["ICER"]))
  data.frame(scenario = label, r = r, theta_on = th_on, theta_off = r * th_on,
             weighted_error = sum(((v_curve - df_calib_masri$km) / df_calib_masri$se)^2),
             S30_on = v_on_off["s_on"], S30_off = v_on_off["s_off"], on_off_ratio = v_on_off["ratio"],
             surv_42 = v_curve[7], surv_54 = v_curve[9],
             LYs_st = res["LYs_st"], LYs_si = res["LYs_si"], inc_Dcosts = res["inc_Dcosts"], inc_Dqalys = res["inc_Dqalys"],
             ICER = res["ICER"], list_price_100k = res["list_price_100k"], list_price_150k = res["list_price_150k"], row.names = NULL)
}

cat("\n--- SCENARIO: Mortality after stopping treatment ---\n")
r_070  <- solve_r(0.70)
r_0625 <- solve_r(log(0.705) / log(0.571))
off_results <- rbind(
  run_case("Base case (one theta)", 1),
  run_case("On/off target: Cox HR 0.70", r_070),
  run_case("On/off target: 30-month ratio 0.625", r_0625)
)
print(off_results, digits = 6, row.names = FALSE)
write.csv(off_results, file.path(dir_out, "Scenario_off_treatment_V.1.5.csv"), row.names = FALSE)
cat("\nSaved:", file.path(dir_out, "Scenario_off_treatment_V.1.5.csv"), "\n")
