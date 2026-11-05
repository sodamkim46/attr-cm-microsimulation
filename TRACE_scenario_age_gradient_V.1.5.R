# Scenario analysis: excess mortality that depends less on age (Silencer vs Stabilizer, base case and subgroups)
# Base case: one theta for all ages. External validation (Fumagalli et al., Cardiol Ther 2026) shows that model survival falls
# with age more steeply than observed.
# Scenario: theta differs by age at baseline (<80, 80-85, >85 years).
#   1. Fit one theta to each Fumagalli age stratum (months 30 and 42), in that stratum's own population.
#   2. Keep those ratios and multiply them by one constant c, chosen so that the model population still fits the Masri curve.
#   3. Run the base case and subgroups, each patient keeping the theta of their baseline age group.

# Step 1. Run the base case model
source("TRACE_microsim_V.1.5.R")

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

### Step 2. Population and model run, split by baseline age ###
draw_population <- function(age, sd, v_start_probs) {
  set.seed(seed + 1)
  list(m_male   = cbind(age = rnorm(num_i, age, sd)),
       m_female = cbind(age = rnorm(num_i, age, sd)),
       v_start  = sample(v_states_names, size = num_i, replace = TRUE, prob = v_start_probs))
}

# Run one arm and sex, one age group at a time, and combine by group size
run_by_age <- function(pop, sex, arm, v_theta, n_cycles) {
  is_si <- arm == "silencer"
  m_features <- if (sex == "male") pop$m_male else pop$m_female
  v_group <- cut(m_features[, "age"], breaks = c(-Inf, 80, 85, Inf), labels = FALSE, right = FALSE)
  l_res <- lapply(1:3, function(g) {
    idx <- which(v_group == g)
    if (length(idx) == 0) return(NULL)
    res <- run_microSimV_hosp_disc(
      v_starting_states = pop$v_start[idx], num_i = length(idx), num_cycles = n_cycles,
      m_indi_features = m_features[idx, , drop = FALSE], v_states_names = v_states_names,
      v_states_costs_arm = if (is_si) v_states_costs_si else v_states_costs_st, v_states_costs_soc = v_states_costs_soc,
      v_states_utilities = v_states_utilities, l_trans_probs_arm = if (is_si) l_trans_probs_si_tv else l_trans_probs_st_tv,
      l_trans_probs_soc = l_trans_probs_soc_tv, p_disc = if (is_si) p_disc_si else p_disc_st, arm_name = arm,
      params_gompertz_sex = params_gompertz[[sex]], theta_arm = v_theta[g], theta_soc = v_theta[g],
      v_hosp_probs = if (is_si) v_hosp_probs_si else v_hosp_probs_st, v_hosp_costs = v_hosp_costs, v_hosp_disutils = v_hosp_disutils,
      discount_rate_costs = discount_rate_costs, discount_rate_QALYs = discount_rate_QALYs,
      cycle_length = cycle_length, starting_seed = seed, l_future = l_future_costs[[sex]]
    )
    list(n = length(idx), surv = colMeans(res$m_States != "Death")[-1], Dcosts = res$mean_Dcosts, Dqalys = res$mean_Dqalys,
         lys = res$mean_lys, DOnTx = res$mean_DOnTx)
  })
  l_res <- Filter(Negate(is.null), l_res)
  v_w <- sapply(l_res, function(x) x$n) / sum(sapply(l_res, function(x) x$n))
  combine <- function(key) Reduce(`+`, Map(function(x, w) x[[key]] * w, l_res, v_w))
  list(surv = combine("surv"), Dcosts = combine("Dcosts"), Dqalys = combine("Dqalys"), lys = combine("lys"), DOnTx = combine("DOnTx"))
}

calc_surv <- function(pop, male, v_theta, n_cycles) {
  male * run_by_age(pop, "male", "stabilizer", v_theta, n_cycles)$surv +
    (1 - male) * run_by_age(pop, "female", "stabilizer", v_theta, n_cycles)$surv
}

### Step 3. Theta by age from the Fumagalli strata ###
# Table 1 and Figure 2: age, SD, male, NYHA I-IV, and survival with 95% CI at months 30 and 42
l_strata <- list(
  list(label = "<80",   age = 72.3, sd = 6.12, male = 0.902, nyha = c(0.182, 0.560, 0.249, 0.009, 0),
       km = c(0.770, 0.709), lower = c(0.734, 0.663), upper = c(0.802, 0.750)),
  list(label = "80-85", age = 82.1, sd = 1.62, male = 0.875, nyha = c(0.099, 0.523, 0.370, 0.008, 0),
       km = c(0.662, 0.572), lower = c(0.599, 0.495), upper = c(0.717, 0.642)),
  list(label = ">85",   age = 88.0, sd = 2.22, male = 0.848, nyha = c(0.093, 0.528, 0.346, 0.033, 0),
       km = c(0.563, 0.429), lower = c(0.463, 0.301), upper = c(0.651, 0.551))
)
fit_stratum <- function(s) {
  pop <- draw_population(s$age, s$sd, s$nyha)
  v_se <- (s$upper - s$lower) / (2 * 1.96)
  optimize(function(th) {
    v_s <- calc_surv(pop, s$male, rep(th, 3), 7)[c(5, 7)] # cycles 5 and 7 = months 30 and 42
    sum(((v_s - s$km) / v_se)^2)
  }, interval = c(0, 3), tol = 0.002)$minimum
}
cat("\n--- SCENARIO: Theta by age ---\n")
v_theta_fuma <- sapply(l_strata, fit_stratum)
names(v_theta_fuma) <- sapply(l_strata, `[[`, "label")
cat("Theta fitted to each Fumagalli stratum:", round(v_theta_fuma, 4), "\n")

### Step 4. Constant c so that the model population still fits the Masri curve ###
pop_base <- draw_population(mean_age, sd_age, p_starting)
calc_error_c <- function(c_scale) {
  v_s <- calc_surv(pop_base, prop_male, c_scale * v_theta_fuma, nrow(df_calib_masri))
  sum(((v_s - df_calib_masri$km) / df_calib_masri$se)^2)
}
opt_c <- optimize(calc_error_c, interval = c(0.2, 3), tol = 0.002)
v_theta_age <- opt_c$minimum * v_theta_fuma
cat("c =", round(opt_c$minimum, 4), "| weighted squared error:", round(opt_c$objective, 3), "| theta by age:", round(v_theta_age, 4), "\n")

### Step 5. Results: base case and subgroups, one theta vs theta by age ###
run_population <- function(label, age, sd, male, v_start_probs, v_theta) {
  pop <- draw_population(age, sd, v_start_probs)
  run_mix <- function(arm) {
    m <- run_by_age(pop, "male", arm, v_theta, num_cycles); f <- run_by_age(pop, "female", arm, v_theta, num_cycles)
    lapply(setNames(names(m), names(m)), function(k) male * m[[k]] + (1 - male) * f[[k]])
  }
  st <- run_mix("stabilizer"); si <- run_mix("silencer")
  inc_c <- si$Dcosts - st$Dcosts
  inc_q <- si$Dqalys - st$Dqalys
  c_si_net <- c_si * (1 - c_discount_si)
  thr <- function(wtp) (wtp * inc_q - (inc_c - c_si_net * si$DOnTx)) / si$DOnTx / cycle_length / (1 - c_discount_si)
  data.frame(population = label, theta_lt80 = v_theta[1], theta_80_85 = v_theta[2], theta_gt85 = v_theta[3],
             surv_30 = st$surv[5], surv_42 = st$surv[7], LYs_st = st$lys, LYs_si = si$lys, Dcosts_st = st$Dcosts, Dcosts_si = si$Dcosts,
             inc_Dcosts = inc_c, inc_Dqalys = inc_q, ICER = inc_c / inc_q, list_price_100k = thr(1e5), list_price_150k = thr(1.5e5))
}
l_pops <- list(list("Base case", mean_age, sd_age, prop_male, p_starting),
               list("Wild-type", 79, 8.1, 0.91, c(0.126, 0.536, 0.326, 0.012, 0)),
               list("Hereditary", 71, 8.15, 0.71, c(0.174, 0.44, 0.385, 0, 0)))
age_results <- do.call(rbind, lapply(l_pops, function(p) {
  rbind(cbind(case = "One theta", run_population(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], rep(calibrated_theta, 3))),
        cbind(case = "Theta by age", run_population(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], v_theta_age)))
}))
print(age_results, digits = 6, row.names = FALSE)

### Step 6. Check: survival by Fumagalli stratum with theta by age, and fit to Masri ###
age_check <- do.call(rbind, lapply(seq_along(l_strata), function(i) {
  s <- l_strata[[i]]
  v_s <- unname(calc_surv(draw_population(s$age, s$sd, s$nyha), s$male, rep(v_theta_age[i], 3), 7)[c(5, 7)])
  data.frame(stratum = s$label, month = c(30, 42), km = s$km, lower = s$lower, upper = s$upper, theta = v_theta_age[i], model = v_s,
             in_ci = v_s >= s$lower & v_s <= s$upper)
}))
print(age_check, digits = 4, row.names = FALSE)
df_masri_age <- data.frame(month = df_calib_masri$month, km = df_calib_masri$km, lower = df_calib_masri$lower, upper = df_calib_masri$upper,
                           model = calc_surv(pop_base, prop_male, v_theta_age, nrow(df_calib_masri)))
print(df_masri_age, digits = 4, row.names = FALSE)

write.csv(age_results,  file.path(dir_out, "Scenario_age_gradient_V.1.5.csv"), row.names = FALSE)
write.csv(age_check,    file.path(dir_out, "Scenario_age_gradient_fumagalli_V.1.5.csv"), row.names = FALSE)
write.csv(df_masri_age, file.path(dir_out, "Scenario_age_gradient_masri_V.1.5.csv"), row.names = FALSE)
cat("\nSaved:", file.path(dir_out, c("Scenario_age_gradient_V.1.5.csv", "Scenario_age_gradient_fumagalli_V.1.5.csv", "Scenario_age_gradient_masri_V.1.5.csv")), sep = "\n  ")
