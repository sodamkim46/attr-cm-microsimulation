# Scenario analysis: societal perspective (Silencer vs Stabilizer, base case and subgroups)
# Productivity and caregiver costs depend on NYHA class and age, so they are added after each run from the simulated states.
# Clinical outcomes and QALYs do not change. All are in addition to the base case costs.
# Societal A: ICER 2024 values as used in its Excel model, patient productivity at all ages + caregiver time.
#   Patient = share not working + overall work impairment, both applied to every patient, x salary x 2021-2024 factor.
# Societal B: patient productivity only below age 65 + caregiver time at all ages.
#   Patient = (share not working - share not working in NYHA I) + (1 - share not working) x work impairment, the formula behind the
#   costs in Cavusoglu et al., Balkan Med J 2022 (Table 5 inputs reproduce its Table 6 totals). NYHA I is the reference, so its share
#   not working is not attributed to disease, and impairment applies only to working patients.
#   No 2021-2024 factor, because the salary is already in January 2024 USD.

# Step 1. Run the base case model
source("TRACE_microsim_V.1.5.R")

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

### Step 2. Inputs, per cycle, Q2 2026 USD ###
# ICER 2024 Excel model, Scenario Analyses sheet: patient + caregiver per cycle, Q1 2024 USD
v_prod_icer <- c(26218.31734, 38085.36387, 62355.75339, 78627.98802) * infl_factor
# ICER Table E2.9 (Lahoz et al. 2021): caregiver hours per week x $28.49 x 52 weeks, annual
v_caregiver <- c(17481.464, 26814.788, 38370.332, 38370.332) / 2 * infl_factor
# Cavusoglu Table 5 by NYHA class; ICER salary $59,384. NYHA II not working is 41.4%; ICER Table E2.8 prints 41.7% and $29,880
v_not_working <- c(0.343, 0.414, 0.700, 0.958)
v_impairment  <- c(0.208, 0.364, 0.661, 0.916)
v_patient_u65 <- ((v_not_working - v_not_working[1]) + (1 - v_not_working) * v_impairment) * 59384 / 2 * infl_factor
age_retire <- 65

### Step 3. Per-patient mean discounted cost from the simulated states ###
# Same discounting and half-cycle correction as the engine
calc_addon <- function(res, v_age0, type) {
  v_w <- calc_discount_wts(discount_rate_costs, num_cycles, cycle_length)
  v_cycle_cost <- sapply(0:num_cycles, function(t) {
    v_state <- match(res$m_States[, t + 1], v_nyha) # NA = dead
    v_age <- v_age0 + t * cycle_length
    v_alive <- !is.na(v_state)
    v_c <- rep(0, length(v_state))
    if (type == "societal_A") v_c[v_alive] <- v_prod_icer[v_state[v_alive]]
    if (type == "societal_B") v_c[v_alive] <- v_caregiver[v_state[v_alive]] + v_patient_u65[v_state[v_alive]] * (v_age[v_alive] < age_retire)
    mean(v_c)
  })
  v_d <- v_cycle_cost * v_w
  sum((v_d[-length(v_d)] + v_d[-1]) / 2)
}

### Step 4. Run each population ###
run_population <- function(label, age, sd, male, v_start_probs) {
  set.seed(seed + 1)
  m_m <- cbind(age = rnorm(num_i, age, sd)); m_f <- cbind(age = rnorm(num_i, age, sd))
  v_start <- sample(v_states_names, size = num_i, replace = TRUE, prob = v_start_probs)
  summarize <- function(m_features, sex, arm) {
    is_si <- arm == "silencer"
    res <- run_microSimV_hosp_disc(
      v_starting_states = v_start, num_i = num_i, num_cycles = num_cycles, m_indi_features = m_features, v_states_names = v_states_names,
      v_states_costs_arm = if (is_si) v_states_costs_si else v_states_costs_st, v_states_costs_soc = v_states_costs_soc,
      v_states_utilities = v_states_utilities, l_trans_probs_arm = if (is_si) l_trans_probs_si_tv else l_trans_probs_st_tv,
      l_trans_probs_soc = l_trans_probs_soc_tv, p_disc = if (is_si) p_disc_si else p_disc_st, arm_name = arm,
      params_gompertz_sex = params_gompertz[[sex]], theta_arm = calibrated_theta, theta_soc = calibrated_theta,
      v_hosp_probs = if (is_si) v_hosp_probs_si else v_hosp_probs_st, v_hosp_costs = v_hosp_costs, v_hosp_disutils = v_hosp_disutils,
      discount_rate_costs = discount_rate_costs, discount_rate_QALYs = discount_rate_QALYs, cycle_length = cycle_length,
      starting_seed = seed, l_future = l_future_costs[[sex]]
    )
    out <- c(Dcosts = res$mean_Dcosts, Dqalys = res$mean_Dqalys, DOnTx = res$mean_DOnTx,
             societal_A = calc_addon(res, m_features[, "age"], "societal_A"), societal_B = calc_addon(res, m_features[, "age"], "societal_B"))
    rm(res); invisible(gc())
    out
  }
  wavg <- function(x, y) male * x + (1 - male) * y
  st <- wavg(summarize(m_m, "male", "stabilizer"), summarize(m_f, "female", "stabilizer"))
  si <- wavg(summarize(m_m, "male", "silencer"), summarize(m_f, "female", "silencer"))

  c_si_net <- c_si * (1 - c_discount_si)
  inc_q <- si["Dqalys"] - st["Dqalys"]
  make_row <- function(case, add_st, add_si) {
    inc_c <- si["Dcosts"] + add_si - st["Dcosts"] - add_st
    thr <- function(wtp) (wtp * inc_q - (inc_c - c_si_net * si["DOnTx"])) / si["DOnTx"] / cycle_length
    data.frame(population = label, case = case, added_cost_st = add_st, added_cost_si = add_si, inc_Dcosts = inc_c,
               inc_Dqalys = inc_q, ICER = inc_c / inc_q, list_price_100k = thr(1e5), list_price_150k = thr(1.5e5), row.names = NULL)
  }
  rbind(make_row("Base case: healthcare sector", 0, 0),
        make_row("Societal A: ICER productivity values", st["societal_A"], si["societal_A"]),
        make_row("Societal B: productivity below 65 + caregiver", st["societal_B"], si["societal_B"]))
}

cat("\n--- SCENARIO: Societal perspective ---\n")
societal_results <- rbind(
  run_population("Base case", mean_age, sd_age, prop_male, p_starting),
  run_population("Wild-type", 79, 8.1, 0.91, c(0.126, 0.536, 0.326, 0.012, 0)),
  run_population("Hereditary", 71, 8.15, 0.71, c(0.174, 0.44, 0.385, 0, 0))
)
print(societal_results, digits = 6, row.names = FALSE)
stopifnot(abs(societal_results$ICER[1] - ICER_overall) < 1e-6 * ICER_overall) # same base case as the model
write.csv(societal_results, file.path(dir_out, "Scenario_societal_V.1.5.csv"), row.names = FALSE)
cat("\nSaved:", file.path(dir_out, "Scenario_societal_V.1.5.csv"), "\n")
