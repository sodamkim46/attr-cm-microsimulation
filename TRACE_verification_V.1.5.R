# Verification of the V1.5 engine: an expected-value version of the model, written separately from the engine,
# extreme-value tests, Monte Carlo error across seeds, and whether the model reproduces the mortality effect of treatment.
# Results do not feed the base case.

source("TRACE_microsim_V.1.5.R") # base case results; it starts with rm(list = ls()), so it must come first

dir_val <- "Validation_output"
dir.create(dir_val, showWarnings = FALSE)

### Block 1. Expected-value model (no random numbers) ###
# Same rules as the engine, tracked as probabilities: age from a normal quantile grid; NYHA I-III by treatment status
# (on, one and two cycles since stopping, off); NYHA IV (no drug, off); death.
# theta_off applies to discontinued patients after the silencer residual effect, and to NYHA IV.
get_P <- function(l) t(sapply(1:4, function(i) unlist(l[paste0("p_", i, "to", 1:4)])))

run_trace <- function(age, sd, v_nyha, sex, arm, theta_on = calibrated_theta, theta_off = theta_on, p_disc_arm = NULL,
                      n_cycles = num_cycles, residual = TRUE, n_grid = 2000) {
  is_si <- arm == "silencer"
  if (is.null(p_disc_arm)) p_disc_arm <- if (is_si) p_disc_si else p_disc_st
  c_drug <- if (is_si) c_si * (1 - c_discount_si) else c_st * (1 - c_discount_st)
  v_age0 <- age + sd * qnorm((seq_len(n_grid) - 0.5) / n_grid)
  gomp <- params_gompertz[[sex]]
  get_hr <- function(th) 1 + th * (rr_HF * c(1, rr_2v1, rr_3v1, rr_4v1) - 1)

  m_on <- m_c1 <- m_c2 <- m_off <- matrix(0, n_grid, 3)
  m_on[, 1:3] <- rep(v_nyha[1:3], each = n_grid)
  v_n4 <- rep(v_nyha[4], n_grid)
  get_occ <- function() c(colMeans(m_on + m_c1 + m_c2 + m_off), mean(v_n4))

  v_u <- v_states_utilities[1:4]; v_c <- v_states_costs_soc[1:4]
  v_ph <- v_hosp_probs_st[1:4]; v_ch <- v_hosp_costs[1:4]; v_dh <- v_hosp_disutils[1:4]
  v_ly <- v_q <- v_cost <- v_ontx <- v_hosp <- numeric(n_cycles + 1)
  l_fut <- l_future_costs[[sex]]
  calc_unrelated <- function(v_a) approx(l_fut$v_age, l_fut$v_unrelated, xout = v_a, rule = 2)$y * cycle_length
  v_alive <- rep(1, n_grid) # by age in the grid
  v_surv <- numeric(n_cycles)
  occ <- get_occ()
  v_ly[1] <- sum(occ) * cycle_length; v_q[1] <- sum(occ * v_u) * cycle_length
  v_ontx[1] <- mean(rowSums(m_on)); v_cost[1] <- sum(occ * v_c) + c_drug * v_ontx[1] + mean(v_alive * calc_unrelated(v_age0))

  for (t in 1:n_cycles) {
    k <- min(t, 5)
    P_on <- get_P(l_trans_probs_st_tv[[k]]); P_soc <- get_P(l_trans_probs_soc_tv[[k]])
    # Discontinuation at the start of the cycle; the counters move on
    m_new_on <- m_on * (1 - p_disc_arm); m_new_c1 <- m_on * p_disc_arm; m_new_c2 <- m_c1; m_new_off <- m_c2 + m_off
    # Death by state
    v_age <- v_age0 + (t - 1) * cycle_length
    v_H <- (gomp$rate / gomp$shape) * (exp(gomp$shape * (v_age + cycle_length)) - exp(gomp$shape * v_age))
    m_pd_on  <- sapply(get_hr(theta_on),  function(h) 1 - exp(-v_H * h))
    m_pd_off <- sapply(get_hr(theta_off), function(h) 1 - exp(-v_H * h))
    move <- function(m, P, m_pd) {
      out <- matrix(0, n_grid, 4)
      for (i in 1:3) for (j in 1:4) out[, j] <- out[, j] + m[, i] * (1 - m_pd[, i]) * P[i, j]
      out
    }
    keep_on <- is_si && residual
    r_on  <- move(m_new_on, P_on, m_pd_on)
    r_c1  <- move(m_new_c1, if (keep_on) P_on else P_soc, if (keep_on) m_pd_on else m_pd_off)
    r_c2  <- move(m_new_c2, if (keep_on) P_on else P_soc, if (keep_on) m_pd_on else m_pd_off)
    r_off <- move(m_new_off, P_soc, m_pd_off)
    v_n4 <- v_n4 * (1 - m_pd_off[, 4]) + r_on[, 4] + r_c1[, 4] + r_c2[, 4] + r_off[, 4]
    m_on <- r_on[, 1:3]; m_c1 <- r_c1[, 1:3]; m_c2 <- r_c2[, 1:3]; m_off <- r_off[, 1:3]
    # Payoffs at the end of the cycle
    occ <- get_occ()
    v_surv[t] <- sum(occ)
    v_ly[t + 1] <- sum(occ) * cycle_length
    v_q[t + 1] <- sum(occ * v_u) * cycle_length - sum(occ * v_ph * v_dh)
    v_ontx[t + 1] <- mean(rowSums(m_on))
    v_hosp[t + 1] <- sum(occ * v_ph)
    v_alive_new <- rowSums(m_on + m_c1 + m_c2 + m_off) + v_n4
    v_age_end <- v_age + cycle_length
    v_future <- mean(v_alive_new * calc_unrelated(v_age_end)) +
      mean((v_alive - v_alive_new) * (l_fut$v_death[1] + l_fut$v_death[2] * pmin(v_age_end, max(l_fut$v_age))))
    v_alive <- v_alive_new
    v_cost[t + 1] <- sum(occ * v_c) + c_drug * v_ontx[t + 1] + sum(occ * v_ph * v_ch) + v_future
  }
  v_w <- calc_discount_wts(discount_rate_costs, n_cycles, cycle_length)
  lt <- function(v, w) sum((v[-(n_cycles + 1)] * w[-(n_cycles + 1)] + v[-1] * w[-1]) / 2)
  list(Dcosts = lt(v_cost, v_w), Dqalys = lt(v_q, v_w), LYs = lt(v_ly, rep(1, n_cycles + 1)), DOnTx = lt(v_ontx, v_w),
       surv = v_surv, hosp = v_hosp, ly_by_cycle = v_ly)
}

# Sex-weighted
run_trace_mix <- function(age, sd, male, v_nyha, arm, ...) {
  m <- run_trace(age, sd, v_nyha, "male", arm, ...); f <- run_trace(age, sd, v_nyha, "female", arm, ...)
  lapply(setNames(names(m), names(m)), function(x) male * m[[x]] + (1 - male) * f[[x]])
}
summarize_trace <- function(label, age, sd, male, v_nyha, ...) {
  st <- run_trace_mix(age, sd, male, v_nyha, "stabilizer", ...); si <- run_trace_mix(age, sd, male, v_nyha, "silencer", ...)
  inc_c <- si$Dcosts - st$Dcosts; inc_q <- si$Dqalys - st$Dqalys
  v_thr <- sapply(c(1e5, 1.5e5), function(w) (w * inc_q + st$Dcosts - (si$Dcosts - c_si * (1 - c_discount_si) * si$DOnTx)) / si$DOnTx / cycle_length)
  data.frame(case = label, Dcosts_st = st$Dcosts, Dcosts_si = si$Dcosts, Dqalys_st = st$Dqalys, Dqalys_si = si$Dqalys,
             LYs_st = st$LYs, LYs_si = si$LYs, inc_Dcosts = inc_c, inc_Dqalys = inc_q, ICER = inc_c / inc_q,
             list_price_100k = v_thr[1], list_price_150k = v_thr[2], S42 = st$surv[7], S54 = st$surv[9])
}

v_nyha_base  <- p_starting[1:4]
v_nyha_wt    <- c(0.126, 0.536, 0.326, 0.012)
v_nyha_hered <- c(0.174, 0.44, 0.385, 0) / 0.999 # sample() rescales these to sum to 1 in the engine

df_compare <- rbind(
  summarize_trace("Base case, expected value", mean_age, sd_age, prop_male, v_nyha_base),
  summarize_trace("Wild-type, expected value", 79, 8.1, 0.91, v_nyha_wt),
  summarize_trace("Hereditary, expected value", 71, 8.15, 0.71, v_nyha_hered)
)
cat("\nExpected-value model (compare with the engine summary above)\n")
print(df_compare, row.names = FALSE, digits = 6)

### Block 2. Extreme-value tests ###
run_m <- function(arm, ...) run_trace(mean_age, sd_age, v_nyha_base, "male", arm, ...)
set.seed(seed + 2)
v_age_check <- rnorm(2e5, mean_age, sd_age)
df_extreme <- data.frame(
  test = c("theta 0, no discontinuation: 10-year survival (male), expected value", "Gompertz general-population survival, same ages",
           "Equal discontinuation, no residual effect: incremental QALYs (expect 0)",
           "No discontinuation in either arm: incremental QALYs (expect 0)",
           "Equal discontinuation, with residual effect: incremental QALYs (expect > 0)",
           "Life-years at theta 3 minus base theta (expect < 0)", "Alive at the end of the horizon (expect about 0)"),
  value = c(run_m("stabilizer", theta_on = 0, p_disc_arm = 0)$surv[20], calc_expected_surv(v_age_check, params_gompertz$male, 10),
            run_m("silencer", p_disc_arm = 0.01, residual = FALSE)$Dqalys - run_m("stabilizer", p_disc_arm = 0.01)$Dqalys,
            run_m("silencer", p_disc_arm = 0)$Dqalys - run_m("stabilizer", p_disc_arm = 0)$Dqalys,
            run_m("silencer", p_disc_arm = 0.01)$Dqalys - run_m("stabilizer", p_disc_arm = 0.01)$Dqalys,
            run_m("stabilizer", theta_on = 3)$LYs - run_m("stabilizer")$LYs, run_m("stabilizer")$surv[num_cycles])
)
print(df_extreme, row.names = FALSE, digits = 4)
# Every transition row sums to 1 (the engine only warns), and no individual ends up in an undefined state
for (l in c(l_trans_probs_st_tv, l_trans_probs_soc_tv)) stopifnot(all(abs(rowSums(get_P(l)) - 1) < 1e-9))
stopifnot(!anyNA(res_st_male$m_States), !anyNA(res_si_female$m_States))

### Block 3. Monte Carlo error: base case at other seeds, theta fixed ###
run_seed <- function(s) {
  set.seed(s + 1)
  m_m <- cbind(age = rnorm(num_i, mean_age, sd_age)); m_f <- cbind(age = rnorm(num_i, mean_age, sd_age))
  v_st <- sample(v_states_names, size = num_i, replace = TRUE, prob = p_starting)
  run_one <- function(m_features, sex, arm) {
    is_si <- arm == "silencer"
    run_microSimV_hosp_disc(
      v_starting_states = v_st, num_i = num_i, num_cycles = num_cycles, m_indi_features = m_features, v_states_names = v_states_names,
      v_states_costs_arm = if (is_si) v_states_costs_si else v_states_costs_st, v_states_costs_soc = v_states_costs_soc,
      v_states_utilities = v_states_utilities, l_trans_probs_arm = l_trans_probs_st_tv, l_trans_probs_soc = l_trans_probs_soc_tv,
      p_disc = if (is_si) p_disc_si else p_disc_st, arm_name = arm, params_gompertz_sex = params_gompertz[[sex]],
      theta_arm = calibrated_theta, theta_soc = calibrated_theta, v_hosp_probs = v_hosp_probs_st, v_hosp_costs = v_hosp_costs,
      v_hosp_disutils = v_hosp_disutils, discount_rate_costs = discount_rate_costs, discount_rate_QALYs = discount_rate_QALYs,
      cycle_length = cycle_length, starting_seed = s, l_future = l_future_costs[[sex]]
    )[c("mean_Dcosts", "mean_Dqalys", "mean_lys")]
  }
  st_m <- run_one(m_m, "male", "stabilizer"); si_m <- run_one(m_m, "male", "silencer")
  st_f <- run_one(m_f, "female", "stabilizer"); si_f <- run_one(m_f, "female", "silencer")
  wavg <- function(x, y) prop_male * x + prop_female * y
  inc_c <- wavg(si_m$mean_Dcosts - st_m$mean_Dcosts, si_f$mean_Dcosts - st_f$mean_Dcosts)
  inc_q <- wavg(si_m$mean_Dqalys - st_m$mean_Dqalys, si_f$mean_Dqalys - st_f$mean_Dqalys)
  data.frame(seed = s, LYs_st = wavg(st_m$mean_lys, st_f$mean_lys), inc_Dcosts = inc_c, inc_Dqalys = inc_q, ICER = inc_c / inc_q)
}
df_seeds <- do.call(rbind, lapply(c(1, 22, 333, 4444, 55555), run_seed))
df_seeds <- rbind(data.frame(seed = seed, LYs_st = mean_lys_st_overall, inc_Dcosts = mean_Dcosts_si_overall - mean_Dcosts_st_overall,
                             inc_Dqalys = mean_Dqalys_si_overall - mean_Dqalys_st_overall, ICER = ICER_overall), df_seeds)
cat("\nBase case ICER by seed; mean", round(mean(df_seeds$ICER)), "SD", round(sd(df_seeds$ICER)), "\n")
print(df_seeds, row.names = FALSE, digits = 6)

### Block 4. Face validity: cardiovascular hospitalizations per person-year alive, stabilizer ###
# ATTR-ACT reports 0.48 per year with tafamidis and 0.70 with placebo (Maurer et al., NEJM 2018, Figure 2C)
tr_st <- run_trace_mix(mean_age, sd_age, prop_male, v_nyha_base, "stabilizer")
calc_rate <- function(n) sum(tr_st$hosp[2:(n + 1)]) / sum((tr_st$ly_by_cycle[1:n] + tr_st$ly_by_cycle[2:(n + 1)]) / 2)
df_hosp <- data.frame(period = c("Months 0-30", "Lifetime"), events_per_person_year = c(calc_rate(5), calc_rate(num_cycles)))
print(df_hosp, row.names = FALSE, digits = 3)

### Block 5. Does the model reproduce the mortality effect of treatment? ###
# All on treatment vs all stopped at cycle 1, 30-month cumulative hazard ratio.
# ATTR-ACT: 70.5% vs 57.1% alive at 30 months (ratio 0.625) and Cox HR 0.70 (Maurer et al., NEJM 2018).
calc_on_off <- function(label, age, sd, male, v_nyha, th_on = calibrated_theta, th_off = th_on) {
  on  <- run_trace_mix(age, sd, male, v_nyha, "stabilizer", theta_on = th_on, theta_off = th_off, p_disc_arm = 0)
  off <- run_trace_mix(age, sd, male, v_nyha, "stabilizer", theta_on = th_on, theta_off = th_off, p_disc_arm = 1)
  data.frame(case = label, S30_on = on$surv[5], S30_off = off$surv[5], cum_hazard_ratio = log(on$surv[5]) / log(off$surv[5]),
             LYs_on = on$LYs, LYs_off = off$LYs, Dqalys_on = on$Dqalys, Dqalys_off = off$Dqalys)
}
v_nyha_aa <- c(0.084, 0.596, 0.32, 0)
df_on_off <- rbind(calc_on_off("Model population", mean_age, sd_age, prop_male, v_nyha_base),
                   calc_on_off("ATTR-ACT population", 74.5, 7.2, 0.913, v_nyha_aa))

### Block 6. Scenario: stopping treatment raises mortality beyond the NYHA mix ###
# theta_off = r x theta_on. r is set so that the 30-month on/off ratio equals the target; theta_on is refitted to the Masri curve each time.
fit_theta_on <- function(r) {
  optimize(function(th) {
    v_s <- run_trace_mix(mean_age, sd_age, prop_male, v_nyha_base, "stabilizer", theta_on = th, theta_off = r * th, n_cycles = 9)$surv
    sum(((v_s - df_calib_masri$km) / df_calib_masri$se)^2)
  }, interval = c(0, 3), tol = 1e-4)$minimum
}
solve_r <- function(target) {
  uniroot(function(r) {
    th <- fit_theta_on(r)
    calc_on_off("x", mean_age, sd_age, prop_male, v_nyha_base, th, r * th)$cum_hazard_ratio - target
  }, interval = c(1, 12), tol = 0.002)$root
}
theta_trace <- fit_theta_on(1)
v_r <- c(HR070 = solve_r(0.70), HR0625 = solve_r(log(0.705) / log(0.571)))
v_th <- sapply(v_r, fit_theta_on)
df_on_off <- rbind(df_on_off,
                   calc_on_off("Target HR 0.70", mean_age, sd_age, prop_male, v_nyha_base, v_th[1], v_r[1] * v_th[1]),
                   calc_on_off("Target ratio 0.625", mean_age, sd_age, prop_male, v_nyha_base, v_th[2], v_r[2] * v_th[2]))
cat("\nOn vs off treatment\n")
print(df_on_off, row.names = FALSE, digits = 4)

df_scenario <- rbind(
  summarize_trace("Base structure, base theta", mean_age, sd_age, prop_male, v_nyha_base),
  summarize_trace("Base structure, theta fitted in the expected-value model", mean_age, sd_age, prop_male, v_nyha_base, theta_on = theta_trace),
  summarize_trace("On/off target HR 0.70", mean_age, sd_age, prop_male, v_nyha_base, theta_on = v_th[1], theta_off = v_r[1] * v_th[1]),
  summarize_trace("On/off target ratio 0.625", mean_age, sd_age, prop_male, v_nyha_base, theta_on = v_th[2], theta_off = v_r[2] * v_th[2])
)
df_scenario$theta_on  <- c(calibrated_theta, theta_trace, v_th)
df_scenario$theta_off <- c(calibrated_theta, theta_trace, v_r * v_th)
cat("\nScenario results (expected-value model)\n")
print(df_scenario, row.names = FALSE, digits = 6)

### Block 7. Save the results ###
write.csv(df_compare,  file.path(dir_val, "Verification_expected_value_V.1.5.csv"), row.names = FALSE)
write.csv(df_extreme,  file.path(dir_val, "Verification_extreme_values_V.1.5.csv"), row.names = FALSE)
write.csv(df_seeds,    file.path(dir_val, "Verification_seeds_V.1.5.csv"), row.names = FALSE)
write.csv(df_hosp,     file.path(dir_val, "Verification_hospitalization_rate_V.1.5.csv"), row.names = FALSE)
write.csv(df_on_off,   file.path(dir_val, "Verification_on_off_treatment_V.1.5.csv"), row.names = FALSE)
write.csv(df_scenario, file.path(dir_val, "Verification_off_treatment_scenario_V.1.5.csv"), row.names = FALSE)
