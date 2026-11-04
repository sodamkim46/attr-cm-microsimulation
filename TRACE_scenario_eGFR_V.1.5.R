# TRACE scenario: eGFR pathway on top of the V1.5 base case (acoramidis vs vutrisiran).
# Each patient gets an eGFR trajectory; a >20% decline over one year sets a permanent flag that multiplies
# the death hazard by 1.48 (Ioannou et al. 2025, adjusted HR). The real-world survival used for theta (Masri)
# already includes deaths after eGFR decline, so one offset HR, applied to every patient in both arms like theta,
# is calibrated so that the acoramidis arm's 42-month survival equals the base case. Arm differences then come
# only from the eGFR slopes; clinical efficacy (NYHA transitions) stays equal, as in the base case.

source("TRACE_microsim_V.1.5.R") # base-case inputs, calibrated theta, base-case results

dir_out_egfr <- "Scenario_eGFR_output"
dir.create(dir_out_egfr, showWarnings = FALSE)

### Block 1. Scenario inputs ###

## Baseline eGFR (mL/min/1.73m2): ATTRibute-CM baseline, truncated at the trial's exclusion threshold
egfr_E0_mean <- 63.3
egfr_E0_sd   <- 18.3
egfr_E0_min  <- 15

## Chronic eGFR slopes (mL/min/1.73m2 per year)
# Placebo and acoramidis: ATTRibute-CM linear spline mixed model (Testani et al. 2026)
egfr_beta_placebo    <- -3.48
egfr_beta_acoramidis <- -1.01

# Vutrisiran: no slope reported; back-calculated from the HELIOS-B >=40% eGFR decline proportions
# (placebo 21.2%, vutrisiran 12.7%, overall population; Sheikh et al. 2026), assuming baseline eGFR 60,
# 3 years of follow-up, a normal 3-year change, and the same spread in both arms
egfr_hb_base  <- 60
egfr_hb_years <- 3
egfr_hb_cut   <- -0.4 * egfr_hb_base                                                   # -24
egfr_hb_sd_3y <- (egfr_hb_cut - egfr_beta_placebo * egfr_hb_years) / qnorm(0.212)      # SD of the 3-year change
egfr_beta_vutrisiran <- (egfr_hb_cut - qnorm(0.127) * egfr_hb_sd_3y) / egfr_hb_years

## Between-patient SD of the chronic slope (u_i ~ N(0, sigma^2)); sensitivity 4 and 10
egfr_sigma <- 6

## Measurement error of a single eGFR test. The decline flag is judged on measured values (true + error), as in
## Ioannou et al., while the trajectory itself follows the true slope. Same value used to derive sigma
## (sigma = sqrt(SD_obs^2 - 2 * SD_meas^2)); sensitivity 0 and 7
egfr_sd_meas <- 5

## Acute eGFR dip with acoramidis at day 28 (Testani et al. 2026: 8.5 +/- 0.48). Base case 0: the dip is
## reversible, hemodynamic, and associated with better outcomes, so it is not treated as loss of kidney function
egfr_dip     <- 0
egfr_dip_day <- 28

## Decline flag: >20% fall over one year -> death hazard x HR from then on (Ioannou et al. 2025).
## Base 1.48: adjusted for increases in NT-proBNP and diuretic intensification (markers of disease progression),
## which the model already captures through NYHA transitions. Sensitivity 1.71 (unadjusted, 95% CI 1.43-2.04).
egfr_decline_cut <- 0.20
hr_egfr_decline  <- 1.48

seed_egfr <- 5678 # separate stream for eGFR draws, so the base-case random numbers are unchanged

## Acoramidis drug cost (the stabilizer arm in this scenario), same net discount as the base case
c_ac <- c_acoramidis_WAC * cycle_length
v_states_costs_ac <- v_states_costs_soc + c(rep(c_ac * (1 - c_discount_st), 4), 0)

cat("\n--- eGFR scenario inputs ---\n")
cat("HELIOS-B 3-year change SD:", round(egfr_hb_sd_3y, 2), "-> implied sigma:", round(egfr_hb_sd_3y / egfr_hb_years, 2), "\n")
cat("Slopes (placebo / acoramidis / vutrisiran):", egfr_beta_placebo, "/", egfr_beta_acoramidis, "/", round(egfr_beta_vutrisiran, 2), "\n")
print(rbind(stabilizer_base = v_states_costs_st, acoramidis = v_states_costs_ac))

### Block 2. Patient eGFR trajectories and the decline flag ###

# Patient-level draws, from their own seed: baseline eGFR (floored at egfr_E0_min), slope deviation u_i,
# and one measurement error per test (column 1 = baseline test, column k + 1 = test at the end of year k)
draw_egfr_patients <- function(num_i, sigma = egfr_sigma, sd_meas = egfr_sd_meas, max_years = 40, seed = seed_egfr) {
  set.seed(seed)
  list(E0 = pmax(rnorm(num_i, egfr_E0_mean, egfr_E0_sd), egfr_E0_min),
       u  = rnorm(num_i, 0, sigma),
       e  = matrix(rnorm(num_i * (max_years + 1), 0, sd_meas), nrow = num_i))
}

# Starting eGFR state. ref = measured value at the reference test: day 28 for acoramidis (the dip applies
# from day 28), baseline for other arms
init_egfr <- function(l_draws, beta_arm, dip = 0, dip_arm = FALSE) {
  v_egfr     <- l_draws$E0 - dip * dip_arm
  v_ref_true <- if (dip_arm) v_egfr + (beta_arm + l_draws$u) * egfr_dip_day / 365 else l_draws$E0
  list(egfr = v_egfr, ref = v_ref_true + l_draws$e[, 1], flag = rep(0, length(v_egfr)),
       e = l_draws$e, k = 1, # k = tests done so far
       cap = l_draws$E0)     # baseline before any dip: upper limit when cap_rise = TRUE
}

# Advance true eGFR by one cycle along each patient's current slope (floored at 0; with cap_rise, also no rise above
# baseline); at the end of each year, test (true + measurement error), flag a >20% fall since the last test
# (permanent), and keep this test as the new reference
update_egfr <- function(l_egfr, t, v_slope, cycle_length, cap_rise = FALSE) {
  l_egfr$egfr <- pmax(l_egfr$egfr + v_slope * cycle_length, 0)
  if (cap_rise) l_egfr$egfr <- pmin(l_egfr$egfr, l_egfr$cap)
  if ((t * cycle_length) %% 1 == 0) {
    l_egfr$k <- l_egfr$k + 1
    v_measured <- l_egfr$egfr + l_egfr$e[, l_egfr$k]
    l_egfr$flag[v_measured < (1 - egfr_decline_cut) * l_egfr$ref] <- 1
    l_egfr$ref <- v_measured
  }
  l_egfr
}

# Check: cumulative % flagged by year when everyone stays alive and on treatment (eGFR process only)
flag_by_year <- function(beta_arm, dip = 0, dip_arm = FALSE, years = 10, sigma = egfr_sigma, sd_meas = egfr_sd_meas,
                         cap_rise = FALSE) {
  l_draws <- draw_egfr_patients(num_i, sigma, sd_meas)
  l_egfr  <- init_egfr(l_draws, beta_arm, dip, dip_arm)
  v_slope <- beta_arm + l_draws$u
  v_out   <- numeric(years)
  for (t in 1:(years / cycle_length)) {
    l_egfr <- update_egfr(l_egfr, t, v_slope, cycle_length, cap_rise)
    if ((t * cycle_length) %% 1 == 0) v_out[t * cycle_length] <- mean(l_egfr$flag)
  }
  v_out
}

df_flag <- data.frame(
  year                  = 1:10,
  placebo               = flag_by_year(egfr_beta_placebo),
  vutrisiran            = flag_by_year(egfr_beta_vutrisiran),
  acoramidis            = flag_by_year(egfr_beta_acoramidis, dip = egfr_dip, dip_arm = TRUE),
  acoramidis_dip8.5     = flag_by_year(egfr_beta_acoramidis, dip = 8.5,      dip_arm = TRUE),
  placebo_no_meas_error = flag_by_year(egfr_beta_placebo, sd_meas = 0),
  placebo_no_rise       = flag_by_year(egfr_beta_placebo, cap_rise = TRUE),
  vutrisiran_no_rise    = flag_by_year(egfr_beta_vutrisiran, cap_rise = TRUE),
  acoramidis_no_rise    = flag_by_year(egfr_beta_acoramidis, dip = egfr_dip, dip_arm = TRUE, cap_rise = TRUE)
)
cat("\n--- Cumulative % with eGFR decline flag (no death, no discontinuation) ---\n")
print(cbind(year = df_flag$year, round(100 * df_flag[, -1], 1)))
write.csv(df_flag, file.path(dir_out_egfr, "eGFR_flag_by_year.csv"), row.names = FALSE)

# Normal slopes give some patients a rising eGFR for life (implausible at these ages). Capping eGFR at baseline
# (sensitivity analysis) removes them but makes every arm's mean slope steeper than in the trials:
u_chk <- draw_egfr_patients(num_i)$u
cat("Share with a rising eGFR (placebo / vutrisiran / acoramidis):",
    round(100 * sapply(c(egfr_beta_placebo, egfr_beta_vutrisiran, egfr_beta_acoramidis), function(b) mean(b + u_chk > 0)), 1),
    "%\nMean slope with eGFR capped at baseline:",
    round(sapply(c(egfr_beta_placebo, egfr_beta_vutrisiran, egfr_beta_acoramidis), function(b) mean(pmin(b + u_chk, 0))), 2), "\n")

# Check against Ioannou et al. 2025 (measured 1-year change: median -5, IQR -12 to 1; 24.0% with >20% decline).
# The cohort's mean change was about -5/year, so the check uses beta = -5 rather than the placebo slope.
l_draws_chk <- draw_egfr_patients(num_i)
v_chg_meas  <- (l_draws_chk$E0 - 5 + l_draws_chk$u + l_draws_chk$e[, 2]) - (l_draws_chk$E0 + l_draws_chk$e[, 1])
cat("Measured 1-year change with beta = -5: median", round(median(v_chg_meas), 1),
    "IQR", round(quantile(v_chg_meas, c(0.25, 0.75)), 1),
    "| >20% decline:", round(100 * mean(v_chg_meas < -0.2 * (l_draws_chk$E0 + l_draws_chk$e[, 1])), 1), "%\n")

# Example trajectories: 30 patients per arm (same patients in every arm)
library(ggplot2)
l_draws_plot <- draw_egfr_patients(30)
df_traj <- do.call(rbind, lapply(c("placebo", "vutrisiran", "acoramidis"), function(arm) {
  beta <- c(placebo = egfr_beta_placebo, vutrisiran = egfr_beta_vutrisiran, acoramidis = egfr_beta_acoramidis)[[arm]]
  v_years <- seq(0, 10, by = cycle_length)
  data.frame(arm = arm, id = rep(1:30, each = length(v_years)), year = v_years,
             egfr = pmax(rep(l_draws_plot$E0, each = length(v_years)) +
                           rep(beta + l_draws_plot$u, each = length(v_years)) * v_years, 0))
}))
df_traj$arm <- factor(df_traj$arm, levels = c("placebo", "vutrisiran", "acoramidis"))
p_traj <- ggplot(df_traj, aes(year, egfr, group = id)) +
  geom_line(alpha = 0.4, colour = "#2a78d6") +
  facet_wrap(~arm) +
  labs(x = "Years", y = "eGFR (mL/min/1.73m2)", title = "Example eGFR trajectories (30 patients, same patients in each arm)") +
  theme_minimal()
ggsave(file.path(dir_out_egfr, "eGFR_trajectories.png"), p_traj, width = 9, height = 3.5, dpi = 150)

### Block 3. Microsimulation engine with the eGFR pathway ###

# Same steps and random-number order as run_microSimV_hosp_disc() (V1.5), plus:
#  - death hazard multiplier for every patient: hr_offset x hr_decline^flag
#  - eGFR moves along the arm slope on treatment and the placebo slope after (effective) discontinuation;
#    a dip, if any, is recovered at discontinuation; yearly tests update the decline flag
# eGFR draws come in pre-drawn (l_egfr_draws), so the engine's own random numbers are unchanged.
# With egfr_on = FALSE the run is identical to run_microSimV_hosp_disc().
run_microSimV_egfr <- function(
    v_starting_states,
    num_i,
    num_cycles,
    m_indi_features,
    v_states_names,
    v_states_costs_arm,
    v_states_costs_soc,
    v_states_utilities,
    l_trans_probs_arm,
    l_trans_probs_soc,
    p_disc,
    arm_name,
    params_gompertz_sex,
    theta_arm = 1,
    theta_soc = 1,
    v_hosp_probs,
    v_hosp_costs,
    v_hosp_disutils,
    discount_rate_costs,
    discount_rate_QALYs,
    cycle_length = 0.5,
    starting_seed = 1,
    waning_years = 1,
    disc_max_cycle = Inf,
    egfr_on = TRUE,          # FALSE = V1.5 engine
    l_egfr_draws,            # draw_egfr_patients(num_i, ...)
    beta_arm,                # chronic eGFR slope on treatment
    beta_soc = egfr_beta_placebo, # chronic eGFR slope after discontinuation
    dip = 0,                 # acute dip (acoramidis)
    hr_offset = 1,           # offset for decline deaths already in the real-world survival (calibrated; all patients)
    hr_decline = hr_egfr_decline,
    cap_rise = FALSE) {      # TRUE = eGFR cannot rise above baseline (sensitivity analysis)

  m_States <- m_Costs <- m_Effs <- m_LYs <- m_OnTx <- m_Flag <- matrix(
    nrow = num_i,
    ncol = num_cycles + 1,
    dimnames = list(paste("ind",   1:num_i,    sep = "_"),
                    paste("cycle", 0:num_cycles, sep = "_"))
  )
  v_discontinued <- as.numeric(v_starting_states == "NYHA4")
  v_cycles_since_disc <- rep(0, num_i)

  # eGFR state (patients starting in NYHA4 are untreated: no dip)
  l_egfr <- init_egfr(l_egfr_draws, beta_arm, dip * (1 - v_discontinued), dip_arm = dip > 0)
  v_dip_left <- dip * (1 - v_discontinued) # dip still in place (recovered at discontinuation)

  set.seed(starting_seed)
  m_States[, 1] <- v_starting_states

  m_Costs[, 1] <- calc_costsV_disc(
    v_occupied_state   = m_States[, 1],
    v_states_costs_arm = v_states_costs_arm,
    v_states_costs_soc = v_states_costs_soc,
    v_discontinued     = v_discontinued
  )
  m_Effs[, 1] <- calc_effsV(
    v_occupied_state   = m_States[, 1],
    v_states_utilities = v_states_utilities,
    cycle_length       = cycle_length
  )
  m_LYs[, 1]  <- (m_States[, 1] != "Death") * cycle_length
  m_OnTx[, 1] <- (m_States[, 1] != "Death") * (1 - v_discontinued)
  m_Flag[, 1] <- l_egfr$flag

  for (t in 1:num_cycles) {

    current_cycle_index <- min(t, length(l_trans_probs_arm))
    l_trans_probs_arm_t <- l_trans_probs_arm[[current_cycle_index]]
    l_trans_probs_soc_t <- l_trans_probs_soc[[current_cycle_index]]

    u_disc  <- runif(num_i)
    u_trans <- runif(num_i)
    u_hosp  <- runif(num_i)

    # 1. Discontinuation
    eligible_to_disc <- (v_discontinued == 0) & (m_States[, t] != "Death") & (t <= disc_max_cycle)
    v_discontinued[eligible_to_disc & (u_disc < p_disc)] <- 1
    v_cycles_since_disc[v_discontinued == 1] <- v_cycles_since_disc[v_discontinued == 1] + 1

    v_effective_disc <- v_discontinued
    if (arm_name == "silencer") {
      waning_period_cycles <- waning_years / cycle_length
      indices_in_waning <- which(v_discontinued == 1 & v_cycles_since_disc <= waning_period_cycles)
      if (length(indices_in_waning) > 0) {
        v_effective_disc[indices_in_waning] <- 0
      }
    }

    # 2. Death hazard multiplier from the eGFR pathway
    v_hr_death <- if (egfr_on) hr_offset * hr_decline^l_egfr$flag else 1

    # 3. Transition probabilities
    m_trans_probs <- update_probsV_disc(
      v_states_names      = v_states_names,
      v_occupied_state    = m_States[, t],
      l_trans_probs_arm   = l_trans_probs_arm_t,
      l_trans_probs_soc   = l_trans_probs_soc_t,
      v_discontinued      = v_effective_disc,
      m_indi_features     = m_indi_features,
      theta_arm           = theta_arm,
      theta_soc           = theta_soc,
      params_gompertz_sex = params_gompertz_sex,
      cycle_length        = cycle_length,
      v_hr_death          = v_hr_death
    )

    # 4. Next state
    m_States[, t + 1] <- sampleV(
      m_trans_probs  = m_trans_probs,
      v_states_names = v_states_names,
      v_rand_values  = u_trans
    )

    # 5. Hospitalization in the new state
    l_hosp_results <- hospV(
      v_occupied_state = m_States[, t + 1],
      v_hosp_probs     = v_hosp_probs,
      v_hosp_costs     = v_hosp_costs,
      v_hosp_disutils  = v_hosp_disutils,
      v_rand_values    = u_hosp
    )

    # 6. eGFR over this cycle: arm slope on treatment, placebo slope after discontinuation; dip recovered
    #    at discontinuation; yearly test at the end of each year (flag applies from the next cycle)
    if (egfr_on) {
      v_recover <- (v_effective_disc == 1) & (v_dip_left > 0)
      l_egfr$egfr[v_recover] <- l_egfr$egfr[v_recover] + v_dip_left[v_recover]
      v_dip_left[v_recover] <- 0
      v_slope <- ifelse(v_effective_disc == 0, beta_arm, beta_soc) + l_egfr_draws$u
      l_egfr  <- update_egfr(l_egfr, t, v_slope, cycle_length, cap_rise)
    }

    # 7. Entering NYHA4 ends treatment
    v_discontinued[m_States[, t + 1] == "NYHA4"] <- 1

    # 8. Age by one cycle
    m_indi_features[, "age"] <- m_indi_features[, "age"] + cycle_length

    # 9. Payoffs for cycle t + 1
    m_Costs[, t + 1] <- calc_costsV_disc(
      v_occupied_state   = m_States[, t + 1],
      v_states_costs_arm = v_states_costs_arm,
      v_states_costs_soc = v_states_costs_soc,
      v_discontinued     = v_discontinued
    ) + l_hosp_results$costs
    m_Effs[, t + 1] <- calc_effsV(
      v_occupied_state   = m_States[, t + 1],
      v_states_utilities = v_states_utilities,
      cycle_length       = cycle_length
    ) - l_hosp_results$disutils
    m_LYs[, t + 1]  <- (m_States[, t + 1] != "Death") * cycle_length
    m_OnTx[, t + 1] <- (m_States[, t + 1] != "Death") * (1 - v_discontinued)
    m_Flag[, t + 1] <- l_egfr$flag
  }

  # Totals with the same half-cycle correction and discounting as V1.5
  v_c_dsc_wts <- calc_discount_wts(discount_rate_costs, num_cycles, cycle_length)
  v_e_dsc_wts <- calc_discount_wts(discount_rate_QALYs, num_cycles, cycle_length)
  lt_sum <- function(m) rowSums((m[, -ncol(m)] + m[, -1]) / 2)

  v_total_costs  <- lt_sum(m_Costs)
  v_total_qalys  <- lt_sum(m_Effs)
  v_total_lys    <- lt_sum(m_LYs)
  v_total_Dcosts <- lt_sum(t(t(m_Costs) * v_c_dsc_wts))
  v_total_Dqalys <- lt_sum(t(t(m_Effs)  * v_e_dsc_wts))
  v_total_Dlys   <- lt_sum(t(t(m_LYs)   * v_e_dsc_wts))
  mean_DOnTx     <- mean(lt_sum(t(t(m_OnTx) * v_c_dsc_wts)))

  v_alive <- colSums(m_States != "Death")
  list(
    m_States       = m_States,
    m_Costs        = m_Costs,
    m_Effs         = m_Effs,
    m_LYs          = m_LYs,
    v_total_costs  = v_total_costs,
    v_total_qalys  = v_total_qalys,
    v_total_lys    = v_total_lys,
    v_total_Dcosts = v_total_Dcosts,
    v_total_Dqalys = v_total_Dqalys,
    v_total_Dlys   = v_total_Dlys,
    mean_costs     = mean(v_total_costs),
    mean_qalys     = mean(v_total_qalys),
    mean_lys       = mean(v_total_lys),
    mean_Dcosts    = mean(v_total_Dcosts),
    mean_Dqalys    = mean(v_total_Dqalys),
    mean_Dlys      = mean(v_total_Dlys),
    mean_DOnTx     = mean_DOnTx,
    v_surv         = v_alive / num_i,                                                   # survival by cycle
    v_flag_alive   = colSums(m_Flag * (m_States != "Death")) / pmax(v_alive, 1)        # % flagged among the living
  )
}

# Regression test: with the eGFR pathway off, or on with no effect (both HRs = 1), every output must equal the
# V1.5 base-case run for the same arm and sex (res_st_male / res_si_male from TRACE_microsim_V.1.5.R)
run_egfr_arm <- function(m_features, sex, arm, n_cycles = num_cycles, costs_arm = NULL, ...) {
  is_si <- arm == "silencer"
  run_microSimV_egfr(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = n_cycles,
    m_indi_features     = m_features,
    v_states_names      = v_states_names,
    v_states_costs_arm  = if (!is.null(costs_arm)) costs_arm else if (is_si) v_states_costs_si else v_states_costs_st,
    v_states_costs_soc  = v_states_costs_soc,
    v_states_utilities  = v_states_utilities,
    l_trans_probs_arm   = if (is_si) l_trans_probs_si_tv else l_trans_probs_st_tv,
    l_trans_probs_soc   = l_trans_probs_soc_tv,
    p_disc              = if (is_si) p_disc_si else p_disc_st,
    arm_name            = arm,
    params_gompertz_sex = params_gompertz[[sex]],
    theta_arm           = calibrated_theta,
    theta_soc           = calibrated_theta,
    v_hosp_probs        = if (is_si) v_hosp_probs_si else v_hosp_probs_st,
    v_hosp_costs        = v_hosp_costs,
    v_hosp_disutils     = v_hosp_disutils,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed,
    ...
  )
}

same_as_base <- function(res, res_base) {
  all(identical(res$m_States, res_base$m_States),
      identical(res$v_total_Dcosts, res_base$v_total_Dcosts),
      identical(res$v_total_Dqalys, res_base$v_total_Dqalys),
      identical(res$v_total_Dlys,   res_base$v_total_Dlys),
      identical(res$mean_DOnTx,     res_base$mean_DOnTx))
}

l_egfr_draws_male <- draw_egfr_patients(num_i)
cat("\n--- Block 3 regression test (males) ---\n")
for (arm in c("stabilizer", "silencer")) {
  res_base <- if (arm == "stabilizer") res_st_male else res_si_male
  beta     <- if (arm == "stabilizer") egfr_beta_acoramidis else egfr_beta_vutrisiran
  res_off  <- run_egfr_arm(m_indi_features_male, "male", arm, egfr_on = FALSE, l_egfr_draws = l_egfr_draws_male, beta_arm = beta)
  res_nohr <- run_egfr_arm(m_indi_features_male, "male", arm, egfr_on = TRUE,  l_egfr_draws = l_egfr_draws_male, beta_arm = beta, hr_decline = 1)
  res_on   <- run_egfr_arm(m_indi_features_male, "male", arm, egfr_on = TRUE,  l_egfr_draws = l_egfr_draws_male, beta_arm = beta)
  cat(arm, "| eGFR off identical:", same_as_base(res_off, res_base),
      "| eGFR on, HR = 1 identical:", same_as_base(res_nohr, res_base), "\n")
  cat("   eGFR on (no offset): 42-month survival", round(res_on$v_surv[8], 4),
      "vs base", round(mean(res_base$m_States[, 8] != "Death"), 4),
      "| % flagged among alive at 42 months:", round(100 * res_on$v_flag_alive[8], 1), "\n")
}

### Block 4. Offset HR: acoramidis 42-month survival with the eGFR pathway = base case; same offset for both arms ###

# eGFR draws for each sex-specific cohort (own seeds)
get_egfr_draws <- function(sigma = egfr_sigma, sd_meas = egfr_sd_meas) {
  list(male   = draw_egfr_patients(num_i, sigma, sd_meas, seed = seed_egfr),
       female = draw_egfr_patients(num_i, sigma, sd_meas, seed = seed_egfr + 1))
}

# Sex-weighted survival by cycle
wsurv <- function(res_m, res_f) prop_male * res_m$v_surv + prop_female * res_f$v_surv
surv_by_cycle <- function(res) colMeans(res$m_States != "Death")

# Base-case survival by cycle and arm; the target is the stabilizer arm at the 42-month anchor point
# (42 months = calib_cycles), the arm theta is calibrated on
v_surv_base <- list(
  stabilizer = prop_male * surv_by_cycle(res_st_male) + prop_female * surv_by_cycle(res_st_female),
  silencer   = prop_male * surv_by_cycle(res_si_male) + prop_female * surv_by_cycle(res_si_female)
)

# Run both sexes for one arm with the eGFR pathway on
run_egfr_both <- function(arm, beta_arm, dip, hr_offset, l_draws, n_cycles = num_cycles, costs_arm = NULL, ...) {
  list(male   = run_egfr_arm(m_indi_features_male,   "male",   arm, n_cycles = n_cycles, costs_arm = costs_arm,
                             l_egfr_draws = l_draws$male,   beta_arm = beta_arm, dip = dip, hr_offset = hr_offset, ...),
       female = run_egfr_arm(m_indi_features_female, "female", arm, n_cycles = n_cycles, costs_arm = costs_arm,
                             l_egfr_draws = l_draws$female, beta_arm = beta_arm, dip = dip, hr_offset = hr_offset, ...))
}

# Root of (acoramidis 42-month survival with eGFR) - (base-case stabilizer 42-month survival);
# only the first calib_cycles are simulated
calibrate_hr_offset <- function(dip, l_draws, beta_arm = egfr_beta_acoramidis, ...) {
  f <- function(hr) {
    l_res <- run_egfr_both("stabilizer", beta_arm, dip, hr, l_draws, n_cycles = calib_cycles, ...)
    wsurv(l_res$male, l_res$female)[calib_cycles + 1] - v_surv_base$stabilizer[calib_cycles + 1]
  }
  uniroot(f, interval = c(0.5, 1.2), tol = 1e-5)$root
}

l_draws_base <- get_egfr_draws()
cat("\n--- Block 4: calibrating the offset HR (acoramidis 42-month survival = base case) ---\n")
hr_offset <- calibrate_hr_offset(egfr_dip, l_draws_base)
cat("Offset HR (both arms, all patients):", round(hr_offset, 4), "\n")

# Check: survival by cycle up to 42 months, eGFR runs vs base case (vutrisiran is not forced to match)
l_chk_ac <- run_egfr_both("stabilizer", egfr_beta_acoramidis, egfr_dip, hr_offset, l_draws_base, n_cycles = calib_cycles)
l_chk_vu <- run_egfr_both("silencer",   egfr_beta_vutrisiran, 0,        hr_offset, l_draws_base, n_cycles = calib_cycles)
df_calib_check <- data.frame(
  month          = (0:calib_cycles) * cycle_length * 12,
  base_st        = v_surv_base$stabilizer[1:(calib_cycles + 1)],
  egfr_ac        = wsurv(l_chk_ac$male, l_chk_ac$female),
  base_si        = v_surv_base$silencer[1:(calib_cycles + 1)],
  egfr_vu        = wsurv(l_chk_vu$male, l_chk_vu$female),
  flag_alive_ac  = prop_male * l_chk_ac$male$v_flag_alive + prop_female * l_chk_ac$female$v_flag_alive,
  flag_alive_vu  = prop_male * l_chk_vu$male$v_flag_alive + prop_female * l_chk_vu$female$v_flag_alive
)
print(round(df_calib_check, 4))
write.csv(df_calib_check, file.path(dir_out_egfr, "eGFR_calibration_check.csv"), row.names = FALSE)

### Block 5. Lifetime results in three steps (vutrisiran vs acoramidis) ###
# (i)   base case: stabilizer = mean tafamidis/acoramidis price, no eGFR pathway
# (ii)  acoramidis price, no eGFR pathway            -> (ii) - (i)  = price effect
# (iii) acoramidis price, eGFR pathway on (calibrated) -> (iii) - (ii) = eGFR pathway effect

summarize_step <- function(label, l_st, l_si) {
  w <- function(l, x) prop_male * l$male[[x]] + prop_female * l$female[[x]]
  out <- data.frame(
    step = label,
    cost_st = w(l_st, "mean_Dcosts"), qaly_st = w(l_st, "mean_Dqalys"), ly_st = w(l_st, "mean_Dlys"), ly_undisc_st = w(l_st, "mean_lys"),
    cost_si = w(l_si, "mean_Dcosts"), qaly_si = w(l_si, "mean_Dqalys"), ly_si = w(l_si, "mean_Dlys"), ly_undisc_si = w(l_si, "mean_lys")
  )
  out$inc_cost <- out$cost_si - out$cost_st
  out$inc_qaly <- out$qaly_si - out$qaly_st
  out$inc_ly   <- out$ly_si   - out$ly_st
  out$icer     <- out$inc_cost / out$inc_qaly
  out
}

cat("\n--- Block 5: lifetime runs ---\n")
l_base_st <- list(male = res_st_male, female = res_st_female)
l_base_si <- list(male = res_si_male, female = res_si_female)

# (ii) acoramidis price only (survival and QALYs as in the base case)
l_ac_noegfr <- list(
  male   = run_egfr_arm(m_indi_features_male,   "male",   "stabilizer", costs_arm = v_states_costs_ac, egfr_on = FALSE,
                        l_egfr_draws = l_draws_base$male,   beta_arm = egfr_beta_acoramidis),
  female = run_egfr_arm(m_indi_features_female, "female", "stabilizer", costs_arm = v_states_costs_ac, egfr_on = FALSE,
                        l_egfr_draws = l_draws_base$female, beta_arm = egfr_beta_acoramidis))

# (iii) eGFR pathway on, calibrated offset HR
l_ac_egfr <- run_egfr_both("stabilizer", egfr_beta_acoramidis, egfr_dip, hr_offset, l_draws_base, costs_arm = v_states_costs_ac)
l_vu_egfr <- run_egfr_both("silencer",   egfr_beta_vutrisiran, 0,        hr_offset, l_draws_base)

df_steps <- rbind(
  summarize_step("(i) base case",                l_base_st,   l_base_si),
  summarize_step("(ii) acoramidis price",        l_ac_noegfr, l_base_si),
  summarize_step("(iii) + eGFR pathway",         l_ac_egfr,   l_vu_egfr)
)
print(format(df_steps[, c("step", "cost_st", "cost_si", "qaly_st", "qaly_si", "ly_undisc_st", "ly_undisc_si", "inc_cost", "inc_qaly", "icer")],
             big.mark = ",", digits = 6), row.names = FALSE)
write.csv(df_steps, file.path(dir_out_egfr, "eGFR_scenario_results.csv"), row.names = FALSE)

# Survival and % flagged among the living, lifetime
v_years <- (0:num_cycles) * cycle_length
df_surv_egfr <- rbind(
  data.frame(year = v_years, arm = "Acoramidis", model = "eGFR pathway", surv = wsurv(l_ac_egfr$male, l_ac_egfr$female)),
  data.frame(year = v_years, arm = "Vutrisiran", model = "eGFR pathway", surv = wsurv(l_vu_egfr$male, l_vu_egfr$female)),
  data.frame(year = v_years, arm = "Acoramidis", model = "Base case",    surv = v_surv_base$stabilizer),
  data.frame(year = v_years, arm = "Vutrisiran", model = "Base case",    surv = v_surv_base$silencer)
)
write.csv(df_surv_egfr, file.path(dir_out_egfr, "eGFR_scenario_survival.csv"), row.names = FALSE)
cat("Survival at 5 / 10 / 15 years:\n")
print(reshape(df_surv_egfr[df_surv_egfr$year %in% c(5, 10, 15), ], idvar = c("arm", "model"), timevar = "year", direction = "wide"), row.names = FALSE)

p_surv_egfr <- ggplot(df_surv_egfr, aes(year, surv, colour = arm, linetype = model)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = calib_months / 12, colour = "grey60", linetype = "dotted") +
  scale_colour_manual(values = c(Acoramidis = "#2a78d6", Vutrisiran = "#eb6834")) +
  coord_cartesian(xlim = c(0, 20)) +
  labs(x = "Years", y = "Survival", colour = NULL, linetype = NULL,
       title = "Survival with and without the eGFR pathway (dotted line: 42-month calibration point)") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(dir_out_egfr, "eGFR_scenario_survival.png"), p_surv_egfr, width = 8, height = 5, dpi = 150)

### Block 6. Sensitivity analyses (offset HR recalibrated in every case) ###

run_egfr_case <- function(label, sigma = egfr_sigma, sd_meas = egfr_sd_meas, beta_vu = egfr_beta_vutrisiran,
                          dip = egfr_dip, hr_decline = hr_egfr_decline, cap_rise = FALSE) {
  l_draws <- get_egfr_draws(sigma, sd_meas)
  hr_off  <- calibrate_hr_offset(dip, l_draws, hr_decline = hr_decline, cap_rise = cap_rise)
  l_ac <- run_egfr_both("stabilizer", egfr_beta_acoramidis, dip, hr_off, l_draws, costs_arm = v_states_costs_ac,
                        hr_decline = hr_decline, cap_rise = cap_rise)
  l_vu <- run_egfr_both("silencer",   beta_vu,              0,   hr_off, l_draws, hr_decline = hr_decline, cap_rise = cap_rise)
  cat("done:", label, "| offset HR", round(hr_off, 4), "\n")
  cbind(summarize_step(label, l_ac, l_vu), hr_offset = hr_off)
}

cat("\n--- Block 6: sensitivity analyses ---\n")
df_sa_egfr <- rbind(
  cbind(df_steps[2, ], hr_offset = NA),                                   # (ii) no eGFR pathway (reference)
  cbind(df_steps[3, ], hr_offset = hr_offset),                            # (iii) eGFR base
  run_egfr_case("sigma 4",                          sigma = 4),
  run_egfr_case("sigma 10",                         sigma = 10),
  run_egfr_case("measurement error SD 0",           sd_meas = 0),
  run_egfr_case("measurement error SD 7",           sd_meas = 7),
  run_egfr_case("vutrisiran slope -3.48 (placebo)", beta_vu = egfr_beta_placebo),
  run_egfr_case("vutrisiran slope -1.01 (= acoramidis)", beta_vu = egfr_beta_acoramidis),
  run_egfr_case("acoramidis dip 8.5",               dip = 8.5),
  run_egfr_case("decline HR 1.71 (unadjusted)",     hr_decline = 1.71),
  run_egfr_case("no eGFR rise above baseline",      cap_rise = TRUE)
)
# eGFR pathway effect = change from (ii), where clinical efficacy and prices are the same
df_sa_egfr$d_inc_qaly_vs_ii <- df_sa_egfr$inc_qaly - df_sa_egfr$inc_qaly[1]
df_sa_egfr$d_icer_vs_ii     <- df_sa_egfr$icer     - df_sa_egfr$icer[1]
print(format(df_sa_egfr[, c("step", "hr_offset", "inc_cost", "inc_qaly", "inc_ly", "icer", "d_inc_qaly_vs_ii", "d_icer_vs_ii")],
             big.mark = ",", digits = 5), row.names = FALSE)
write.csv(df_sa_egfr, file.path(dir_out_egfr, "eGFR_scenario_sensitivity.csv"), row.names = FALSE)
