# TRACE: individual-level microsimulation of TTR silencer (vutrisiran) vs TTR stabilizers (tafamidis/acoramidis)
# in ATTR-CM. Health states: NYHA I-IV and death, 6-month cycles, lifetime horizon, US healthcare sector perspective.
# Features: time-varying NYHA transitions, treatment discontinuation with a residual effect for the silencer,
# cardiovascular hospitalization events, age/sex-specific Gompertz background mortality, and an excess-hazard
# factor (theta) calibrated to real-world survival (THAOS). Outputs: costs, QALYs, life-years, ICER,
# threshold price, subgroups, and a cost breakdown.

rm(list = ls())

### Model functions ###

# Transition probability matrix (one row per individual) for one set of NYHA transition probabilities.
# Death probability per state = background Gompertz hazard for the individual's age x state-specific HR.
get_m_probs <- function(
    l_trans_probs,
    v_occupied_state,
    cycle_length,
    v_states_names,
    m_indi_features,
    params_gompertz_sex,
    theta = 1) { # excess-hazard factor: 1 = full NYHA HRs, 0 = general population mortality

  num_i <- nrow(m_indi_features)
  m_probs <- matrix(0,
                    nrow = num_i,
                    ncol = length(v_states_names),
                    dimnames = list(NULL, v_states_names))

  # NYHA transition probabilities (conditional on survival)
  p_1to1 <- l_trans_probs$p_1to1; p_1to2 <- l_trans_probs$p_1to2; p_1to3 <- l_trans_probs$p_1to3; p_1to4 <- l_trans_probs$p_1to4
  p_2to1 <- l_trans_probs$p_2to1; p_2to2 <- l_trans_probs$p_2to2; p_2to3 <- l_trans_probs$p_2to3; p_2to4 <- l_trans_probs$p_2to4
  p_3to1 <- l_trans_probs$p_3to1; p_3to2 <- l_trans_probs$p_3to2; p_3to3 <- l_trans_probs$p_3to3; p_3to4 <- l_trans_probs$p_3to4
  p_4to1 <- l_trans_probs$p_4to1; p_4to2 <- l_trans_probs$p_4to2; p_4to3 <- l_trans_probs$p_4to3; p_4to4 <- l_trans_probs$p_4to4

  # Mortality hazard ratios
  rr_HF  <- l_trans_probs$rr_HF
  rr_2v1 <- l_trans_probs$rr_2v1
  rr_3v1 <- l_trans_probs$rr_3v1
  rr_4v1 <- l_trans_probs$rr_4v1

  # 1. Background cumulative hazard over the cycle, Gompertz h(a) = rate * exp(shape * a):
  #    H(a, a + L) = (rate / shape) * (exp(shape * (a + L)) - exp(shape * a)): The cumulative hazard $H$ is calculated by integrating the hazard function $h(t)$ over that specific time interval.
  v_age <- m_indi_features[, "age"]
  alpha <- params_gompertz_sex$rate
  beta  <- params_gompertz_sex$shape
  if (abs(beta) < 1e-9) {cumulative_hazard_base <- rep(alpha * cycle_length, length(v_age))}
  else {cumulative_hazard_base <- (alpha / beta) * (exp(beta * (v_age + cycle_length)) - exp(beta * v_age))}

  # 2. State-specific hazards with excess-hazard scaling: HR_k = 1 + theta * (rr_HF * HR_k - 1), theta >= 0. theta: excess-hazard theta
  #    Keeps the NYHA ordering and never puts any state below general population mortality.
  v_hr <- 1 + theta * (rr_HF * c(1, rr_2v1, rr_3v1, rr_4v1) - 1)
  cum_h_1 <- cumulative_hazard_base * v_hr[1]
  cum_h_2 <- cumulative_hazard_base * v_hr[2]
  cum_h_3 <- cumulative_hazard_base * v_hr[3]
  cum_h_4 <- cumulative_hazard_base * v_hr[4]

  # Cumulative hazard -> probability of death in the cycle
  p_1toD <- (1 - exp(-cum_h_1)); p_2toD <- (1 - exp(-cum_h_2))
  p_3toD <- (1 - exp(-cum_h_3)); p_4toD <- (1 - exp(-cum_h_4))
  p_1toD[p_1toD > 1] <- 1; p_2toD[p_2toD > 1] <- 1
  p_3toD[p_3toD > 1] <- 1; p_4toD[p_4toD > 1] <- 1

  # 3. Fill each individual's row by current state; NYHA transitions are scaled by (1 - p_death); 'which' extracts vectors of indiviuals in specific states
  idx_nyha1 <- which(v_occupied_state == "NYHA1"); idx_nyha2 <- which(v_occupied_state == "NYHA2")
  idx_nyha3 <- which(v_occupied_state == "NYHA3"); idx_nyha4 <- which(v_occupied_state == "NYHA4")
  idx_death <- which(v_occupied_state == "Death")

  if (length(idx_nyha1) > 0) {
    scaler <- 1 - p_1toD[idx_nyha1] # scaler is a vector, a set of age-specific mortality
    m_probs[idx_nyha1, "NYHA1"] <- p_1to1 * scaler; m_probs[idx_nyha1, "NYHA2"] <- p_1to2 * scaler
    m_probs[idx_nyha1, "NYHA3"] <- p_1to3 * scaler; m_probs[idx_nyha1, "NYHA4"] <- p_1to4 * scaler
    m_probs[idx_nyha1, "Death"] <- p_1toD[idx_nyha1]
  }
  if (length(idx_nyha2) > 0) {
    scaler <- 1 - p_2toD[idx_nyha2]
    m_probs[idx_nyha2, "NYHA1"] <- p_2to1 * scaler; m_probs[idx_nyha2, "NYHA2"] <- p_2to2 * scaler
    m_probs[idx_nyha2, "NYHA3"] <- p_2to3 * scaler; m_probs[idx_nyha2, "NYHA4"] <- p_2to4 * scaler
    m_probs[idx_nyha2, "Death"] <- p_2toD[idx_nyha2]
  }
  if (length(idx_nyha3) > 0) {
    scaler <- 1 - p_3toD[idx_nyha3]
    m_probs[idx_nyha3, "NYHA1"] <- p_3to1 * scaler; m_probs[idx_nyha3, "NYHA2"] <- p_3to2 * scaler
    m_probs[idx_nyha3, "NYHA3"] <- p_3to3 * scaler; m_probs[idx_nyha3, "NYHA4"] <- p_3to4 * scaler
    m_probs[idx_nyha3, "Death"] <- p_3toD[idx_nyha3]
  }
  if (length(idx_nyha4) > 0) {
    scaler <- 1 - p_4toD[idx_nyha4]
    m_probs[idx_nyha4, "NYHA1"] <- p_4to1 * scaler; m_probs[idx_nyha4, "NYHA2"] <- p_4to2 * scaler
    m_probs[idx_nyha4, "NYHA3"] <- p_4to3 * scaler; m_probs[idx_nyha4, "NYHA4"] <- p_4to4 * scaler
    m_probs[idx_nyha4, "Death"] <- p_4toD[idx_nyha4]
  }
  if (length(idx_death) > 0) {
    m_probs[idx_death, "Death"] <- 1
  }

  return(m_probs)
}

# Transition matrix for the cycle: on-treatment probabilities for patients on treatment,
# SoC (placebo) probabilities for patients who have discontinued
update_probsV_disc <- function(
    v_states_names,
    v_occupied_state,
    l_trans_probs_arm,   # probabilities on treatment
    l_trans_probs_soc,   # probabilities after discontinuation (SoC)
    v_discontinued,      # 1 if discontinued, 0 otherwise
    m_indi_features,     # individual ages (for mortality)
    params_gompertz_sex, # sex-specific Gompertz parameters
    theta_arm = 1,       # excess-hazard factor, on treatment
    theta_soc = 1,       # excess-hazard factor, discontinued
    cycle_length = 0.5) {

  m_probs_arm <- get_m_probs(l_trans_probs_arm, v_occupied_state, cycle_length, v_states_names, m_indi_features, params_gompertz_sex, theta = theta_arm)
  m_probs_soc <- get_m_probs(l_trans_probs_soc, v_occupied_state, cycle_length, v_states_names, m_indi_features, params_gompertz_sex, theta = theta_soc)
  m_probs <- m_probs_arm * (1 - v_discontinued) + m_probs_soc * v_discontinued

  if (!all(abs(rowSums(m_probs) - 1) < 1e-8)) {
    warning("Transition probabilities do not sum to 1")
  }

  return(m_probs)
}

# Sample each individual's next health state from their row of transition probabilities
sampleV <- function(m_trans_probs,
  v_states_names,
  v_rand_values = runif(n = nrow(m_trans_probs))) { # uniforms supplied by the caller (common random numbers)

  # Row-wise cumulative probabilities
  m_upper_tri <- upper.tri(x = diag(ncol(m_trans_probs)), diag = TRUE)
  m_cum_probs <- m_trans_probs %*% m_upper_tri
  colnames(m_cum_probs) <- v_states_names

  if (any(m_cum_probs[, ncol(m_cum_probs)] > 1.0000001)) {
    stop("Error in multinomial sampling: probabilities do not sum to 1")
  }

  # The next state is the first state whose cumulative probability exceeds the uniform draw
  m_rand_values <- matrix(
    data  = rep(x = v_rand_values, each = length(v_states_names)),
    nrow  = nrow(m_trans_probs),
    ncol  = length(v_states_names),
    byrow = TRUE
  )
  v_transitions <- rowSums(m_rand_values > m_cum_probs)
  v_health_states <- v_states_names[1 + v_transitions]
  return(v_health_states)
}

# Per-cycle costs by health state: arm costs (state cost + drug) on treatment, SoC costs (no drug) once discontinued
calc_costsV_disc <- function(
    v_occupied_state,
    v_states_costs_arm, # state costs including drug
    v_states_costs_soc, # state costs excluding drug
    v_discontinued) {

  calc_costs_internal <- function(states, base_costs) {
    v_c <- rep(NA, length(states))
    v_c[states == "NYHA1"] <- base_costs["NYHA1"]
    v_c[states == "NYHA2"] <- base_costs["NYHA2"]
    v_c[states == "NYHA3"] <- base_costs["NYHA3"]
    v_c[states == "NYHA4"] <- base_costs["NYHA4"]
    v_c[states == "Death"] <- base_costs["Death"]
    return(v_c)
  }

  v_costs_arm <- calc_costs_internal(v_occupied_state, v_states_costs_arm)
  v_costs_soc <- calc_costs_internal(v_occupied_state, v_states_costs_soc)
  v_final_costs <- v_costs_arm * (1 - v_discontinued) + v_costs_soc * v_discontinued #binary mathematical switch

  return(v_final_costs)
}

# Per-cycle QALYs: state utility x cycle length
calc_effsV <- function(
  v_occupied_state,
  v_states_utilities,
  cycle_length = 1) {

  v_state_utility <- rep(NA, length(v_occupied_state))
  v_state_utility[v_occupied_state == "NYHA1"] <- v_states_utilities["NYHA1"]
  v_state_utility[v_occupied_state == "NYHA2"] <- v_states_utilities["NYHA2"]
  v_state_utility[v_occupied_state == "NYHA3"] <- v_states_utilities["NYHA3"]
  v_state_utility[v_occupied_state == "NYHA4"] <- v_states_utilities["NYHA4"]
  v_state_utility[v_occupied_state == "Death"] <- v_states_utilities["Death"]

  QALYs <- v_state_utility * cycle_length
  return(QALYs)
}

# Cardiovascular hospitalization: at most one Bernoulli event per cycle, with a cost and a one-time disutility
hospV <- function(
    v_occupied_state,
    v_hosp_probs,    # hospitalization probability per cycle, by state
    v_hosp_costs,    # cost per hospitalization, by state
    v_hosp_disutils, # QALY decrement per hospitalization, by state
    v_rand_values = runif(length(v_occupied_state))) { # uniforms supplied by the caller (common random numbers)

  v_probs <- v_hosp_probs[v_occupied_state]
  v_probs[is.na(v_probs)] <- 0
  v_hosp_ind <- as.numeric(v_rand_values < v_probs)

  v_c <- v_hosp_ind * v_hosp_costs[v_occupied_state]
  v_d <- v_hosp_ind * v_hosp_disutils[v_occupied_state]
  v_c[is.na(v_c)] <- 0
  v_d[is.na(v_d)] <- 0

  return(list(hosp_ind = v_hosp_ind, costs = v_c, disutils = v_d))
}

# Discount weights for cycles 0..num_cycles
calc_discount_wts <- function(discount_rate, num_cycles, cycle_length) {
  v_discount_wts <- 1 / (1 + discount_rate) ^ ((0:num_cycles) * cycle_length)
  return(v_discount_wts)
}

#------------------------------------------------------------------------------#

# Run one arm for one sex-specific cohort
run_microSimV_hosp_disc <- function(
    v_starting_states,
    num_i,
    num_cycles,
    m_indi_features,
    v_states_names,
    v_states_costs_arm,  # state costs including drug (on treatment)
    v_states_costs_soc,  # state costs excluding drug (discontinued)
    v_states_utilities,
    l_trans_probs_arm,   # time-varying on-treatment transition probabilities (list by cycle)
    l_trans_probs_soc,   # time-varying SoC transition probabilities (list by cycle)
    p_disc,              # discontinuation probability per cycle
    arm_name,            # "stabilizer" or "silencer" (the silencer keeps a residual effect after discontinuation)
    params_gompertz_sex, # sex-specific Gompertz parameters for background mortality
    theta_arm = 1,       # excess-hazard factor, on treatment
    theta_soc = 1,       # excess-hazard factor, discontinued
    v_hosp_probs,
    v_hosp_costs,
    v_hosp_disutils,
    discount_rate_costs,
    discount_rate_QALYs,
    cycle_length = 0.5,
    starting_seed = 1,
    waning_years = 1) {  # silencer residual effect after discontinuation (years)

  # Matrices of states, costs, QALYs, life-years, and on-treatment status (individuals x cycles)
  m_States <- m_Costs <- m_Effs <- m_LYs <- m_OnTx <- matrix(
    nrow = num_i,
    ncol = num_cycles + 1,
    dimnames = list(paste("ind",   1:num_i,    sep = "_"),
                    paste("cycle", 0:num_cycles, sep = "_"))
  )
  # [DEMO 1] Empty matrices (all NA). Try: dim(m_States); m_States[1:5, 1:6]   (c = continue, Q = quit)
  browser(expr = demo_browser)

  # Discontinuation status (0 = on treatment, 1 = discontinued); patients starting in NYHA4 are untreated
  v_discontinued <- as.numeric(v_starting_states == "NYHA4")
  v_cycles_since_disc <- rep(0, num_i)

  set.seed(starting_seed)
  m_States[, 1] <- v_starting_states

  # Cycle 0 payoffs
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
  m_OnTx[, 1] <- (m_States[, 1] != "Death") * (1 - v_discontinued) # who incurs drug costs
  # [DEMO 2] Cycle 0 filled. 
  browser(expr = demo_browser)

  for (t in 1:num_cycles) {

    # Time-varying transition probabilities for cycles 1-5 (to 30 months); the cycle-5 set is used afterwards
    current_cycle_index <- min(t, length(l_trans_probs_arm))
    l_trans_probs_arm_t <- l_trans_probs_arm[[current_cycle_index]]
    l_trans_probs_soc_t <- l_trans_probs_soc[[current_cycle_index]]

    # Common random numbers: one uniform per individual and process every cycle, whether or not the individual
    # is alive or eligible, so both arms use the same random numbers and differ only through their inputs
    u_disc  <- runif(num_i)
    u_trans <- runif(num_i)
    u_hosp  <- runif(num_i)

    # 1. Discontinuation (alive patients on treatment), at the start of the cycle
    eligible_to_disc <- (v_discontinued == 0) & (m_States[, t] != "Death") 
    v_discontinued[eligible_to_disc & (u_disc < p_disc)] <- 1
    v_cycles_since_disc[v_discontinued == 1] <- v_cycles_since_disc[v_discontinued == 1] + 1

    # Residual effect: discontinued silencer patients keep on-treatment transitions for waning_years
    v_effective_disc <- v_discontinued #decouple the economic discontinuation from the clinical discontinuation 
    if (arm_name == "silencer") {
      waning_period_cycles <- waning_years / cycle_length
      indices_in_waning <- which(v_discontinued == 1 & v_cycles_since_disc <= waning_period_cycles)
      if (length(indices_in_waning) > 0) {
        v_effective_disc[indices_in_waning] <- 0
      }
    }

    # 2. Transition probabilities
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
      cycle_length        = cycle_length
    )

    # 3. Next state
    m_States[, t + 1] <- sampleV(
      m_trans_probs  = m_trans_probs,
      v_states_names = v_states_names,
      v_rand_values  = u_trans
    )

    # 4. Hospitalization in the new state
    l_hosp_results <- hospV(
      v_occupied_state = m_States[, t + 1],
      v_hosp_probs     = v_hosp_probs,
      v_hosp_costs     = v_hosp_costs,
      v_hosp_disutils  = v_hosp_disutils,
      v_rand_values    = u_hosp
    )

    # 5. Entering NYHA4 ends treatment (no drug cost from this cycle; SoC transitions from the next cycle)
    v_discontinued[m_States[, t + 1] == "NYHA4"] <- 1

    # 6. Age by one cycle
    m_indi_features[, "age"] <- m_indi_features[, "age"] + cycle_length

    # 7. Payoffs for cycle t + 1 (drug costs stop immediately at discontinuation)
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
    m_LYs[, t + 1]  <- (m_States[, t + 1] != "Death") * cycle_length #logical coercion. TRUE = 1, FALSE = 0 when multiplied by cycle_length
    m_OnTx[, t + 1] <- (m_States[, t + 1] != "Death") * (1 - v_discontinued)
    # [DEMO 3] Cycle t filled (pauses for t = 1-3). t; round(m_trans_probs[1:5, ], 3); m_States[1:5, 1:(t + 2)]
    #          m_Costs[1:5, 1:(t + 2)]; m_Effs[1:5, 1:(t + 2)]; l_hosp_results$hosp_ind[1:5]; v_discontinued[1:5]
    browser(expr = demo_browser && t <= 3)
  }
  # [DEMO 4] All cycles filled. m_States["ind_1", 1:12]; table(m_States[, num_cycles + 1]); rowSums(m_LYs)[1:5]
  browser(expr = demo_browser)

  # Totals with life-table (half-cycle) correction: average of consecutive cycles, (M[, -last] + M[, -1]) / 2
  v_c_dsc_wts <- calc_discount_wts(discount_rate_costs, num_cycles, cycle_length)
  v_e_dsc_wts <- calc_discount_wts(discount_rate_QALYs, num_cycles, cycle_length)
  lt_sum <- function(m) rowSums((m[, -ncol(m)] + m[, -1]) / 2) # (1 - 80) + (2 - 81 )

  # Undiscounted
  v_total_costs <- lt_sum(m_Costs)
  v_total_qalys <- lt_sum(m_Effs)
  v_total_lys   <- lt_sum(m_LYs)

  # Discounted (weights applied to each cycle before the correction)
  v_total_Dcosts <- lt_sum(t(t(m_Costs) * v_c_dsc_wts)) # transpose, then transpose again
  v_total_Dqalys <- lt_sum(t(t(m_Effs)  * v_e_dsc_wts))
  v_total_Dlys   <- lt_sum(t(t(m_LYs)   * v_e_dsc_wts))
  # Discounted on-treatment cycles: discounted drug cost = per-cycle drug cost x mean_DOnTx
  mean_DOnTx <- mean(lt_sum(t(t(m_OnTx) * v_c_dsc_wts))) # we don't need individual OnTx

  results <- list(
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
    mean_DOnTx     = mean_DOnTx
  )

  return(results)
}

#------------------------------------------------------------------------------#

### Model inputs ###

## General
num_i               <- 1e3   # simulated individuals per sex-specific run (demo; base case uses 1e5)
demo_browser        <- FALSE # demo: TRUE pauses inside run_microSimV_hosp_disc() (set only around one base-case run)
num_cycles          <- 80    # lifetime horizon: 40 years of 6-month cycles
cycle_length        <- 0.5   # years
seed                <- 1234
discount_rate_costs <- 0.03
discount_rate_QALYs <- 0.03

## Population: separate male and female cohorts (sex-specific background mortality), weighted by prop_male
mean_age    <- 78.0
sd_age      <- 5
prop_male   <- 0.875
prop_female <- 1 - prop_male
set.seed(seed)
m_indi_features_male   <- cbind("age" = rnorm(n = num_i, mean = mean_age, sd = sd_age))
m_indi_features_female <- cbind("age" = rnorm(n = num_i, mean = mean_age, sd = sd_age))

## Health states and starting NYHA distribution
v_states_names <- c("NYHA1", "NYHA2", "NYHA3", "NYHA4", "Death")
p_starting <- c(0.134, 0.519, 0.337, 0.01, 0)
v_starting_states <- sample(x = v_states_names, size = num_i, replace = TRUE, prob = p_starting)

## Background mortality: Gompertz fitted by maximum likelihood to the 2022 US Actuarial Life Table,
## ages >= 60 (TRACE_mortality_fitting.R, section 5)
params_gompertz <- list(
  male = list(
    shape = 0.09656038924,
    rate  = 2.842950706e-05
  ),
  female = list(
    shape = 0.1072543865,
    rate  = 8.615557054e-06
  )
)

## Mortality hazard ratios by NYHA class (vs NYHA I)
rr_HF   <- 1    # NYHA I vs general population (NYHA I = general population mortality)
rr_2v1  <- 1.78
rr_3v1  <- 3.51
rr_4v1  <- 5.74

## Time-varying NYHA transition probabilities, one list per 6-month cycle (cycles 1-5);
## the cycle-5 (30-month) set is used for all later cycles. Source: ATTR-ACT (via ICER 2024).

# On-treatment: tafamidis arm of ATTR-ACT
l_trans_probs_st_tv <- list(
  # Cycle 1 (6 months)
  list("p_1to1" = 0.565, "p_1to2" = 0.392, "p_1to3" = 0.043, "p_1to4" = 0,
       "p_2to1" = 0.072, "p_2to2" = 0.751, "p_2to3" = 0.17,  "p_2to4" = 0.007,
       "p_3to1" = 0,     "p_3to2" = 0.29,  "p_3to3" = 0.678, "p_3to4" = 0.032,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 2 (12 months)
  list("p_1to1" = 0.522, "p_1to2" = 0.478, "p_1to3" = 0,     "p_1to4" = 0,
       "p_2to1" = 0.069, "p_2to2" = 0.758, "p_2to3" = 0.166, "p_2to4" = 0.007,
       "p_3to1" = 0.019, "p_3to2" = 0.396, "p_3to3" = 0.566, "p_3to4" = 0.019,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 3 (18 months)
  list("p_1to1" = 0.381, "p_1to2" = 0.476, "p_1to3" = 0.143, "p_1to4" = 0,
       "p_2to1" = 0.096, "p_2to2" = 0.697, "p_2to3" = 0.207, "p_2to4" = 0,
       "p_3to1" = 0.023, "p_3to2" = 0.273, "p_3to3" = 0.681, "p_3to4" = 0.023,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 4 (24 months)
  list("p_1to1" = 0.5,   "p_1to2" = 0.3,   "p_1to3" = 0.15,  "p_1to4" = 0.05,
       "p_2to1" = 0.103, "p_2to2" = 0.675, "p_2to3" = 0.214, "p_2to4" = 0.008,
       "p_3to1" = 0,     "p_3to2" = 0.27,  "p_3to3" = 0.622, "p_3to4" = 0.108,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 5 (30 months), used for all later cycles
  list("p_1to1" = 0.368, "p_1to2" = 0.368, "p_1to3" = 0.211, "p_1to4" = 0.053,
       "p_2to1" = 0.115, "p_2to2" = 0.598, "p_2to3" = 0.287, "p_2to4" = 0,
       "p_3to1" = 0,     "p_3to2" = 0.3,   "p_3to3" = 0.633, "p_3to4" = 0.067,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1)
)

# Silencer: same transitions as stabilizers (equal-efficacy assumption; no head-to-head trial,
# indirect comparisons show no significant difference; as in NICE TA1115)
l_trans_probs_si_tv <- l_trans_probs_st_tv

# After discontinuation (SoC): placebo arm of ATTR-ACT
l_trans_probs_soc_tv <- list(
  # Cycle 1 (6 months)
  list("p_1to1" = 0.538, "p_1to2" = 0.462, "p_1to3" = 0,     "p_1to4" = 0,
       "p_2to1" = 0.062, "p_2to2" = 0.763, "p_2to3" = 0.175, "p_2to4" = 0,
       "p_3to1" = 0.039, "p_3to2" = 0.216, "p_3to3" = 0.706, "p_3to4" = 0.039,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 2 (12 months)
  list("p_1to1" = 0.273, "p_1to2" = 0.545, "p_1to3" = 0.182, "p_1to4" = 0,
       "p_2to1" = 0.07,  "p_2to2" = 0.651, "p_2to3" = 0.267, "p_2to4" = 0.012,
       "p_3to1" = 0,     "p_3to2" = 0.239, "p_3to3" = 0.696, "p_3to4" = 0.065,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 3 (18 months)
  list("p_1to1" = 0.222, "p_1to2" = 0.556, "p_1to3" = 0.111, "p_1to4" = 0.111,
       "p_2to1" = 0.039, "p_2to2" = 0.649, "p_2to3" = 0.299, "p_2to4" = 0.013,
       "p_3to1" = 0.054, "p_3to2" = 0.243, "p_3to3" = 0.676, "p_3to4" = 0.027,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 4 (24 months)
  list("p_1to1" = 0.125, "p_1to2" = 0.75,  "p_1to3" = 0.125, "p_1to4" = 0,
       "p_2to1" = 0.016, "p_2to2" = 0.5,   "p_2to3" = 0.452, "p_2to4" = 0.032,
       "p_3to1" = 0,     "p_3to2" = 0.28,  "p_3to3" = 0.6,   "p_3to4" = 0.12,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 5 (30 months), used for all later cycles
  list("p_1to1" = 0.125, "p_1to2" = 0.75,  "p_1to3" = 0.125, "p_1to4" = 0,
       "p_2to1" = 0.016, "p_2to2" = 0.5,   "p_2to3" = 0.452, "p_2to4" = 0.032,
       "p_3to1" = 0,     "p_3to2" = 0.28,  "p_3to3" = 0.6,   "p_3to4" = 0.12,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1)
)

## Discontinuation: AE-related discontinuation in the pivotal trials (same definition; both enrolled 2019+),
## converted to a 6-month probability assuming a constant hazard: p_cycle = 1 - (1 - P_T)^(1 / cycles in T)
p_disc_2.5y_st <- 0.093 # ATTRibute-CM (acoramidis): 39/421 over 30 months; applied to all stabilizers
p_disc_3y_si   <- 0.03  # HELIOS-B (vutrisiran): 10/326 over up to 36 months

p_disc_st <- 1 - (1 - p_disc_2.5y_st)^(1/5) # 1.93% per cycle
p_disc_si <- 1 - (1 - p_disc_3y_si)^(1/6)   # 0.51% per cycle

## Drug costs (annual, 2026 USD)
c_tafamidis_WAC  <- 286653.60 # WAC (Whole Acquisition Cost), Micromedex RED BOOK
c_acoramidis_WAC <- 257976.79 # WAC, Micromedex RED BOOK
c_vutrisiran_ASP <- 477404.00 # ASP (Average Sales Price), CMS ASP Pricing File
c_discount_st  <- 0.275 # WAC-to-net discount for stabilizers (assumption, ICER 2024; OWSA 20-40%)
c_discount_si  <- 0     # ASP already reflects net price
c_st_annual    <- mean(c(c_tafamidis_WAC, c_acoramidis_WAC)) # stabilizer = mean of tafamidis and acoramidis WAC
c_si_annual    <- c_vutrisiran_ASP
c_st           <- c_st_annual * cycle_length # per cycle
c_si           <- c_si_annual * cycle_length # per cycle

## Health-state costs per cycle (ICER 2024; Wang et al. 2023), Q1 2024 USD inflated to Q2 2026 USD
## with the BEA PCE health care services price index
infl_factor    <- 124.229 / 116.774
c_NYHA1_base   <- 2911  * infl_factor
c_NYHA2_base   <- 4129  * infl_factor
c_NYHA3_base   <- 6194  * infl_factor
c_NYHA4_base   <- 10208 * infl_factor
c_D            <- 0

# State costs on treatment = state cost + net drug cost
c_NYHA1_st     <- c_NYHA1_base + (c_st * (1 - c_discount_st))
c_NYHA2_st     <- c_NYHA2_base + (c_st * (1 - c_discount_st))
c_NYHA3_st     <- c_NYHA3_base + (c_st * (1 - c_discount_st))
c_NYHA4_st     <- c_NYHA4_base + (c_st * (1 - c_discount_st))

c_NYHA1_si     <- c_NYHA1_base + (c_si * (1 - c_discount_si))
c_NYHA2_si     <- c_NYHA2_base + (c_si * (1 - c_discount_si))
c_NYHA3_si     <- c_NYHA3_base + (c_si * (1 - c_discount_si))
c_NYHA4_si     <- c_NYHA4_base + (c_si * (1 - c_discount_si))

v_states_costs_st  <- c("NYHA1" = c_NYHA1_st,   "NYHA2" = c_NYHA2_st,   "NYHA3" = c_NYHA3_st,   "NYHA4" = c_NYHA4_st,   "Death" = c_D)
v_states_costs_si  <- c("NYHA1" = c_NYHA1_si,   "NYHA2" = c_NYHA2_si,   "NYHA3" = c_NYHA3_si,   "NYHA4" = c_NYHA4_si,   "Death" = c_D)
v_states_costs_soc <- c("NYHA1" = c_NYHA1_base, "NYHA2" = c_NYHA2_base, "NYHA3" = c_NYHA3_base, "NYHA4" = c_NYHA4_base, "Death" = c_D)

## Utilities (ICER 2024; ATTR-ACT EQ-5D, adjusted to US norms)
u_NYHA1   <- 0.82
u_NYHA2   <- 0.729
u_NYHA3   <- 0.633
u_NYHA4   <- 0.333
u_Death   <- 0
v_states_utilities <- c("NYHA1" = u_NYHA1, "NYHA2" = u_NYHA2, "NYHA3" = u_NYHA3, "NYHA4" = u_NYHA4, "Death" = u_Death)

## Cardiovascular hospitalization (ICER 2024): probability per cycle by NYHA class (ATTR-ACT pooled across arms;
## HAS 2020), same in both arms, so treatment affects hospitalization only through NYHA class.
## Cost per admission inflated as above; disutility per admission (Griffiths et al. 2017).
v_hosp_probs_st <- c("NYHA1" = 0.1683, "NYHA2" = 0.3107, "NYHA3" = 0.698, "NYHA4" = 0.8627, "Death" = 0)
v_hosp_probs_si <- v_hosp_probs_st
v_hosp_costs    <- c("NYHA1" = 30584.15689, "NYHA2" = 17400.10754, "NYHA3" = 17694.56047, "NYHA4" = 21041.54226, "Death" = 0) * infl_factor
v_hosp_disutils <- c("NYHA1" = 0.023, "NYHA2" = 0.01, "NYHA3" = 0.027, "NYHA4" = 0.07, "Death" = 0)

### Mortality calibration (excess-hazard factor theta) ###
# Target: THAOS tafamidis-treated registry cohort (Garcia-Pavia et al., J Card Fail 2025;31(3):525-533):
# n = 587, median age 77.7, 91.8% male, 91.8% wild-type; survival 84.4% at 30 months, 76.8% at 42 months.
# Relative survival transport: registry survival / registry expected (general population) survival = relative survival; 
# target = model-population expected survival x relative survival. theta is calibrated on the
# stabilizer arm (registry patients were tafamidis-treated) and applied to both arms.

cat("\n--- Calibrating theta to 42-month survival ---\n")

# Sex-weighted model survival at calib_cycles minus the target (root found with uniroot)
get_survival_diff <- function(
    theta_val,
    target_survival,
    calib_cycles,
    arm_params,    # arm-specific inputs
    common_params  # inputs shared by both arms
) {

  run_sex <- function(m_features, gompertz_sex) {
    sim_results <- run_microSimV_hosp_disc(
      v_starting_states   = common_params$v_starting_states,
      num_i               = common_params$num_i,
      num_cycles          = calib_cycles,
      m_indi_features     = m_features,
      v_states_names      = common_params$v_states_names,
      v_states_costs_arm  = arm_params$v_states_costs_arm,
      v_states_costs_soc  = common_params$v_states_costs_soc,
      v_states_utilities  = common_params$v_states_utilities,
      l_trans_probs_arm   = arm_params$l_trans_probs_arm,
      l_trans_probs_soc   = common_params$l_trans_probs_soc,
      p_disc              = arm_params$p_disc,
      arm_name            = arm_params$arm_name,
      params_gompertz_sex = gompertz_sex,
      theta_arm           = theta_val,
      theta_soc           = theta_val,
      v_hosp_probs        = arm_params$v_hosp_probs,
      v_hosp_costs        = common_params$v_hosp_costs,
      v_hosp_disutils     = common_params$v_hosp_disutils,
      discount_rate_costs = common_params$discount_rate_costs,
      discount_rate_QALYs = common_params$discount_rate_QALYs,
      cycle_length        = common_params$cycle_length,
      starting_seed       = common_params$seed
    )
    sum(sim_results$m_States[, calib_cycles + 1] != "Death") / common_params$num_i
  }

  surv_male   <- run_sex(common_params$m_indi_features_male,   common_params$params_gompertz$male)
  surv_female <- run_sex(common_params$m_indi_features_female, common_params$params_gompertz$female)
  survival_at_target <- surv_male * common_params$prop_male + surv_female * (1 - common_params$prop_male)

  return(survival_at_target - target_survival)
}

common_params <- list(
  v_starting_states      = v_starting_states,
  num_i                  = num_i,
  m_indi_features_male   = m_indi_features_male,
  m_indi_features_female = m_indi_features_female,
  prop_male              = prop_male,
  v_states_names         = v_states_names,
  v_states_utilities     = v_states_utilities,
  params_gompertz        = params_gompertz,
  l_trans_probs_soc      = l_trans_probs_soc_tv,
  v_states_costs_soc     = v_states_costs_soc,
  v_hosp_costs           = v_hosp_costs,
  v_hosp_disutils        = v_hosp_disutils,
  discount_rate_costs    = discount_rate_costs,
  discount_rate_QALYs    = discount_rate_QALYs,
  cycle_length           = cycle_length,
  seed                   = seed
)

st_params_calib <- list(
  v_states_costs_arm  = v_states_costs_st,
  l_trans_probs_arm   = l_trans_probs_st_tv,
  p_disc              = p_disc_st,
  v_hosp_probs        = v_hosp_probs_st,
  arm_name            = "stabilizer"
)

calib_months       <- 42
calib_cycles       <- calib_months / 12 / cycle_length # 7 cycles
ref_surv_obs       <- 0.768  # THAOS 42-month survival, tafamidis-treated
ref_mean_age       <- 77.7   # THAOS median age
ref_sd_age         <- sd_age # SD not reported; same spread as the model population
ref_prop_male      <- 0.918

# Expected (general population) survival over t years from Gompertz, averaged over individuals
calc_expected_surv <- function(v_age, gompertz_sex, t) {
  mean(exp(-(gompertz_sex$rate / gompertz_sex$shape) *
             (exp(gompertz_sex$shape * (v_age + t)) - exp(gompertz_sex$shape * v_age))))
}

set.seed(seed)
v_age_ref_male   <- rnorm(n = num_i, mean = ref_mean_age, sd = ref_sd_age)
v_age_ref_female <- rnorm(n = num_i, mean = ref_mean_age, sd = ref_sd_age)
exp_surv_ref <- ref_prop_male * calc_expected_surv(v_age_ref_male, params_gompertz$male, calib_months / 12) +
  (1 - ref_prop_male) * calc_expected_surv(v_age_ref_female, params_gompertz$female, calib_months / 12)
exp_surv_model <- prop_male * calc_expected_surv(m_indi_features_male[, "age"], params_gompertz$male, calib_months / 12) +
  (1 - prop_male) * calc_expected_surv(m_indi_features_female[, "age"], params_gompertz$female, calib_months / 12)

rel_surv_ref      <- ref_surv_obs / exp_surv_ref # relative survival: THAOS survival / expected survival with THAOS cohort
target_surv_calib <- exp_surv_model * rel_surv_ref # calibration target (age, sex adjusted): expected survival with our model cohort * relative survival

cat("Registry expected (general population) survival at", calib_months, "months:", round(exp_surv_ref, 4), "\n")
cat("Registry relative survival:", round(rel_surv_ref, 4), "\n")
cat("Model population expected (general population) survival:", round(exp_surv_model, 4), "\n")
cat("Target", calib_months, "-month survival (stabilizer arm, sex-weighted):", target_surv_calib, "\n")
# theta = 0 gives general population survival, so the root lies in [0, 3] whenever relative survival <= 1
# to find the exact theta, make the different between the simulated survival and the terget 0.
theta_calibration_result <- uniroot(get_survival_diff, interval = c(0, 3), target_survival = target_surv_calib, calib_cycles = calib_cycles, arm_params = st_params_calib, common_params = common_params, tol = 1e-4)
calibrated_theta <- theta_calibration_result$root
cat("-> Calibrated theta (applied to both arms):", calibrated_theta, "\n")
cat("-> Calibrated HRs vs general population (NYHA1-4):", round(1 + calibrated_theta * (rr_HF * c(1, rr_2v1, rr_3v1, rr_4v1) - 1), 3), "\n")

### Base case: run both arms for each sex (same seed in every run: common random numbers) ###

run_arm_base <- function(m_features, sex, arm) {
  is_si <- arm == "silencer"
  run_microSimV_hosp_disc(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_features,
    v_states_names      = v_states_names,
    v_states_costs_arm  = if (is_si) v_states_costs_si else v_states_costs_st,
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
    starting_seed       = seed
  )
}

cat("\n--- Running base case ---\n")
demo_browser  <- TRUE  # demo: pause only in this run (calibration runs above do not stop)
res_st_male   <- run_arm_base(m_indi_features_male,   "male",   "stabilizer") # Stabilizer arm, males
demo_browser  <- FALSE
res_si_male   <- run_arm_base(m_indi_features_male,   "male",   "silencer")   # Silencer arm, males
res_st_female <- run_arm_base(m_indi_features_female, "female", "stabilizer") # Stabilizer arm, females
res_si_female <- run_arm_base(m_indi_features_female, "female", "silencer")   # Silencer arm, females

### Results ###

report_results <- function(label, dc_st, dc_si, dq_st, dq_si, ly_st, ly_si, start_age) {
  cat("\n---", label, "(Silencer vs Stabilizer) ---\n")
  cat("Mean Discounted Costs (Stabilizer):", dc_st, "\n")
  cat("Mean Discounted Costs (Silencer):",   dc_si, "\n")
  cat("Mean Discounted QALYs (Stabilizer):", dq_st, "\n")
  cat("Mean Discounted QALYs (Silencer):",   dq_si, "\n")
  cat("Mean Life Years (Stabilizer):", ly_st, "\n")
  cat("Mean Life Years (Silencer):",   ly_si, "\n")
  cat("Life Expectancy (Age at Death) (Stabilizer):", start_age + ly_st, "\n")
  cat("Life Expectancy (Age at Death) (Silencer):",   start_age + ly_si, "\n")
  icer <- (dc_si - dc_st) / (dq_si - dq_st)
  cat("ICER (Silencer vs Stabilizer):", icer, "\n")
  icer
}

ICER_male <- report_results("MALE COHORT", res_st_male$mean_Dcosts, res_si_male$mean_Dcosts,
                            res_st_male$mean_Dqalys, res_si_male$mean_Dqalys,
                            res_st_male$mean_lys, res_si_male$mean_lys, mean(m_indi_features_male[, "age"]))
ICER_female <- report_results("FEMALE COHORT", res_st_female$mean_Dcosts, res_si_female$mean_Dcosts,
                              res_st_female$mean_Dqalys, res_si_female$mean_Dqalys,
                              res_st_female$mean_lys, res_si_female$mean_lys, mean(m_indi_features_female[, "age"]))

# Overall population: sex-weighted
mean_Dcosts_st_overall <- (res_st_male$mean_Dcosts * prop_male) + (res_st_female$mean_Dcosts * prop_female)
mean_Dcosts_si_overall <- (res_si_male$mean_Dcosts * prop_male) + (res_si_female$mean_Dcosts * prop_female)
mean_Dqalys_st_overall <- (res_st_male$mean_Dqalys * prop_male) + (res_st_female$mean_Dqalys * prop_female)
mean_Dqalys_si_overall <- (res_si_male$mean_Dqalys * prop_male) + (res_si_female$mean_Dqalys * prop_female)
mean_lys_st_overall    <- (res_st_male$mean_lys * prop_male) + (res_st_female$mean_lys * prop_female)
mean_lys_si_overall    <- (res_si_male$mean_lys * prop_male) + (res_si_female$mean_lys * prop_female)
mean_start_age_overall <- (mean(m_indi_features_male[, "age"]) * prop_male) + (mean(m_indi_features_female[, "age"]) * prop_female)

ICER_overall <- report_results("OVERALL POPULATION", mean_Dcosts_st_overall, mean_Dcosts_si_overall,
                               mean_Dqalys_st_overall, mean_Dqalys_si_overall,
                               mean_lys_st_overall, mean_lys_si_overall, mean_start_age_overall)

### Threshold price analysis ###
# Drug cost is linear in price: Dcosts_si = non-drug costs + per-cycle net drug cost x DOnTx_si.
# Solve for the silencer annual price at which the ICER vs stabilizers equals the WTP threshold
# (stabilizer price fixed).
calc_threshold_price <- function(res_st_m, res_si_m, res_st_f, res_si_f, p_male, v_wtp,
                                 c_si_cycle = c_si, discount = c_discount_si, cyc_len = cycle_length) {
  p_female  <- 1 - p_male
  dc_st     <- res_st_m$mean_Dcosts * p_male + res_st_f$mean_Dcosts * p_female
  dc_si     <- res_si_m$mean_Dcosts * p_male + res_si_f$mean_Dcosts * p_female
  dq_st     <- res_st_m$mean_Dqalys * p_male + res_st_f$mean_Dqalys * p_female
  dq_si     <- res_si_m$mean_Dqalys * p_male + res_si_f$mean_Dqalys * p_female
  ontx_si   <- res_si_m$mean_DOnTx  * p_male + res_si_f$mean_DOnTx  * p_female

  dc_si_nondrug  <- dc_si - c_si_cycle * (1 - discount) * ontx_si # non-drug cost
  net_cycle_cost <- (v_wtp * (dq_si - dq_st) + dc_st - dc_si_nondrug) / ontx_si # icer = ((new silencer price + non-drug cost) - stabilizer cost)/diff in cost = WTP
  data.frame(
    wtp               = v_wtp,
    annual_net_price  = net_cycle_cost / cyc_len,
    annual_list_price = net_cycle_cost / cyc_len / (1 - discount) 
  ) # how much the manufacturer would need to discount their drug to hit a specific value benchmark
}

v_wtp_thresholds <- c(100000, 150000)

cat("\n--- THRESHOLD PRICE ANALYSIS (Base case, Silencer annual price) ---\n")
threshold_base <- calc_threshold_price(res_st_male, res_si_male, res_st_female, res_si_female, prop_male, v_wtp_thresholds)
print(threshold_base)

### Cost breakdown helpers (used for the base case and subgroups) ###
# Discounted costs split into drug costs (net per-cycle drug cost x discounted on-treatment cycles) and
# non-intervention costs: health-state costs from discounted NYHA state occupancy (same discounting and
# half-cycle correction as the model) and hospitalization costs (remainder).
v_nyha  <- v_states_names[1:4]
v_c_dsc <- calc_discount_wts(discount_rate_costs, num_cycles, cycle_length) # pre-calculated vector of discount weights for every cycle
c_drug_st_cycle <- c_st * (1 - c_discount_st)
c_drug_si_cycle <- c_si * (1 - c_discount_si)

# Mean discounted (or undiscounted) cycles in each NYHA state, half-cycle corrected
calc_state_cycles <- function(res, wts = v_c_dsc) {
  sapply(v_nyha, function(k) { # instead of for loop
    occ <- colMeans(res$m_States == k) * wts # population proportions per cycle, then apply discount
    sum((occ[-length(occ)] + occ[-1]) / 2) # apply life-table correction
  })
}

summarize_arm <- function(res, drug_cycle) {
  cyc   <- calc_state_cycles(res)
  drug  <- drug_cycle * res$mean_DOnTx #price * mean on-treatment cycles
  state <- sum(v_states_costs_soc[v_nyha] * cyc) # SOC prices * time spent in each state
  c(cyc, cyc_undisc = sum(calc_state_cycles(res, rep(1, num_cycles + 1))), ontx = res$mean_DOnTx,
    drug = drug, state = state, hosp = res$mean_Dcosts - drug - state, total = res$mean_Dcosts)
} # cyc_undisc: undiscounted life yrs

### Subgroup analysis ###
# Subgroup-specific age, sex mix, and starting NYHA distribution; base-case theta kept fixed.
run_subgroup <- function(label, sg_mean_age, sg_sd_age, sg_prop_male, sg_p_starting) {
  set.seed(seed)
  m_male   <- cbind("age" = rnorm(n = num_i, mean = sg_mean_age, sd = sg_sd_age))
  m_female <- cbind("age" = rnorm(n = num_i, mean = sg_mean_age, sd = sg_sd_age))
  v_start  <- sample(x = v_states_names, size = num_i, replace = TRUE, prob = sg_p_starting)

  run_arm <- function(m_features, sex, arm) {
    is_si <- arm == "silencer"
    run_microSimV_hosp_disc(
      v_starting_states   = v_start,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_features,
      v_states_names      = v_states_names,
      v_states_costs_arm  = if (is_si) v_states_costs_si else v_states_costs_st,
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
      starting_seed       = seed # used a fixed global seed -> random number generator resets to the identical state each time 
                                 #-> synchronizing their stocahstic events and minimizaing model variance
    )
  }

  r_st_m <- run_arm(m_male,   "male",   "stabilizer") #module
  r_si_m <- run_arm(m_male,   "male",   "silencer") #module
  r_st_f <- run_arm(m_female, "female", "stabilizer") #module
  r_si_f <- run_arm(m_female, "female", "silencer") #module

  w <- c(sg_prop_male, 1 - sg_prop_male)
  wavg <- function(x_m, x_f) sum(c(x_m, x_f) * w) # demographically weighted average
  out <- data.frame(
    subgroup    = label,
    Dcosts_st   = wavg(r_st_m$mean_Dcosts, r_st_f$mean_Dcosts),
    Dcosts_si   = wavg(r_si_m$mean_Dcosts, r_si_f$mean_Dcosts),
    Dqalys_st   = wavg(r_st_m$mean_Dqalys, r_st_f$mean_Dqalys),
    Dqalys_si   = wavg(r_si_m$mean_Dqalys, r_si_f$mean_Dqalys),
    LYs_st      = wavg(r_st_m$mean_lys,    r_st_f$mean_lys),
    LYs_si      = wavg(r_si_m$mean_lys,    r_si_f$mean_lys)
  )
  # Life expectancy (mean age at death) = mean starting age + undiscounted life-years
  start_age <- wavg(mean(m_male[, "age"]), mean(m_female[, "age"]))
  out$LE_st <- start_age + out$LYs_st
  out$LE_si <- start_age + out$LYs_si
  out$ICER  <- (out$Dcosts_si - out$Dcosts_st) / (out$Dqalys_si - out$Dqalys_st)

  thr <- calc_threshold_price(r_st_m, r_si_m, r_st_f, r_si_f, sg_prop_male, v_wtp_thresholds)
  out$list_price_100k <- thr$annual_list_price[thr$wtp == 100000]
  out$list_price_150k <- thr$annual_list_price[thr$wtp == 150000]

  # Cost breakdown (sex-weighted with the subgroup's sex mix)
  b_st <- summarize_arm(r_st_m, c_drug_st_cycle) * w[1] + summarize_arm(r_st_f, c_drug_st_cycle) * w[2]
  b_si <- summarize_arm(r_si_m, c_drug_si_cycle) * w[1] + summarize_arm(r_si_f, c_drug_si_cycle) * w[2]
  for (k in c("drug", "state", "hosp")) {
    out[[paste0(k, "_st")]] <- b_st[[k]]
    out[[paste0(k, "_si")]] <- b_si[[k]]
  }
  out
}

cat("\n--- SUBGROUP ANALYSIS (Silencer vs Stabilizer) ---\n")
subgroup_results <- rbind(
  run_subgroup("Wild-type (wtATTR)", 79, 5, 0.91, c(0.126, 0.536, 0.326, 0.012, 0)),
  run_subgroup("Hereditary (vATTR)", 71, 5, 0.71, c(0.174, 0.44, 0.385, 0, 0))
)
print(subgroup_results)

### Summary table (base case + subgroups) ###
base_row <- data.frame(
  subgroup        = "Base case",
  Dcosts_st       = mean_Dcosts_st_overall,
  Dcosts_si       = mean_Dcosts_si_overall,
  Dqalys_st       = mean_Dqalys_st_overall,
  Dqalys_si       = mean_Dqalys_si_overall,
  LYs_st          = mean_lys_st_overall,
  LYs_si          = mean_lys_si_overall,
  LE_st           = mean_start_age_overall + mean_lys_st_overall,
  LE_si           = mean_start_age_overall + mean_lys_si_overall,
  ICER            = ICER_overall,
  list_price_100k = threshold_base$annual_list_price[threshold_base$wtp == 100000],
  list_price_150k = threshold_base$annual_list_price[threshold_base$wtp == 150000]
)
summary_results <- rbind(base_row, subgroup_results[, names(base_row)])
summary_results$inc_Dcosts <- summary_results$Dcosts_si - summary_results$Dcosts_st
summary_results$inc_Dqalys <- summary_results$Dqalys_si - summary_results$Dqalys_st
summary_results <- summary_results[, c("subgroup", "Dcosts_st", "Dcosts_si", "inc_Dcosts", "Dqalys_st", "Dqalys_si", "inc_Dqalys",
                                       "LYs_st", "LYs_si", "LE_st", "LE_si", "ICER", "list_price_100k", "list_price_150k")] # subset & reorder

cat("\n--- SUMMARY RESULTS (Silencer vs Stabilizer) ---\n")
print(summary_results, digits = 6, row.names = FALSE)

### Cost breakdown and life-extension decomposition (base case) ###
dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)
wavg_vec <- function(x_m, x_f) x_m * prop_male + x_f * prop_female
arm_st <- wavg_vec(summarize_arm(res_st_male, c_drug_st_cycle), summarize_arm(res_st_female, c_drug_st_cycle))
arm_si <- wavg_vec(summarize_arm(res_si_male, c_drug_si_cycle), summarize_arm(res_si_female, c_drug_si_cycle))

inc_total <- arm_si["total"] - arm_st["total"]
cost_breakdown <- data.frame(
  component  = c("Drug", "Health-state (supportive care)", "Hospitalization", "Non-intervention subtotal", "Total"),
  stabilizer = c(arm_st["drug"], arm_st["state"], arm_st["hosp"], arm_st["state"] + arm_st["hosp"], arm_st["total"]),
  silencer   = c(arm_si["drug"], arm_si["state"], arm_si["hosp"], arm_si["state"] + arm_si["hosp"], arm_si["total"]),
  row.names = NULL)
cost_breakdown$incremental   <- cost_breakdown$silencer - cost_breakdown$stabilizer
cost_breakdown$pct_of_inc    <- cost_breakdown$incremental / inc_total # where the cost difference is coming from
cost_breakdown$pct_of_arm_st <- cost_breakdown$stabilizer / arm_st["total"]
cost_breakdown$pct_of_arm_si <- cost_breakdown$silencer / arm_si["total"]

# Decomposition of the incremental cost (exact accounting identities on sex-weighted, discounted quantities).
# L = discounted life-cycles, phi = share of life on treatment, pi_k = share of life in NYHA k,
# c_bar_k = non-intervention cost per cycle in NYHA k (state cost + expected hospitalization cost).
# The two survival terms are the costs of the added life (impact of life extension on costs).
L_st <- sum(arm_st[v_nyha]); L_si <- sum(arm_si[v_nyha])
phi_st <- arm_st["ontx"] / L_st; phi_si <- arm_si["ontx"] / L_si
pi_st <- arm_st[v_nyha] / L_st; pi_si <- arm_si[v_nyha] / L_si
c_bar <- v_states_costs_soc[v_nyha] + v_hosp_probs_st[v_nyha] * v_hosp_costs[v_nyha]
ni_expected_st <- sum(c_bar * arm_st[v_nyha]); ni_expected_si <- sum(c_bar * arm_si[v_nyha]) # expected non-intervention cost
ni_residual <- (arm_si["state"] + arm_si["hosp"] - ni_expected_si) - (arm_st["state"] + arm_st["hosp"] - ni_expected_st) # to account for stochastic uncertainty

cost_decomposition <- data.frame( # pharmacy cost = drug price * ontx % * life years 
  component = c("Drug: price difference (at stabilizer time on treatment)",
                "Drug: persistence (larger share of life on treatment)",
                "Drug: survival (treatment during added life)",
                "Non-intervention: NYHA mix (time in better classes)",
                "Non-intervention: survival (care during added life)",
                "Non-intervention: hospitalization residual (realized vs expected)",
                "Total incremental cost"),
  value = c((c_drug_si_cycle - c_drug_st_cycle) * phi_st * L_st, # same L, ontx, but attributed to price
            c_drug_si_cycle * (phi_si - phi_st) * L_st, # same price, L, but attributed to ontx
            c_drug_si_cycle * phi_si * (L_si - L_st), # same price, ontx, but attributed to L - drug penalty 
            sum(c_bar * (pi_si - pi_st)) * L_st, # same L, attributed to different proportion of NYHA (thanks to better health)
            sum(c_bar * pi_si) * (L_si - L_st), # same NYHA proportion, but attributed to L (punishment for extended life) - medical cost penalty
            ni_residual,
            inc_total),
  row.names = NULL)
cost_decomposition$pct_of_inc <- cost_decomposition$value / inc_total
stopifnot(abs(sum(cost_decomposition$value[1:6]) - inc_total) < 1e-6 * abs(inc_total))
life_ext_cost <- sum(cost_decomposition$value[c(3, 5)])

cat("\n--- COST BREAKDOWN (discounted, per patient; base case) ---\n")
print(cost_breakdown, digits = 6, row.names = FALSE)
cat("Discounted years on treatment: stabilizer", round(arm_st["ontx"] * cycle_length, 3), "| silencer", round(arm_si["ontx"] * cycle_length, 3), "\n")
cat("Discounted life-years: stabilizer", round(L_st * cycle_length, 3), "| silencer", round(L_si * cycle_length, 3), "\n")
cat("\n--- INCREMENTAL COST DECOMPOSITION (life extension = survival terms) ---\n")
print(cost_decomposition, digits = 6, row.names = FALSE)
cat("Impact of life extension on costs:", round(life_ext_cost), "(", round(100 * life_ext_cost / inc_total, 2), "% of incremental cost )\n")
ICER_drug_only <- cost_breakdown$incremental[1] / (mean_Dqalys_si_overall - mean_Dqalys_st_overall) # to exclude survival penalty
cat("ICER excluding non-intervention costs (drug costs only):", round(ICER_drug_only), "\n")

arm_detail <- data.frame(quantity = c(paste0("Discounted years in ", v_nyha), "Undiscounted life-years", "Discounted years on treatment", "Discounted life-years"),
                         stabilizer = c(arm_st[v_nyha] * cycle_length, arm_st["cyc_undisc"] * cycle_length, arm_st["ontx"] * cycle_length, L_st * cycle_length),
                         silencer   = c(arm_si[v_nyha] * cycle_length, arm_si["cyc_undisc"] * cycle_length, arm_si["ontx"] * cycle_length, L_si * cycle_length),
                         row.names = NULL)
arm_detail$incremental <- arm_detail$silencer - arm_detail$stabilizer
write.csv(cost_breakdown,     file.path(dir_out, "Cost_breakdown_V.1.5_test.csv"), row.names = FALSE)
write.csv(cost_decomposition, file.path(dir_out, "Cost_decomposition_V.1.5_test.csv"), row.names = FALSE)
write.csv(arm_detail,         file.path(dir_out, "Time_in_state_V.1.5_test.csv"), row.names = FALSE)

### Drug vs non-intervention costs: base case and subgroups ###
make_breakdown <- function(label, drug, state, hosp) { # each argument: c(stabilizer, silencer)
  nonint <- state + hosp
  total  <- drug + nonint
  out <- data.frame(population = label,
                    component  = c("Drug", "Health-state (supportive care)", "Hospitalization", "Non-intervention subtotal", "Total"),
                    stabilizer = c(drug[1], state[1], hosp[1], nonint[1], total[1]),
                    silencer   = c(drug[2], state[2], hosp[2], nonint[2], total[2]))
  out$incremental   <- out$silencer - out$stabilizer
  out$pct_of_arm_st <- out$stabilizer / total[1]
  out$pct_of_arm_si <- out$silencer / total[2]
  out
}
cost_breakdown_pop <- rbind(
  make_breakdown("Base case", c(arm_st[["drug"]], arm_si[["drug"]]),
                 c(arm_st[["state"]], arm_si[["state"]]), c(arm_st[["hosp"]], arm_si[["hosp"]])),
  do.call(rbind, lapply(seq_len(nrow(subgroup_results)), function(i) {
    sg <- subgroup_results[i, ]
    make_breakdown(sg$subgroup, c(sg$drug_st, sg$drug_si), c(sg$state_st, sg$state_si), c(sg$hosp_st, sg$hosp_si))
  }))
)
# Check: totals equal the discounted costs from the model
stopifnot(all(abs(cost_breakdown_pop$stabilizer[cost_breakdown_pop$component == "Total"] - summary_results$Dcosts_st) < 1e-6),
          all(abs(cost_breakdown_pop$silencer[cost_breakdown_pop$component == "Total"]   - summary_results$Dcosts_si) < 1e-6))
cat("\n--- DRUG VS NON-INTERVENTION COSTS (discounted, per patient) ---\n")
print(cost_breakdown_pop, digits = 6, row.names = FALSE)
write.csv(cost_breakdown_pop, file.path(dir_out, "Cost_breakdown_by_population_V.1.5_test.csv"), row.names = FALSE)
