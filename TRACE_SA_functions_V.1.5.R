# Helper functions for the sensitivity and scenario analyses
# Run TRACE_microsim_V.1.5.R first: these functions reuse its model functions, inputs, and calibrated theta

### Base case inputs ###
# Every SA script starts from this list and changes only the inputs it needs
base_inputs <- list(
  p_disc_3y_si      = p_disc_3y_si,
  p_disc_2.5y_st    = p_disc_2.5y_st,
  waning_years      = 1,
  c_si_annual       = c_si_annual,
  c_st_annual       = c_st_annual,
  c_discount_st     = c_discount_st,
  c_discount_si     = c_discount_si,
  discount_rate     = discount_rate_costs,
  u_NYHA1           = u_NYHA1,
  u_NYHA2           = u_NYHA2,
  u_NYHA3           = u_NYHA3,
  u_NYHA4           = u_NYHA4,
  state_cost_mult   = 1,
  hosp_prob_mult    = 1,
  hosp_rate_mult    = 1, # scales the rate, not the probability: p' = 1 - (1 - p)^mult
  hosp_cost_mult    = 1,
  hosp_disutil_mult = 1,
  theta             = calibrated_theta,
  p_disc_cycle_st   = NULL, # per-cycle (6-month) discontinuation; if set, used instead of the trial probability above
  p_disc_cycle_si   = NULL,
  # PSA only: values by NYHA class (I-IV); if set, used instead of the multipliers above
  c_states          = NULL, # health state cost per cycle (inflated)
  c_hosp            = NULL, # cost per hospitalization (inflated)
  u_hosp            = NULL, # disutility per hospitalization
  # PSA only: multipliers on the transition matrices (ICER 2024 approach)
  trans_mult_tx     = 1,    # on-treatment matrix (same in both arms, equal efficacy)
  trans_mult_soc    = 1     # SoC matrix, used after discontinuation
)

### Scale transition probabilities (ICER 2024 PSA approach) ###
# Moves to another NYHA class are multiplied by mult, and staying in the same class = 1 - sum of the moves
# If the moves add up to more than 1, they are rescaled to sum to 1 (staying = 0)
scale_trans_probs <- function(l_tv, mult) {
  if (mult == 1) return(l_tv)
  lapply(l_tv, function(l) { # one list per cycle (cycles 1-5)
    for (i in 1:4) { # i = current NYHA class
      off <- setdiff(1:4, i) # except my current status
      p_off <- unlist(l[paste0("p_", i, "to", off)]) * mult
      if (sum(p_off) > 1) p_off <- p_off / sum(p_off) # force to be 1
      l[paste0("p_", i, "to", off)] <- as.list(p_off)
      l[[paste0("p_", i, "to", i)]] <- 1 - sum(p_off) # probability of staying
    }
    l
  })
}

### Run one case: both arms x both sexes, then sex-weighted results ###
# inp: list of inputs (start from base_inputs)
# pop: optional smaller population, list(v_start, m_male, m_female); the PSA uses this
run_owsa_case <- function(inp, pop = NULL) {
  if (is.null(pop)) pop <- list(v_start = v_starting_states, m_male = m_indi_features_male, m_female = m_indi_features_female)
  n_pop <- length(pop$v_start)

  # 1. Discontinuation per cycle: from the trial probability, unless a per-cycle value is given
  p_disc_st_case <- if (is.null(inp$p_disc_cycle_st)) 1 - (1 - inp$p_disc_2.5y_st)^(1/5) else inp$p_disc_cycle_st
  p_disc_si_case <- if (is.null(inp$p_disc_cycle_si)) 1 - (1 - inp$p_disc_3y_si)^(1/6) else inp$p_disc_cycle_si

  # 2. Costs: health state costs (SoC), then add the net drug cost for each arm
  c_states_case <- if (is.null(inp$c_states)) c(c_NYHA1_base, c_NYHA2_base, c_NYHA3_base, c_NYHA4_base) * inp$state_cost_mult else inp$c_states
  v_costs_soc_case <- setNames(c(c_states_case, 0), v_states_names)
  drug_st <- inp$c_st_annual * cycle_length * (1 - inp$c_discount_st)
  drug_si <- inp$c_si_annual * cycle_length * (1 - inp$c_discount_si)
  v_costs_st_case <- v_costs_soc_case + c(rep(drug_st, 4), 0)
  v_costs_si_case <- v_costs_soc_case + c(rep(drug_si, 4), 0)

  # 3. Utilities and hospitalization inputs
  v_utils_case <- c("NYHA1" = inp$u_NYHA1, "NYHA2" = inp$u_NYHA2, "NYHA3" = inp$u_NYHA3, "NYHA4" = inp$u_NYHA4, "Death" = 0)
  v_hosp_probs_case    <- pmin((1 - (1 - v_hosp_probs_si)^inp$hosp_rate_mult) * inp$hosp_prob_mult, 1) # same in both arms
  v_hosp_costs_case    <- if (is.null(inp$c_hosp)) v_hosp_costs * inp$hosp_cost_mult else setNames(c(inp$c_hosp, 0), v_states_names)
  v_hosp_disutils_case <- if (is.null(inp$u_hosp)) v_hosp_disutils * inp$hosp_disutil_mult else setNames(c(inp$u_hosp, 0), v_states_names)

  # 4. Transition probabilities
  l_tx_case  <- scale_trans_probs(l_trans_probs_st_tv,  inp$trans_mult_tx) # silencer = stabilizer (equal efficacy)
  l_soc_case <- scale_trans_probs(l_trans_probs_soc_tv, inp$trans_mult_soc)

  # 5. Run the model (same function and seed as the base case)
  run_arm <- function(m_features, sex, arm) {
    is_si <- arm == "silencer"
    run_microSimV_hosp_disc(
      v_starting_states   = pop$v_start,
      num_i               = n_pop,
      num_cycles          = num_cycles,
      m_indi_features     = m_features,
      v_states_names      = v_states_names,
      v_states_costs_arm  = if (is_si) v_costs_si_case else v_costs_st_case,
      v_states_costs_soc  = v_costs_soc_case,
      v_states_utilities  = v_utils_case,
      l_trans_probs_arm   = l_tx_case,
      l_trans_probs_soc   = l_soc_case,
      p_disc              = if (is_si) p_disc_si_case else p_disc_st_case,
      arm_name            = arm,
      params_gompertz_sex = params_gompertz[[sex]],
      theta_arm           = inp$theta,
      theta_soc           = inp$theta,
      v_hosp_probs        = v_hosp_probs_case,
      v_hosp_costs        = v_hosp_costs_case,
      v_hosp_disutils     = v_hosp_disutils_case,
      discount_rate_costs = inp$discount_rate,
      discount_rate_QALYs = inp$discount_rate,
      cycle_length        = cycle_length,
      starting_seed       = seed,
      waning_years        = inp$waning_years
    )
  }

  r_st_m <- run_arm(pop$m_male,   "male",   "stabilizer")
  r_si_m <- run_arm(pop$m_male,   "male",   "silencer")
  r_st_f <- run_arm(pop$m_female, "female", "stabilizer")
  r_si_f <- run_arm(pop$m_female, "female", "silencer")

  # 6. Sex-weighted results
  wavg <- function(x_m, x_f) x_m * prop_male + x_f * prop_female
  Dcosts_st <- wavg(r_st_m$mean_Dcosts, r_st_f$mean_Dcosts); Dcosts_si <- wavg(r_si_m$mean_Dcosts, r_si_f$mean_Dcosts)
  Dqalys_st <- wavg(r_st_m$mean_Dqalys, r_st_f$mean_Dqalys); Dqalys_si <- wavg(r_si_m$mean_Dqalys, r_si_f$mean_Dqalys)
  inc_c <- Dcosts_si - Dcosts_st
  inc_q <- Dqalys_si - Dqalys_st
  # 7. Threshold price: drug cost is linear in price, so solve for the silencer list price where ICER = WTP
  #    (stabilizer price fixed)
  ontx_si <- wavg(r_si_m$mean_DOnTx, r_si_f$mean_DOnTx)
  inc_c_nondrug <- inc_c - drug_si * ontx_si
  thr_list <- function(wtp) (wtp * inc_q - inc_c_nondrug) / ontx_si / cycle_length / (1 - inp$c_discount_si)
  c(Dcosts_st = Dcosts_st, Dcosts_si = Dcosts_si, Dqalys_st = Dqalys_st, Dqalys_si = Dqalys_si,
    inc_Dcosts = inc_c, inc_Dqalys = inc_q, ICER = inc_c / inc_q,
    list_price_100k = thr_list(100000), list_price_150k = thr_list(150000))
}

### Re-calibrate theta for a different 42-month survival target (used in the OWSA) ###
# Same method as the base case calibration
# If the target is at or above general population survival, theta = 0
calibrate_theta <- function(surv_obs) {
  target <- exp_surv_model * (surv_obs / exp_surv_ref)
  f_at_0 <- get_survival_diff(0, target, calib_cycles, st_params_calib, common_params)
  if (f_at_0 <= 0) {
    cat("  [Note] Survival target", surv_obs, "is at/above general population survival; theta floored at 0\n")
    return(0)
  }
  uniroot(get_survival_diff, interval = c(0, 3), target_survival = target, calib_cycles = calib_cycles,
          arm_params = st_params_calib, common_params = common_params, tol = 1e-4)$root
}
