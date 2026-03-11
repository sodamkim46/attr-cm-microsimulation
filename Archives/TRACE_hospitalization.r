#
# Vectorized Hospitalization Function
# This function determines if an individual is hospitalized based on their health state.
# It returns a list containing the binary indicator, costs, and disutilities.
hospV <- function(
    v_occupied_state, 
    v_hosp_probs,      #Named vector of hospitalization probabilities per state
    v_hosp_costs,      #Named vector of costs per hospitalization event per state
    v_hosp_disutils) { #Named vector of utility decrements (disutility) per event
  
  n_ind <- length(v_occupied_state)
  
  # 1. Map probabilities to each individual based on their current state
  # Handles cases where the state might not have a defined prob (e.g., "Death")
  v_probs <- v_hosp_probs[v_occupied_state]
  v_probs[is.na(v_probs)] <- 0
  
  # 2. Stochastic determination using rbinom (Vectorized)
  v_hosp_ind <- rbinom(n = n_ind, size = 1, prob = v_probs)
  
  # 3. Calculate event-driven Costs and Disutilities
  # Only incurred if v_hosp_ind == 1
  v_c <- v_hosp_ind * v_hosp_costs[v_occupied_state]
  v_d <- v_hosp_ind * v_hosp_disutils[v_occupied_state]
  
  # Ensure NAs (from states like Death) are treated as 0
  v_c[is.na(v_c)] <- 0
  v_d[is.na(v_d)] <- 0
  
  return(list(hosp_ind = v_hosp_ind, costs = v_c, disutils = v_d))
}

# (Functions update_probsV_st, update_probsV_si, sampleV, calc_costsV, 
# calc_effsV, calc_discount_wts remain UNCHANGED as requested)

# ------------------------------------------------------------------------------

## Run Microsimulation function (Modified for Hospitalization)
run_microSimV_hosp<- function(
    v_starting_states,
    num_i,
    num_cycles,
    m_indi_features,
    v_states_names,
    v_states_costs,
    v_cost_coeffs,
    v_states_utilities,
    v_util_coeffs,
    v_util_t_decs,
    l_trans_probs,
    # [MODIFIED] Added hospitalization parameters
    v_hosp_probs,     
    v_hosp_costs,     
    v_hosp_disutils,  
    discount_rate_costs,
    discount_rate_QALYs,
    cycle_length = 0.5,
    starting_seed = 1) {
  
  # Create matrices to capture states, costs, and QALYs
  # [MODIFIED] Added m_Hosp to track hospitalization events over time
  m_States <- m_Costs <- m_Effs <- m_Hosp <- matrix(
    nrow = num_i,
    ncol = num_cycles + 1,
    dimnames = list(paste("ind",   1:num_i,    sep ="_"),
                    paste("cycle", 0:num_cycles, sep ="_"))
  )
  
  # Set the seed for reproducibility
  set.seed(starting_seed)
  
  # Initialize parameter tracking time in current state
  v_time_in_state <- rep(1, times = num_i)
  
  # Get the initial health state
  m_States[, 1] <- v_starting_states
  
  # [MODIFIED] Initialize hospitalization record for cycle 0
  m_Hosp[, 1] <- 0 
  
  # Calculate the costs incurred in their starting health state
  m_Costs[, 1]  <- calc_costsV(
    v_occupied_state = m_States[, 1],
    v_states_costs   = v_states_costs,
    m_indi_features  = m_indi_features,
    v_cost_coeffs    = v_cost_coeffs
  )
  
  # Calculate the QALYs accrued in their starting health state
  m_Effs[, 1]   <- calc_effsV(
    v_occupied_state   = m_States[, 1],
    v_states_utilities = v_states_utilities,
    m_indi_features    = m_indi_features,
    v_util_coeffs      = v_util_coeffs,
    v_util_t_decs      = v_util_t_decs,
    v_time_in_state    = v_time_in_state,
    cycle_length       = cycle_length
  )
  
  # Loop through each cycle 't'
  for (t in 1:num_cycles) {
    # Update the transition probabilities at cycle 't'
    m_trans_probs     <- update_probsV_st(
      v_states_names   = v_states_names,
      v_occupied_state = m_States[, t],
      l_trans_probs    = l_trans_probs,
      v_time_in_state  = v_time_in_state
    )
    
    # Sample the health state at 't + 1'
    m_States[, t + 1] <- sampleV(
      m_trans_probs  = m_trans_probs,
      v_states_names = v_states_names
    )
    
    # [MODIFIED] Simulate Hospitalization events
    # Based on the state the individual just transitioned into (t+1)
    l_hosp_results <- hospV(
      v_occupied_state = m_States[, t + 1],
      v_hosp_probs     = v_hosp_probs,
      v_hosp_costs     = v_hosp_costs,
      v_hosp_disutils  = v_hosp_disutils
    )
    m_Hosp[, t + 1] <- l_hosp_results$hosp_ind
    
    # Keep track of time in state at 't + 1'
    stayed                   <- m_States[, t] == m_States[, t + 1] 
    v_time_in_state[stayed]  <- v_time_in_state[stayed] + 1        
    v_time_in_state[!stayed] <- 1                                  
    
    # Keep track of time in the model (aging)
    m_indi_features[, "age"] <- m_indi_features[, "age"] + 1
    
    # Calculate the costs incurred in their 't + 1' health state
    # [MODIFIED] Total Cost = State-based Cost + Hospitalization Event Cost
    m_Costs[, t + 1]  <- calc_costsV(
      v_occupied_state = m_States[, t + 1],
      v_states_costs   = v_states_costs,
      m_indi_features  = m_indi_features,
      v_cost_coeffs    = v_cost_coeffs
    ) + l_hosp_results$costs
    
    # Calculate the QALYs accrued in their 't + 1' health state
    # [MODIFIED] Total QALY = State-based Utility - Hospitalization Disutility
    # Note: This assumes disutility is an 'event-based' lump sum deduction.
    m_Effs[, t + 1]   <- calc_effsV(
      v_occupied_state   = m_States[, t + 1],
      v_states_utilities = v_states_utilities,
      m_indi_features    = m_indi_features,
      v_util_coeffs      = v_util_coeffs,
      v_util_t_decs      = v_util_t_decs,
      v_time_in_state    = v_time_in_state,
      cycle_length       = cycle_length
    ) - l_hosp_results$disutils
    
  }