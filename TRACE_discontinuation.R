# Vectorised Healthy-Sick-Dead microsimulation model with Discontinuation Logic

# clear R's session memory (Global Environment)
rm(list = ls())

### Defining model functions ###

# Helper function to calculate transition matrix for a given set of parameters
# This extracts the logic from the original update_probsV so it can be reused for different arms
get_m_probs <- function(l_trans_probs, v_occupied_state, v_time_in_state, cycle_length, v_states_names) {
  with(
    data = l_trans_probs,
    expr = {
      
      # update probabilities of death after first converting them to rates and applying the rate ratio
      p_Dto1  <- 0
      p_Dto2  <- 0
      p_Dto3  <- 0
      p_Dto4  <- 0
      p_DtoD  <- 1
      
      p_1toD  <- rr_HF*p_HtoD
      
      r_1toD  <- -log(1- p_1toD)        #rate of death in NYHA1
      r_2toD  <-  rr_2v1*r_1toD         #rate of death in NYHA2
      r_3toD  <-  rr_3v1*r_1toD         #rate of death in NYHA2
      r_4toD  <-  rr_4v1*r_1toD         #rate of death in NYHA2
      
      p_2toD  <-  1 - exp(-r_2toD*cycle_length)      #probability to die in NYHA2
      
      # calculate p_S1D conditional on current state and duration of being sick
      p_3toD  <-  1 - exp(-r_3toD*cycle_length*(1 + v_time_in_state[v_occupied_state == "NYHA3"] * rp_sv))      #probability to die in NYHA3
      p_4toD  <-  1 - exp(-r_4toD*cycle_length*(1 + v_time_in_state[v_occupied_state == "NYHA4"] * rp_sv))      #probability to die in NYHA4
      
      p_1toD  <- rep(p_1toD, length(which(v_occupied_state == "NYHA1")))
      p_2toD  <- rep(p_2toD, length(which(v_occupied_state == "NYHA2")))
      
      p_Dto1  <- rep(0, length(v_occupied_state[v_occupied_state == "Death"]))
      p_Dto2  <- rep(0, length(v_occupied_state[v_occupied_state == "Death"]))
      p_Dto3  <- rep(0, length(v_occupied_state[v_occupied_state == "Death"]))
      p_Dto4  <- rep(0, length(v_occupied_state[v_occupied_state == "Death"]))
      p_DtoD  <- rep(1, length(v_occupied_state[v_occupied_state == "Death"]))
      
      # Create a state transition probabilities matrix
      m_probs <- matrix(
        nrow = length(v_time_in_state), # a row for each individual
        ncol = length(v_states_names),  # a column for each state
        dimnames = list(
          v_occupied_state,             # name each row based on the occupied state
          v_states_names                # give each column one of the states names
        )
      )
      
      # update m_probs with the appropriate probabilities
      m_probs[v_occupied_state == "NYHA1", ] <- cbind(p_1to1*(1-p_1toD), p_1to2*(1-p_1toD), p_1to3*(1-p_1toD), p_1to4*(1-p_1toD), p_1toD) # transition probabilities when in NYHA1
      m_probs[v_occupied_state == "NYHA2", ] <- cbind(p_2to1*(1-p_2toD), p_2to2*(1-p_2toD), p_2to3*(1-p_2toD), p_2to4*(1-p_2toD), p_2toD) # transition probabilities when in NYHA2
      m_probs[v_occupied_state == "NYHA3", ] <- cbind(p_3to1*(1-p_3toD), p_3to2*(1-p_3toD), p_3to3*(1-p_3toD), p_3to4*(1-p_3toD), p_3toD) # transition probabilities when in NYHA3
      m_probs[v_occupied_state == "NYHA4", ] <- cbind(p_4to1*(1-p_4toD), p_4to2*(1-p_4toD), p_4to3*(1-p_4toD), p_4to4*(1-p_4toD), p_4toD) # transition probabilities when in NYHA4
      m_probs[v_occupied_state == "Death", ] <- cbind(p_Dto1, p_Dto2, p_Dto3, p_Dto4, p_DtoD)                                            # transition probabilities when dead
      
      return(m_probs)
    }
  )
}

## Update Transition Probability function (Modified for Discontinuation)
update_probsV_disc <- function(
    v_states_names,
    v_occupied_state,
    l_trans_probs_arm, # Probabilities for those on treatment
    l_trans_probs_soc, # Probabilities for those discontinued (SoC)
    v_discontinued,    # Binary vector: 1 if discontinued, 0 otherwise
    v_time_in_state,
    cycle_length = 0.5) {
  
  # 1. Calculate probabilities assuming everyone is on treatment
  m_probs_arm <- get_m_probs(l_trans_probs_arm, v_occupied_state, v_time_in_state, cycle_length, v_states_names)
  
  # 2. Calculate probabilities assuming everyone is on SoC
  m_probs_soc <- get_m_probs(l_trans_probs_soc, v_occupied_state, v_time_in_state, cycle_length, v_states_names)
  
  # 3. Combine based on discontinuation status
  # v_discontinued is a vector of length N. We multiply to broadcast across columns.
  m_probs <- m_probs_arm * (1 - v_discontinued) + m_probs_soc * v_discontinued
  
  # Sanity check
  # We check if row sums are close to 1 (using a small tolerance for floating point errors)
  if(!all(abs(rowSums(m_probs) - 1) < 1e-8)) {
    print("Probabilities do not sum to 1")
  }
  
  return(m_probs)
}

## Sample Health States function (Unchanged)
sampleV <- function(m_trans_probs, v_states_names) {
  m_upper_tri <- upper.tri(x = diag(ncol(m_trans_probs)), diag = TRUE)
  m_cum_probs <- m_trans_probs %*% m_upper_tri
  colnames(m_cum_probs) <- v_states_names
  if (any(m_cum_probs[, ncol(m_cum_probs)] > 1.0000001)) { # Slight tolerance
    stop("Error in multinomial sampling: probabilities do not sum to 1")
  }
  v_rand_values <- runif(n = nrow(m_trans_probs))
  m_rand_values <- matrix(data = rep(x = v_rand_values, each = length(v_states_names)), nrow = nrow(m_trans_probs), ncol = length(v_states_names), byrow = TRUE)
  m_transitions <- m_rand_values > m_cum_probs
  v_transitions <- rowSums(m_transitions)
  v_health_states <- v_states_names[1 + v_transitions]
  return(v_health_states)
}

## Calculate Costs function (Modified for Discontinuation)
calc_costsV_disc <- function (
    v_occupied_state,
    v_states_costs_arm, # Costs including drug
    v_states_costs_soc, # Costs excluding drug (SoC)
    m_indi_features,
    v_cost_coeffs,
    v_discontinued) {
  
  # calculate individual-specific costs based on costs regression coefficients
  v_indi_costs <- m_indi_features %*% v_cost_coeffs
  
  # Helper to calculate costs for a specific cost vector
  calc_costs_internal <- function(states, base_costs, indi_costs) {
    v_c <- rep(NA, length(states))
    v_c[states == "NYHA1"] <- base_costs["NYHA1"]
    v_c[states == "NYHA2"] <- base_costs["NYHA2"]
    v_c[states == "NYHA3"] <- base_costs["NYHA3"] + indi_costs[states == "NYHA3"]
    v_c[states == "NYHA4"] <- base_costs["NYHA4"] + indi_costs[states == "NYHA4"]
    v_c[states == "Death"] <- base_costs["Death"]
    return(v_c)
  }
  
  # Calculate costs assuming on treatment
  v_costs_arm <- calc_costs_internal(v_occupied_state, v_states_costs_arm, v_indi_costs)
  
  # Calculate costs assuming SoC (discontinued)
  v_costs_soc <- calc_costs_internal(v_occupied_state, v_states_costs_soc, v_indi_costs)
  
  # Combine
  v_final_costs <- v_costs_arm * (1 - v_discontinued) + v_costs_soc * v_discontinued
  
  return(v_final_costs)
}

## Calculate Health Outcomes function (Unchanged)
calc_effsV <- function (v_occupied_state, v_states_utilities, m_indi_features, v_util_coeffs, v_util_t_decs, v_time_in_state, cycle_length = 1) {
  v_ind_decrement <- (m_indi_features %*% v_util_coeffs)[,1]
  time_decrement <- rep(0, length(v_occupied_state))
  time_decrement[v_occupied_state == "NYHA3"] <- v_util_t_decs["NYHA3 or 4"] * v_time_in_state[v_occupied_state == "NYHA3"]
  time_decrement[v_occupied_state == "NYHA4"] <- v_util_t_decs["NYHA3 or 4"] * v_time_in_state[v_occupied_state == "NYHA4"]
  decrement <- v_ind_decrement + time_decrement
  v_state_utility <- rep(NA, length(v_occupied_state))
  v_state_utility[v_occupied_state == "NYHA1"]  <- v_states_utilities["NYHA1"]
  v_state_utility[v_occupied_state == "NYHA2"]  <- v_states_utilities["NYHA2"]
  v_state_utility[v_occupied_state == "NYHA3"]  <- v_states_utilities["NYHA3"]  + decrement[v_occupied_state == "NYHA3"]
  v_state_utility[v_occupied_state == "NYHA4"]  <- v_states_utilities["NYHA4"]  + decrement[v_occupied_state == "NYHA4"]
  v_state_utility[v_occupied_state == "Death"]  <- v_states_utilities["Death"]
  QALYs <-  v_state_utility * cycle_length
  return(QALYs)
}

## Hospitalization Function
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

## Calculate Discount Weights function (Unchanged)
calc_discount_wts <- function(discount_rate, num_cycles, cycle_length) {
  v_discount_wts <- 1 / (1 + discount_rate) ^ ((0:num_cycles) * cycle_length)
  return(v_discount_wts)
}

#------------------------------------------------------------------------------#

## Run Microsimulation function (Modified for Discontinuation)
run_microSimV_hosp_disc <- function(
    v_starting_states,
    num_i,
    num_cycles,
    m_indi_features,
    v_states_names,
    v_states_costs_arm, # Costs for the treatment arm
    v_states_costs_soc, # Costs for SoC (discontinued)
    v_cost_coeffs,
    v_states_utilities,
    v_util_coeffs,
    v_util_t_decs,
    l_trans_probs_arm,  # Probs for treatment arm
    l_trans_probs_soc,  # Probs for SoC
    p_disc,             # Probability of discontinuation per cycle
    v_hosp_probs,     
    v_hosp_costs,     
    v_hosp_disutils,  
    discount_rate_costs,
    discount_rate_QALYs,
    cycle_length = 0.5,
    starting_seed = 1) {
  
  # Create matrices to capture states, costs, QALYs, and Hospitalizations
  m_States <- m_Costs <- m_Effs <- m_Hosp <- matrix(
    nrow = num_i,
    ncol = num_cycles + 1,
    dimnames = list(paste("ind",   1:num_i,    sep ="_"),
                    paste("cycle", 0:num_cycles, sep ="_"))
  )
  
  # Track discontinuation status (0 = on treatment, 1 = discontinued)
  v_discontinued <- rep(0, num_i)
  
  set.seed(starting_seed)
  v_time_in_state <- rep(1, times = num_i)
  m_States[, 1] <- v_starting_states
  m_Hosp[, 1] <- 0 
  
  # Initial costs (Cycle 0) - assume everyone starts on treatment
  m_Costs[, 1]  <- calc_costsV_disc(
    v_occupied_state   = m_States[, 1],
    v_states_costs_arm = v_states_costs_arm,
    v_states_costs_soc = v_states_costs_soc,
    m_indi_features    = m_indi_features,
    v_cost_coeffs      = v_cost_coeffs,
    v_discontinued     = v_discontinued
  )
  
  m_Effs[, 1]   <- calc_effsV(
    v_occupied_state   = m_States[, 1],
    v_states_utilities = v_states_utilities,
    m_indi_features    = m_indi_features,
    v_util_coeffs      = v_util_coeffs,
    v_util_t_decs      = v_util_t_decs,
    v_time_in_state    = v_time_in_state,
    cycle_length       = cycle_length
  )
  
  for (t in 1:num_cycles) {
    
    # 1. Check for Discontinuation
    # Only those currently on treatment (v_discontinued == 0) and alive can discontinue
    # We assume discontinuation happens at the start of the interval, affecting the transition and costs for this cycle
    eligible_to_disc <- (v_discontinued == 0) & (m_States[, t] != "Death")
    if (any(eligible_to_disc)) {
      n_eligible <- sum(eligible_to_disc)
      # Sample discontinuation events
      new_disc <- rbinom(n = n_eligible, size = 1, prob = p_disc)
      v_discontinued[eligible_to_disc] <- new_disc
    }
    
    # 2. Update transition probabilities
    # Uses v_discontinued to select between arm probs and SoC probs
    m_trans_probs <- update_probsV_disc(
      v_states_names    = v_states_names,
      v_occupied_state  = m_States[, t],
      l_trans_probs_arm = l_trans_probs_arm,
      l_trans_probs_soc = l_trans_probs_soc,
      v_discontinued    = v_discontinued,
      v_time_in_state   = v_time_in_state,
      cycle_length      = cycle_length
    )
    
    # 3. Sample next state
    m_States[, t + 1] <- sampleV(
      m_trans_probs  = m_trans_probs,
      v_states_names = v_states_names
    )
    
    # 4. Simulate Hospitalization
    l_hosp_results <- hospV(
      v_occupied_state = m_States[, t + 1],
      v_hosp_probs     = v_hosp_probs,
      v_hosp_costs     = v_hosp_costs,
      v_hosp_disutils  = v_hosp_disutils
    )
    m_Hosp[, t + 1] <- l_hosp_results$hosp_ind
    
    # 5. Update time in state and age
    stayed                   <- m_States[, t] == m_States[, t + 1] 
    v_time_in_state[stayed]  <- v_time_in_state[stayed] + 1        
    v_time_in_state[!stayed] <- 1                                  
    m_indi_features[, "age"] <- m_indi_features[, "age"] + 1
    
    # 6. Calculate Costs
    # Uses v_discontinued to remove drug costs if applicable
    m_Costs[, t + 1]  <- calc_costsV_disc(
      v_occupied_state   = m_States[, t + 1],
      v_states_costs_arm = v_states_costs_arm,
      v_states_costs_soc = v_states_costs_soc,
      m_indi_features    = m_indi_features,
      v_cost_coeffs      = v_cost_coeffs,
      v_discontinued     = v_discontinued
    ) + l_hosp_results$costs
    
    # 7. Calculate QALYs
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
  
  # Discounting and Summary
  v_c_dsc_wts <- calc_discount_wts(discount_rate_costs, num_cycles, cycle_length)
  v_e_dsc_wts <- calc_discount_wts(discount_rate_QALYs, num_cycles, cycle_length)
  
  v_total_costs <- rowSums(m_Costs)
  v_total_qalys <- rowSums(m_Effs)
  mean_costs    <- mean(v_total_costs)
  mean_qalys    <- mean(v_total_qalys)
  
  v_total_Dcosts <- m_Costs %*% v_c_dsc_wts
  v_total_Dqalys <- m_Effs  %*% v_e_dsc_wts
  mean_Dcosts    <- mean(v_total_Dcosts)
  mean_Dqalys    <- mean(v_total_Dqalys)
  
  results <- list(
    m_States       = m_States,
    m_Costs        = m_Costs,
    m_Effs         = m_Effs,
    v_total_costs  = v_total_costs,
    v_total_qalys  = v_total_qalys,
    v_total_Dcosts = v_total_Dcosts,
    v_total_Dqalys = v_total_Dqalys,
    mean_costs     = mean_costs,
    mean_qalys     = mean_qalys,
    mean_Dcosts    = mean_Dcosts,
    mean_Dqalys    = mean_Dqalys
  )
  
  return(results)
}

#------------------------------------------------------------------------------#

### Defining model parameters ###

# Define model inputs
## General parameters
num_i               <- 1e5               # number of simulated individuals
num_cycles          <- 50                # time horizon
cycle_length        <- 0.5               # length of cycle (in years)
seed                <- 1234              
wtp                 <- 150000            
discount_rate_costs <- 0.03              
discount_rate_QALYs <- 0.03              

## Population characteristics
mean_age            <- 77                
sd_age              <- 5                 
set.seed(seed)                           
m_indi_features     <- cbind("age" = rnorm(n = num_i, mean = mean_age, sd = sd_age))

## Health states
v_states_names <- c("NYHA1","NYHA2", "NYHA3", "NYHA4", "Death")
p_starting <- c(0.108, 0.72, 0.172, 0, 0)
v_starting_states <- sample(x = v_states_names, size = num_i, replace = TRUE, prob = p_starting)

## Transition probabilities (per cycle)

# Common Mortality Parameters
p_HtoD  <- 0.00299                
rr_HF   <- 3.17                   
rr_2v1  <- 1.78                   
rr_3v1  <- 3.51                   
rr_4v1  <- 5.74                   
rp_sv   <- 0.2                       

# 1. Stabilizers
l_trans_probs_st <- list(
  "p_1to1" = 0.565, "p_1to2" = 0.392, "p_1to3" = 0.043, "p_1to4" = 0,
  "p_2to1" = 0.072, "p_2to2" = 0.751, "p_2to3" = 0.17,  "p_2to4" = 0.007,
  "p_3to1" = 0,     "p_3to2" = 0.29,  "p_3to3" = 0.678, "p_3to4" = 0.032,
  "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
  "p_HtoD" = p_HtoD, "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1, "rp_sv" = rp_sv
)

# 2. Silencers
l_trans_probs_si <- list(
  "p_1to1" = 0.570, "p_1to2" = 0.397, "p_1to3" = 0.033, "p_1to4" = 0,
  "p_2to1" = 0.092, "p_2to2" = 0.771, "p_2to3" = 0.13,  "p_2to4" = 0.007,
  "p_3to1" = 0,     "p_3to2" = 0.30,  "p_3to3" = 0.668, "p_3to4" = 0.032,
  "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
  "p_HtoD" = p_HtoD, "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1, "rp_sv" = rp_sv
)

# 3. Standard of Care (SoC) - For Discontinued
l_trans_probs_soc <- list(
  "p_1to1" = 0.538, "p_1to2" = 0.462, "p_1to3" = 0,     "p_1to4" = 0,
  "p_2to1" = 0.062, "p_2to2" = 0.763, "p_2to3" = 0.175, "p_2to4" = 0,
  "p_3to1" = 0.039, "p_3to2" = 0.216, "p_3to3" = 0.706, "p_3to4" = 0.039,
  "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
  "p_HtoD" = p_HtoD, "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1, "rp_sv" = rp_sv
)

## Discontinuation Probabilities
# Convert 2.5-year probability to 6-month probability
# 2.5 years = 5 cycles of 0.5 years.
# P(Disc in 2.5y) = 1 - (1 - p_cycle)^5
# p_cycle = 1 - (1 - P_2.5y)^(1/5)

p_disc_2.5y_st <- 0.212
p_disc_2.5y_si <- 0.093

p_disc_st <- 1 - (1 - p_disc_2.5y_st)^(1/5)
p_disc_si <- 1 - (1 - p_disc_2.5y_si)^(1/5)

## Cost and utility inputs
c_discount     <- 0.28
c_st           <- 112555
c_si           <- 238850

# Base Costs (SoC Costs) - Derived from original file logic
c_NYHA1_base   <- 2911
c_NYHA2_base   <- 4129
c_NYHA3_base   <- 6194
c_NYHA4_base   <- 10208
c_D            <- 0

# Costs with Treatment
c_NYHA1_st     <- c_NYHA1_base + (c_st * (1-c_discount))
c_NYHA2_st     <- c_NYHA2_base + (c_st * (1-c_discount))
c_NYHA3_st     <- c_NYHA3_base + (c_st * (1-c_discount))
c_NYHA4_st     <- c_NYHA4_base + (c_st * (1-c_discount))

c_NYHA1_si     <- c_NYHA1_base + (c_si * (1-c_discount))
c_NYHA2_si     <- c_NYHA2_base + (c_si * (1-c_discount))
c_NYHA3_si     <- c_NYHA3_base + (c_si * (1-c_discount))
c_NYHA4_si     <- c_NYHA4_base + (c_si * (1-c_discount))

c_age_cof <- 11.5
v_cost_coeffs <- c("age" = c_age_cof)

# Utilities (Unchanged)
u_NYHA1   <- 0.82
u_NYHA2   <- 0.729
u_NYHA3   <- 0.633
u_NYHA4   <- 0.333
u_Death   <- 0
u_age_cof <- -0.0018
ru_sv     <- -0.0015

v_util_coeffs <- c("age" = u_age_cof)
v_util_t_decs <- c("NYHA3 or 4" = ru_sv)

# Cost Vectors
v_states_costs_st  <- c("NYHA1" = c_NYHA1_st, "NYHA2" = c_NYHA2_st, "NYHA3" = c_NYHA3_st, "NYHA4" = c_NYHA4_st,"Death" = c_D)
v_states_costs_si  <- c("NYHA1" = c_NYHA1_si, "NYHA2" = c_NYHA2_si, "NYHA3" = c_NYHA3_si, "NYHA4" = c_NYHA4_si,"Death" = c_D)
v_states_costs_soc <- c("NYHA1" = c_NYHA1_base, "NYHA2" = c_NYHA2_base, "NYHA3" = c_NYHA3_base, "NYHA4" = c_NYHA4_base,"Death" = c_D)

# Utility Vectors
v_states_utilities <- c("NYHA1" = u_NYHA1, "NYHA2" = u_NYHA2, "NYHA3" = u_NYHA3, "NYHA4" = u_NYHA4, "Death" = u_Death)

# Hospitalization parameters
v_hosp_probs_st <- c("NYHA1" = 0.1683, "NYHA2" = 0.3107, "NYHA3" = 0.698, "NYHA4" = 0.8627, "Death" = 0)
v_hosp_probs_si <- c("NYHA1" = 0.1481, "NYHA2" = 0.2134, "NYHA3" = 0.61424, "NYHA4" = 0.7591, "Death" = 0)

v_hosp_costs    <- c("NYHA1" = 30584.15689, "NYHA2" = 17400.10754, "NYHA3" = 17694.56047, "NYHA4" = 21041.54226, "Death" = 0)
v_hosp_disutils <- c("NYHA1" = 0.023, "NYHA2" = 0.01, "NYHA3" = 0.027, "NYHA4" = 0.07, "Death" = 0)

### Running the simulation ###

# Run the simulation:
## For being treated with stabilizers
res_st <- run_microSimV_hosp_disc(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_names      = v_states_names,
  v_states_costs_arm  = v_states_costs_st,
  v_states_costs_soc  = v_states_costs_soc,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs_arm   = l_trans_probs_st,
  l_trans_probs_soc   = l_trans_probs_soc,
  p_disc              = p_disc_st,
  v_hosp_probs        = v_hosp_probs_st,
  v_hosp_costs        = v_hosp_costs,
  v_hosp_disutils     = v_hosp_disutils,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)

## For being treated with silencers
res_si <- run_microSimV_hosp_disc(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_names      = v_states_names,
  v_states_costs_arm  = v_states_costs_si,
  v_states_costs_soc  = v_states_costs_soc,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs_arm   = l_trans_probs_si,
  l_trans_probs_soc   = l_trans_probs_soc,
  p_disc              = p_disc_si,
  v_hosp_probs        = v_hosp_probs_si,
  v_hosp_costs        = v_hosp_costs,
  v_hosp_disutils     = v_hosp_disutils,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)

nb_st   <- res_st$mean_Dqalys   * wtp - res_st$mean_Dcosts
nb_st
nb_si   <- res_si$mean_Dqalys   * wtp - res_si$mean_Dcosts
nb_si
which.max(c("TTR stabilizers" = nb_st, "TTR silencers" = nb_si))

ICER <- (res_st$mean_Dcosts - res_si$mean_Dcosts)/(res_st$mean_Dqalys - res_si$mean_Dqalys)
ICER