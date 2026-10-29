# Vectorised Healthy-Sick-Dead microsimulation model with Discontinuation Logic

# clear R's session memory (Global Environment)
rm(list = ls())

### Defining model functions ###

# Helper function to calculate transition matrix for a given set of parameters
# This extracts the logic from the original update_probsV so it can be reused for different arms
# [MODIFIED] to handle age/sex-dependent background mortality
get_m_probs <- function(
    l_trans_probs,
    v_occupied_state,
    cycle_length,
    v_states_names,
    m_indi_features,
    params_gompertz_sex) {
  
  num_i <- nrow(m_indi_features)
  m_probs <- matrix(0,
                    nrow = num_i,
                    ncol = length(v_states_names),
                    dimnames = list(NULL, v_states_names))
  
  # Extract transition probabilities for non-death states
  p_1to1 <- l_trans_probs$p_1to1; p_1to2 <- l_trans_probs$p_1to2; p_1to3 <- l_trans_probs$p_1to3; p_1to4 <- l_trans_probs$p_1to4
  p_2to1 <- l_trans_probs$p_2to1; p_2to2 <- l_trans_probs$p_2to2; p_2to3 <- l_trans_probs$p_2to3; p_2to4 <- l_trans_probs$p_2to4
  p_3to1 <- l_trans_probs$p_3to1; p_3to2 <- l_trans_probs$p_3to2; p_3to3 <- l_trans_probs$p_3to3; p_3to4 <- l_trans_probs$p_3to4
  p_4to1 <- l_trans_probs$p_4to1; p_4to2 <- l_trans_probs$p_4to2; p_4to3 <- l_trans_probs$p_4to3; p_4to4 <- l_trans_probs$p_4to4
  
  # Extract hazard ratios for mortality
  rr_HF  <- l_trans_probs$rr_HF
  rr_2v1 <- l_trans_probs$rr_2v1
  rr_3v1 <- l_trans_probs$rr_3v1
  rr_4v1 <- l_trans_probs$rr_4v1
  
  # 1. Calculate age- and sex-specific background mortality HAZARD (annual) from Gompertz
  v_age <- m_indi_features[, "age"]
  h_H <- params_gompertz_sex$rate * exp(params_gompertz_sex$shape * v_age)
  
  # 2. Calculate state-specific mortality PROBABILITIES for the cycle
  # Apply hazard ratios to the background hazard to get state-specific annual hazards
  h_1 <- h_H * rr_HF
  h_2 <- h_1 * rr_2v1
  h_3 <- h_1 * rr_3v1
  h_4 <- h_1 * rr_4v1
  
  # Convert annual hazards to probabilities for the current cycle length
  p_1toD <- 1 - exp(-h_1 * cycle_length)
  p_2toD <- 1 - exp(-h_2 * cycle_length)
  p_3toD <- 1 - exp(-h_3 * cycle_length) # rp_sv logic removed
  p_4toD <- 1 - exp(-h_4 * cycle_length) # rp_sv logic removed
  
  # Ensure probabilities are capped at 1
  p_1toD[p_1toD > 1] <- 1; p_2toD[p_2toD > 1] <- 1
  p_3toD[p_3toD > 1] <- 1; p_4toD[p_4toD > 1] <- 1
  
  # 3. Populate the transition matrix based on current state
  idx_nyha1 <- which(v_occupied_state == "NYHA1"); idx_nyha2 <- which(v_occupied_state == "NYHA2")
  idx_nyha3 <- which(v_occupied_state == "NYHA3"); idx_nyha4 <- which(v_occupied_state == "NYHA4")
  idx_death <- which(v_occupied_state == "Death")
  
  # Probabilities to non-death states are scaled by (1 - p_death)
  if (length(idx_nyha1) > 0) {
    scaler <- 1 - p_1toD[idx_nyha1]
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

## Update Transition Probability function (Modified for Discontinuation)
update_probsV_disc <- function(
    v_states_names,
    v_occupied_state,
    l_trans_probs_arm, # Probabilities for those on treatment
    l_trans_probs_soc, # Probabilities for those discontinued (SoC)
    v_discontinued,    # Binary vector: 1 if discontinued, 0 otherwise
    m_indi_features,   # For age/sex to calculate mortality
    params_gompertz_sex, # Sex-specific gompertz parameters
    cycle_length = 0.5) {
  
  # 1. Calculate probabilities assuming everyone is on treatment
  m_probs_arm <- get_m_probs(l_trans_probs_arm, v_occupied_state, cycle_length, v_states_names, m_indi_features, params_gompertz_sex)
  
  # 2. Calculate probabilities assuming everyone is on SoC
  m_probs_soc <- get_m_probs(l_trans_probs_soc, v_occupied_state, cycle_length, v_states_names, m_indi_features, params_gompertz_sex)
  
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
### This function identifies the health state each individual will transition
### to in the next model cycle

sampleV <- function(m_trans_probs, 
  v_states_names) {
  
  # create an upper triangular matrix of ones
  m_upper_tri <- upper.tri(
    x = diag(ncol(m_trans_probs)), 
    diag = TRUE
  )

  # create matrix with row-wise cumulative transition probabilities
  m_cum_probs <- m_trans_probs %*% m_upper_tri
  colnames(m_cum_probs) <- v_states_names

  # ensure that the maximum cumulative probabilities are equal to 1
  if (any(m_cum_probs[, ncol(m_cum_probs)] > 1.0000001)) { # Slight tolerance
    stop("Error in multinomial sampling: probabilities do not sum to 1")
  }

  # sample random values from Uniform standard distribution for each individual
  v_rand_values <- runif(n = nrow(m_trans_probs))

  # repeat each sampled value to have as many copies as the number of states
  m_rand_values <- matrix(
    data = rep(
      x = v_rand_values, 
      each = length(v_states_names)), 
      nrow = nrow(m_trans_probs), 
      ncol = length(v_states_names), 
      byrow = TRUE
    )
  
  # identify transitions, compare random samples to cumulative probabilities
  m_transitions <- m_rand_values > m_cum_probs # transitions from first state

  # sum transitions to identify health state in next cycle
  v_transitions <- rowSums(m_transitions)

  # sum transitions to identify health state in next cycle
  v_health_states <- v_states_names[1 + v_transitions]
  return(v_health_states)
}

## Calculate Costs function (Modified for Discontinuation)
### This function estimates the costs at every cycle based on the health state
### occupied by each individuals at cycle 't' and relevant individuals features
calc_costsV_disc <- function (
    v_occupied_state,
    v_states_costs_arm, # Costs including drug
    v_states_costs_soc, # Costs excluding drug (SoC)
    v_discontinued) {
  
  # Helper to calculate costs for a specific cost vector
  calc_costs_internal <- function(states, base_costs) {
    v_c <- rep(NA, length(states))
    v_c[states == "NYHA1"] <- base_costs["NYHA1"]
    v_c[states == "NYHA2"] <- base_costs["NYHA2"]
    v_c[states == "NYHA3"] <- base_costs["NYHA3"]
    v_c[states == "NYHA4"] <- base_costs["NYHA4"]
    v_c[states == "Death"] <- base_costs["Death"]
    return(v_c)
  }
  
  # Calculate costs assuming on treatment
  v_costs_arm <- calc_costs_internal(v_occupied_state, v_states_costs_arm)
  
  # Calculate costs assuming SoC (discontinued)
  v_costs_soc <- calc_costs_internal(v_occupied_state, v_states_costs_soc)
  
  # Combine
  v_final_costs <- v_costs_arm * (1 - v_discontinued) + v_costs_soc * v_discontinued
  
  return(v_final_costs)
}

## Calculate Health Outcomes function (Unchanged)
### This function estimates the Quality Adjusted Life Years (QALYs) at every
### cycle based on the health state occupied by each individuals at cycle 't',
### time spent in the states and the cycle_length (measured in years)
calc_effsV <- function (
  v_occupied_state,
  v_states_utilities,
  cycle_length = 1) {
  
  # estimate utilities based on occupied state 
  v_state_utility <- rep(NA, length(v_occupied_state))
  v_state_utility[v_occupied_state == "NYHA1"]  <- v_states_utilities["NYHA1"]
  v_state_utility[v_occupied_state == "NYHA2"]  <- v_states_utilities["NYHA2"]
  v_state_utility[v_occupied_state == "NYHA3"]  <- v_states_utilities["NYHA3"]
  v_state_utility[v_occupied_state == "NYHA4"]  <- v_states_utilities["NYHA4"]
  v_state_utility[v_occupied_state == "Death"]  <- v_states_utilities["Death"]
  
  # calculate Quality Adjusted Life Years (QALYs)
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
    v_states_utilities,
    l_trans_probs_arm,  # Probs for treatment arm
    l_trans_probs_soc,  # Probs for SoC
    p_disc,             # Probability of discontinuation per cycle
    params_gompertz_sex,# Sex-specific Gompertz parameters for background mortality
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
  m_States[, 1] <- v_starting_states
  m_Hosp[, 1] <- 0 
  
  # Initial costs (Cycle 0) - assume everyone starts on treatment
  m_Costs[, 1]  <- calc_costsV_disc(
    v_occupied_state   = m_States[, 1],
    v_states_costs_arm = v_states_costs_arm,
    v_states_costs_soc = v_states_costs_soc,
    v_discontinued     = v_discontinued
  )
  
  m_Effs[, 1]   <- calc_effsV(
    v_occupied_state   = m_States[, 1],
    v_states_utilities = v_states_utilities,
    cycle_length       = cycle_length
  )
  
  for (t in 1:num_cycles) {
    
    # NEW: Select the appropriate transition probabilities for the current cycle
    # Use time-varying probabilities for the first 5 cycles (up to 30 months)
    # and the 30-month probabilities thereafter.
    # The index is capped at the number of time-varying probability sets.
    current_cycle_index <- min(t, length(l_trans_probs_arm))
    l_trans_probs_soc_t <- l_trans_probs_soc[[current_cycle_index]] # NEW: Select time-dependent SoC probs
    l_trans_probs_arm_t <- l_trans_probs_arm[[current_cycle_index]]
    
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
      l_trans_probs_arm = l_trans_probs_arm_t, # Use time-dependent probabilities
      l_trans_probs_soc = l_trans_probs_soc_t, # NEW: Use time-dependent SoC probabilities
      v_discontinued    = v_discontinued,
      m_indi_features   = m_indi_features,
      params_gompertz_sex = params_gompertz_sex,
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
    
# 4a. Update discontinuation status for patients entering NYHA4
    # Patients who move to the NYHA4 state are considered to have discontinued treatment.
    # This is a model assumption. Once discontinued, they remain so.
    # This will affect their costs for the current cycle (step 6) and transition probabilities for the next cycle.
    v_discontinued[m_States[, t + 1] == "NYHA4"] <- 1

    # 5. Update time in state and age
    # Age individuals by the cycle length
    m_indi_features[, "age"] <- m_indi_features[, "age"] + cycle_length
    
    # 6. Calculate Costs
    # Uses v_discontinued to remove drug costs if applicable
    m_Costs[, t + 1]  <- calc_costsV_disc(
      v_occupied_state   = m_States[, t + 1],
      v_states_costs_arm = v_states_costs_arm,
      v_states_costs_soc = v_states_costs_soc,
      v_discontinued     = v_discontinued
    ) + l_hosp_results$costs
    
    # 7. Calculate QALYs
    m_Effs[, t + 1]   <- calc_effsV(
      v_occupied_state   = m_States[, t + 1],
      v_states_utilities = v_states_utilities,
      cycle_length       = cycle_length
    ) - l_hosp_results$disutils
  }
  
  # Discounting and Summary
  v_c_dsc_wts <- calc_discount_wts(discount_rate_costs, num_cycles, cycle_length)
  v_e_dsc_wts <- calc_discount_wts(discount_rate_QALYs, num_cycles, cycle_length)
  
  # Apply Life-Table Correction to Undiscounted Costs/QALYs
  # Logic: (Current + Next) / 2. Vectorized for matrices: (M[, -last] + M[, -1]) / 2
  m_Costs_LT <- (m_Costs[, -ncol(m_Costs)] + m_Costs[, -1]) / 2
  m_Effs_LT  <- (m_Effs[, -ncol(m_Effs)]  + m_Effs[, -1]) / 2
  
  v_total_costs <- rowSums(m_Costs_LT)
  v_total_qalys <- rowSums(m_Effs_LT)
  mean_costs    <- mean(v_total_costs)
  mean_qalys    <- mean(v_total_qalys)
  
  # Apply Life-Table Correction to Discounted Costs/QALYs
  # First, apply discount weights to the raw matrices (broadcast column-wise)
  m_Costs_D <- t(t(m_Costs) * v_c_dsc_wts)
  m_Effs_D  <- t(t(m_Effs)  * v_e_dsc_wts)
  
  # Then apply correction to the discounted values and sum
  v_total_Dcosts <- rowSums((m_Costs_D[, -ncol(m_Costs_D)] + m_Costs_D[, -1]) / 2)
  v_total_Dqalys <- rowSums((m_Effs_D[, -ncol(m_Effs_D)]  + m_Effs_D[, -1]) / 2)
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
# Create separate populations for males and females to account for different mortality
m_indi_features_male <- cbind(
  "age" = rnorm(n = num_i, mean = mean_age, sd = sd_age)
)
m_indi_features_female <- cbind(
  "age" = rnorm(n = num_i, mean = mean_age, sd = sd_age)
)

## Health states
v_states_names <- c("NYHA1","NYHA2", "NYHA3", "NYHA4", "Death")
p_starting <- c(0.108, 0.72, 0.172, 0, 0)
v_starting_states <- sample(x = v_states_names, size = num_i, replace = TRUE, prob = p_starting)

## Transition and Mortality Parameters

# --- Background Mortality Parameters (from workingTRACE_mortality_fitting.R) --- #
# !!! IMPORTANT !!!
# Please replace these placeholder values with the actual values from your fitting script output.
params_gompertz <- list(
  male = list(
    shape = 0.1014027, 
    rate  = 1.934313e-05    
  ),
  female = list(
    shape = 0.1089245, 
    rate  = 7.472652e-06     
  )
)

# Hazard Ratios for mortality (relative to background mortality)
rr_HF   <- 3.17                   
rr_2v1  <- 1.78                   
rr_3v1  <- 3.51                   
rr_4v1  <- 5.74                   
# NOTE: p_HtoD and rp_sv are removed as they are replaced by the Gompertz model

# NEW: Time-varying transition probabilities for treatment arms
# A list of lists, where each inner list corresponds to a 6-month cycle.
# Probabilities for cycles > 5 (after 30 months) will use the 5th set of probabilities.

# 1. Stabilizers (Tafamidis) - Time-Varying
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
  # Cycle 5 (30 months) - This will be used for all subsequent cycles
  list("p_1to1" = 0.368, "p_1to2" = 0.368, "p_1to3" = 0.211, "p_1to4" = 0.053,
       "p_2to1" = 0.115, "p_2to2" = 0.598, "p_2to3" = 0.287, "p_2to4" = 0,
       "p_3to1" = 0,     "p_3to2" = 0.3,   "p_3to3" = 0.633, "p_3to4" = 0.067,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1)
)

# 2. Silencers - Time-Varying (ARBITRARY PLACEHOLDER VALUES)
# NOTE: These values are illustrative, created to be slightly more effective than stabilizers.
l_trans_probs_si_tv <- list(
  # Cycle 1 (6 months)
  list("p_1to1" = 0.570, "p_1to2" = 0.397, "p_1to3" = 0.033, "p_1to4" = 0,
       "p_2to1" = 0.092, "p_2to2" = 0.771, "p_2to3" = 0.13,  "p_2to4" = 0.007,
       "p_3to1" = 0,     "p_3to2" = 0.30,  "p_3to3" = 0.668, "p_3to4" = 0.032,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 2 (12 months)
  list("p_1to1" = 0.54, "p_1to2" = 0.46, "p_1to3" = 0,     "p_1to4" = 0,
       "p_2to1" = 0.08, "p_2to2" = 0.77, "p_2to3" = 0.14,  "p_2to4" = 0.01,
       "p_3to1" = 0.02, "p_3to2" = 0.40, "p_3to3" = 0.56,  "p_3to4" = 0.02,
       "p_4to1" = 0,    "p_4to2" = 0,    "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 3 (18 months)
  list("p_1to1" = 0.40, "p_1to2" = 0.47, "p_1to3" = 0.13,  "p_1to4" = 0,
       "p_2to1" = 0.10, "p_2to2" = 0.71, "p_2to3" = 0.19,  "p_2to4" = 0,
       "p_3to1" = 0.03, "p_3to2" = 0.28, "p_3to3" = 0.67,  "p_3to4" = 0.02,
       "p_4to1" = 0,    "p_4to2" = 0,    "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 4 (24 months)
  list("p_1to1" = 0.52, "p_1to2" = 0.29, "p_1to3" = 0.14,  "p_1to4" = 0.05,
       "p_2to1" = 0.11, "p_2to2" = 0.69, "p_2to3" = 0.19,  "p_2to4" = 0.01,
       "p_3to1" = 0.01, "p_3to2" = 0.28, "p_3to3" = 0.61,  "p_3to4" = 0.10,
       "p_4to1" = 0,    "p_4to2" = 0,    "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1),
  # Cycle 5 (30 months)
  list("p_1to1" = 0.38, "p_1to2" = 0.36, "p_1to3" = 0.21,  "p_1to4" = 0.05,
       "p_2to1" = 0.12, "p_2to2" = 0.61, "p_2to3" = 0.27,  "p_2to4" = 0,
       "p_3to1" = 0.01, "p_3to2" = 0.31, "p_3to3" = 0.62,  "p_3to4" = 0.06,
       "p_4to1" = 0,    "p_4to2" = 0,    "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1)
)

# 3. Standard of Care (SoC) - Time-Varying
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
  # Cycle 5 (30 months) - This will be used for all subsequent cycles
  list("p_1to1" = 0.125, "p_1to2" = 0.75,  "p_1to3" = 0.125, "p_1to4" = 0,
       "p_2to1" = 0.016, "p_2to2" = 0.5,   "p_2to3" = 0.452, "p_2to4" = 0.032,
       "p_3to1" = 0,     "p_3to2" = 0.28,  "p_3to3" = 0.6,   "p_3to4" = 0.12,
       "p_4to1" = 0,     "p_4to2" = 0,     "p_4to3" = 0,     "p_4to4" = 1,
       "rr_HF" = rr_HF, "rr_2v1" = rr_2v1, "rr_3v1" = rr_3v1, "rr_4v1" = rr_4v1)
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

# Utilities (Unchanged)
u_NYHA1   <- 0.82
u_NYHA2   <- 0.729
u_NYHA3   <- 0.633
u_NYHA4   <- 0.333
u_Death   <- 0

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

cat("\n--- Running Simulation for MALE Cohort ---\n")

## Stabilizers vs SoC for MALES
res_st_male <- run_microSimV_hosp_disc(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features_male,
  v_states_names      = v_states_names,
  v_states_costs_arm  = v_states_costs_st,
  v_states_costs_soc  = v_states_costs_soc,
  v_states_utilities  = v_states_utilities,
  l_trans_probs_arm   = l_trans_probs_st_tv,
  l_trans_probs_soc   = l_trans_probs_soc_tv,
  p_disc              = p_disc_st,
  params_gompertz_sex = params_gompertz$male,
  v_hosp_probs        = v_hosp_probs_st,
  v_hosp_costs        = v_hosp_costs,
  v_hosp_disutils     = v_hosp_disutils,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed # Use same seed for comparability
)

## Silencers vs SoC for MALES
res_si_male <- run_microSimV_hosp_disc(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features_male,
  v_states_names      = v_states_names,
  v_states_costs_arm  = v_states_costs_si,
  v_states_costs_soc  = v_states_costs_soc,
  v_states_utilities  = v_states_utilities,
  l_trans_probs_arm   = l_trans_probs_si_tv,
  l_trans_probs_soc   = l_trans_probs_soc_tv,
  p_disc              = p_disc_si,
  params_gompertz_sex = params_gompertz$male,
  v_hosp_probs        = v_hosp_probs_si,
  v_hosp_costs        = v_hosp_costs,
  v_hosp_disutils     = v_hosp_disutils,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed # Use same seed for comparability
)

cat("\n--- Running Simulation for FEMALE Cohort ---\n")

## Stabilizers vs SoC for FEMALES
res_st_female <- run_microSimV_hosp_disc(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features_female,
  v_states_names      = v_states_names,
  v_states_costs_arm  = v_states_costs_st,
  v_states_costs_soc  = v_states_costs_soc,
  v_states_utilities  = v_states_utilities,
  l_trans_probs_arm   = l_trans_probs_st_tv,
  l_trans_probs_soc   = l_trans_probs_soc_tv,
  p_disc              = p_disc_st,
  params_gompertz_sex = params_gompertz$female,
  v_hosp_probs        = v_hosp_probs_st,
  v_hosp_costs        = v_hosp_costs,
  v_hosp_disutils     = v_hosp_disutils,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed # Use same seed for comparability
)

## Silencers vs SoC for FEMALES
res_si_female <- run_microSimV_hosp_disc(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features_female,
  v_states_names      = v_states_names,
  v_states_costs_arm  = v_states_costs_si,
  v_states_costs_soc  = v_states_costs_soc,
  v_states_utilities  = v_states_utilities,
  l_trans_probs_arm   = l_trans_probs_si_tv,
  l_trans_probs_soc   = l_trans_probs_soc_tv,
  p_disc              = p_disc_si,
  params_gompertz_sex = params_gompertz$female,
  v_hosp_probs        = v_hosp_probs_si,
  v_hosp_costs        = v_hosp_costs,
  v_hosp_disutils     = v_hosp_disutils,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed # Use same seed for comparability
)

### Analyzing the results ###

cat("\n--- MALE COHORT RESULTS ---\n")
nb_st_male <- res_st_male$mean_Dqalys * wtp - res_st_male$mean_Dcosts
nb_si_male <- res_si_male$mean_Dqalys * wtp - res_si_male$mean_Dcosts
cat("Net Monetary Benefit (Stabilizers):", nb_st_male, "\n")
cat("Net Monetary Benefit (Silencers):", nb_si_male, "\n")

ICER_male <- (res_st_male$mean_Dcosts - res_si_male$mean_Dcosts) / (res_st_male$mean_Dqalys - res_si_male$mean_Dqalys)
cat("ICER (Stabilizers vs Silencers):", ICER_male, "\n")

cat("\n--- FEMALE COHORT RESULTS ---\n")
nb_st_female <- res_st_female$mean_Dqalys * wtp - res_st_female$mean_Dcosts
nb_si_female <- res_si_female$mean_Dqalys * wtp - res_si_female$mean_Dcosts
cat("Net Monetary Benefit (Stabilizers):", nb_st_female, "\n")
cat("Net Monetary Benefit (Silencers):", nb_si_female, "\n")

ICER_female <- (res_st_female$mean_Dcosts - res_si_female$mean_Dcosts) / (res_st_female$mean_Dqalys - res_si_female$mean_Dqalys)
cat("ICER (Stabilizers vs Silencers):", ICER_female, "\n")
#------------------------------------------------------------------------------#

## @knitr micro_run_microSimV_diagram

if (requireNamespace("DiagrammeR", quietly = TRUE)) {
  DiagrammeR::grViz("
  digraph flowchart {
    node [fontname = 'Helvetica', shape = box, style=filled, fillcolor='grey', fontsize=16]
    edge [fontname = 'Helvetica']

    input [shape = box, label = 'Inputs: v_starting_states, num_i, num_cycles, m_indi_features, v_states_names, v_states_costs, v_cost_coeffs, v_states_utilities, v_util_coeffs, v_util_t_decs, l_trans_probs, cycle_length, starting_seed', style=filled, fillcolor='yellow', width=3.5]
    initialize_matrices [label = 'Initialize m_States, m_Costs, m_Effs', style=filled, fillcolor='palegreen']
    calc_initial_costs [label = 'Set seed, get starting health state, and calculate initial costs and QALYs', style=filled, fillcolor='palegreen']
    loop_cycles [shape = diamond, style=filled, fillcolor='skyblue', fontsize=24, fontname='Helvetica-Bold', label = 'For each cycle t']
    update_probs [label = 'Update transition probabilities for cycle (t)', style=filled, fillcolor='palegreen']
    sample_state [label = 'Sample health state for cycle (t + 1)', style=filled, fillcolor='palegreen']
    calculate_cycle_costs [label = 'Calculate payoffs (costs and QALYs) for cycle (t + 1)', style=filled, fillcolor='palegreen']
    update_time [label = 'Update time in state for cycle (t + 1)', style=filled, fillcolor='palegreen']
    update_age [label = 'Advance age if alive for cycle (t + 1)', style=filled, fillcolor='palegreen']
    check_cycles [label = 'Was this the last cycle?', shape = diamond, style=filled, fillcolor='skyblue', fontsize=24, fontname='Helvetica-Bold']
    summarize_results [label = 'Discount and summarize costs and QALYs', style=filled, fillcolor='palegreen']
    return_results [shape = ellipse, label = 'Return results', style=filled, fillcolor='yellow']

    input -> initialize_matrices
    initialize_matrices -> calc_initial_costs
    calc_initial_costs -> loop_cycles
    loop_cycles -> update_probs
    update_probs -> sample_state
    sample_state ->  update_time -> update_age -> calculate_cycle_costs
    calculate_cycle_costs -> check_cycles
    check_cycles -> loop_cycles [label = 'No\nNext cycle']
    check_cycles -> summarize_results [label = 'Yes']
    summarize_results -> return_results
  }
  ")
} else {
  message("DiagrammeR package not installed. Flowchart skipped.")
}
