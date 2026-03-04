# Vectorised Healthy-Sick-Dead microsimulation model

# clear R's session memory (Global Environment)
rm(list = ls())



### Defining model functions ###

# Define model functions

## Update Transition Probability function
### This function updates the transition probabilities at every cycle based on
### the health state occupied by each individual at cycle 't' and the time spent
### in the states
update_probsV <- function(
    v_states_names,
    v_occupied_state,
    l_trans_probs,
    v_time_in_state) {
  
  with(
    data = l_trans_probs,
    expr = {

      # update probabilities of death after first converting them to rates and applying the rate ratio
      p_1to1  <- 0.565                      
      p_2to1  <- 0.072
      p_3to1  <- 0
      p_4to1  <- 0
      p_1to2  <- 0.392                      
      p_2to2  <- 0.751
      p_3to2  <- 0.29
      p_4to2  <- 0
      p_1to3  <- 0.043                      
      p_2to3  <- 0.17
      p_3to3  <- 0.678
      p_4to3  <- 0
      p_1to4  <- 0                      
      p_2to4  <- 0.007
      p_3to4  <- 0.032
      p_4to4  <- 1
      p_Dto1  <- 0
      p_Dto2  <- 0
      p_Dto3  <- 0
      p_Dto4  <- 0
      p_DtoD  <- 1
      p_HtoD  <- 0.00299                # all cause mortality, H = healthy
      rr_HF   <- 3.17                   # HF vs healthy (HR)
      p_1toD  <- rr_HF*p_HtoD
      rr_2v1  <- 1.78                   # I vs II (HR)
      rr_3v1  <- 3.51                   # I vs III (HR)
      rr_4v1  <- 5.74                   # I vs IV (HR)
      r_1toD  <- -log(1- p_1toD)        #rate of death in NYHA1
      r_2toD  <-  rr_2v1*r_1toD         #rate of death in NYHA2
      p_2toD  <-  1 - exp(-r_2toD*0.5)      #probability to die in NYHA2
      r_3toD  <-  rr_3v1*r_1toD         #rate of death in NYHA2
      r_4toD  <-  rr_4v1*r_1toD         #rate of death in NYHA2
      rp_sv   <-  0.2                       # increase in mortality rate with every additional year being NYHA 3 or 4 

            # calculate p_S1D conditional on current state and duration of being sick
      p_3toD  <-  1 - exp(-r_3toD*0.5*(1 + v_time_in_state[v_occupied_state == "NYHA3"] * rp_sv))      #probability to die in NYHA3
      p_4toD  <-  1 - exp(-r_4toD*0.5*(1 + v_time_in_state[v_occupied_state == "NYHA4"] * rp_sv))      #probability to die in NYHA4
      
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

      # sanity check
      ifelse(
        test = rowSums(m_probs) == 1,               # check if the transition probabilities add up to 1
        yes = return(m_probs),                      # return the transition probabilities
        no = print("Probabilities do not sum to 1") # or produce an error
      )
    }
  )
}

## Sample Health States function
### This function identifies the health state each individual will transition
### to in the next model cycle
sampleV <- function(
    m_trans_probs,
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
  if (any(m_cum_probs[, ncol(m_cum_probs)] > 1.000000)) {
    stop("Error in multinomial sampling: probabilities do not sum to 1")
  }

  # sample random values from Uniform standard distribution for each individual
  v_rand_values <- runif(n = nrow(m_trans_probs))

  # repeat each sampled value to have as many copies as the number of states
  m_rand_values <- matrix(
    data  = rep(
      x = v_rand_values,
      each = length(v_states_names)
    ),
    nrow  = nrow(m_trans_probs),
    ncol  = length(v_states_names),
    byrow = TRUE
  )

  # identify transitions, compare random samples to cumulative probabilities
  m_transitions <- m_rand_values > m_cum_probs # transitions from first state

  # sum transitions to identify health state in next cycle
  v_transitions <- rowSums(m_transitions)

  # identify health state to which each individual is transitioning
  v_health_states <- v_states_names[1 + v_transitions]

  return(v_health_states)
}

## Calculate Costs function
### This function estimates the costs at every cycle based on the health state
### occupied by each individuals at cycle 't' and relevant individuals features
calc_costsV <- function (
    v_occupied_state,
    v_states_costs,
    m_indi_features,
    v_cost_coeffs) {

  # calculate individual-specific costs based on costs regression coefficients
  v_indi_costs <- m_indi_features %*% v_cost_coeffs

  # estimate costs based on occupied state
  v_state_costs                           <- rep(NA, length(v_occupied_state))
  v_state_costs[v_occupied_state == "NYHA1"]  <- v_states_costs["NYHA1"]                                   # update the cost if in NYHA1
  v_state_costs[v_occupied_state == "NYHA2"]  <- v_states_costs["NYHA2"]                                   # update the cost if in NYHA2
  v_state_costs[v_occupied_state == "NYHA3"]  <- v_states_costs["NYHA3"] + v_indi_costs[v_occupied_state == "NYHA3"] # update the cost if in NYHA3
  v_state_costs[v_occupied_state == "NYHA4"]  <- v_states_costs["NYHA4"] + v_indi_costs[v_occupied_state == "NYHA4"] # update the cost if in NYHA4
  v_state_costs[v_occupied_state == "Death"]  <- v_states_costs["Death"]                                   # update the cost if Dead
 
  return(v_state_costs)                                                                                    # return the costs
}

## Calculate Health Outcomes function
### This function estimates the Quality Adjusted Life Years (QALYs) at every
### cycle based on the health state occupied by each individuals at cycle 't',
### time spent in the states and the cycle_length (measured in years)
calc_effsV <- function (
    v_occupied_state,
    v_states_utilities,
    m_indi_features,
    v_util_coeffs,
    v_util_t_decs,
    v_time_in_state,
    cycle_length = 0.5) {

  # calculate individual-specific utility decrements based on utilities regression coefficients
  v_ind_decrement <- (m_indi_features %*% v_util_coeffs)[,1]

  # calculate time-dependent state-specific utility decrements
  time_decrement <- rep(0, length(v_occupied_state))
  time_decrement[v_occupied_state == "NYHA3"] <- v_util_t_decs["NYHA3 or 4"] * v_time_in_state[v_occupied_state == "NYHA3"]
  time_decrement[v_occupied_state == "NYHA4"] <- v_util_t_decs["NYHA3 or 4"] * v_time_in_state[v_occupied_state == "NYHA4"]

  # estimate total decrements
  decrement <- v_ind_decrement + time_decrement

  # estimate utilities based on occupied state
  v_state_utility                           <- rep(NA, length(v_occupied_state))
  v_state_utility[v_occupied_state == "NYHA1"]  <- v_states_utilities["NYHA1"]                                            # update the utility if in NYHA1
  v_state_utility[v_occupied_state == "NYHA2"]  <- v_states_utilities["NYHA2"]                                            # update the utility if in NYHA2
  v_state_utility[v_occupied_state == "NYHA3"]  <- v_states_utilities["NYHA3"]  + decrement[v_occupied_state == "NYHA3"]  # update the utility if in NYHA3
  v_state_utility[v_occupied_state == "NYHA4"]  <- v_states_utilities["NYHA4"]  + decrement[v_occupied_state == "NYHA4"]  # update the utility if in NYHA4
  v_state_utility[v_occupied_state == "Death"]  <- v_states_utilities["Death"]                                            # update the utility if dead

  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  v_state_utility * cycle_length                                                                    # calculate the QALYs during cycle `t`

  return(QALYs)                                                                                               # return the QALYs
}

## Calculate Discount Weights function
### This function estimates the discount weights to be applied to the outputs of
### each cycle in order to scale their values back to the present.
calc_discount_wts <- function(
    discount_rate,
    num_cycles,
    cycle_length) {

  # calculate discount weights based on the number & length (in years) of cycles
  v_discount_wts <- 1 / (1 + discount_rate) ^ ((0:num_cycles) * cycle_length)

  return(v_discount_wts)
}

#------------------------------------------------------------------------------#

## @knitr micro_run_microSimV_foo

## Run Microsimulation function
### This function runs the microsimulation function of the Healthy-Sick-Dead
### model vectorised and extended to handle non-homogeneous populations.
run_microSimV <- function(
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
    discount_rate_costs,
    discount_rate_QALYs,
    cycle_length = 0.5,
    starting_seed = 1) {

  # create matrices to capture states' names, associated costs and QALYs
  m_States <- m_Costs <- m_Effs <-  matrix(
    nrow = num_i,
    ncol = num_cycles + 1,
    dimnames = list(paste("ind",   1:num_i,    sep ="_"),
                    paste("cycle", 0:num_cycles, sep ="_"))
  )

  # set the seed for every individual for the random number generator
  set.seed(starting_seed)

  # initialize parameter tracking time in current state
  v_time_in_state <- rep(1, times = num_i)

  # get the initial health state
  m_States[, 1] <- v_starting_states

  # calculate the costs incurred in their starting health state
  m_Costs[, 1]  <- calc_costsV(
    v_occupied_state = m_States[, 1],
    v_states_costs   = v_states_costs,
    m_indi_features  = m_indi_features,
    v_cost_coeffs    = v_cost_coeffs
  )

  # calculate the QALYs accrued in their starting health state
  m_Effs[, 1]   <- calc_effsV(
    v_occupied_state   = m_States[, 1],
    v_states_utilities = v_states_utilities,
    m_indi_features    = m_indi_features,
    v_util_coeffs      = v_util_coeffs,
    v_util_t_decs      = v_util_t_decs,
    v_time_in_state    = v_time_in_state,
    cycle_length       = cycle_length
  )

  # for each 't' of the 'num_cycles' cycles:
  for (t in 1:num_cycles) {
    # update the transition probabilities at cycle 't'
    m_trans_probs     <- update_probsV(
      v_states_names   = v_states_names,
      v_occupied_state = m_States[, t],
      l_trans_probs    = l_trans_probs,
      v_time_in_state  = v_time_in_state
    )

    # sample the health state at 't + 1'
    m_States[, t + 1] <- sampleV(
      m_trans_probs  = m_trans_probs,
      v_states_names = v_states_names
    )

    # keep track of time in state at 't + 1'
    stayed                   <- m_States[, t] == m_States[, t + 1] # check if remains in current state at 't + 1'
    v_time_in_state[stayed]  <- v_time_in_state[stayed] + 1        # increment time spent in state
    v_time_in_state[!stayed] <- 1                                  # reset time once transitioned

    # keep track of time in the model
    m_indi_features[, "age"] <- m_indi_features[, "age"] + 1

    # calculate the costs incurred in their 't + 1' health state
    m_Costs[, t + 1]  <- calc_costsV(
      v_occupied_state = m_States[, t + 1],
      v_states_costs   = v_states_costs,
      m_indi_features  = m_indi_features,
      v_cost_coeffs    = v_cost_coeffs
    )

    # calculate the QALYs accrued in their 't + 1' health state
    m_Effs[, t + 1]   <- calc_effsV(
      v_occupied_state   = m_States[, t + 1],
      v_states_utilities = v_states_utilities,
      m_indi_features    = m_indi_features,
      v_util_coeffs      = v_util_coeffs,
      v_util_t_decs      = v_util_t_decs,
      v_time_in_state    = v_time_in_state,
      cycle_length       = cycle_length
    )

  } # close the loop for the cycles 't'

  # Calculate discount weights for both outcomes:
  v_c_dsc_wts <- calc_discount_wts(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
  v_e_dsc_wts <- calc_discount_wts(
    discount_rate = discount_rate_QALYs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
  # Compute costs and QALYs:
  v_total_costs <- rowSums(m_Costs)         # calculate total costs per individual
  v_total_qalys <- rowSums(m_Effs)          # calculate total QALYs per individual
  mean_costs    <- mean(v_total_costs)      # calculate average costs
  mean_qalys    <- mean(v_total_qalys)      # calculate average QALYs

  # Compute discounted costs and QALYs:
  v_total_Dcosts <- m_Costs %*% v_c_dsc_wts # calculate total discounted costs per individual
  v_total_Dqalys <- m_Effs  %*% v_e_dsc_wts # calculate total discounted QALYs per individual
  mean_Dcosts    <- mean(v_total_Dcosts)    # calculate average discounted costs
  mean_Dqalys    <- mean(v_total_Dqalys)    # calculate average discounted QALYs

  # store the results in a list:
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

  # return the results
  return(results)
}

#------------------------------------------------------------------------------#

## @knitr micro_run_microSimV_demo

### Defining model parameters ###

# Define model inputs
## General parameters
num_i               <- 1e6               # number of simulated individuals
num_cycles          <- 50                # time horizon if each cycle is a year long
cycle_length        <- 0.5               # length of cycle (in years)
seed                <- 1234              # random number generator state
wtp                 <- 150000            # Willingness to pay for each QALY ($)
discount_rate_costs <- 0.03              # annual discount rate for costs
discount_rate_QALYs <- 0.03              # annual discount rate for health outcomes

## Population characteristics/features
mean_age            <- 77                # mean age in the simulated population
sd_age              <- 5                 # standard deviation of the age in the simulated population
set.seed(seed)                           # set a seed to ensure reproducible samples
m_indi_features     <- cbind(            # simulate individuals characteristics
  "age" = rnorm(                         # get random samples for 'age' from a normal distribution
    n = num_i,
    mean = mean_age,
    sd = sd_age
  )
)

## Health states
v_states_names <- c("NYHA1","NYHA2", "NYHA3", "NYHA4", "Death")       # the model states: Healthy (H), Sick (S1), Dead (D)
p_starting <- c(0.108, 0.72, 0.172, 0, 0)
v_starting_states <- sample(
  x       = v_states_names,
  size    = num_i,
  replace = TRUE,
  prob    = p_starting
  )     # Assign individuals to starting states based on the defined distribution

## Transition probabilities (per cycle)
# Stabilizers
p_1to1  <- 0.565                      
p_2to1  <- 0.072
p_3to1  <- 0
p_4to1  <- 0
p_1to2  <- 0.392                      
p_2to2  <- 0.751
p_3to2  <- 0.29
p_4to2  <- 0
p_1to3  <- 0.043                      
p_2to3  <- 0.17
p_3to3  <- 0.678
p_4to3  <- 0
p_1to4  <- 0                      
p_2to4  <- 0.007
p_3to4  <- 0.032
p_4to4  <- 1
p_Dto1  <- 0
p_Dto2  <- 0
p_Dto3  <- 0
p_Dto4  <- 0
p_DtoD  <- 1
p_HtoD  <- 0.00299                # all cause mortality, H = healthy
rr_HF   <- 3.17                   # HF vs healthy (HR)
p_1toD  <- rr_HF*p_HtoD
rr_2v1  <- 1.78                   # I vs II (HR)
rr_3v1  <- 3.51                   # I vs III (HR)
rr_4v1  <- 5.74                   # I vs IV (HR)
r_1toD  <- -log(1- p_1toD)        #rate of death in NYHA1
r_2toD  <-  rr_2v1*r_1toD         #rate of death in NYHA2
r_3toD  <-  rr_3v1*r_1toD         #rate of death in NYHA2
r_4toD  <-  rr_4v1*r_1toD         #rate of death in NYHA2
p_2toD  <-  1 - exp(-r_2toD*0.5)      #probability to die in NYHA2
p_3toD  <-  1 - exp(-r_3toD*0.5)      #probability to die in NYHA3
p_4toD  <-  1 - exp(-r_4toD*0.5)      #probability to die in NYHA4
rp_sv   <-  0.2                       # increase in mortality rate with every additional year being NYHA 3 or 4 

l_trans_probs <- list(                # pack the transition probabilities and rates in a list
  "p_1to1" = p_1to1*(1-p_1toD),
  "p_1to2" = p_1to2*(1-p_1toD),
  "p_1to3" = p_1to3*(1-p_1toD),
  "p_1to4" = p_1to4*(1-p_1toD),
  "p_1toD" = p_1toD,
  "p_2to1" = p_2to1*(1-p_2toD),
  "p_2to2" = p_2to2*(1-p_2toD),
  "p_2to3" = p_2to3*(1-p_2toD),
  "p_2to4" = p_2to4*(1-p_2toD),
  "p_2toD" = p_2toD,
  "p_3to1" = p_3to1*(1-p_3toD),
  "p_3to2" = p_3to2*(1-p_3toD),
  "p_3to3" = p_3to3*(1-p_3toD),
  "p_3to4" = p_3to4*(1-p_3toD),
  "p_3toD" = p_3toD,
  "p_4to1" = p_4to1*(1-p_4toD),
  "p_4to2" = p_4to2*(1-p_4toD),
  "p_4to3" = p_4to3*(1-p_4toD),
  "p_4to4" = p_4to4*(1-p_4toD),
  "p_4toD" = p_4toD,
  "rp_sv"  = rp_sv
)



## Cost and utility inputs
c_discount     <- 0.28
c_st           <- 112555
c_si           <- 238850
c_NYHA1_st     <- 2911  + (c_st * (1-c_discount))   # cost of remaining one cycle NYHA1 under stabilizers
c_NYHA2_st     <- 4129  + (c_st * (1-c_discount))   # cost of remaining one cycle NYHA2 under stabilizers
c_NYHA3_st     <- 6194  + (c_st * (1-c_discount))   # cost of remaining one cycle NYHA3 under stabilizers
c_NYHA4_st     <- 10208 + (c_st * (1-c_discount))   # cost of remaining one cycle NYHA4 under stabilizers
c_NYHA1_si     <- 2911  + (c_si * (1-c_discount))    # cost of remaining one cycle NYHA1 under silencers
c_NYHA2_si     <- 4129  + (c_si * (1-c_discount))    # cost of remaining one cycle NYHA2 under silencers
c_NYHA3_si     <- 6194  + (c_si * (1-c_discount))    # cost of remaining one cycle NYHA3 under silencers
c_NYHA4_si     <- 10208 + (c_si * (1-c_discount))    # cost of remaining one cycle NYHA4 under silencers
c_D       <- 0            # cost associated with being dead
c_age_cof <- 11.5         # cost age coefficient
v_cost_coeffs <- c(       # pack the cost regression coefficients in a vector
  "age" = c_age_cof
)


u_NYHA1   <- 0.82         # utility when NYHA1
u_NYHA2   <- 0.729        # utility when NYHA1
u_NYHA3   <- 0.633        # utility when NYHA1
u_NYHA4   <- 0.333        # utility when NYHA1
u_Death   <- 0            # utility when dead
u_age_cof <- -0.0018      # utility age coefficient
ru_sv     <- -0.0015      # change in utility of individuals with every additional year being NYHA 3 or 4

v_util_coeffs <- c(       # pack the utility regression coefficients in a vector
  "age" = u_age_cof
)
v_util_t_decs <- c(       # pack the state-specific utility decrements in a vector
  "NYHA3 or 4" = ru_sv
)

### Payoffs - stabilizers
v_states_costs_st <- c("NYHA1" = c_NYHA1_st, "NYHA2" = c_NYHA2_st, "NYHA3" = c_NYHA3_st, "NYHA4" = c_NYHA4_st,"Death" = c_D)
v_states_utilities_st <- c("NYHA1" = u_NYHA1, "NYHA2" = u_NYHA2, "NYHA3" = u_NYHA3, "NYHA4" = u_NYHA4, "Death" = u_Death)

### Payoffs - silencers
v_states_costs_si <- c("NYHA1" = c_NYHA1_si, "NYHA2" = c_NYHA2_si, "NYHA3" = c_NYHA3_si, "NYHA4" = c_NYHA4_si,"Death" = c_D)
v_states_utilities_si <- c("NYHA1" = u_NYHA1, "NYHA2" = u_NYHA2, "NYHA3" = u_NYHA3, "NYHA4" = u_NYHA4, "Death" = u_Death)


### Running the simulation ###

# Run the simulation:
## For being treated with stabilizers
res_st <- run_microSimV(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_names      = v_states_names,
  v_states_costs      = v_states_costs_st,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities_st,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)
## For being treated with silencers
res_si <- run_microSimV(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_names      = v_states_names,
  v_states_costs      = v_states_costs_si,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities_si,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
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
#------------------------------------------------------------------------------#

## @knitr micro_run_microSimV_diagram

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
