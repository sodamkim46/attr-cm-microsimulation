# Healthy-Sick-Dead microsimulation model demo

# clear R's session memory (Global Environment)
rm(list = ls())



### Defining model functions ###

# Define model functions

## Update Transition Probability function
### This function updates the transition probabilities at every cycle based on
### the health state occupied by individual 'i' at cycle 't'
update_probs1 <- function(
    occupied_state, m_trans_probs) {
  
  v_probs <- rep(NA, num_states)                                 # create vector of state transition probabilities
  
  # update v_probs with the appropriate probabilities
  v_probs[occupied_state == "NYHA1"] <- m_trans_probs["NYHA1",]
  v_probs[occupied_state == "NYHA2"] <- m_trans_probs["NYHA2",]
  v_probs[occupied_state == "NYHA3"] <- m_trans_probs["NYHA3",]
  v_probs[occupied_state == "NYHA4"] <- m_trans_probs["NYHA4",]
  v_probs[occupied_state == "Death"] <- m_trans_probs["Death",]
  
  # sanity check
  ifelse(
    test = abs(sum(v_probs) - 1) < 1e-12,        # check if the transition probabilities add up to 1
    yes = return(v_probs),                        # return the transition probabilities           
    no = stop("Probabilities do not sum to 1")   # or produce an error    
  )
  
}

## Calculate Costs function
### This function estimates the costs at every cycle based on the health state
### occupied by individual 'i' at cycle 't'
calc_costs1 <- function (
    occupied_state,
    v_states_costs) {
  
  state_costs <- NA # If we see NA, that means we didn't find occupied state
  state_costs[occupied_state == "NYHA1"]  <- v_states_costs["NYHA1"]  # update the cost if NYHA1
  state_costs[occupied_state == "NYHA2"] <- v_states_costs["NYHA2"] # update the cost if NYHA2
  state_costs[occupied_state == "NYHA3"]  <- v_states_costs["NYHA3"]  # update the cost if NYHA3
  state_costs[occupied_state == "NYHA4"]  <- v_states_costs["NYHA4"]  # update the cost if NYHA4
  state_costs[occupied_state == "Death"]  <- v_states_costs["Death"]  # update the cost if NYHA3
  
  return(state_costs)                                         # return the costs
}

## Calculate Health Outcomes function
### This function estimates the Quality Adjusted Life Years (QALYs) at every
### cycle based on the health state occupied by individual 'i' at cycle 't' and
### the cycle_length (measured in years)
calc_effs1 <- function (
    occupied_state,
    v_states_utilities,
    cycle_length = 1) {
  
  state_utility <-NA
  state_utility[occupied_state == "NYHA1"]  <- v_states_utilities["NYHA1"]  # update the utility if NYHA1
  state_utility[occupied_state == "NYHA2"]  <- v_states_utilities["NYHA2"]  # update the utility if NYHA2
  state_utility[occupied_state == "NYHA3"]  <- v_states_utilities["NYHA3"]  # update the utility if NYHA3
  state_utility[occupied_state == "NYHA4"]  <- v_states_utilities["NYHA4"]  # update the utility if NYHA4
  state_utility[occupied_state == "Death"]  <- v_states_utilities["Death"]  # update the utility if dead
  
  QALYs <-  state_utility * cycle_length                            # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                     # return the QALYs
}

#------------------------------------------------------------------------------#

## @knitr micro_run_microSim1_foo

# ATTR-CM microsimulation model - `run_microSim1`

## Run Microsimulation function
### This function runs the microsimulation function of the NYHA1-4
### model
run_microSim1 <- function(
    v_starting_states,
    num_i,
    num_cycles,
    v_states_names,
    v_states_costs,
    v_states_utilities,
    m_trans_probs,
    cycle_length = 0.5,
    starting_seed = 1) {

  # create matrices to capture states' names, associated costs and QALYs
  m_States <- m_Costs <- m_Effs <-  matrix(nrow = num_i, ncol = num_cycles + 1)

  for (i in 1:num_i) {                     # for each 'i' of the 'num_i' simulated individual:
    set.seed(starting_seed + i)            # set the seed for every individual for the random number generator

    # Step 1:
    m_States[i, 1] <- v_starting_states[i] # indicate the initial health state
    m_Costs[i, 1]  <- calc_costs1(         # calculate the costs incurred in their starting health state
      occupied_state = m_States[i, 1],
      v_states_costs = v_states_costs
    )
    m_Effs[i, 1]   <- calc_effs1(          # calculate the QALYs accrued in their starting health state
      occupied_state = m_States[i, 1],
      v_states_utilities = v_states_utilities,
      cycle_length = cycle_length
    )

    for (t in 1:num_cycles) {              # for each 't' of the 'num_cycles' cycles:
      # Step 2:
      v_trans_probs <- update_probs1(      # update the transition probabilities at cycle 't'
        occupied_state = m_States[i, t],
        m_trans_probs = m_trans_probs
      )

      # Step 3:
      m_States[i, t + 1] <- sample(        # sample the health state at 't + 1'
        x = v_states_names,
        prob = v_trans_probs,
        size = 1
      )

      # Step 4:
      m_Costs[i, t + 1]  <- calc_costs1(    # calculate the costs incurred in their 't + 1' health state
        occupied_state = m_States[i, t + 1],
        v_states_costs = v_states_costs
      )
      m_Effs[i, t + 1]   <- calc_effs1(     # calculate the QALYs accrued in their 't + 1' health state
        occupied_state = m_States[i, t + 1],
        v_states_utilities = v_states_utilities,
        cycle_length = cycle_length
      )

    } # close the loop for the cycles 't'

    if(i/100 == round(i/100,0)) {          # display the progress of the simulation
      cat('\r', paste(i/num_i * 100, "% done", sep = " "))
    }

  } # close the loop for the individuals 'i'

  # Step 5:
  v_total_costs <- rowSums(m_Costs)           # calculate total cost per individual
  v_total_qalys <- rowSums(m_Effs)            # calculate total QALYs per individual
  mean_costs    <- mean(v_total_costs)        # calculate average cost
  mean_qalys    <- mean(v_total_qalys)        # calculate average QALYs

  # store the results in a list:
  results <- list(
    m_States = m_States,
    m_Costs = m_Costs,
    m_Effs = m_Effs,
    v_total_costs = v_total_costs,
    v_total_qalys = v_total_qalys,
    mean_costs = mean_costs,
    mean_qalys = mean_qalys
  )

  # return the results
  return(results)
}

#------------------------------------------------------------------------------#

## @knitr micro_run_microSim1_demo

### Defining model parameters ###

# Define model inputs
## General parameters
num_i      <- 1e5                     # number of simulated individuals
num_cycles <- 50                      # time horizon

## Health states
v_states_names <- c("NYHA1","NYHA2", "NYHA3", "NYHA4", "Death")   # the model states
num_states <- length(v_states_names) # the number of health states
p_starting <- c(0.108, 0.72, 0.172, 0, 0)
v_starting_states <- sample(v_states_names, size = num_i, replace = TRUE, prob = p_starting)   # individuals start in the state they are assigned by starting probability

## Transition probabilities (per cycle)
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

m_trans_probs <- matrix(                      # create a transition probability matrix
  data = c(
    p_1to1*(1-p_1toD), p_1to2*(1-p_1toD), p_1to3*(1-p_1toD), p_1to4*(1-p_1toD), p_1toD,            # transition probabilities when NYHA1
    p_2to1*(1-p_2toD), p_2to2*(1-p_2toD), p_2to3*(1-p_2toD), p_2to4*(1-p_2toD), p_2toD,            # transition probabilities when NYHA2
    p_3to1*(1-p_3toD), p_3to2*(1-p_3toD), p_3to3*(1-p_3toD), p_3to4*(1-p_3toD), p_3toD,            # transition probabilities when NYHA3
    p_4to1*(1-p_4toD), p_4to2*(1-p_4toD), p_4to3*(1-p_4toD), p_4to4*(1-p_4toD), p_4toD,            # transition probabilities when NYHA4
    p_Dto1, p_Dto2, p_Dto3, p_Dto4, p_DtoD
  ),
  nrow = num_states,
  byrow = TRUE,
  dimnames = list(v_states_names, v_states_names)
)

## Cost inputs
c_NYHA1        <- 2911                               # cost of remaining one cycle NYHA1
c_NYHA2        <- 4129                               # cost of remaining one cycle NYHA2
c_NYHA3        <- 6194                               # cost of remaining one cycle NYHA3
c_NYHA4        <- 10208                              # cost of remaining one cycle NYHA4
v_states_costs <- c("NYHA1" = c_NYHA1, "NYHA2" = c_NYHA2, "NYHA3" = c_NYHA3, "NYHA4" = c_NYHA4,"Death" = 0) # named costs vector

## Utility inputs
u_NYHA1   <- 0.82         # utility when NYHA1
u_NYHA2   <- 0.729        # utility when NYHA1
u_NYHA3   <- 0.633        # utility when NYHA1
u_NYHA4   <- 0.333        # utility when NYHA1
u_Death   <- 0            # utility when dead
v_states_utilities <- c("NYHA1" = u_NYHA1, "NYHA2" = u_NYHA2, "NYHA3" = u_NYHA3, "NYHA4" = u_NYHA4, "Death" = u_Death) # named utilities vector


### Running the simulation ###

# Run the simulation:
microsim_results <- run_microSim1(
  v_starting_states = v_starting_states,
  num_i = num_i,
  num_cycles = num_cycles,
  v_states_names = v_states_names,
  v_states_costs = v_states_costs,
  v_states_utilities = v_states_utilities,
  m_trans_probs = m_trans_probs,
  cycle_length = 0.5,
  starting_seed = 1
)

# View the results:
str(microsim_results)

microsim_results$v_total_costs[1:10]
microsim_results$mean_costs

microsim_results$v_total_qalys[1:10]
microsim_results$mean_qalys

#------------------------------------------------------------------------------#
