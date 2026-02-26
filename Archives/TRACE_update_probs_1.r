# ATTR-CM NYHA class-based microsimulation model - `update_probs1`

## Update Transition Probability function
### This function updates the transition probabilities at every cycle based on
### the health state occupied by individual 'i' at cycle 't'

# define the update_probs1 function:
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
    yes = return(v_probs),                       # return the transition probabiltiies           
    no = stop("Probabilities do not sum to 1")   # or produce an error    
  )
    
}
#------------------------------------------------------------------------------#

## @knitr micro_update_probs1_demo

# define function inputs:
v_states_names <- c("NYHA1","NYHA2", "NYHA3", "NYHA4", "Death")   # the model states
num_states <- length(v_states_names) # the number of health states
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
p_HtoD  <- 0.039962712            # all cause mortality, H = healthy
rr_HF   <- 3.17                   # HF vs healthy (HR)
p_1toD  <- rr_HF*p_HtoD
rr_2v1  <- 1.78                   # I vs II (HR)
rr_3v1  <- 3.51                   # I vs III (HR)
rr_4v1  <- 5.74                   # I vs IV (HR)
r_1toD  <- -log(1- p_1toD)        #rate of death in NYHA1
r_2toD  <-  rr_2v1*r_1toD         #rate of death in NYHA2
r_3toD  <-  rr_3v1*r_1toD         #rate of death in NYHA2
r_4toD  <-  rr_4v1*r_1toD         #rate of death in NYHA2
p_2toD  <-  1 - exp(-r_2toD)      #probability to die in NYHA2
p_3toD  <-  1 - exp(-r_3toD)      #probability to die in NYHA3
p_4toD  <-  1 - exp(-r_4toD)      #probability to die in NYHA4

m_trans_probs <- matrix(                      # create a transition probability matrix
  data = c(
    p_1to1, p_1to2, p_1to3, p_1to4, p_1toD,            # transition probabilities when NYHA1
    p_2to1, p_2to2, p_2to3, p_2to4, p_2toD,            # transition probabilities when NYHA2
    p_3to1, p_3to2, p_3to3, p_3to4, p_3toD,            # transition probabilities when NYHA3
    p_4to1, p_4to2, p_4to3, p_4to4, p_4toD,            # transition probabilities when NYHA4
    p_Dto1, p_Dto2, p_Dto3, p_Dto4, p_DtoD
  ),
  nrow = num_states,
  byrow = TRUE,
  dimnames = list(v_states_names, v_states_names)
)

# run the update_probs1 function:
update_probs1(
  occupied_state = "NYHA1",
  m_trans_probs = m_trans_probs
)
update_probs1(
  occupied_state = "NYHA2",
  m_trans_probs = m_trans_probs
)
update_probs1(
  occupied_state = "NYHA3",
  m_trans_probs = m_trans_probs
)
pdate_probs1(
  occupied_state = "NYHA4",
  m_trans_probs = m_trans_probs
)
pdate_probs1(
  occupied_state = "Death",
  m_trans_probs = m_trans_probs
)
