## @knitr micro_calc_costs1_foo

# ATTR-CM microsimulation model - `calc_costs1`

## Calculate Costs function
### This function estimates the costs at every cycle based on the health state
### occupied by individual 'i' at cycle 't'

# define the calc_costs1 function:
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

#------------------------------------------------------------------------------#

## @knitr micro_calc_costs1_demo

# define function inputs:
c_NYHA1        <- 2911                               # cost of remaining one cycle NYHA1
c_NYHA2        <- 4129                               # cost of remaining one cycle NYHA2
c_NYHA3        <- 6194                               # cost of remaining one cycle NYHA3
c_NYHA4        <- 10208                              # cost of remaining one cycle NYHA4
v_states_costs <- c("NYHA1" = c_NYHA1, "NYHA2" = c_NYHA2, "NYHA3" = c_NYHA3, "NYHA4" = c_NYHA4,"Death" = 0) # named costs vector

# run the calc_costs1 function:
calc_costs1(
  occupied_state = "NYHA1",
  v_states_costs = v_states_costs
)
calc_costs1(
  occupied_state = "NYHA2",
  v_states_costs = v_states_costs
)
calc_costs1(
  occupied_state = "NYHA3",
  v_states_costs = v_states_costs
)
calc_costs1(
  occupied_state = "NYHA4",
  v_states_costs = v_states_costs
)
#------------------------------------------------------------------------------#

## @knitr micro_calc_costs1_diagram

DiagrammeR::grViz("
  digraph flowchart {
    graph [ratio=1.0]
    node [fontname = 'Helvetica', fontsize=20, shape = box, style=filled, fillcolor='grey']
    edge [fontname = 'Helvetica']

    input [label = 'Inputs: occupied_state, v_states_costs', style=filled, fillcolor='yellow']
    initialize [label = 'Initialize state_costs', style=filled, fillcolor='palegreen']
    check_h [shape = diamond, style=filled, fillcolor='skyblue', label = 'occupied_state == \\\"H\\\"']
    update_h [label = 'state_costs <- v_states_costs[\\\"H\\\"]', style=filled, fillcolor='palegreen']
    check_s1 [shape = diamond, style=filled, fillcolor='skyblue', label = 'occupied_state == \\\"S1\\\"']
    update_s1 [label = 'state_costs <- v_states_costs[\\\"S1\\\"]', style=filled, fillcolor='palegreen']
    check_d [shape = diamond, style=filled, fillcolor='skyblue', label = 'occupied_state == \\\"D\\\"']
    update_d [label = 'state_costs <- v_states_costs[\\\"D\\\"]', style=filled, fillcolor='palegreen']
    return_state [shape = ellipse, label = 'Return state_costs', style=filled, fillcolor='yellow']

    input -> initialize
    initialize -> check_h
    check_h -> update_h [label = 'TRUE']
    check_h -> check_s1 [label = 'FALSE']
    update_h -> check_s1
    check_s1 -> update_s1 [label = 'TRUE']
    check_s1 -> check_d [label = 'FALSE']
    update_s1 -> check_d
    check_d -> update_d [label = 'TRUE']
    check_d -> return_state [label = 'FALSE']
    update_d -> return_state
  }
")
