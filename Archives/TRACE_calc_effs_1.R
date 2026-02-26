## @knitr micro_calc_effs1_foo

# ATTR-CM microsimulation model - `calc_effs1`

## Calculate Health Outcomes function
### This function estimates the Quality Adjusted Life Years (QALYs) at every
### cycle based on the health state occupied by individual 'i' at cycle 't' and
### the cycle_length (measured in years)

# define the calc_effs1 function:
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

## @knitr micro_calc_effs1_demo

# define function inputs:
u_NYHA1   <- 0.82         # utility when NYHA1
u_NYHA2   <- 0.729        # utility when NYHA1
u_NYHA3   <- 0.633        # utility when NYHA1
u_NYHA4   <- 0.333        # utility when NYHA1
u_Death   <- 0            # utility when dead
v_states_utilities <- c("NYHA1" = u_NYHA1, "NYHA2" = u_NYHA2, "NYHA3" = u_NYHA3, "NYHA4" = u_NYHA4, "Death" = u_Death) # named utilities vector

# run the calc_effs1 function:
calc_effs1(
  occupied_state = "NYHA1",
  v_states_utilities = v_states_utilities,
  cycle_length = 0.5                         # 1 year
)
calc_effs1(
  occupied_state = "NYHA2",
  v_states_utilities = v_states_utilities,
  cycle_length = 0.5
)
calc_effs1(
  occupied_state = "NYHA3",
  v_states_utilities = v_states_utilities,
  cycle_length = 0.5                       # 6 months
)
calc_effs1(
  occupied_state = "NYHA4",
  v_states_utilities = v_states_utilities,
  cycle_length = 0.5                       # 6 months
)

#------------------------------------------------------------------------------#

## @knitr micro_calc_effs1_diagram

DiagrammeR::grViz("
  digraph flowchart {
    graph [ratio=1.0]
    node [fontname = 'Helvetica', fontsize=20, shape = box, style=filled, fillcolor='grey']
    edge [fontname = 'Helvetica']

    input [label = 'Inputs: occupied_state, v_states_utilities, cycle_length', style=filled, fillcolor='yellow']
    initialize [label = 'Initialize state_utility', style=filled, fillcolor='palegreen']
    check_h [shape = diamond, style=filled, fillcolor='skyblue', label = 'occupied_state == \\\"H\\\"']
    update_h [label = 'state_utility <- v_states_utilities[\\\"H\\\"]', style=filled, fillcolor='palegreen']
    check_s1 [shape = diamond, style=filled, fillcolor='skyblue', label = 'occupied_state == \\\"S1\\\"']
    update_s1 [label = 'state_utility <- v_states_utilities[\\\"S1\\\"]', style=filled, fillcolor='palegreen']
    check_d [shape = diamond, style=filled, fillcolor='skyblue', label = 'occupied_state == \\\"D\\\"']
    update_d [label = 'state_utility <- v_states_utilities[\\\"D\\\"]', style=filled, fillcolor='palegreen']
    calculate_qalys [label = 'QALYs <- state_utility x cycle_length', style=filled, fillcolor='palegreen']
    return_qalys [shape = ellipse, label = 'Return QALYs', style=filled, fillcolor='yellow']

    input -> initialize
    initialize -> check_h
    check_h -> update_h [label = 'TRUE']
    check_h -> check_s1 [label = 'FALSE']
    update_h -> check_s1
    check_s1 -> update_s1 [label = 'TRUE']
    check_s1 -> check_d [label = 'FALSE']
    update_s1 -> check_d
    check_d -> update_d [label = 'TRUE']
    check_d -> calculate_qalys [label = 'FALSE']
    update_d -> calculate_qalys
    calculate_qalys -> return_qalys
  }
")
