# Scenario analysis: health-state cost inputs and future medical costs (Silencer vs Stabilizer, base case and subgroups)
# Base case: Wang et al. 2023 annual all-cause costs (visits + pharmacy), no future medical costs because all-cause costs include unrelated care.
# Scenarios: the health-state costs in the ICER 2024 model (each annual visit cost divided by the annual visits), with future medical costs
# (Jiao & Basu 2021) and without them (as in the ICER 2024 model).

# Step 1. Run the base case model and load the SA helper functions
source("TRACE_microsim_V.1.5.R")
source("TRACE_SA_functions_V.1.5.R") # base_inputs, run_owsa_case()

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

### Step 2. Populations: same draws as run_subgroup() ###
make_pop <- function(age, sd, male, v_start_probs) {
  set.seed(seed + 1)
  m_male   <- cbind("age" = rnorm(n = num_i, mean = age, sd = sd))
  m_female <- cbind("age" = rnorm(n = num_i, mean = age, sd = sd))
  v_start  <- sample(x = v_states_names, size = num_i, replace = TRUE, prob = v_start_probs)
  list(v_start = v_start, m_male = m_male, m_female = m_female, prop_male = male)
}

### Step 3. Run each case ###
inp_icer <- base_inputs
inp_icer$c_states <- v_costs_icer
inp_future <- inp_icer
inp_future$future_cost_mult <- 1

run_population <- function(label, pop = NULL) {
  run_case <- function(case, inp) {
    res <- run_owsa_case(inp, pop)
    data.frame(population = label, case = case, t(res[c("Dcosts_st", "Dcosts_si", "inc_Dcosts", "inc_Dqalys", "ICER", "list_price_100k", "list_price_150k")]))
  }
  rbind(run_case("Base case: Wang annual costs, no future medical costs", base_inputs),
        run_case("ICER 2024 values + future medical costs", inp_future),
        run_case("ICER 2024 values, no future medical costs (as ICER 2024)", inp_icer))
}

cat("\n--- SCENARIO: Health-state cost inputs and future medical costs ---\n")
cost_results <- rbind(
  run_population("Base case"),
  run_population("Wild-type", make_pop(79, 8.1, 0.91, c(0.126, 0.536, 0.326, 0.012, 0))),
  run_population("Hereditary", make_pop(71, 8.15, 0.71, c(0.174, 0.44, 0.385, 0, 0)))
)
print(cost_results, digits = 6, row.names = FALSE)
stopifnot(abs(cost_results$ICER[1] - ICER_overall) < 1e-6 * ICER_overall) # same base case as the model
write.csv(cost_results, file.path(dir_out, "Scenario_costs_V.1.5.csv"), row.names = FALSE)
cat("\nSaved:", file.path(dir_out, "Scenario_costs_V.1.5.csv"), "\n")
