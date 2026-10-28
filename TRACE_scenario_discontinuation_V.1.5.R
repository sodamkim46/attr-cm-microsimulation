# Scenario analysis: discontinuation rates (Silencer vs Stabilizer, base case population)
# Base case: AE-related discontinuation from the pivotal trials
#   Stabilizer: ATTRibute-CM, 9.3% over 30 months
#   Silencer:   HELIOS-B, 3% over 36 months
# Scenarios: ATTR-ACT and real-world persistence data
# Everything else stays at base case (incl. theta and 1-year silencer waning)

# Step 1. Run the base case model and load the SA helper functions
source("TRACE_microsim_V.1.5.R")
source("TRACE_SA_functions_V.1.5.R") # base_inputs, run_owsa_case()

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

### Step 2. Discontinuation probabilities per cycle (6 months) ###
# 12-month probability -> per-cycle probability
p12_to_cycle <- function(p12) 1 - (1 - p12)^cycle_length

p_cycle_st_trial <- 1 - (1 - p_disc_2.5y_st)^(1/5)
p_cycle_si_trial <- 1 - (1 - p_disc_3y_si)^(1/6)
p_cycle_st_attract <- 1 - (1 - 0.212)^(1/5) # ATTR-ACT tafamidis pooled, 21.2% over 30 months

# Each row: label, stabilizer per-cycle prob, silencer per-cycle prob
disc_scenarios <- list(
  list("Base case (trial-based: ST ATTRibute-CM 9.3%/30mo, SI HELIOS-B 3%/36mo)", p_cycle_st_trial,    p_cycle_si_trial),
  list("ST ATTR-ACT (21.2%/30mo), SI HELIOS-B",                                  p_cycle_st_attract,  p_cycle_si_trial),
  list("Real-world US: ST 7%/yr (Medicare+PAP), SI 12.6%/yr (Komodo, net of restarts)", p12_to_cycle(0.07), p12_to_cycle(0.126)),
  list("Real-world US: ST 7%/yr (Medicare+PAP), SI 12%/yr (Optum)",              p12_to_cycle(0.07),  p12_to_cycle(0.12)),
  list("Real-world upper: ST 22%/yr (German claims, incl. deaths), SI 20%/yr (Komodo, incl. restarts)", p12_to_cycle(0.22), p12_to_cycle(0.20)),
  # Change one arm at a time, to see which arm drives the result
  list("Real-world ST only (7%/yr), SI trial-based",                             p12_to_cycle(0.07),  p_cycle_si_trial),
  list("Real-world SI only (12.6%/yr), ST trial-based",                          p_cycle_st_trial,    p12_to_cycle(0.126))
)

### Step 3. Run each scenario ###
cat("\n--- SCENARIO: Real-world discontinuation ---\n")
disc_results <- do.call(rbind, lapply(disc_scenarios, function(sc) {
  inp <- base_inputs
  inp$p_disc_cycle_st <- sc[[2]]
  inp$p_disc_cycle_si <- sc[[3]]
  r <- run_owsa_case(inp)
  data.frame(scenario         = sc[[1]],
             annual_disc_st   = 1 - (1 - sc[[2]])^(1 / cycle_length), # back to annual, for the table
             annual_disc_si   = 1 - (1 - sc[[3]])^(1 / cycle_length),
             inc_Dcosts       = r["inc_Dcosts"],
             inc_Dqalys       = r["inc_Dqalys"],
             ICER             = r["ICER"],
             list_price_100k  = r["list_price_100k"],
             list_price_150k  = r["list_price_150k"],
             row.names = NULL)
}))
print(disc_results, digits = 6, row.names = FALSE)
write.csv(disc_results, file.path(dir_out, "Scenario_discontinuation_V.1.5.csv"), row.names = FALSE)
cat("\nSaved:", file.path(dir_out, "Scenario_discontinuation_V.1.5.csv"), "\n")
