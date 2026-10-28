# Scenario analysis: hospitalization rates from a more recent trial (Silencer vs Stabilizer, base case population)
# Base case: CV hospitalization probabilities by NYHA class from ATTR-ACT (ICER 2024 Table E2.3; HAS 2020),
# same in both arms
# ATTRibute-CM (enrolled 2019+) had lower CV hospitalization rates than ATTR-ACT (2013-2018)
# No trial reports rates by NYHA class, so all NYHA-specific rates are scaled by the ratio of overall annual rates (k):
#   p' = 1 - (1 - p)^k
# Everything else stays at base case (incl. theta)

# Step 1. Run the base case model and load the SA helper functions
source("TRACE_microsim_V.1.5.R")
source("TRACE_SA_functions_V.1.5.R") # base_inputs, run_owsa_case()

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

### Step 2. Rate ratio k ###
# Annual CV hospitalization rates (per patient-year)
cvh_attract_pbo <- 0.70  # ATTR-ACT placebo (Maurer et al., NEJM 2018)
cvh_attract_taf <- 0.48  # ATTR-ACT tafamidis pooled
cvh_attrib_pbo  <- 0.450 # ATTRibute-CM placebo (EMA EPAR, imputed)
cvh_attrib_aco  <- 0.224 # ATTRibute-CM acoramidis

k_placebo <- cvh_attrib_pbo / cvh_attract_pbo # option 1: placebo arms only
k_pooled  <- ((409 * cvh_attrib_aco + 202 * cvh_attrib_pbo) / 611) /   # option 2: both arms, weighted by n
             ((264 * cvh_attract_taf + 177 * cvh_attract_pbo) / 441)    # ATTRibute-CM mITT 409 vs 202; ATTR-ACT 264 vs 177

# Each row: label, rate multiplier k
hosp_scenarios <- list(
  list("Base case (ATTR-ACT, ICER/HAS)",                                1),
  list(sprintf("Scaled to ATTRibute-CM: placebo ratio (k = %.3f)", k_placebo), k_placebo),
  list(sprintf("Scaled to ATTRibute-CM: pooled-arm ratio (k = %.3f)", k_pooled), k_pooled)
)

### Step 3. Run each scenario ###
cat("\n--- SCENARIO: Hospitalization rates scaled to ATTRibute-CM ---\n")
hosp_results <- do.call(rbind, lapply(hosp_scenarios, function(sc) {
  inp <- base_inputs
  inp$hosp_rate_mult <- sc[[2]]
  r <- run_owsa_case(inp)
  data.frame(scenario        = sc[[1]],
             rate_mult       = sc[[2]],
             inc_Dcosts      = r["inc_Dcosts"],
             inc_Dqalys      = r["inc_Dqalys"],
             ICER            = r["ICER"],
             list_price_100k = r["list_price_100k"],
             list_price_150k = r["list_price_150k"],
             row.names = NULL)
}))
print(hosp_results, digits = 6, row.names = FALSE)
write.csv(hosp_results, file.path(dir_out, "Scenario_hospitalization_V.1.5.csv"), row.names = FALSE)
cat("\nSaved:", file.path(dir_out, "Scenario_hospitalization_V.1.5.csv"), "\n")
