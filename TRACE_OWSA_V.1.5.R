# One-way sensitivity analysis (OWSA): Silencer vs Stabilizer, base case population
# Change one input at a time to its low and high value, keep everything else at base case
# Theta stays fixed, except for the survival target row (theta is re-calibrated there)

# Step 1. Run the base case model (functions, inputs, calibrated theta)
source("TRACE_microsim_V.1.5.R")
library(ggplot2)

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

source("TRACE_SA_functions_V.1.5.R") # base_inputs, run_owsa_case(), calibrate_theta()

### Step 2. Input ranges ###
# Each row: label, input name, low value, high value
owsa_params <- list(
  # Discontinuation per cycle (6 months)
  list("Silencer discontinuation per cycle (0.51% HELIOS-B; 0.25-6.5%)",     "p_disc_cycle_si", p_disc_si / 2, 1 - (1 - 0.126)^0.5), # high: Komodo 12-month, net of restarts
  list("Stabilizer discontinuation per cycle (1.93% ATTRibute-CM; 0.96-4.65%)", "p_disc_cycle_st", p_disc_st / 2, 1 - (1 - 0.212)^(1/5)), # high: ATTR-ACT
  list("Silencer waning period (1 yr; 0-2 yr)",                 "waning_years",   0,      2),
  # Drug prices: one input per arm
  # Stabilizer WAC is a list price, so only the WAC discount varies
  # Vutrisiran ASP is already close to net, so only the payment level varies
  list("Silencer annual price ($477K ASP; ASP to ASP+6%)",       "c_si_annual",    c_si_annual,       c_si_annual * 1.06), # high: ASP + 6% = Medicare Part B payment limit
  list("Stabilizer WAC discount (27.5%; 20-40%)",               "c_discount_st",  0.20,   0.40), # low: 20% = Part D Manufacturer Discount Program, catastrophic phase
  list("Discount rate, costs & QALYs (3%; 0-5%)",               "discount_rate",  0,      0.05),
  list("Utility NYHA I (0.820; ±10%)",                          "u_NYHA1",        u_NYHA1 * 0.9, u_NYHA1 * 1.1),
  list("Utility NYHA II (0.729; ±10%)",                         "u_NYHA2",        u_NYHA2 * 0.9, u_NYHA2 * 1.1),
  list("Utility NYHA III (0.633; ±10%)",                        "u_NYHA3",        u_NYHA3 * 0.9, u_NYHA3 * 1.1),
  list("Utility NYHA IV (0.333; ±10%)",                         "u_NYHA4",        u_NYHA4 * 0.9, u_NYHA4 * 1.1),
  list("Health state costs (±25%)",                             "state_cost_mult",   0.75, 1.25),
  list("Hospitalization probabilities (±25%)",                  "hosp_prob_mult",    0.75, 1.25),
  list("Hospitalization costs (±25%)",                          "hosp_cost_mult",    0.75, 1.25),
  list("Hospitalization disutilities (±25%)",                   "hosp_disutil_mult", 0.75, 1.25),
  list("42-month survival target, THAOS 95% CI (0.768; 0.709-0.817)", "surv_target", 0.709, 0.817)
)

### Step 3. Run OWSA ###
cat("\n--- Running OWSA ---\n")
base_owsa <- run_owsa_case(base_inputs)
cat("Base case ICER:", base_owsa["ICER"], "\n")

owsa_rows <- lapply(owsa_params, function(p) {
  res <- lapply(c(low = p[[3]], high = p[[4]]), function(val) {
    inp <- base_inputs
    if (p[[2]] == "surv_target") { # survival target: re-calibrate theta
      inp$theta <- calibrate_theta(val)
      cat(sprintf("  Survival target %.3f -> theta = %.4f\n", val, inp$theta))
    } else inp[[p[[2]]]] <- val # all other inputs: just replace the value
    run_owsa_case(inp)
  })
  cat(sprintf("  %-62s low: %12.0f  high: %12.0f\n", p[[1]], res$low["ICER"], res$high["ICER"]))
  data.frame(parameter = p[[1]], low_value = p[[3]], high_value = p[[4]],
             ICER_low = res$low["ICER"], ICER_high = res$high["ICER"],
             inc_Dqalys_low = res$low["inc_Dqalys"], inc_Dqalys_high = res$high["inc_Dqalys"],
             inc_Dcosts_low = res$low["inc_Dcosts"], inc_Dcosts_high = res$high["inc_Dcosts"],
             list_price_150k_low = res$low["list_price_150k"], list_price_150k_high = res$high["list_price_150k"],
             row.names = NULL)
})
owsa_results <- do.call(rbind, owsa_rows)
# Sort by ICER range (widest bar goes on top of the tornado)
owsa_results$range <- abs(owsa_results$ICER_high - owsa_results$ICER_low)
owsa_results <- owsa_results[order(-owsa_results$range), ]

cat("\n--- OWSA RESULTS (sorted by ICER range) ---\n")
print(owsa_results[, c("parameter", "ICER_low", "ICER_high", "range")], row.names = FALSE)
write.csv(owsa_results, file.path(dir_out, "OWSA_results_V.1.5.csv"), row.names = FALSE)

### Step 4. Scenario: silencer waning period 0-3 years ###
cat("\n--- SCENARIO: Silencer waning period ---\n")
waning_scenarios <- do.call(rbind, lapply(0:3, function(w) {
  inp <- base_inputs; inp$waning_years <- w
  r <- run_owsa_case(inp)
  data.frame(waning_years = w, inc_Dcosts = r["inc_Dcosts"], inc_Dqalys = r["inc_Dqalys"], ICER = r["ICER"], row.names = NULL)
}))
print(waning_scenarios, row.names = FALSE)
write.csv(waning_scenarios, file.path(dir_out, "OWSA_waning_scenarios_V.1.5.csv"), row.names = FALSE)

### Step 5. Tornado diagram ###
source("TRACE_OWSA_tornado_V.1.5.R") # saves OWSA_tornado_V.1.5.png
cat("\nSaved:", file.path(dir_out, c("OWSA_results_V.1.5.csv", "OWSA_waning_scenarios_V.1.5.csv", "OWSA_tornado_V.1.5.png")), sep = "\n  ")
