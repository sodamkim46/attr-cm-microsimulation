# Probabilistic sensitivity analysis (PSA): Silencer vs Stabilizer, base case population
# Distributions follow the ICER 2024 ATTR-CM model (Inputs sheet of the ICER Excel model), except drug prices
#   Utilities                 Beta, SE 0.02 / 0.01 / 0.01 / 0.06 (NYHA I-IV)
#   Hospital disutilities     Beta, SE = 10% of the mean
#   State and hospital costs  Gamma, ICER SEs (CV 16.3% / 11.0% / 13.2% / 29.8% for NYHA I-IV)
#   Transitions               Normal(1, 0.1) multiplier on moves between NYHA classes (staying = 1 - sum),
#                             one for the on-treatment matrix and one for the matrix after discontinuation
#   Stabilizer WAC discount   PERT(20%, 27.5%, 40%)
#   Silencer price            Uniform(ASP, ASP + 6%)
# Fixed at base case (same as ICER): mortality HRs, theta, hospitalization probabilities, discontinuation,
# population, discount rate

# Step 1. Run the base case model and load the SA helper functions
source("TRACE_microsim_V.1.5.R")
source("TRACE_SA_functions_V.1.5.R") # base_inputs, run_owsa_case()
library(ggplot2)
library(parallel)

dir_out <- "OWSA_output"
dir.create(dir_out, showWarnings = FALSE)

### Step 2. PSA settings ### - 40 M
n_sim     <- 1000  # number of PSA iterations (ICER 2024 also used 1,000) - parameter uncertainty - extracting parameters
n_psa_i   <- 20000 # patients per sex-specific run: first 20,000 of the base case population, to save time  
n_cores   <- 4
psa_seed  <- 20261020
v_wtp_psa <- c(100000, 150000, 200000)

# Same population and same model seed in every iteration,
# so the iterations differ only by the sampled parameters
pop_psa <- list(v_start  = v_starting_states[1:n_psa_i],
                m_male   = m_indi_features_male[1:n_psa_i, , drop = FALSE],
                m_female = m_indi_features_female[1:n_psa_i, , drop = FALSE])

### Step 3. Distribution functions ###
beta_mom <- function(n, mu, se) { # Beta from mean and SE (method of moments)
  a <- mu * (mu * (1 - mu) / se^2 - 1)
  rbeta(n, a, a * (1 - mu) / mu)
}
gamma_cv <- function(n, mu, cv) rgamma(n, shape = 1 / cv^2, scale = mu * cv^2) # Gamma from mean and CV (coefficient of variation) = SE/mean
rpert <- function(n, min, mode, max) { # Beta-PERT from min, mode, max
  a <- 1 + 4 * (mode - min) / (max - min)
  b <- 1 + 4 * (max - mode) / (max - min)
  min + (max - min) * rbeta(n, a, b)
}

# ICER 2024 state cost means and SEs (Q1 2024 USD)
se_state_icer <- c(473.77027, 453.26261, 820.38952, 3043.51831)
mu_state_icer <- c(2910.89329, 4129.41891, 6193.80504, 10208.27074)
cv_cost       <- se_state_icer / mu_state_icer # CV = SE / mean; ICER uses the same CVs for hospitalization costs

### Step 4. Sample parameters (one row per iteration) ###
set.seed(psa_seed)
c_states_base <- c(c_NYHA1_base, c_NYHA2_base, c_NYHA3_base, c_NYHA4_base)
c_hosp_base   <- v_hosp_costs[1:4]
u_hosp_base   <- v_hosp_disutils[1:4]

psa_params <- data.frame(
  c_discount_st  = rpert(n_sim, 0.20, c_discount_st, 0.40),
  c_si_annual    = runif(n_sim, c_si_annual, c_si_annual * 1.06),
  u_NYHA1        = beta_mom(n_sim, u_NYHA1, 0.02),
  u_NYHA2        = beta_mom(n_sim, u_NYHA2, 0.01),
  u_NYHA3        = beta_mom(n_sim, u_NYHA3, 0.01),
  u_NYHA4        = beta_mom(n_sim, u_NYHA4, 0.06),
  u_hosp1        = beta_mom(n_sim, u_hosp_base[1], u_hosp_base[1] / 10),
  u_hosp2        = beta_mom(n_sim, u_hosp_base[2], u_hosp_base[2] / 10),
  u_hosp3        = beta_mom(n_sim, u_hosp_base[3], u_hosp_base[3] / 10),
  u_hosp4        = beta_mom(n_sim, u_hosp_base[4], u_hosp_base[4] / 10),
  c_state1       = gamma_cv(n_sim, c_states_base[1], cv_cost[1]),
  c_state2       = gamma_cv(n_sim, c_states_base[2], cv_cost[2]),
  c_state3       = gamma_cv(n_sim, c_states_base[3], cv_cost[3]),
  c_state4       = gamma_cv(n_sim, c_states_base[4], cv_cost[4]),
  c_hosp1        = gamma_cv(n_sim, c_hosp_base[1], cv_cost[1]),
  c_hosp2        = gamma_cv(n_sim, c_hosp_base[2], cv_cost[2]),
  c_hosp3        = gamma_cv(n_sim, c_hosp_base[3], cv_cost[3]),
  c_hosp4        = gamma_cv(n_sim, c_hosp_base[4], cv_cost[4]),
  trans_mult_tx  = rnorm(n_sim, 1, 0.1), 
  # as transition counts were not available, same approach as the ICER model. X capture the uncertainty in each individual transition or the correlation
  trans_mult_soc = rnorm(n_sim, 1, 0.1)
)

# Put one row of sampled values into the input list
psa_inputs <- function(par) {
  inp <- base_inputs
  inp$c_discount_st  <- par$c_discount_st
  inp$c_si_annual    <- par$c_si_annual
  inp$u_NYHA1 <- par$u_NYHA1; inp$u_NYHA2 <- par$u_NYHA2; inp$u_NYHA3 <- par$u_NYHA3; inp$u_NYHA4 <- par$u_NYHA4
  inp$u_hosp         <- c(par$u_hosp1, par$u_hosp2, par$u_hosp3, par$u_hosp4)
  inp$c_states       <- c(par$c_state1, par$c_state2, par$c_state3, par$c_state4)
  inp$c_hosp         <- c(par$c_hosp1, par$c_hosp2, par$c_hosp3, par$c_hosp4)
  inp$trans_mult_tx  <- par$trans_mult_tx # maintain the equal efficacy assumption
  inp$trans_mult_soc <- par$trans_mult_soc
  inp
}

### Step 5. Check: base case on the PSA population ###
# Compare with the full population to see the Monte Carlo error from using 20,000 patients
cat("\n--- PSA: base case on the PSA population (n =", n_psa_i, "per sex) ---\n")
det_psa_pop <- run_owsa_case(base_inputs, pop = pop_psa) # one-row for base case
print(round(det_psa_pop, 4))

### Step 6. Run PSA (parallel) ###
cat("\n--- Running PSA:", n_sim, "iterations on", n_cores, "cores ---\n")
t_start <- Sys.time()
psa_list <- mclapply(seq_len(n_sim), function(k) {
  run_owsa_case(psa_inputs(psa_params[k, ]), pop = pop_psa)
}, mc.cores = n_cores)
cat("PSA run time:", format(Sys.time() - t_start), "\n")

# Stop if any iteration failed
failed <- which(!vapply(psa_list, is.numeric, logical(1)))
if (length(failed) > 0) stop("PSA iterations failed: ", paste(failed, collapse = ", "))
psa_results <- cbind(iteration = seq_len(n_sim), psa_params, do.call(rbind, psa_list)) # 1000 rows
write.csv(psa_results, file.path(dir_out, "PSA_results_V.1.5.csv"), row.names = FALSE)

### Step 7. Summary: mean and 95% interval ###
q95 <- function(x) quantile(x, c(0.025, 0.975))
summ_row <- function(label, x) data.frame(outcome = label, mean = mean(x), lower_95 = q95(x)[1], upper_95 = q95(x)[2], row.names = NULL)
psa_summary <- rbind(
  summ_row("Discounted costs, stabilizer", psa_results$Dcosts_st),
  summ_row("Discounted costs, silencer",   psa_results$Dcosts_si),
  summ_row("Incremental costs",            psa_results$inc_Dcosts),
  summ_row("Discounted QALYs, stabilizer", psa_results$Dqalys_st),
  summ_row("Discounted QALYs, silencer",   psa_results$Dqalys_si),
  summ_row("Incremental QALYs",            psa_results$inc_Dqalys),
  summ_row("ICER (per iteration)",         psa_results$ICER), # just for reference, DO NOT USE this VALUE to interpret the results
  summ_row("Threshold annual price, $100K/QALY", psa_results$list_price_100k),
  summ_row("Threshold annual price, $150K/QALY", psa_results$list_price_150k)
)
icer_psa <- mean(psa_results$inc_Dcosts) / mean(psa_results$inc_Dqalys) # ICER = mean incremental cost / mean incremental QALYs (ratio of means not  mean of ind ICERs)
prob_ce  <- sapply(v_wtp_psa, function(w) mean(w * psa_results$inc_Dqalys - psa_results$inc_Dcosts > 0)) # cost-effective if NMB > 0
# Net monertary benefit (NMB) = WTP * delta (QALY) - delta (cost)

cat("\n--- PSA RESULTS (Silencer vs Stabilizer) ---\n")
print(psa_summary, digits = 6, row.names = FALSE)
cat("ICER (ratio of mean incremental costs to mean incremental QALYs):", round(icer_psa), "\n")
cat("Deterministic ICER on the PSA population:", round(det_psa_pop["ICER"]), "| full population:", round(ICER_overall), "\n")
cat("Share of iterations with incremental QALYs > 0:", mean(psa_results$inc_Dqalys > 0), "\n")
cat("Probability cost-effective at", paste0("$", v_wtp_psa / 1000, "K"), ":", prob_ce, "\n")

psa_summary <- rbind(psa_summary,
  data.frame(outcome = "ICER (ratio of means)", mean = icer_psa, lower_95 = NA, upper_95 = NA),
  data.frame(outcome = paste0("Probability cost-effective at $", v_wtp_psa / 1000, "K/QALY"), mean = prob_ce, lower_95 = NA, upper_95 = NA))
write.csv(psa_summary, file.path(dir_out, "PSA_summary_V.1.5.csv"), row.names = FALSE)

### Step 8. Cost-effectiveness acceptability curve (CEAC) ###
# WTP goes up to $30M per QALY because the ICER is far above the usual thresholds
v_wtp_ceac <- c(seq(0, 1e6, by = 25000), seq(1.25e6, 30e6, by = 250000))
ceac <- data.frame(wtp = v_wtp_ceac,
                   p_ce = sapply(v_wtp_ceac, function(w) mean(w * psa_results$inc_Dqalys - psa_results$inc_Dcosts > 0)))
write.csv(ceac, file.path(dir_out, "PSA_CEAC_V.1.5.csv"), row.names = FALSE)

# Plot style
theme_trace <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0"),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        axis.text = element_text(colour = "#52514e"), plot.title = element_text(colour = "#0b0b0b"),
        plot.subtitle = element_text(colour = "#52514e"))

### Step 9. Cost-effectiveness plane ###
# Dots = PSA iterations, diamond = mean, dashed line = $150,000/QALY
p_plane <- ggplot(psa_results, aes(x = inc_Dqalys, y = inc_Dcosts / 1e6)) +
  geom_hline(yintercept = 0, colour = "#a3a29e", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "#a3a29e", linewidth = 0.4) +
  geom_abline(slope = 150000 / 1e6, intercept = 0, colour = "#eb6834", linewidth = 0.6, linetype = "dashed") +
  geom_point(colour = "#2a78d6", alpha = 0.35, size = 1.6, stroke = 0) +
  geom_point(aes(x = mean(inc_Dqalys), y = mean(inc_Dcosts) / 1e6), colour = "#0b0b0b", size = 3, shape = 18) +
  annotate("text", x = Inf, y = 150000 / 1e6 * max(psa_results$inc_Dqalys), label = "$150,000/QALY",
           hjust = 1.05, vjust = -0.6, colour = "#eb6834", size = 3.3) +
  scale_x_continuous(limits = c(min(0, min(psa_results$inc_Dqalys)), NA)) +
  scale_y_continuous(limits = c(min(0, min(psa_results$inc_Dcosts / 1e6)), NA)) +
  labs(title = "Cost-effectiveness plane: Silencer vs Stabilizer",
       subtitle = sprintf("%d PSA iterations; diamond = mean (ICER $%.2fM per QALY gained)", n_sim, icer_psa / 1e6),
       x = "Incremental discounted QALYs", y = "Incremental discounted costs ($ millions)") +
  theme_trace
ggsave(file.path(dir_out, "PSA_CE_plane_V.1.5.png"), p_plane, width = 8, height = 6, dpi = 200)

### Step 10. CEAC plot ###
p_ceac <- ggplot(ceac, aes(x = wtp / 1e6, y = p_ce)) +
  geom_vline(xintercept = 0.15, colour = "#eb6834", linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = 0.15, y = 0.95, label = "$150,000/QALY", hjust = -0.08, colour = "#eb6834", size = 3.3) +
  geom_line(colour = "#2a78d6", linewidth = 0.9) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), limits = c(0, 1)) +
  labs(title = "Cost-effectiveness acceptability curve: Silencer vs Stabilizer",
       subtitle = sprintf("Probability cost-effective: %.0f%% at $100K, %.0f%% at $150K per QALY", prob_ce[1] * 100, prob_ce[2] * 100),
       x = "Willingness to pay ($ millions per QALY gained)", y = "Probability silencer is cost-effective") +
  theme_trace
ggsave(file.path(dir_out, "PSA_CEAC_V.1.5.png"), p_ceac, width = 8, height = 5, dpi = 200)

cat("\nSaved:", file.path(dir_out, c("PSA_results_V.1.5.csv", "PSA_summary_V.1.5.csv", "PSA_CEAC_V.1.5.csv",
                                    "PSA_CE_plane_V.1.5.png", "PSA_CEAC_V.1.5.png")), sep = "\n  ")
