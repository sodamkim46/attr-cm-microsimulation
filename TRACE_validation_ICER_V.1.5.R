# TRACE V1.5 cross-validation against the ICER ATTR-CM model (Excel, 06dec2024): tafamidis vs placebo.
# The V1.5 model functions are loaded unchanged; only the inputs are replaced with the ICER model inputs.
# Known structural differences (kept and quantified, not removed): Gompertz vs life-table background mortality,
# and the ICER mortality adjustment factor applied to the cycle probability (Excel) vs the hazard (V1.5).

rm(list = ls())

### 1. Load the V1.5 model functions only ###
# Evaluates only the top-level `name <- function(...)` definitions, so the V1.5 base case (inputs, calibration,
# runs, output files) is not executed and no V1.5 input enters this script
source_functions_only <- function(file) {
  exprs <- parse(file)
  is_fun_def <- vapply(exprs, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
  }, logical(1))
  for (e in exprs[is_fun_def]) eval(e, envir = globalenv())
  invisible(vapply(exprs[is_fun_def], function(e) as.character(e[[2]]), character(1)))
}

v_loaded_funs <- source_functions_only("TRACE_microsim_V.1.5.R")
print(v_loaded_funs)

# Check: the workspace holds only functions (no V1.5 inputs such as num_i or calibrated_theta)
v_objs <- setdiff(ls(), "v_loaded_funs")
stopifnot(all(vapply(v_objs, function(x) is.function(get(x)), logical(1))))
cat("Workspace check passed:", length(v_objs), "objects, all functions\n")

### 2. Read the ICER model results and traces (Excel) ###
# Read straight from the workbook to avoid copying errors. Structure checked in the workbook:
# - Mortality adjustment factor on the cycle death probability: 2.25 in both arms for cycles 1-3
#   (Mortality Calibration!AG5:AG7), then m_ta = 1.2 (tafamidis) and m_pta = 2.75 (placebo) from cycle 4
# - Discontinued tafamidis patients use the placebo transitions (p_*_tadis_* = p_*_pta_*) but keep the
#   tafamidis mortality factor
# - Discontinuation (dis_ta) applies in cycles 1-5 only (trial period, 30 months); from cycle 6 the formula has no dis_ta
# - Drug costs apply to on-treatment NYHA I-III only (NYHA IV on treatment has no drug cost; NYHA IV is absorbing)
# - Totals sum rows 5-106 with SUMIF(">0"): row 5 is cycle 0 counted as a full cycle (start occupancy x 0.5 years),
#   then rows 6+ are half-cycle (trapezoid) averages. V1.5 has no separate cycle-0 term, so both are kept below
library(readxl)

path_icer <- file.path("..", "..", "ATTR-CM Model _ with acoramidis price 06dec2024 (drug cost input updated, PSA, OWSA updated).xlsm")
v_icer_sheets <- c(tafamidis = "Tafamidis", placebo = "Placebo (Tafamidis)")

read_block <- function(sheet, range, col_names, ...) { # ...: passed to read_excel (e.g., col_types)
  df <- suppressMessages(read_excel(path_icer, sheet = sheet, range = range, col_names = FALSE, ...))
  setNames(as.data.frame(df), col_names)
}

# Trace rows 5-105 = cycles 0-100 (50 years); cells past age 100 are errors in Excel, read as NA
read_icer_trace <- function(sheet) {
  # K:N on treatment NYHA I-IV, O:R discontinued NYHA I-IV, S dead (start of each cycle)
  occ <- read_block(sheet, "K5:S105", c(paste0("tx_", 1:4), paste0("dc_", 1:4), "dead"))
  # AY:BF discounted payoffs per row (half-cycle corrected)
  pay <- read_block(sheet, "AY5:BF105", c("drug", "hosp", "nondrug", "total", "ly", "qaly", "evly", "nyha12"))
  tr <- cbind(cycle = 0:100, occ, pay)
  tr[is.na(tr)] <- 0
  tr$dead[tr$cycle > 0 & rowSums(occ, na.rm = TRUE) == 0] <- 1 # everyone dead after the life table ends
  tr$alive <- 1 - tr$dead
  tr
}

l_icer_trace <- lapply(v_icer_sheets, read_icer_trace)

# Discounted totals: as reported (all rows) and without the cycle-0 row (comparable to V1.5)
v_outcomes <- c("drug", "hosp", "nondrug", "total", "ly", "qaly", "nyha12")
icer_totals <- function(tr, drop_cycle0) {
  keep <- if (drop_cycle0) tr$cycle > 0 else rep(TRUE, nrow(tr))
  colSums(tr[keep, v_outcomes])
}
m_icer_reported <- sapply(l_icer_trace, icer_totals, drop_cycle0 = FALSE) # rows = outcomes, columns = arms
m_icer_no_c0    <- sapply(l_icer_trace, icer_totals, drop_cycle0 = TRUE)
calc_icer <- function(m) (m["total", "tafamidis"] - m["total", "placebo"]) / (m["qaly", "tafamidis"] - m["qaly", "placebo"])

# Check: totals rebuilt from the traces match the Results sheet
icer_results <- read_block("Results", "C6:J7", c("drug", "hosp", "nondrug", "total", "ly", "qaly", "evly", "nyha12"))
stopifnot(isTRUE(all.equal(unname(m_icer_reported["total", ]), icer_results$total)),
          isTRUE(all.equal(unname(m_icer_reported["qaly", ]),  icer_results$qaly)))
icer_icer_reported <- read_block("Results", "E15:E15", "icer")$icer

cat("\nICER model, discounted totals (as reported):\n");   print(round(m_icer_reported, 3))
cat("\nICER model, discounted totals (without cycle-0 row):\n"); print(round(m_icer_no_c0, 3))
stopifnot(isTRUE(all.equal(unname(calc_icer(m_icer_reported)), icer_icer_reported)))
cat("\nICER (as reported):", round(icer_icer_reported), "\n")
cat("ICER (without cycle-0 row):", round(calc_icer(m_icer_no_c0)), "\n")
cat("30-month survival (cycle 5):", round(sapply(l_icer_trace, function(tr) tr$alive[tr$cycle == 5]), 4), "\n")

### 3. ICER model inputs in the V1.5 format ###
# Inputs sheet: column B = input name, column F = live value; read once, then look up by name
icer_in_tab <- suppressWarnings(read_block("Inputs", "B1:F800", c("name", "desc", "value"),
                                           col_types = c("text", "text", "skip", "skip", "numeric")))
icer_in <- setNames(icer_in_tab$value, icer_in_tab$name)
get_in <- function(name) {
  x <- unname(icer_in[name])
  if (is.na(x)) stop("ICER input not found: ", name)
  x
}

## General
icer_num_i        <- 1e5
icer_seed         <- 1234
icer_cycle_length <- read_block("Settings", "J5:J5", "x")$x   # 0.5 years
icer_num_cycles   <- 100                                     # 50 years (trace rows 6-105)
icer_disc_costs   <- read_block("Settings", "D15:D15", "x")$x
icer_disc_qalys   <- read_block("Settings", "D17:D17", "x")$x

## Population: one starting age (cohort model), separate male and female runs weighted by sex
icer_v_states_names <- c("NYHA1", "NYHA2", "NYHA3", "NYHA4", "Death") # state names used inside the V1.5 functions
icer_prop_male  <- 1 - get_in("female")
icer_m_features <- cbind("age" = rep(get_in("age"), icer_num_i))
icer_p_start    <- c(get_in("start_I"), get_in("start_II"), get_in("start_III"), get_in("start_IV"), 0)
set.seed(icer_seed)
icer_v_start <- sample(x = icer_v_states_names, size = icer_num_i, replace = TRUE, prob = icer_p_start)

## Mortality
# V1.5: HR_k = 1 + theta * (rr_HF * rr_k - 1). With theta = 1 this is rr_HF * rr_k, so rr_HF carries the
# ICER mortality adjustment factor (applied to the hazard in V1.5, to the cycle probability in Excel)
icer_theta    <- 1
v_mort_first3 <- read_block("Mortality Calibration", "AG5:AG7", "x")$x # cycles 1-3, both arms
icer_mort_factor <- list(                                               # cycles 1-5 (cycle 5 value used afterwards)
  tafamidis = c(v_mort_first3, get_in("m_ta"),  get_in("m_ta")),
  placebo   = c(v_mort_first3, get_in("m_pta"), get_in("m_pta"))
)
# Background mortality: the V1.5 Gompertz parameters, kept on purpose (Excel looks up the life table directly)
icer_params_gompertz <- list(
  male   = list(shape = 0.09656038924, rate = 2.842950706e-05),
  female = list(shape = 0.1072543865,  rate = 8.615557054e-06)
)

## NYHA transitions: one list per 6-month cycle (cycles 1-5), in the V1.5 format
# Excel names: p_<from>2<to>_<arm>_<month>, e.g. p_II2III_ta_6 (arm: ta, pta, tadis = tafamidis after discontinuation)
v_roman  <- c("I", "II", "III", "IV")
v_months <- c(6, 12, 18, 24, 30)
make_icer_trans <- function(arm, v_factor) {
  lapply(seq_along(v_months), function(k) {
    l <- list()
    for (i in 1:4) for (j in 1:4) {
      l[[paste0("p_", i, "to", j)]] <- get_in(paste0("p_", v_roman[i], "2", v_roman[j], "_", arm, "_", v_months[k]))
    }
    c(l, rr_HF = v_factor[k], rr_2v1 = get_in("m_II v I"), rr_3v1 = get_in("m_III v I"), rr_4v1 = get_in("m_IV v I"))
  })
}
l_icer_trans <- list(
  ta    = make_icer_trans("ta",    icer_mort_factor$tafamidis), # tafamidis, on treatment
  tadis = make_icer_trans("tadis", icer_mort_factor$tafamidis), # tafamidis, discontinued (keeps tafamidis mortality)
  pta   = make_icer_trans("pta",   icer_mort_factor$placebo)    # placebo
)
# Check: each row of every matrix sums to 1
for (l in unlist(l_icer_trans, recursive = FALSE)) {
  for (i in 1:4) stopifnot(abs(sum(unlist(l[paste0("p_", i, "to", 1:4)])) - 1) < 1e-9)
}

## Costs per cycle: health state costs (supportive care + productivity, 0 in the health care perspective) + net drug cost
get_by_class <- function(prefix, suffix = "") sapply(v_roman, function(r) get_in(paste0(prefix, r, suffix)))
icer_v_costs_soc <- setNames(c(get_by_class("c_suppcare") + get_by_class("c_prod"), 0), icer_v_states_names)
icer_v_costs_ta  <- icer_v_costs_soc + c(rep(get_in("c_ta"), 4), 0) # NYHA IV has no drug cost: treatment stops on entry

## Utilities and hospitalization (same in both arms in Excel; checked)
stopifnot(all(get_by_class("u_", "_ta") == get_by_class("u_", "_pta")),
          all(get_by_class("p_hosp", "_ta") == get_by_class("p_hosp", "_pta")))
icer_v_utils         <- setNames(c(get_by_class("u_", "_ta"), 0), icer_v_states_names)
icer_v_hosp_probs    <- setNames(c(get_by_class("p_hosp", "_ta"), 0), icer_v_states_names)
icer_v_hosp_costs    <- setNames(c(get_by_class("c_hosp"), 0), icer_v_states_names)
icer_v_hosp_disutils <- setNames(c(-get_by_class("u_dishosp"), 0), icer_v_states_names) # Excel stores negatives and adds them

## Discontinuation: per-cycle probability, cycles 1-5 only (trial period)
icer_p_disc         <- get_in("dis_ta")
icer_disc_max_cycle <- 5

cat("\nICER inputs loaded. Mortality factor by cycle (1-5): tafamidis", icer_mort_factor$tafamidis,
    "| placebo", icer_mort_factor$placebo, "\n")
print(rbind(cost_soc = icer_v_costs_soc, cost_ta = icer_v_costs_ta, utility = icer_v_utils,
            hosp_prob = icer_v_hosp_probs, hosp_cost = icer_v_hosp_costs, hosp_disutil = icer_v_hosp_disutils))

### 4. Run the V1.5 model with the ICER inputs: 2 arms x 2 sexes ###
# Same seed in every run, so both arms use the same random numbers (common random numbers)
run_icer_arm <- function(arm, sex, gompertz = icer_params_gompertz) { # gompertz: background mortality (Block 7 swaps it)
  is_ta <- arm == "tafamidis"
  run_microSimV_hosp_disc(
    v_starting_states   = icer_v_start,
    num_i               = icer_num_i,
    num_cycles          = icer_num_cycles,
    m_indi_features     = icer_m_features,
    v_states_names      = icer_v_states_names,
    v_states_costs_arm  = if (is_ta) icer_v_costs_ta else icer_v_costs_soc,
    v_states_costs_soc  = icer_v_costs_soc,
    v_states_utilities  = icer_v_utils,
    l_trans_probs_arm   = if (is_ta) l_icer_trans$ta else l_icer_trans$pta,
    l_trans_probs_soc   = if (is_ta) l_icer_trans$tadis else l_icer_trans$pta, # used after discontinuation
    p_disc              = if (is_ta) icer_p_disc else 0,
    arm_name            = arm, # not "silencer": no residual effect after discontinuation
    params_gompertz_sex = gompertz[[sex]],
    theta_arm           = icer_theta,
    theta_soc           = icer_theta,
    v_hosp_probs        = icer_v_hosp_probs,
    v_hosp_costs        = icer_v_hosp_costs,
    v_hosp_disutils     = icer_v_hosp_disutils,
    discount_rate_costs = icer_disc_costs,
    discount_rate_QALYs = icer_disc_qalys,
    cycle_length        = icer_cycle_length,
    starting_seed       = icer_seed,
    disc_max_cycle      = icer_disc_max_cycle
  )
}

# Keep only what the comparison needs (the full individual-level matrices are large)
v_nyha <- icer_v_states_names[1:4] # global name used inside calc_state_cycles() (V1.5)
icer_v_c_dsc <- calc_discount_wts(icer_disc_costs, icer_num_cycles, icer_cycle_length)
icer_v_e_dsc <- calc_discount_wts(icer_disc_qalys, icer_num_cycles, icer_cycle_length)
summarize_run <- function(res) list(
  totals   = c(total = res$mean_Dcosts, ly = res$mean_Dlys, qaly = res$mean_Dqalys, ontx = res$mean_DOnTx),
  cycles_c = calc_state_cycles(res, wts = icer_v_c_dsc), # discounted cycles in NYHA I-IV (cost discount rate)
  cycles_e = calc_state_cycles(res, wts = icer_v_e_dsc), # discounted cycles in NYHA I-IV (outcome discount rate)
  occ      = sapply(icer_v_states_names, function(s) colMeans(res$m_States == s)) # proportion by state, cycles 0-100
)

# Both arms x both sexes, then sex-weighted: every element (vector or matrix) weighted the same way
wavg <- function(x_m, x_f) x_m * icer_prop_male + x_f * (1 - icer_prop_male)
run_icer_all <- function(gompertz = icer_params_gompertz) {
  l_icer_runs <- list()
  for (arm in c("tafamidis", "placebo")) {
    for (sex in c("male", "female")) {
      cat("Running", arm, sex, "...\n")
      l_icer_runs[[arm]][[sex]] <- summarize_run(run_icer_arm(arm, sex, gompertz))
    }
  }
  lapply(l_icer_runs, function(r) Map(wavg, r$male, r$female))
}
l_trace_res <- run_icer_all()

cat("\n30-month survival (cycle 5), V1.5 vs ICER Excel:\n")
print(round(rbind(V1.5  = sapply(l_trace_res,  function(r)  1 - r$occ[6, "Death"]), # row 6 = cycle 5
                  Excel = sapply(l_icer_trace, function(tr) tr$alive[tr$cycle == 5])), 4))

### 5. Comparison table: V1.5 vs ICER Excel (discounted) ###
# V1.5 totals split like the Excel results: drug = net drug cost per cycle x discounted on-treatment cycles,
# non-drug = health state costs x discounted cycles in each NYHA class, hospital = the remainder
v15_totals <- function(r, drug_cycle) {
  drug    <- drug_cycle * r$totals[["ontx"]]
  nondrug <- sum(icer_v_costs_soc[v_nyha] * r$cycles_c)
  c(drug = drug, hosp = r$totals[["total"]] - drug - nondrug, nondrug = nondrug, total = r$totals[["total"]],
    ly = r$totals[["ly"]], qaly = r$totals[["qaly"]],
    nyha12 = sum(r$cycles_e[c("NYHA1", "NYHA2")]) * icer_cycle_length) # years in NYHA I-II
}
m_v15 <- cbind(tafamidis = v15_totals(l_trace_res$tafamidis, get_in("c_ta")),
               placebo   = v15_totals(l_trace_res$placebo,   0))

# One row per arm x outcome, then incremental (tafamidis - placebo) and the ICER
# excel_reported = Results sheet; excel = without the cycle-0 row (the comparison basis for V1.5)
make_rows <- function(arm, f) {
  data.frame(arm = arm, outcome = v_outcomes,
             excel_reported = f(m_icer_reported), excel = f(m_icer_no_c0), v15 = f(m_v15))
}
inc <- function(m) m[v_outcomes, "tafamidis"] - m[v_outcomes, "placebo"]
df_compare <- rbind(
  make_rows("tafamidis",   function(m) m[v_outcomes, "tafamidis"]),
  make_rows("placebo",     function(m) m[v_outcomes, "placebo"]),
  make_rows("incremental", inc),
  data.frame(arm = "incremental", outcome = "icer", excel_reported = calc_icer(m_icer_reported),
             excel = calc_icer(m_icer_no_c0), v15 = calc_icer(m_v15))
)
df_compare$diff     <- df_compare$v15 - df_compare$excel
df_compare$pct_diff <- 100 * df_compare$diff / df_compare$excel
rownames(df_compare) <- NULL

dir_val <- "Validation_output"
dir.create(dir_val, showWarnings = FALSE)
write.csv(df_compare, file.path(dir_val, "ICER_validation_comparison_V.1.5.csv"), row.names = FALSE)

cat("\n--- V1.5 vs ICER Excel (discounted; % diff vs Excel without the cycle-0 row) ---\n")
df_print <- df_compare
df_print[, 3:6] <- round(df_print[, 3:6], 3)
df_print$pct_diff <- round(df_print$pct_diff, 1)
print(df_print)

### 6. Curves: survival and NYHA class over time, V1.5 vs ICER Excel ###
library(ggplot2)

# Long data: one row per arm x model x cycle x state (proportion of the starting cohort)
v_arm_labels <- c(tafamidis = "Tafamidis", placebo = "Placebo")
v_years      <- (0:icer_num_cycles) * icer_cycle_length
df_curves <- do.call(rbind, lapply(names(v_arm_labels), function(arm) {
  tr <- l_icer_trace[[arm]]
  occ_excel <- cbind(sapply(1:4, function(k) tr[[paste0("tx_", k)]] + tr[[paste0("dc_", k)]]), tr$dead) # on tx + discontinued
  colnames(occ_excel) <- icer_v_states_names
  l_occ <- list("TRACE V1.5" = l_trace_res[[arm]]$occ, "ICER model (Excel)" = occ_excel)
  do.call(rbind, lapply(names(l_occ), function(model) {
    data.frame(arm = v_arm_labels[[arm]], model = model, years = v_years,
               state = rep(icer_v_states_names, each = length(v_years)), prop = as.vector(l_occ[[model]]))
  }))
}))
df_surv <- subset(df_curves, state == "Death")
df_surv$surv <- 1 - df_surv$prop

# Survival at 30 months, 5 and 10 years
df_surv_milestones <- reshape(subset(df_surv, years %in% c(2.5, 5, 10))[, c("arm", "model", "years", "surv")],
                              idvar = c("arm", "years"), timevar = "model", direction = "wide")
names(df_surv_milestones) <- c("arm", "years", "v15", "excel")
df_surv_milestones$diff <- df_surv_milestones$v15 - df_surv_milestones$excel
write.csv(df_surv_milestones, file.path(dir_val, "ICER_validation_survival_V.1.5.csv"), row.names = FALSE)
cat("\nSurvival milestones:\n"); print(df_surv_milestones, digits = 3, row.names = FALSE)

# Shared style (same as the poster figures): arm = colour, model = line (V1.5) vs points (Excel, yearly)
col_arm   <- c(Tafamidis = "#2a78d6", Placebo = "#eb6834")
max_years <- 20
plot_v15_vs_excel <- function(df, y) {
  ggplot() +
    geom_line(data = subset(df, model == "TRACE V1.5" & years <= max_years),
              aes(x = years, y = .data[[y]], colour = arm, linetype = model), linewidth = 1.1) +
    geom_point(data = subset(df, model == "ICER model (Excel)" & years <= max_years & years %% 1 == 0),
               aes(x = years, y = .data[[y]], colour = arm, shape = model), size = 2.6, stroke = 1, fill = "#fcfcfb") +
    scale_colour_manual(values = col_arm, name = NULL) +
    scale_linetype_manual(values = c("TRACE V1.5" = "solid"), name = NULL) +
    scale_shape_manual(values = c("ICER model (Excel)" = 21), name = NULL) +
    scale_x_continuous(breaks = seq(0, max_years, 5), expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    guides(linetype = guide_legend(order = 1, override.aes = list(colour = "#52514e")),
           shape    = guide_legend(order = 2, override.aes = list(colour = "#52514e"))) +
    theme_minimal(base_size = 14) +
    theme(plot.title.position = "plot", plot.caption.position = "plot",
          plot.title = element_text(colour = "#0b0b0b"),
          plot.subtitle = element_text(colour = "#52514e", size = 11),
          plot.caption = element_text(colour = "#52514e", hjust = 0, size = 9.5),
          panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0"),
          axis.text = element_text(colour = "#52514e"), axis.title = element_text(colour = "#52514e"),
          strip.text = element_text(colour = "#0b0b0b", face = "bold"),
          legend.position = "top", legend.justification = "left",
          plot.background = element_rect(fill = "#fcfcfb", colour = NA))
}
caption_diff <- paste("Known structural differences: background mortality (V1.5 Gompertz fitted to the 2022 SSA life table vs",
                      "the Excel life table looked up at the rounded-down age),\nmortality factor applied to the hazard (V1.5)",
                      "vs the cycle probability (Excel). Excel points shown yearly.")

# Figure 1: overall survival
p_surv <- plot_v15_vs_excel(df_surv, "surv") +
  labs(title = "Overall survival: TRACE V1.5 engine vs ICER model, same inputs",
       subtitle = "Tafamidis vs placebo, ICER ATTR-CM model inputs; V1.5 sex-weighted (90.2% male)",
       x = "Years since model start", y = "Alive", caption = caption_diff)
ggsave(file.path(dir_val, "ICER_validation_survival_V.1.5.png"), p_surv, width = 10, height = 6.5, dpi = 200)

# Figure 2: proportion of the starting cohort in each NYHA class
df_nyha <- subset(df_curves, state != "Death")
df_nyha$state <- factor(df_nyha$state, levels = v_nyha, labels = paste("NYHA", v_roman))
p_nyha <- plot_v15_vs_excel(df_nyha, "prop") +
  facet_wrap(~ state, nrow = 2) +
  labs(title = "NYHA class over time: TRACE V1.5 engine vs ICER model, same inputs",
       subtitle = "Proportion of the starting cohort in each class (on treatment + discontinued)",
       x = "Years since model start", y = "Share of starting cohort", caption = caption_diff)
ggsave(file.path(dir_val, "ICER_validation_NYHA_V.1.5.png"), p_nyha, width = 10, height = 7.5, dpi = 200)

### 7. Decomposition: how much of the gap is the life table data? ###
# Refit the V1.5 background mortality with the same method (left-truncated maximum likelihood Gompertz, ages >= 60;
# fit_gompertz_lt() in TRACE_mortality_fitting.R) to the Excel life table (Mortality sheet, sex-specific columns),
# rerun V1.5, and compare again. What remains is the calculation (structural) differences.
source_functions_only("TRACE_mortality_fitting.R") # loads fit_gompertz_lt() without running the fitting script
df_excel_lt <- read_block("Mortality", "C12:D112", c("qx_male", "qx_female")) # ages 0-100
icer_params_gompertz_excel_lt <- list(
  male   = fit_gompertz_lt(0:100, df_excel_lt$qx_male,   age_min = 60),
  female = fit_gompertz_lt(0:100, df_excel_lt$qx_female, age_min = 60)
)
cat("\nGompertz fitted to the Excel life table:\n"); str(icer_params_gompertz_excel_lt)

l_trace_res_excel_lt <- run_icer_all(icer_params_gompertz_excel_lt)
m_v15_excel_lt <- cbind(tafamidis = v15_totals(l_trace_res_excel_lt$tafamidis, get_in("c_ta")),
                        placebo   = v15_totals(l_trace_res_excel_lt$placebo,   0))

# Key outcomes, % difference vs Excel (without the cycle-0 row)
v_key <- c("total", "ly", "qaly")
key_vals <- function(m, l_res) {
  surv_at <- function(arm, years) 1 - l_res[[arm]]$occ[years / icer_cycle_length + 1, "Death"]
  c(m[v_key, "tafamidis"], m[v_key, "placebo"], inc(m)[v_key], calc_icer(m),
    surv_at("tafamidis", 2.5), surv_at("tafamidis", 5), surv_at("placebo", 2.5), surv_at("placebo", 5))
}
excel_surv_at <- function(arm, years) l_icer_trace[[arm]]$alive[l_icer_trace[[arm]]$cycle == years / icer_cycle_length]
df_decomp <- data.frame(
  outcome = c(paste0("tafamidis_", v_key), paste0("placebo_", v_key), paste0("incremental_", v_key), "icer",
              "tafamidis_surv_2.5y", "tafamidis_surv_5y", "placebo_surv_2.5y", "placebo_surv_5y"),
  excel            = c(c(m_icer_no_c0[v_key, "tafamidis"], m_icer_no_c0[v_key, "placebo"], inc(m_icer_no_c0)[v_key],
                         calc_icer(m_icer_no_c0)),
                       excel_surv_at("tafamidis", 2.5), excel_surv_at("tafamidis", 5),
                       excel_surv_at("placebo", 2.5),   excel_surv_at("placebo", 5)),
  v15_ssa_lt       = key_vals(m_v15, l_trace_res),                   # V1.5 base: Gompertz fitted to the 2022 SSA table
  v15_excel_lt     = key_vals(m_v15_excel_lt, l_trace_res_excel_lt)  # Gompertz fitted to the Excel life table
)
df_decomp$pct_ssa_lt   <- 100 * (df_decomp$v15_ssa_lt   / df_decomp$excel - 1)
df_decomp$pct_excel_lt <- 100 * (df_decomp$v15_excel_lt / df_decomp$excel - 1)
write.csv(df_decomp, file.path(dir_val, "ICER_validation_decomposition_V.1.5.csv"), row.names = FALSE)

cat("\n--- Decomposition: % difference vs Excel, before and after swapping in the Excel life table ---\n")
df_print <- df_decomp
df_print[, 2:4] <- signif(df_print[, 2:4], 5)
df_print[, 5:6] <- round(df_print[, 5:6], 1)
print(df_print)
