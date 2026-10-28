# Poster figure: background survival from age 60, left-truncated Gompertz vs 2022 US Actuarial Life Table
# Standalone (no model run). Gompertz parameters are those in TRACE_microsim_V.1.5.R (params_gompertz),
# fitted by maximum likelihood to ages >= 60 in TRACE_mortality_fitting.R (section 5).

library(ggplot2)

dir_out <- "Poster_figures"
dir.create(dir_out, showWarnings = FALSE)

df_raw_lifetable <- read.csv("2022 Actuarial Life Table.csv")

params_gompertz <- list(
  male   = list(shape = 0.09656038924, rate = 2.842950706e-05),
  female = list(shape = 0.1072543865,  rate = 8.615557054e-06)
)
age_start <- 60
age_end   <- 110

# Gompertz survival from age a to age a + t (same formula as surv_gompertz() in TRACE_mortality_fitting.R)
surv_gompertz <- function(a, g, t) exp(-(g$rate / g$shape) * (exp(g$shape * (a + t)) - exp(g$shape * a)))

# Life table survival from age a to age a + t (log-linear interpolation between integer ages)
surv_lifetable <- function(a, qx, t) {
  l <- approx(df_raw_lifetable$Exact.age, log(cumprod(c(1, 1 - qx))[seq_along(qx)]), xout = c(a, a + t))$y
  exp(l[2] - l[1])
}

qx_sex <- list(male = df_raw_lifetable$Male.Death.probability, female = df_raw_lifetable$Female.Death.probability)
sex_labels <- c(male = "Male", female = "Female")

### Curves conditional on survival to age 60 ###
v_age_fine <- seq(age_start, age_end, by = 0.1)
v_age_int  <- age_start:age_end
df_gompertz <- do.call(rbind, lapply(names(params_gompertz), function(s)
  data.frame(sex = sex_labels[[s]], age = v_age_fine, surv = surv_gompertz(age_start, params_gompertz[[s]], v_age_fine - age_start))))
df_lt <- do.call(rbind, lapply(names(qx_sex), function(s) {
  lx <- cumprod(c(1, 1 - qx_sex[[s]]))[seq_along(qx_sex[[s]])] # survivors at the start of each integer age
  keep <- df_raw_lifetable$Exact.age %in% v_age_int
  data.frame(sex = sex_labels[[s]], age = df_raw_lifetable$Exact.age[keep], surv = lx[keep] / lx[df_raw_lifetable$Exact.age == age_start])
}))
df_gompertz$sex <- factor(df_gompertz$sex, levels = sex_labels)
df_lt$sex       <- factor(df_lt$sex, levels = sex_labels)

# Fit quality: maximum absolute difference in survival from age 60, ages 60-100
fit_err <- sapply(names(params_gompertz), function(s) {
  d <- df_lt[df_lt$sex == sex_labels[[s]] & df_lt$age <= 100, ]
  max(abs(d$surv - surv_gompertz(age_start, params_gompertz[[s]], d$age - age_start)))
})
cat("Max absolute difference in survival from age 60 (ages 60-100): male", round(fit_err["male"], 4),
    "| female", round(fit_err["female"], 4), "\n")

### Survival table for the poster text: 1-, 5-, 10-year survival from ages 70 / 78 / 85 ###
tab_surv <- do.call(rbind, lapply(names(params_gompertz), function(s) do.call(rbind, lapply(c(70, 78, 85), function(a)
  do.call(rbind, lapply(c(1, 5, 10), function(t)
    data.frame(sex = sex_labels[[s]], start_age = a, years = t,
               life_table = surv_lifetable(a, qx_sex[[s]], t),
               gompertz   = surv_gompertz(a, params_gompertz[[s]], t))))))))
tab_surv$difference <- tab_surv$gompertz - tab_surv$life_table
print(tab_surv, digits = 3, row.names = FALSE)
write.csv(tab_surv, file.path(dir_out, "Survival_Gompertz_vs_lifetable_table_V.1.5.csv"), row.names = FALSE)

# Check against TRACE_mortality_fitting.R section 5: 3.5-year survival
for (a in c(71, 77.5, 80)) {
  cat(sprintf("Age %.1f, 3.5-year survival | Male: LT %.3f Gompertz %.3f | Female: LT %.3f Gompertz %.3f\n", a,
              surv_lifetable(a, qx_sex$male, 3.5),   surv_gompertz(a, params_gompertz$male, 3.5),
              surv_lifetable(a, qx_sex$female, 3.5), surv_gompertz(a, params_gompertz$female, 3.5)))
}

### Figure ###
col_sex <- c("Male" = "#2a78d6", "Female" = "#eb6834")
# Direct labels: male below-left of its curve (age 85), female above-right of its curve (age 90)
df_lab <- data.frame(sex = factor(c("Male", "Female"), levels = sex_labels), age = c(85, 90),
                     surv = c(surv_gompertz(age_start, params_gompertz$male, 25) - 0.08,
                              surv_gompertz(age_start, params_gompertz$female, 30) + 0.08),
                     hjust = c(1, 0))

p_surv <- ggplot() +
  geom_line(data = df_gompertz, aes(x = age, y = surv, colour = sex, linetype = "Gompertz (fitted)"), linewidth = 1.1) +
  geom_point(data = df_lt, aes(x = age, y = surv, colour = sex, shape = "2022 US life table"), size = 2.6, stroke = 1, fill = "#fcfcfb") +
  geom_text(data = df_lab, aes(x = age, y = surv, label = sex, colour = sex, hjust = hjust), size = 5.5, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = col_sex, guide = "none") +
  scale_linetype_manual(values = c("Gompertz (fitted)" = "solid"), name = NULL) +
  scale_shape_manual(values = c("2022 US life table" = 21), name = NULL) +
  scale_x_continuous(breaks = seq(60, 110, 10), limits = c(age_start, age_end), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  guides(linetype = guide_legend(order = 1, override.aes = list(colour = "#52514e")),
         shape    = guide_legend(order = 2, override.aes = list(colour = "#52514e"))) +
  labs(title = "Background mortality: Gompertz fit vs 2022 US life table",
       subtitle = sprintf("Survival from age 60 by sex; maximum absolute difference, ages 60-100: %.1f (male) and %.1f (female) percentage points",
                          fit_err["male"] * 100, fit_err["female"] * 100),
       caption = "Life table: SSA 2022 period life table. Gompertz: maximum likelihood fit to ages 60+, by sex.",
       x = "Age (years)", y = "Survival from age 60") +
  theme_minimal(base_size = 16) +
  theme(plot.title.position = "plot", plot.caption.position = "plot",
        plot.caption = element_text(colour = "#52514e", hjust = 0, size = 11),
        plot.subtitle = element_text(colour = "#52514e", size = 12),
        plot.title = element_text(colour = "#0b0b0b"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0"),
        axis.text = element_text(colour = "#52514e"),
        legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1), legend.spacing.y = unit(0, "pt"),
        legend.background = element_blank(), legend.box.background = element_rect(fill = "#fcfcfb", colour = NA),
        legend.key = element_blank(),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

ggsave(file.path(dir_out, "Survival_Gompertz_vs_lifetable_V.1.5.png"), p_surv, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(dir_out, "Survival_Gompertz_vs_lifetable_V.1.5.pdf"), p_surv, width = 10, height = 6.5, device = pdf) # base PDF device (vector); Cairo needs XQuartz
cat("\nSaved:", file.path(dir_out, c("Survival_Gompertz_vs_lifetable_V.1.5.png", "Survival_Gompertz_vs_lifetable_V.1.5.pdf",
                                    "Survival_Gompertz_vs_lifetable_table_V.1.5.csv")), sep = "\n  ")
