# Tornado diagram for the OWSA
# Called at the end of TRACE_OWSA_V.1.5.R, or can be run alone from the saved CSV
# When the silencer is dominated (more cost, fewer QALYs), the ICER is negative and the bar would look favorable,
# so those bars are drawn to the right edge (least favorable) and labeled "Dominated"
library(ggplot2)

# If run alone, read the saved results
if (!exists("dir_out")) dir_out <- "OWSA_output"
if (!exists("owsa_results")) owsa_results <- read.csv(file.path(dir_out, "OWSA_results_V.1.5.csv"))
if (!exists("base_owsa")) base_owsa <- c(ICER = 16062431.1251219) # base case ICER from TRACE_microsim_V.1.5.R

# One row per input x bound (low / high)
df_tornado <- rbind(
  data.frame(parameter = owsa_results$parameter, bound = "Low input value",  ICER = owsa_results$ICER_low,  inc_q = owsa_results$inc_Dqalys_low),
  data.frame(parameter = owsa_results$parameter, bound = "High input value", ICER = owsa_results$ICER_high, inc_q = owsa_results$inc_Dqalys_high)
)
# Dominated bars go to the right edge
df_tornado$dominated <- df_tornado$inc_q < 0 & df_tornado$ICER < 0
x_max <- max(df_tornado$ICER[!df_tornado$dominated]) * 1.15
df_tornado$ICER_plot <- ifelse(df_tornado$dominated, x_max, df_tornado$ICER)
# Widest range on top
df_tornado$parameter <- factor(df_tornado$parameter, levels = rev(owsa_results$parameter))
df_tornado$bound <- factor(df_tornado$bound, levels = c("Low input value", "High input value"))
base_icer <- base_owsa["ICER"]

# Each bar starts at the base case ICER and ends at the low / high ICER
p_tornado <- ggplot(df_tornado) +
  geom_segment(aes(x = base_icer / 1e6, xend = ICER_plot / 1e6, y = parameter, yend = parameter, colour = bound),
               linewidth = 5) +
  geom_text(data = df_tornado[df_tornado$dominated, ], aes(x = ICER_plot / 1e6, y = parameter),
            label = "Dominated", hjust = 1.05, colour = "#fcfcfb", size = 3.2, fontface = "bold") +
  geom_vline(xintercept = base_icer / 1e6, colour = "#52514e", linewidth = 0.5) +
  scale_colour_manual(values = c("Low input value" = "#2a78d6", "High input value" = "#eb6834"), name = NULL) +
  labs(title = "One-way sensitivity analysis: Silencer vs Stabilizer",
       subtitle = sprintf("Base case ICER: $%.2fM per QALY gained", base_icer / 1e6),
       caption = "Dominated: silencer more costly with fewer QALYs (ICER not meaningful); bar drawn to the axis edge.",
       x = "ICER ($ millions per QALY gained)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title.position = "plot", plot.caption.position = "plot", plot.caption = element_text(colour = "#52514e", hjust = 0),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e5e0"),
        legend.position = "top", legend.justification = "left",
        plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        axis.text = element_text(colour = "#52514e"), plot.title = element_text(colour = "#0b0b0b"))
ggsave(file.path(dir_out, "OWSA_tornado_V.1.5.png"), p_tornado, width = 10, height = 6.5, dpi = 200)
