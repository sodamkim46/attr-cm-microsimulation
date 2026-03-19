library(pacman)
p_load(readxl, tidyverse, flexsurv)

# Read in the SSA actuarial life table.
# This file gives us, for each exact age x:
#   - prob  = the probability of dying between age x and x+1
#   - lives = the number of people alive at exact age x out of the life table radix
#
# In other words, the life table is already a grouped survival dataset:
# it tells us how many people are alive at the start of each age interval
# and what proportion die during that interval.
X2022_Actuarial_Life_Table <- read_excel("2022 Actuarial Life Table.xlsx")
# View(X2022_Actuarial_Life_Table)

lt <- X2022_Actuarial_Life_Table %>%
  mutate(
    # Rename variables to make downstream code easier to read.
    age = `Exact age`,
    prob = `Male Death probability`,
    lives = `Male Number of lives`,
    
    # lives_next is the number alive at the next exact age.
    # So if 91,000 are alive at age 50 and 90,200 are alive at age 51,
    # then about 800 died in the interval [50, 51).
    lives_next = lead(lives),
    
    # wt is being used as a frequency weight.
    #
    # Intuition:
    # We do NOT have one row per person.
    # Instead, we have one row per age interval.
    #
    # The quantity lives - lives_next is approximately the number of deaths
    # that occurred during that one-year interval.
    #
    # By passing this as the weight, we are telling flexsurv that this interval
    # represents many identical death events occurring somewhere between age and age+1.
    #
    # This is the key trick:
    # rather than simulating a synthetic individual-level dataset,
    # we use the life table as grouped survival data and let weights stand in
    # for repeated identical observations.
    wt = if_else(!is.na(lives_next), lives - lives_next, lives * prob),
    
    # age_next defines the end of the interval.
    # So each row now represents the interval [age, age+1).
    age_next = age + 1
  ) %>%
  # Restrict to ages 50-100.
  #
  # This is important because a pure Gompertz model is generally intended
  # to describe adult mortality, where the hazard rises approximately exponentially.
  # It is usually not a good model for childhood mortality or very early adulthood.
  filter(age >= 50, age <= 100)

# Fit a Gompertz survival model using flexsurv.
#
# Why does this work?
# -------------------
# Surv(age, age_next, type = "interval2") tells flexsurv that the event time
# is not known exactly, but is known to lie somewhere in the interval [age, age+1).
#
# That matches the life table structure very naturally:
# we usually do NOT know the exact death time for the grouped individuals,
# only that they died during that one-year age band.
#
# So each row is being treated like:
#   "wt deaths occurred sometime between age and age+1"
#
# This is conceptually similar to a counting-process / start-stop dataset,
# except that here the row represents grouped interval-censored events rather
# than one person's follow-up interval.
#
# In other words:
#   - start = age
#   - stop  = age + 1
#   - event occurred somewhere in that interval
#   - weight = how many such deaths that row represents
#
# So we are leveraging the fact that the life table already contains
# the age-specific event structure needed for survival modeling.
fit <- flexsurvreg(
  Surv(age, age_next, type = "interval2") ~ 1,
  data = lt,
  dist = "gompertz",
  weights = wt
)

# Get predicted survival from the fitted Gompertz model at the start of each interval.
#
# Sx = S(x) = probability of surviving to exact age x under the fitted model.
# For example, at x = 50, this is the model-implied probability of surviving to age 50.
Sx <- summary(fit, type = "survival", t = lt$age, tidy = TRUE)$est

# Get predicted survival one year later, at age x+1.
#
# Sx1 = S(x+1) = probability of surviving to exact age x+1 under the fitted model.
Sx1 <- summary(fit, type = "survival", t = lt$age + 1, tidy = TRUE)$est

# Convert survival into one-year death probabilities.
#
# Why does this formula work?
# ---------------------------
# The conditional probability of surviving from age x to age x+1 is:
#
#   P(T > x+1 | T > x) = S(x+1) / S(x)
#
# Therefore, the conditional probability of dying between age x and x+1 is:
#
#   qx = 1 - S(x+1) / S(x)
#
# This is directly analogous to the life table's annual death probability column.
# So this gives us model-based predicted qx values that are directly comparable
# to the observed life-table probabilities.
qx_pred <- 1 - (Sx1 / Sx)

# Attach the model-implied annual death probabilities back to the life table.
#
# Now each age interval has:
#   - prob      = observed SSA one-year death probability
#   - pred_prob = Gompertz model one-year death probability
lt <- lt %>%
  mutate(pred_prob = qx_pred)

# Create survival estimates from observed and predicted 1-year death probabilities.
#
# Why can we reconstruct survival from annual death probabilities?
# ---------------------------------------------------------------
# In a life table, survival across ages is built recursively.
#
# If q50 is the probability of dying between 50 and 51,
# then the probability of surviving that interval is:
#
#   1 - q50
#
# To survive from age 50 to age 52, you must survive both [50, 51) and [51, 52):
#
#   (1 - q50) * (1 - q51)
#
# More generally, conditional survival from the starting age is the cumulative
# product of the interval-specific survival probabilities.
#
# So cumprod(1 - prob) gives us the survival curve implied by the life table.
# Likewise, cumprod(1 - pred_prob) gives the survival curve implied by the
# Gompertz model's predicted annual death probabilities.
#
# This is why the conversion from qx to survival is so straightforward:
# life tables and survival curves are just two different ways of representing
# the same underlying mortality process.
lt_surv <- lt %>%
  arrange(age) %>%
  mutate(
    surv_obs  = cumprod(1 - prob),
    surv_pred = cumprod(1 - pred_prob)
  ) %>%
  select(age, surv_obs, surv_pred) %>%
  pivot_longer(
    cols = c(surv_obs, surv_pred),
    names_to = "source",
    values_to = "survival"
  ) %>%
  mutate(
    source = recode(
      source,
      surv_obs = "SSA life table",
      surv_pred = "Gompertz model"
    )
  )

# Overlay survival curves.
#
# This plot makes comparison easier because:
#   - the life table curve shows the empirical/grouped survival pattern
#   - the Gompertz curve shows the smoothed parametric pattern implied by the model
#
# If the fit is good, the two curves should lie close together.
# If they diverge, that suggests the Gompertz form is not capturing the age pattern
# of mortality well over this range.
ggplot(lt_surv, aes(x = age, y = survival, color = source, linetype = source)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Male survival: SSA life table vs Gompertz model",
    subtitle = "Derived from annual death probabilities, ages 50-100",
    x = "Age",
    y = "Survival probability",
    color = "",
    linetype = ""
  ) +
  theme_minimal(base_size = 13)