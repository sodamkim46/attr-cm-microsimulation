# --- Code Snippet for Verifying Mortality Rates ---

# 1. Define key parameters used in the main script
cycle_length <- 0.5 # 6 months
params_gompertz_male <- list(
  shape = 0.1014027, 
  rate  = 1.934313e-05    
)
rr_HF   <- 3.17
rr_2v1  <- 1.78
rr_3v1  <- 3.51
rr_4v1  <- 5.74

# 2. Define information for a hypothetical patient
initial_age <- 77 # Starting age
num_cycles_to_check <- 50 # Number of cycles to check (25 years)

# 3. Create a data frame to store the results
mortality_over_time <- data.frame(
  Cycle = 0:num_cycles_to_check,
  Age = initial_age + (0:num_cycles_to_check) * cycle_length,
  p_NYHA1_Death = numeric(num_cycles_to_check + 1),
  p_NYHA2_Death = numeric(num_cycles_to_check + 1),
  p_NYHA3_Death = numeric(num_cycles_to_check + 1),
  p_NYHA4_Death = numeric(num_cycles_to_check + 1)
)

# 4. Calculate mortality rates for each time point by iterating through cycles
for (i in 1:nrow(mortality_over_time)) {
  current_age <- mortality_over_time$Age[i]
  
  # Calculate annual background mortality hazard for the current age using the Gompertz distribution
  h_H <- params_gompertz_male$shape * exp(params_gompertz_male$rate * current_age)
  
  # Calculate state-specific annual mortality hazards
  h_1 <- h_H * rr_HF
  h_2 <- h_1 * rr_2v1
  h_3 <- h_1 * rr_3v1
  h_4 <- h_1 * rr_4v1
  
  # Convert annual hazards to 6-month (cycle_length) mortality probabilities
  p_1toD <- 1 - exp(-h_1 * cycle_length)
  p_2toD <- 1 - exp(-h_2 * cycle_length)
  p_3toD <- 1 - exp(-h_3 * cycle_length)
  p_4toD <- 1 - exp(-h_4 * cycle_length)
  
  # Store the results in the data frame
  mortality_over_time$p_NYHA1_Death[i] <- p_1toD
  mortality_over_time$p_NYHA2_Death[i] <- p_2toD
  mortality_over_time$p_NYHA3_Death[i] <- p_3toD
  mortality_over_time$p_NYHA4_Death[i] <- p_4toD
}

# 5. Print the results (first 5 years and last 5 years)
print("Change in 6-month mortality probability for a hypothetical patient (male, starting at age 77):")
print(head(mortality_over_time, 11)) # Data for the first 5 years of the simulation
print(tail(mortality_over_time, 11)) # Data for the last 5 years of the simulation
