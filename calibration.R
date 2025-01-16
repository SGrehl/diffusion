rm(list=ls())

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(magick)

source("functions_simulation.R")
source("functions_calibration.R")

load("historic/historic_data.RData")

historic_report <- report_generate(historic_data)

# run a single model

model_sim <- list(
  seed                     = 1234,
  world_size               = 13,
  world_agents_per_cell    = 6,
  agent_contact            = "geometric",   
  agent_contact_parameter  = 0.5,
  agent_contact_contagion  = 0.5
)

sim_data <- run_simulation_old(model_sim)
sim_report <- report_generate(sim_data)

measure_fit(historic_report, sim_report)



# If agent_contact_parameter is a single number (e.g. p in a geometric distribution),
# you might explore a range from 0.1 to 1 in steps of 0.1
param_range_contact   <- seq(0.1, 1.0, by = 0.3)

# If agent_contact_contagion is also a probability,
# perhaps you do something similar:
param_range_contagion <- seq(0.1, 1.0, by = 0.3)

base_model <- list(
  seed                     = 1234,
  world_size               = 13,
  world_agents_per_cell    = 6,
  agent_contact            = "geometric",   
  agent_contact_parameter  = NA_real_,  # placeholder
  agent_contact_contagion  = NA_real_
)

calibration_result <- calibrate_model(
  historic_report         = historic_report,
  model                   = base_model,
  param_range_contact     = param_range_contact,
  param_range_contagion   = param_range_contagion,
  n_sims                  = 5,       # for example, 5 replicates per param set
  max_steps               = 25
)

# View the combination(s) with the best fit
calibration_result$best

# See full grid
calibration_result$results %>% arrange(fit) %>% head()
