rm(list=ls())

source("functions_simulation.R")
source("functions_calibration.R")

load("data/historic_data.RData")

historic_report <- report_generate(historic_data)

# run a single model

model_sim <- list(
  seed                     = 1234,
  world_size               = 13,
  world_agents_per_cell    = 6,
  agent_distance_function  = "geometric",   
  agent_distance_parameter = 0.5,
  agent_adoption_p         = 0.5
)

sim_data <- run_simulation(model_sim)
sim_report <- report_generate(sim_data)

measure_fit(historic_report, sim_report)

# define the base model
base_model <- list(
  seed                     = 1234,
  world_size               = 13,
  world_agents_per_cell    = 6,
  agent_distance_function  = "geometric",   
  agent_distance_parameter = NA_real_,    # placeholder
  agent_adoption_p         = NA_real_     # placeholder
)

# run the calibration
calibration_result <- calibrate_model_mc(
  historic_report          = historic_report,
  model                    = base_model,
  agent_distance_parameter = seq(0.05, 1.0, by = 0.05),
  agent_adoption_p         = seq(0.05, 1.0, by = 0.05),
  n_sims                   = 50,
  max_steps                = 25
)

save(calibration_result, file =  "data/calibration_result.Rdata")

# generate a heatmap based on calibration
calibration_heatmap_pdf(calibration_result)
