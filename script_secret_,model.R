rm(list=ls())

# Load necessary libraries

source("functions_simulation.R")

model <- list(seed                     = 7654,
              world_size               = 13,
              world_agents_per_cell    = 6,
              agent_distance_function  = "geometric",
              agent_distance_parameter = c(0.6),
              agent_adoption_p         = 0.55)

historic_data <- run_simulation(model)

save(historic_data, file = "data/historic_data.RData")

source("functions_animations.R")

animation <- world_animation(historic_data)
