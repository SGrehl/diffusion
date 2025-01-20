rm(list=ls())

source("functions_simulation.R")
source("functions_animations.R")

model <- list(seed                     = 12345,
              world_size               = 13,
              world_agents_per_cell    = 6,
              agent_distance_function  = "geometric",
              agent_distance_parameter = c(0.7),
              agent_adoption_p         = 0.8)


sim_data <- run_simulation(model)

world_animation(sim_data)