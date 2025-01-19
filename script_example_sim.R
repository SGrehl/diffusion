rm(list=ls())

source("functions_simulation.R")
source("functions_animations.R")

model <- list(seed = 12345,
              world_size = 13,
              world_agents_per_cell = 6,
              agent_contact = "geometric",
              agent_contact_parameter = c(0.8),
              agent_contact_contagion = 0.8)


sim_data <- run_simulation(model)

world_animation(sim_data)

