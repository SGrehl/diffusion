rm(list=ls())

source("functions_simulation.R")

model <- list(seed = 12345,
              world_size = 13,
              world_agents_per_cell = 6,
              agent_contact = "geometric",
              agent_contact_parameter = c(0.9),
              agent_contact_contagion = 0.9)


run_simulation(model, graph = TRUE)
