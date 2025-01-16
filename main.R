rm(list=ls())

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(magick)

source("functions_simulation.R")

model <- list(seed = 12300,
              world_size = 13,
              world_agents_per_cell = 6,
              agent_contact = "geometric", # or  weibull with c(0.9, 0.5)
              agent_contact_parameter = c(0.6),
              agent_contact_contagion = 0.75)


run_simulation(model, graph = TRUE)
