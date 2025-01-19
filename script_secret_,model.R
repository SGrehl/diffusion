rm(list=ls())

# Load necessary libraries
#library(magick)
library(ragg) # only used when on a RStudio Server
#library(magick)

source("functions_simulation.R")

model <- list(seed = 12300,
              world_size = 13,
              world_agents_per_cell = 6,
              agent_contact = "geometric",
              agent_contact_parameter = c(0.6),
              agent_contact_contagion = 0.55)

historic_data <- run_simulation(model, save_path = "data/", ragg = TRUE)

save(historic_data, file = "data/historic_data.RData")

#png_files <- list.files(path = "img", pattern = "world_.*\\.png", full.names = TRUE)
#gif <- image_read(png_files) %>%
#  image_animate(fps = 2) %>%  # Adjust FPS (frames per second) for speed
#  image_write("world_evolution.gif")

