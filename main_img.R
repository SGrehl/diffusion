run_simulation(model, save_path= "img/")

png_files <- list.files(path = "img", pattern = "world_.*\\.png", full.names = TRUE)
gif <- image_read(png_files) %>%
  image_animate(fps = 2) %>%  # Adjust FPS (frames per second) for speed
  image_write("world_evolution.gif")
