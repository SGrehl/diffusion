# Create the world
world_init <- function(size = 9,
                       agents_per_cell = 3) {
  # Create the grid
  cells <- expand.grid(x = 1:size, y = 1:size)
  
  # Assign agents to grid cells
  world <- cells %>%
    slice(rep(1:n(), each = agents_per_cell)) %>%
    mutate(id   = row_number(),
           innovation = ifelse(x == ceiling(size/2) & y == ceiling(size/2) & row_number() == ceiling((size*size*agents_per_cell)/2), 1, 0)  # Initialize innovation to 0
    )
  return(world)
}

# Generate a picture of this square world. Cells without agents are filled black, cells with agents are filled according to the share of agents that have innovation = 1.
# The colors are defined by the variable color. The first color is used when the share is 0, the second color is used when the share is 0.5 and the third color is used if the share is 1. 
# If the share is anything in between, the colors reflect this accordingly
world_print <- function(world, 
                        i = NA,
                        palette = "RdYlBu") {
  # Summarize innovation per cell
  world_summary <- world %>%
    group_by(x, y) %>%
    summarise(share_innovation = mean(innovation), .groups = 'drop')
  
  # Create the plot
  plot <- ggplot(world_summary, aes(x = x, y = y, fill = share_innovation)) +
    # Draw a black background for the entire grid
    geom_tile(data = expand.grid(x = seq(min(world$x), max(world$x)), y = seq(min(world$y), max(world$y))), aes(x = x, y = y), fill = "black") +
    # Overlay the tiles for non-empty cells
    geom_tile() +
    scale_x_continuous(breaks = 1:max(world_summary$x)) +
    scale_y_continuous(breaks = 1:max(world_summary$y)) +
    coord_fixed() +
    theme_minimal() +
    scale_fill_distiller(palette = palette, limits = c(0, 1)) +  # Fill color for tiles
    labs(title = paste0("World Innovation Map: Round ",i),
         x = "X Coordinate",
         y = "Y Coordinate")
  
  print(plot)
  flush.console()
}

world_step <- function(world, model) {
  # Select all agents with innovation == 1
  innovators <- world %>% filter(innovation == 1)
  
  # if everybody is already converted don't run the simulation step
  if (nrow(innovators) == nrow(world)) return(world)
  
  # Draw random distance from Weibull distribution if specified in the model
  if (model$agent_contact == "weibull") {
    max_distance <- function(){round(rweibull(1, shape = model$agent_contact_parameter[1], scale = model$agent_contact_parameter[2]))}
  } else if (model$agent_contact == "geometric") {
    max_distance <- function(){round(rgeom(1, prob = model$agent_contact_parameter[1]))}
  } else {
    stop("Unsupported agent_contact type")
  }
  
  get_a_cell_at_manhattan_distance <- function(x, y, d, grid_size) {
    expand.grid(dx = -d:d, dy = -d:d) %>%
      filter(abs(dx) + abs(dy) == d) %>%
      mutate(x = x + dx, y = y + dy) %>%
      filter(x >= 1 & x <= grid_size & y >= 1 & y <= grid_size) %>%
      select(x, y) %>%
      slice_sample(n = 1)
  }
  
  # Iterate through each innovator
  for (i in 1:nrow(innovators)) {
    innovator <- innovators[i, ]
    
    # select a random cell in this range
    cell <- get_a_cell_at_manhattan_distance(innovator$x, innovator$y, max_distance(), model$world_size)
    
    # select a random agent at this cell
    
    # check whether there are cells 
    if (nrow(cell)>0) {
      contact <- world %>%
        filter(x == cell$x, y == cell$y) %>% 
        select(id) %>%
        slice_sample(n = 1)
      
      # If there are candidates, randomly select one to gain innovation
      if (nrow(contact) > 0 && runif(1) <= model$agent_contact_contagion) {
        world <- world %>%
          mutate(innovation = ifelse(id == contact$id, 1, innovation))
      } 
    }
  }
  
  return(world)
}

run_simulation <- function(model, 
                           max_steps = 25,
                           graph = FALSE,
                           save_path = NA){
  set.seed(model$seed)
  worlds <- list()
  
  if (graph) x11()
  
  for (i in 1:max_steps){
    if (i==1) worlds[[i]] <- world_init(model$world_size, model$world_agents_per_cell)
    else    worlds[[i]] <- world_step_old(worlds[[i-1]], model)
    
    if (!is.na(save_path)) {
      png_filename <- paste0(save_path,file.path(sprintf("world_%02d.png", i)))
      png(png_filename, width = 800, height = 800)
      world_print(worlds[[i]], i)
      dev.off()
    } else if (graph) {
      Sys.sleep(1)
      world_print(worlds[[i]], i)
      Sys.sleep(0)  
    }
  } 

  return(worlds)
}

report_generate <- function(worlds){
  report <- tibble()
  
  for (i in seq_along(worlds)) {
    
    world <- worlds[[i]]
    
    share_global <- sum(world$innovation, na.rm = TRUE) / nrow(world)
    
    world_summary <- world %>%
      group_by(x, y) %>%
      summarise(share_innovation = mean(innovation), .groups = 'drop') %>%
      group_by(share_innovation) %>%
      summarise(cells = n(), .groups = 'drop') %>%
      mutate(step = i) %>% 
      pivot_wider(
        names_from = share_innovation, 
        values_from = cells,
        values_fill = 0  # Fill missing values with 0
      ) %>% 
      ungroup() %>% 
      mutate(share_global = share_global)
    
    report <- bind_rows(report, world_summary)
    
  }
  
  report <- report %>%
    mutate(across(everything(), ~ replace_na(., 0))) %>% 
    relocate(step, share_global)
  
  share_cols <- names(report)[!names(report) %in% c("share_global", "step")]  # Exclude original columns
  new_names <- paste0("s_", seq_along(share_cols))              # Generate s_1, s_2, ...
  names(report)[!names(report) %in% c("share_global", "step")] <- new_names
  
  return(report)
}
