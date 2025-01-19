library(tidyverse)

# Create the world
world_init <- function(size = 9,
                       agents_per_cell = 3) {
  # Create the grid
  cells <- expand.grid(x = 1:size, y = 1:size)
  
  # Assign agents to grid cells
  world <- cells %>%
    slice(rep(1:n(), each = agents_per_cell)) %>%
    mutate(id   = row_number(),
           cell_id = x + (y-1) * size,
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
}

get_a_cell_at_manhattan_distance <- function(x, y, d, world_size) {
  # fast
  if (d == 0) return(x+(y-1)*world_size)
  
  # we are sure there are non 
  if (d > (world_size-1)*2) return(-1)
  
  cells_at_distance <- do.call(
    rbind,
    lapply(0:d, function(k) {
      dist_rem <- d - k
      # Four possible points for each k
      rbind(
        c(x + k, y + dist_rem),
        c(x + k, y - dist_rem),
        c(x - k, y + dist_rem),
        c(x - k, y - dist_rem)
      )
    })
  )
  
  # Remove duplicates (especially when k=0 or dist_rem=0)
  cells_at_distance <- unique(cells_at_distance)
  
  # Filter out-of-bounds points
  in_bounds <- cells_at_distance[,1] >= 1 & cells_at_distance[,1] <= world_size &
               cells_at_distance[,2] >= 1 & cells_at_distance[,2] <= world_size
  cells_at_distance <- cells_at_distance[in_bounds, , drop = FALSE]
  
  # If no valid points remain (for example, near edges), return NA or error
  if (nrow(cells_at_distance) == 0) return(-1) 
  
  # Pick one cell at random
  chosen_cell <- cells_at_distance[sample(nrow(cells_at_distance), 1), ]
  
  return(chosen_cell[1]+(chosen_cell[2]-1)*world_size)
}

create_position_lookup <- function(world) {
  # Group by cell_id, store a list of all IDs in that cell
  # 'deframe()' makes a named list with 'cell_id' as the name
  world %>%
    group_by(cell_id) %>%
    summarize(agent_ids = list(id), .groups = "drop") %>%
    deframe()
}

world_step <- function(world, model, position_lookup) {
  
  # Select all agents with innovation == 1
  innovators <- world %>% filter(innovation == 1)
  
  # 2) Filter out innovators whose random_value is above the threshold
  innovators <- innovators %>% 
    mutate(random_value = runif(n())) %>%
    filter(!(innovation == 1 & random_value > model$agent_contact_contagion))
    
  # if there are no innovators or everybody is already converted don't run the simulation step
  if (nrow(innovators) == 0 || nrow(innovators) == nrow(world)) return(world)

  innovators <- innovators %>% 
    mutate(
      distances = model$draw_distance(nrow(innovators))
    ) 
  
  
  # instead of innovators %>% rowwise() %>% mutate(target_cell = get_a_cell_at_manhattan_distance(x, y, distances, model$world_size))
  innovators <- innovators %>%
    mutate(
      target_cell = mapply(
        function(x, y, d) get_a_cell_at_manhattan_distance(x, y, d, model$world_size),
        x, y, distances
      )
    )
  
  # sample for each innovators target_cell a random target
  innovators <- innovators %>%
    mutate(
      target_id = mapply(
        function(tc, my_id) {
          # If invalid cell, return -1
          if (tc == -1) return(-1)
          
          # Possible targets in that cell
          targets_id <- position_lookup[[tc]]
          
          # Exclude myself
          targets_id <- targets_id[targets_id != my_id]
          
          # If no one left in that cell, return -1
          if (length(targets_id) == 0) return(-1)
          
          # Otherwise pick one at random
          sample(targets_id, 1)
        },
        target_cell,  # from innovators
        id            # from innovators
      )
    )
  
  targets = innovators %>% filter(target_id != -1) %>% pull(target_id)
  
  if (length(targets) > 0) world$innovation[world$id %in% targets] <- 1

  return(world)
}

run_simulation <- function(model, 
                           max_steps = 25
                           ){
  set.seed(model$seed)
  worlds <- list()

  for (i in 1:max_steps){
    if (i==1) {
      # initialize the world
      worlds[[i]] <- world_init(model$world_size, model$world_agents_per_cell)
      
      # generate a lookup table for optimizing the simulation speed
      position_lookup <- create_position_lookup(worlds[[i]]) # since agent don't move is sufficient to do this only once at the start of the simulation
      
      # define the draw_distance() function, which probability function is used
      model$draw_distance <- switch(
        model$agent_contact,
        "weibull" = function(n = 1) {
          round(rweibull(n,
                         shape = model$agent_contact_parameter[1],
                         scale = model$agent_contact_parameter[2]))
        },
        "geometric" = function(n = 1) {
          round(rgeom(n, prob = model$agent_contact_parameter[1]))
        },
        stop("Unsupported agent_contact type")
      )
    }
    else worlds[[i]] <- world_step(worlds[[i-1]], model, position_lookup)
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
