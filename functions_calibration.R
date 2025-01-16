# Example measure of fit between two reports
measure_fit <- function(historic_report, sim_report) {
  # Join by step (inner join: only compare matching steps)
  df <- dplyr::inner_join(historic_report, sim_report,
                          by = "step",
                          suffix = c("_hist", "_sim"))
  
  # Identify numeric columns to compare (excluding "step")
  compare_cols <- setdiff(names(df), "step")
  compare_cols_hist <- grep("_hist$", compare_cols, value = TRUE)
  compare_cols_sim  <- grep("_sim$",  compare_cols, value = TRUE)
  
  # We'll assume the columns match in order, i.e. 
  # for every X_hist there is a corresponding X_sim
  # For safety, you could match them by sub-string if needed.
  
  # For each pair, compute squared differences
  diff <- c()
  for (col_hist in compare_cols_hist) {
    col_sim <- sub("_hist$", "_sim", col_hist)  # get the matching sim column
    if (col_sim %in% names(df)) {
      # numeric difference
      v_hist <- df[[col_hist]]
      v_sim  <- df[[col_sim]]
      diff <- c(diff, ((v_hist - v_sim)/((v_hist +v_sim)/2))^2)
      diff[is.nan(diff)] <- 0
    }
  }
  
  # Return mean difference (or sum, or RMSE, etc.)
  mean(diff, na.rm = TRUE)
}


# a more complex function that allows for different weights with respect to certain cols or steps
measure_fit_weighted <- function(historic_report, 
                                 sim_report, 
                                 weight_steps = NULL,    # e.g. c(1,2,3,4,5)
                                 weight_cols  = NULL    # e.g. c("share_global"=2, "s_1"=1, "s_2"=1, ...)
) {
  # Join by step (inner join: only compare matching steps)
  df <- dplyr::inner_join(
    historic_report, sim_report,
    by = "step", suffix = c("_hist", "_sim")
  )
  
  # Identify numeric columns to compare (excluding "step")
  compare_cols <- setdiff(names(df), "step")
  compare_cols_hist <- grep("_hist$", compare_cols, value = TRUE)
  compare_cols_sim  <- grep("_sim$",  compare_cols, value = TRUE)
  
  # We will store our weighted differences in a vector:
  all_diff <- numeric(0)  # empty numeric vector
  
  for (col_hist in compare_cols_hist) {
    # Find the matching simulated column
    col_sim <- sub("_hist$", "_sim", col_hist)
    if (!(col_sim %in% names(df))) next  # skip if missing
    
    # Extract the underlying variable name (without _hist / _sim)
    var_name <- sub("_hist$", "", col_hist)
    
    # Grab vectors from the data frame
    v_hist <- df[[col_hist]]
    v_sim  <- df[[col_sim]]
    
    # We'll compute the unweighted difference
    raw_diff <- ((v_hist - v_sim)/((v_hist + v_sim)/2))^2
    
    # Replace NaN with 0
    raw_diff[is.nan(raw_diff)] <- 0
    
    # ------------------------------------------------------------------
    # Column-based weighting
    # ------------------------------------------------------------------
    # If `weight_cols` is provided, it should be something like:
    #   weight_cols = c("share_global"=2, "s_1"=1, "s_2"=3, ...)
    # If var_name isn't in weight_cols, we'll assume weight = 1
    if (!is.null(weight_cols) && var_name %in% names(weight_cols)) {
      w_col <- weight_cols[[var_name]]
    } else {
      w_col <- 1
    }
    
    # ------------------------------------------------------------------
    # Step-based weighting
    # ------------------------------------------------------------------
    # If `weight_steps` is provided, e.g. weight_steps = c(1,2,3,4,5)
    # We assume the step values go from 1..N, so we can do:
    #   w_step[i] = weight_steps[df$step[i]]
    # If no match, default to 1.
    
    # Create a vector of step weights, same length as raw_diff
    if (!is.null(weight_steps)) {
      # Some caution: step might not start at 1. 
      # If your steps are strictly 1..max, we can do:
      w_steps <- weight_steps[df$step]
      # If step is not guaranteed to be an integer index, you can do 
      # a named lookup, e.g. weight_steps = c("1"=1, "2"=2, ...)
      # w_steps <- weight_steps[as.character(df$step)]
      # fallback for missing names
      w_steps[is.na(w_steps)] <- 1
    } else {
      w_steps <- rep(1, length(raw_diff))
    }
    
    # ------------------------------------------------------------------
    # Combine the weights: multiply raw_diff by (w_col * w_step)
    # ------------------------------------------------------------------
    weighted_diff <- raw_diff * w_col * w_steps
    
    # Accumulate in our all_diff vector
    all_diff <- c(all_diff, weighted_diff)
  }
  
  # Finally, compute your desired summary (mean, sum, RMSE, etc.)
  mean(all_diff, na.rm = TRUE)
}



calibrate_model <- function(
    model,
    param_range_contact,      
    param_range_contagion,
    historic_report, 
    n_sims = 3,
    max_steps = 25
) {
  # We'll store results in a data frame
  results <- tibble()
  
  # Loop over all combinations in the grid
  for (param_contact in param_range_contact) {
    for (param_contagion in param_range_contagion) {
      cat(param_contact, " and ", param_contagion , "\n")
      
      # Update the model with these parameters
      # *If agent_contact_parameter is a scalar, assign param_contact directly
      # *If it's a vector, you might need to do something else
      model$agent_contact_parameter <- param_contact
      model$agent_contact_contagion <- param_contagion
      
      fit_values <- numeric(n_sims)  # holds error score for each replicate
      
      # Repeat the simulation multiple times to average out randomness
      for (rep_i in seq_len(n_sims)) {
        # (Re-)set seed if you want each replicate to be identical or 
        # if you want them to differ, do not reset the seed here.
        # set.seed(...) # optional
        model$seed = model$seed + 1
        
        # Run the simulation
        sim_worlds <- run_simulation(model, max_steps)
        
        # Generate the simulated report
        sim_report <- report_generate(sim_worlds)
        
        # Compute error/fit relative to historic data
        fit_values[rep_i] <- measure_fit(historic_report, sim_report)
      }
      
      # Average fit score across replicates
      avg_fit <- mean(fit_values)
      
      # Store in a results table
      results <- results %>%
        bind_rows(tibble(
          agent_contact_parameter = param_contact,
          agent_contact_contagion = param_contagion,
          fit = avg_fit
        ))
    }
  }
  
  # Identify the best (lowest) fit score
  best_result <- results %>% filter(fit == min(fit))
  
  # Return the full table and the best row(s)
  list(
    best = best_result,
    results = results
  )
}
