#' Run Markov chain Monte Carlo (MCMC) function
reslr_mcmc_fun <- function(data, data_grid,
                           n_iterations = 5000,
                           n_burnin = 1000,
                           n_thin = 4,
                           n_chains = 3,
                           spline_nseg_t = 5,#20,#15,#7,
                           spline_nseg_st = 6,#6,
                           spline_nseg_c = 20
                           ) {

  # Section: i.e. north south id set up-------
  section_tibble <- tibble::tibble(
    section = c( "South","North"),
    section_index = c(1, 2)
  )

  # Input Data -------------
  data <- data %>%
    dplyr::mutate(Age = Age / 1000, Age_err = Age_err / 1000) %>%
    mutate(section = factor(section,
      levels = c( "South","North")
    ))
  data <- left_join(data, section_tibble)


  # Input Data grid--------------
  data_grid <- data_grid %>%
    dplyr::mutate(Age = Age / 1000) %>%
    mutate(section = factor(section,
      levels = c( "South","North")
    ))
  data_grid <- left_join(data_grid, section_tibble)

  
  # Basis functions -----------------------------
  spline_basis_fun_list <- spline_basis_fun(
    data = data,
    data_grid = data_grid,
    spline_nseg_t = spline_nseg_t,
    spline_nseg_st = spline_nseg_st,
    spline_nseg_c = spline_nseg_c
  )
  
  # JAGS file for Common Component----------------------------------------------------
  common_jags_file <- "jags_models/common_model_ni_gam_decomp.jags"
  
  # JAGS data
  jags_data <- list(
    section = data$section_index,
    n_section = length(unique(data$section_index)),
    y = data$RSL,
    y_err = data$RSL_err,
    t = data$Age,
    site = as.factor(data$SiteName),
    n_sites = length(unique(data$SiteName)),
    n_obs = nrow(data),
    B_c = spline_basis_fun_list$B_c,
    n_knots_c = ncol(spline_basis_fun_list$B_c),
    linear_rate = data %>%
      dplyr::group_by(SiteName) %>%
      dplyr::slice(1) %>%
      dplyr::select(linear_rate) %>%
      dplyr::pull(),
    linear_rate_err = data %>%
      dplyr::group_by(SiteName) %>%
      dplyr::slice(1) %>%
      dplyr::select(linear_rate_err) %>%
      dplyr::pull()
  )
  
  # Parameters to save in JAGs
  jags_pars <- c(
    "mu_y",
    "sigma_y",
    "b_c",
    "b_g",
    "h_z_x",
    "g_h_z_x",
    "g_z_x",
    "c",
    "intercept",
    "sigma_beta_h",
    "sigma_beta_c",
    "sigmasq_all"
  )
  
  # Run JAGS------------------------
  common_model_run <- suppressWarnings(R2jags::jags.parallel(
    data = jags_data,
    parameters.to.save = jags_pars,
    model.file = common_jags_file,
    n.iter = n_iterations,
    n.burnin = n_burnin,
    n.thin = n_thin,
    n.chains = n_chains
  ))
  
  # Adding Noisy Input-------------------
  update_input_df <- add_noisy_input(
    data = data,
    data_grid = data_grid,
    model_run = common_model_run,
    jags_data = jags_data,
    spline_nseg_t = spline_nseg_t,
    spline_nseg_st = spline_nseg_st,
    spline_nseg_c = spline_nseg_c
  )
  data <- update_input_df$data
  data_grid <- update_input_df$data_grid

  # JAGS file----------------------------------------------------
  jags_file <- "jags_models/reg_model_ni_gam_decomp.jags"
  
  # JAGS data
  jags_data <- list(
    section = data$section_index,
    section_obs_id = data$section_obs_id,
    n_section = length(unique(data$section_index)),
    b_c_value = common_model_run$BUGSoutput$median$b_c,
    b_c_sd_value = common_model_run$BUGSoutput$sd$b_c,
    h_value = common_model_run$BUGSoutput$median$intercept,
    h_sd_value = common_model_run$BUGSoutput$sd$intercept,
    y = data$RSL,
    y_err = data$RSL_err,
    t = data$Age,
    site = as.factor(data$SiteName),
    n_sites = length(unique(data$SiteName)),
    n_obs = nrow(data),
    B_t = spline_basis_fun_list$B_t,
    B_t_deriv = spline_basis_fun_list$B_t_deriv,
    n_knots_t = ncol(spline_basis_fun_list$B_t),
    B_c = spline_basis_fun_list$B_c,
    n_knots_c = ncol(spline_basis_fun_list$B_c),
    linear_rate = data %>%
      dplyr::group_by(SiteName) %>%
      dplyr::slice(1) %>%
      dplyr::select(linear_rate) %>%
      dplyr::pull(),
    linear_rate_err = data %>%
      dplyr::group_by(SiteName) %>%
      dplyr::slice(1) %>%
      dplyr::select(linear_rate_err) %>%
      dplyr::pull()
  )

  # Parameters to save in JAGs
  jags_pars <- c(
    "mu_y",
    "sigma_y",
    "b_t",
    "h_z_x",
    "g_h_z_x",
    "g_z_x",
    "r",
    "c",
    "b_c",
    "intercept",
    "sigma_beta_h",
    "sigma_beta_r",
    "b_g",
    "sigmasq_all",
    "sigma_beta_r_world"
  )

  # Run JAGS------------------------
  model_run <- suppressWarnings(R2jags::jags.parallel(
    data = jags_data,
    parameters.to.save = jags_pars,
    model.file = jags_file,
    n.iter = n_iterations,
    n.burnin = n_burnin,
    n.thin = n_thin,
    n.chains = n_chains
  ))
  
  # Include Noise-----------------------
  noisy_jags_file <- "jags_models/noisy_model_ni_gam_decomp.jags"

  # JAGS input data
  jags_data <- list(
    section = data$section_index,
    n_section = length(unique(data$section_index)),
    NI_var_term = data$NI_var_term,
    b_t_value = model_run$BUGSoutput$median$b_t,
    b_t_sd_value = model_run$BUGSoutput$sd$b_t,
    b_c_value = common_model_run$BUGSoutput$median$b_c,
    b_c_sd_value = common_model_run$BUGSoutput$sd$b_c,
    h_value = common_model_run$BUGSoutput$median$intercept,
    h_sd_value = common_model_run$BUGSoutput$sd$intercept,
    y = data$RSL,
    y_err = data$RSL_err,
    t = data$Age,
    t_pred = data_grid$Age,
    site = as.factor(data$SiteName),
    site_pred = factor(data_grid$SiteName, levels = levels(as.factor(data$SiteName))),
    n_sites = length(unique(data$SiteName)),
    n_sites_pred = length(unique(data_grid$SiteName)),
    section_pred = as.numeric(data_grid$section_index),
    n_section_pred = length(unique(data_grid$section)),
    n_pred = length(data_grid$Age),
    n_obs = nrow(data),
    B_c = spline_basis_fun_list$B_c,
    B_c_deriv = spline_basis_fun_list$B_c_deriv,
    B_c_pred = spline_basis_fun_list$B_c_pred,
    B_c_pred_deriv = spline_basis_fun_list$B_c_pred_deriv,
    n_knots_c = ncol(spline_basis_fun_list$B_c),
    B_t = spline_basis_fun_list$B_t,
    B_t_deriv = spline_basis_fun_list$B_t_deriv,
    B_t_pred = spline_basis_fun_list$B_t_pred,
    B_t_pred_deriv = spline_basis_fun_list$B_t_pred_deriv,
    n_knots_t = ncol(spline_basis_fun_list$B_t),
    B_st = spline_basis_fun_list$B_st,
    B_st_pred = spline_basis_fun_list$B_st_pred,
    B_st_deriv = spline_basis_fun_list$B_st_deriv,
    B_st_deriv_pred = spline_basis_fun_list$B_st_deriv_pred,
    n_knots_st = ncol(spline_basis_fun_list$B_st),
    linear_rate = data %>%
      dplyr::group_by(SiteName) %>%
      dplyr::slice(1) %>%
      dplyr::select(linear_rate) %>%
      dplyr::pull(),
    linear_rate_err = data %>%
      dplyr::group_by(SiteName) %>%
      dplyr::slice(1) %>%
      dplyr::select(linear_rate_err) %>%
      dplyr::pull()
  )

  # Parameters to save in JAGs
  jags_pars <- c(
    "mu_y",
    "mu_deriv",
    "mu_pred",
    "mu_pred_deriv",
    "sigma_y",
    "sigma_beta_h",
    "sigma_beta_r",
    "sigma_beta_l",
    "sigmasq_all",
    "b_t",
    "b_st",
    "b_g",
    "b_c",
    "intercept",
    "h_z_x",
    "h_z_x_pred",
    "g_h_z_x",
    "g_h_z_x_pred",
    "g_z_x",
    "g_z_x_deriv",
    "g_z_x_pred",
    "g_z_x_pred_deriv",
    "c",
    "diff_r",
    "diff_r_deriv",
    "diff_r_pred",
    "diff_r_pred_deriv",
    "c_deriv",
    "c_pred",
    "c_pred_deriv",
    "r",
    "r_deriv",
    "r_pred",
    "r_pred_deriv",
    "l",
    "l_deriv",
    "l_pred",
    "l_pred_deriv",
    "sigmasq_all_pred",
    "y_pred"
  )
  # Run JAGS--------------
  noisy_model_run_output <-
    suppressWarnings(R2jags::jags.parallel(
      data = jags_data,
      parameters.to.save = jags_pars,
      model.file = noisy_jags_file,
      n.iter = n_iterations,
      n.burnin = n_burnin,
      n.thin = n_thin,
      n.chains = n_chains
    ))
  # Convert back to 1000
  data <- data %>%
    dplyr::mutate(Age = Age * 1000, Age_err = Age_err * 1000)
  data_grid <- data_grid %>%
    dplyr::mutate(Age = Age * 1000)

  # Output with everything-------------
  jags_output <- list(
    noisy_model_run_output = noisy_model_run_output,
    #jags_data = jags_data,
    data = data,
    data_grid = data_grid
  )
  

  return(jags_output)
}
