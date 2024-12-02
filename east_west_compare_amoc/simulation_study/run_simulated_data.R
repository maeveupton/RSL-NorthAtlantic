devtools::load_all()
library(reslr)


sea_level_data$linear_rate <- 0
sea_level_data$linear_rate_err <- 0.0001


multi_16_sites <- reslr_load(
  data = sea_level_data,
  prediction_grid_res = 1)

plot(
  x = multi_16_sites,
  title = "Plot of the raw data",
  xlab = "Year (CE)",
  ylab = "Relative Sea Level (m)",
  plot_tide_gauges = FALSE,
  plot_caption = TRUE
)


res_ni_gam_decomp <- reslr_mcmc(
  input_data = multi_16_sites,
  model_type = "ni_gam_decomp",
  # Update these values
  n_iterations = 2000, # Number of iterations
  n_burnin = 1000, # Number of iterations to discard at the beginning
  n_thin = 2, # Reduces number of output samples to save memory and computation time
  n_chains = 2 # Number of Markov chains
)

print(res_ni_gam_decomp)
summary(res_ni_gam_decomp)

plot(res_ni_gam_decomp,
     plot_type = "model_fit_plot",
     plot_tide_gauge = FALSE
)

plot(res_ni_gam_decomp, plot_type = "regional_plot") +
  geom_line(data = data.frame(Time = time, SeaLevel = common_sea_level), 
            aes(x = Time, y = SeaLevel), 
            color = "black", linetype = "dashed", size = 1.2)


