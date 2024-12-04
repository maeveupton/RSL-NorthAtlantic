# TO DO:
# 1. 8vs8 50 50 split of data --> then check that the common stays the same
# 2. 5 vs 11 split of data --> check common stays same
# 3. 3 vs 13 split of data --> check common stays same
# 4. With real data check that using the full data set for the common and then doing the splits

# Clear workspace
rm(list = ls())
# Code to simulate data for simulation study
setwd("/Users/maeve.upton/Dropbox/PhD_Work_2019_2023/RSL_North_Atlantic_chp3/Codes_1Nov/simulation_study")
# Load required packages
if (!require(splines)) install.packages("splines")
if (!require(MASS)) install.packages("MASS")  # for mvrnorm function
if (!require(ggplot2)) install.packages("ggplot2")
library(splines)
library(MASS)
library(tidyverse)
library(ggplot2)
# R functions-------
source("R/plot_data.R") # Plotting separate east & west data & map plot
source("R/internal_functions.R") # New basis functions for east & west data
source("R/reslr_mcmc_fun.R") # New basis functions for east & west data
source("R/reslr_load_fun.R")

# Set parameters
set.seed(123)                    # For reproducibility
n <- 200                         # Original number of time points
time <- unique(floor(seq(0, 100, length.out = n)))  # Regular time variable over 100 years
num_sites <- 16
n<- length(time)# Number of sites

# Generate random Latitude and Longitude coordinates for each site
latitudes <- runif(num_sites, min = -10, max = 10)  # Example latitude range
longitudes <- runif(num_sites, min = 120, max = 130)  # Example longitude range
site_coords <- data.frame(Site = factor(1:num_sites), Latitude = latitudes, Longitude = longitudes)

# Define B-spline basis for the common trend
xr <- max(time)
xl <- min(time)
dx <- (xr-xl)/20
degree <- 3 
common_knots <- seq(xl-degree*dx, xr+degree*dx, by = dx) # Internal knots for the common spline
# Degree of the spline (cubic)
common_spline_basis <- bs(time, knots = common_knots, degree = degree, Boundary.knots = c(common_knots[1],common_knots[length(common_knots)]))

# Generate synthetic coefficients for the common spline basis function
common_coefficients <- rnorm(ncol(common_spline_basis), mean = 0, sd = 5)

# Calculate the common sea level component based on the common spline
common_sea_level <- as.vector(common_spline_basis %*% common_coefficients)

# Site-specific intercepts
site_intercepts <- rnorm(num_sites, mean = 0, sd = 3)

# # Region-specific splines
num_regions <- 2
region_splines <- list()
#region_specific_spline_knots <- seq(xl, xr, by = 30)  # Fewer knots for site-specific splines
#region_specific_spline_knots <- seq(xl, xr, by = 8)  # Fewer knots for site-specific splines
dx_r = 8
region_specific_spline_knots <- seq(xl-degree*dx_r, xr+degree*dx_r, by = dx_r)
for (region in 1:num_regions) {
  # Generate a B-spline basis for each site with site-specific knots
  region_splines[[region]] <- bs(time, knots = region_specific_spline_knots, degree = degree, intercept = TRUE)
}

# Generate coefficients for site-specific splines
region_specific_coefficients <- matrix(rnorm(num_regions * ncol(region_splines[[1]]), mean = 0, sd = 2),
                                       nrow = num_regions, ncol = ncol(region_splines[[1]]))

# Calculate the regional sea level components
region_sea_level1 <- as.vector(region_splines[[1]] %*% region_specific_coefficients[1,])
region_sea_level2 <- as.vector(region_splines[[2]] %*% region_specific_coefficients[2,])

# Space-time component
site_locations <- seq(1, num_sites, by = 1)           # Simple 1D location for each site
distance_matrix <- as.matrix(dist(site_locations))    # Pairwise distances
spatial_range <- 3                                    # Spatial correlation range
spatial_cov <- exp(-distance_matrix / spatial_range)  # Exponential decay spatial covariance

# Temporal correlation setup
temporal_range <- 10                                  # Temporal correlation range
time_diff_matrix <- as.matrix(dist(time))             # Pairwise time differences
temporal_cov <- exp(-time_diff_matrix / temporal_range)  # Exponential decay temporal covariance

# Full space-time covariance structure
space_time_cov <- kronecker(temporal_cov, spatial_cov)  # Kronecker product for space-time covariance

# Generate space-time effects as multivariate normal across space and time
space_time_effects <- mvrnorm(1, mu = rep(0, num_sites * n), Sigma = space_time_cov)
space_time_effects <- matrix(space_time_effects, nrow = n, ncol = num_sites)

# Create uneven time points for each site
## fix so that you never get one data point at any site
uneven_time_indices <- lapply(1:num_sites, function(x) sort(sample(1:n, size = sample(40:50, 1))))
uneven_times <- lapply(uneven_time_indices, function(idx) time[idx])

# Initialize a data frame to store sea level data for each site
sea_level_data <- data.frame()
# Even split of sites----
#region_index <- c(rep(1,num_sites/2),rep(2,num_sites/2)) # 1st 8 sites to region 1 and 2nd 8 to region 2

# Un-Even split of sites----
region_index <- c(rep(1,num_sites - num_sites/4),rep(2,num_sites/4)) # 1st 12 sites to region 1 and 2nd 4 to region 2

# Generate sea level data for each site with uneven time steps
for (site in 1:num_sites) {
  selected_times <- uneven_times[[site]]
  idx <- uneven_time_indices[[site]]
  
  # Calculate site-specific spline component
  region_specific_spline <- as.vector(region_splines[[region_index[site]]][idx, ] %*% region_specific_coefficients[region_index[site], ])
  
  # Combine intercept, common sea level, site-specific spline, space-time effect, and random noise
  sea_level <- site_intercepts[site] +
    common_sea_level[idx] +
    region_specific_spline +
    space_time_effects[idx, site] +
    rnorm(length(idx), mean = 0, sd = 1.5)
  
  # Append to the data frame, including site coordinates
  sea_level_data <- rbind(sea_level_data, 
                          data.frame(Age = selected_times, 
                                     Site = factor(paste0("site",site)), 
                                     Latitude = site_coords$Latitude[site],
                                     Longitude = site_coords$Longitude[site],
                                     RSL = sea_level,
                                     RSL_err = rep(0.1,length(selected_times)),
                                     Age_err = rep(0.5, length(selected_times))))
}

# Simulated data for 16 sites with uneven time
sim_data <- ggplot(sea_level_data, aes(x = Age, y = RSL)) +
  geom_point(size = 1, alpha = 0.4) +
  #geom_line(data = data.frame(Time = time, SeaLevel = common_sea_level), 
  #          aes(x = Time, y = SeaLevel), 
  #         color = "black", linetype = "dashed", size = 1.2) +
  #geom_line(data = data.frame(Time = time, SeaLevel = region_sea_level1), 
  #          aes(x = Time, y = SeaLevel, colour = "region 1"), 
  #          linetype = "dashed", size = 1.2) +
  #geom_line(data = data.frame(Time = time, SeaLevel = region_sea_level2), 
  #          aes(x = Time, y = SeaLevel, colour = "region 2"), 
  #          linetype = "dashed", size = 1.2) +
  facet_wrap(~Site)+
  labs(title = "Simulated Sea Level Record",
       x = "Time",
       y = "Sea Level (synthetic data)") +
  theme_minimal() +
  scale_color_viridis_d(name = "Site") +
  theme(legend.position = "right")
sim_data
ggsave(sim_data,filename = "fig/simulate_data.pdf", width = 10, height = 7)


# Plot using ggplot2
library(ggplot2)

component_sim_plot <-ggplot(sea_level_data, aes(x = Age, y = RSL)) +
  # Points for sites
  geom_point(size = 1, alpha = 0.4, colour = "black") +
  
  # Line for the common component
  geom_line(data = data.frame(Time = time, SeaLevel = common_sea_level), 
            aes(x = Time, y = SeaLevel, colour = "Common"), 
            linetype = "dashed", size = 1.2) +
  
  # Line for Region 1
  geom_line(data = data.frame(Time = time, SeaLevel = region_sea_level1), 
            aes(x = Time, y = SeaLevel, colour = "Region 1"), 
            linetype = "dashed", size = 1.2) +
  
  # Line for Region 2
  geom_line(data = data.frame(Time = time, SeaLevel = region_sea_level2), 
            aes(x = Time, y = SeaLevel, colour = "Region 2"), 
            linetype = "dashed", size = 1.2) +
  
  # Title and axis labels
  labs(title = "Simulated Sea Level Record with Space-Time Component (Uneven Time Steps)",
       x = "Time",
       y = "Sea Level (synthetic data)") +
  
  # Minimal theme
  theme_minimal() +
  
  # Custom legend
  scale_color_manual(
    name = "Components",
    values = c("Common" = "red", 
               "Region 1" = "blue", 
               "Region 2" = "purple")
  ) +
  
  # Adjust legend position
  theme(legend.position = "right")
component_sim_plot
ggsave(component_sim_plot,filename = "fig/component_simulate_data.pdf", width = 10, height = 7)


#devtools::load_all()
library(reslr)


sea_level_data$linear_rate <- 0
sea_level_data$linear_rate_err <- 0.0001
sea_level_data$Region <- "Testing"
sea_level_data$Age_err <- 0.001
# Input data setting up 
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

# Testing different data splits , ie subsampling different number of sites either side
input_data <- multi_16_sites$data 
input_data_grid <- multi_16_sites$data_grid
# # 50:50 split
# input_data_split <- input_data %>% 
#   mutate(section = ifelse(Longitude >= 125.9, "East","West"))
# input_data_grid_split <- input_data_grid %>% 
#   mutate(section = ifelse(Longitude >= 125.9, "East","West"))

# 12 and 4 split
input_data_split <- input_data %>%
  mutate(section = ifelse(Longitude >= 123, "East","West"))
input_data_grid_split <- input_data_grid %>%
  mutate(section = ifelse(Longitude >= 123, "East","West"))

## Running the extended ni-gam using updated reslr code------------------
res_ni_gam_decomp <-
  reslr_mcmc_fun(
    data = input_data_split,
    data_grid = input_data_grid_split,
    n_iterations = 6000,
    n_burnin = 1000,
    n_thin = 5,
    n_chains = 2,
    spline_nseg_t = 8,
    spline_nseg_st = 6,
    spline_nseg_c = 20
  )

#saveRDS(res_ni_gam_decomp,file = "reslr_outputs/full_5050split.rds")
saveRDS(res_ni_gam_decomp,file = "reslr_outputs/full_12_4split.rds")
res_ni_gam_decomp <- readRDS("reslr_outputs/full_12_4split.rds")
print(res_ni_gam_decomp)
summary(res_ni_gam_decomp)

# Plot outputs--------
# Total model fit
tot_pred_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$mu_pred

tot_post_df_full <- data.frame(pred = colMeans(tot_pred_post_full),
                               lwr = apply(tot_pred_post_full,2,quantile, probs = 0.025),
                               upr = apply(tot_pred_post_full,2,quantile, probs = 0.975),
                               lwr_50 = apply(tot_pred_post_full,2,quantile, probs = 0.25),
                               upr_50 = apply(tot_pred_post_full,2,quantile, probs = 0.75),
                               res_ni_gam_decomp$data_grid)

tot_post_df <- tot_post_df_full 
data <- res_ni_gam_decomp$data 
tot_plot <-
  ggplot()+
  ggplot2::geom_rect(data = data, ggplot2::aes(
    xmin = Age - Age_err, xmax = Age + Age_err,
    ymin = RSL - RSL_err, ymax = RSL + RSL_err,
    fill = "Uncertainty",
  ), alpha = 0.7) +
  ggplot2::geom_point(
    data = data,
    ggplot2::aes(y = RSL, x = Age, colour = "black"), size = 0.3
  ) +
  geom_line(data = tot_post_df, aes(x = Age, y = pred,colour = "mean"))+
  geom_ribbon(data = tot_post_df, aes(x = Age,ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(data = tot_post_df,aes(x = Age,ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c(
                               "Uncertainty" = ggplot2::alpha("grey", 0.3),
                               "CI" = ggplot2::alpha("purple3", 0.2)
                             ),
                             labels = c(
                               CI = "50% & 95% Credible Interval",
                               expression(paste("1-sigma Error")))
  )+
  ggplot2::scale_colour_manual("",
                               values = c(
                                 "black" = "black",
                                 "mean" = "purple3"
                               ),
                               labels = c("Data", "Posterior Fit")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(legend.box = "horizontal",
                 legend.position = "bottom",
                 axis.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10)) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Relative Sea Level (m)",
    colour = ""
  )+
  facet_wrap(~SiteName)
tot_plot
ggsave(tot_plot,filename = "fig/tot_mod_fit.pdf", width = 10, height = 6)

# Rate of change for Common Component
tot_pred_rate_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$mu_pred_deriv
tot_rate_post_df_full <- data.frame(pred = colMeans(tot_pred_rate_post_full),
                                    lwr = apply(tot_pred_rate_post_full,2,quantile, probs = 0.025),
                                    upr = apply(tot_pred_rate_post_full,2,quantile, probs = 0.975),
                                    lwr_50 = apply(tot_pred_rate_post_full,2,quantile, probs = 0.25),
                                    upr_50 = apply(tot_pred_rate_post_full,2,quantile, probs = 0.75),
                                    res_ni_gam_decomp$data_grid)

tot_rate_post_df <- tot_rate_post_df_full 

tot_rate_plot <- ggplot(data = tot_rate_post_df, aes(x = Age, y = pred))+
  geom_line(aes(colour = "mean"))+
  geom_ribbon(aes(ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(aes(ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  ggplot2::scale_fill_manual("",
                             values = c(
                               "CI" = ggplot2::alpha("purple3", 0.2)
                             ),
                             labels = c(
                               CI = "50% & 95% Credible Interval")
  )+
  ggplot2::scale_colour_manual("",
                               values = c(
                                 "mean" = "purple3"
                               ),
                               labels = c("Posterior Fit")
  ) +
  theme_bw()+
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(legend.box = "horizontal",
                 legend.position = "bottom",
                 axis.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10),
                 strip.text.x = ggplot2::element_text(size = 7),
                 strip.background = ggplot2::element_rect(fill = c("white"))) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Rate of Change (mm/year)",
    colour = ""
  )+
  facet_wrap(~SiteName)


ggsave(tot_rate_plot,filename = "fig/tot_mod_fit_rate.pdf", width = 10, height = 6)


# Common Component
c_pred_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$c_pred
c_post_df <- data.frame(pred = colMeans(c_pred_post_full),
                        lwr = apply(c_pred_post_full,2,quantile, probs = 0.025),
                        upr = apply(c_pred_post_full,2,quantile, probs = 0.975),
                        lwr_50 = apply(c_pred_post_full,2,quantile, probs = 0.25),
                        upr_50 = apply(c_pred_post_full,2,quantile, probs = 0.75),
                        res_ni_gam_decomp$data_grid)


c_plot <-
  ggplot(data = c_post_df, aes(x = Age, y = pred))+
  geom_line(aes(colour = "mean"))+
  geom_ribbon(aes(ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(aes(ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c(
                               "CI" = ggplot2::alpha("brown2", 0.2)
                             ),
                             labels = c(
                               CI = "50% & 95% Credible Interval")
  )+
  ggplot2::scale_colour_manual("",
                               values = c(
                                 "mean" = "brown2"
                               ),
                               labels = c("Posterior Fit")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(legend.box = "horizontal",
                 legend.position = "bottom",
                 axis.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10)) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Sea Level (m)",
    colour = ""
  )
ggsave(c_plot,filename = "fig/common.pdf", width = 10, height = 6)
# COmparing output with simulation
common_res_plot <- c_plot +
  geom_line(data = data.frame(Time = time, SeaLevel = common_sea_level),
            aes(x = Time, y = SeaLevel),
            color = "red", linetype = "dashed", size = 1.2)
ggsave(common_res_plot,filename = "common_res_plot_compare.pdf", width = 10, height = 7)


# Rate of change for Common Component
c_pred_rate_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$c_pred_deriv
c_rate_post_df <- data.frame(pred = colMeans(c_pred_rate_post_full),
                             lwr = apply(c_pred_rate_post_full,2,quantile, probs = 0.025),
                             upr = apply(c_pred_rate_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(c_pred_rate_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(c_pred_rate_post_full,2,quantile, probs = 0.75),
                             res_ni_gam_decomp$data_grid)

c_rate_plot <- ggplot(data = c_rate_post_df, aes(x = Age, y = pred))+
  geom_line(aes(colour = "mean"))+
  geom_ribbon(aes(ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(aes(ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  ggplot2::scale_fill_manual("",
                             values = c(
                               "CI" = ggplot2::alpha("brown2", 0.2)
                             ),
                             labels = c(
                               CI = "50% & 95% Credible Interval")
  )+
  ggplot2::scale_colour_manual("",
                               values = c(
                                 "mean" = "brown2"
                               ),
                               labels = c("Posterior Fit")
  ) +
  theme_bw()+
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(legend.box = "horizontal",
                 legend.position = "bottom",
                 axis.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10)) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Rate of Change (mm/year)",
    colour = ""
  )

ggsave(c_rate_plot,filename = "fig/common_rate.pdf", width = 10, height = 6)

# Non Lin Local Component
local_pred_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$l_pred
l_post_df_full <- data.frame(pred = colMeans(local_pred_post_full),
                             lwr = apply(local_pred_post_full,2,quantile, probs = 0.025),
                             upr = apply(local_pred_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(local_pred_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(local_pred_post_full,2,quantile, probs = 0.75),
                             res_ni_gam_decomp$data_grid)
#global_output$data)
l_post_df <- l_post_df_full 
l_plot <- ggplot()+
  geom_line(data = l_post_df, aes(x = Age, y = pred,colour = "mean"))+
  geom_ribbon(data = l_post_df,aes(x = Age,ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(data = l_post_df,aes(x = Age,ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c(
                               "CI" = ggplot2::alpha("#ad4c14", 0.2)
                             ),
                             labels = c(
                               CI = "50% & 95% Credible Interval")
  )+
  ggplot2::scale_colour_manual("",
                               values = c(
                                 "mean" = "#ad4c14"
                               ),
                               labels = c("Posterior Fit")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(legend.box = "horizontal",
                 legend.position = "bottom",
                 axis.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10)) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Sea Level (m)",
    colour = ""
  )+
  facet_wrap(~SiteName)


ggsave(l_plot,filename = "fig/non_lin_loc.pdf", width = 10, height = 6)

# Rate of change of Non Lin Local Component
local_pred_rate_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$l_pred_deriv
l_post_rate_df_full <- data.frame(pred = colMeans(local_pred_rate_post_full),
                                  lwr = apply(local_pred_rate_post_full,2,quantile, probs = 0.025),
                                  upr = apply(local_pred_rate_post_full,2,quantile, probs = 0.975),
                                  lwr_50 = apply(local_pred_rate_post_full,2,quantile, probs = 0.25),
                                  upr_50 = apply(local_pred_rate_post_full,2,quantile, probs = 0.75),
                                  res_ni_gam_decomp$data_grid)
l_post_rate_df <- l_post_rate_df_full
l_rate_plot <- ggplot()+
  geom_line(data = l_post_rate_df, aes(x = Age, y = pred,colour = "mean"))+
  geom_ribbon(data = l_post_rate_df,aes(x = Age,ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(data = l_post_rate_df,aes(x = Age,ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  ggplot2::theme_bw() +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c(
                               "CI" = ggplot2::alpha("#ad4c14", 0.2)
                             ),
                             labels = c(
                               CI = "50% & 95% Credible Interval")
  )+
  ggplot2::scale_colour_manual("",
                               values = c(
                                 "mean" = "#ad4c14"
                               ),
                               labels = c("Posterior Fit")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(legend.box = "horizontal",
                 legend.position = "bottom",
                 axis.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10)) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Rate of Change (mm/year)",
    colour = ""
  )+
  facet_wrap(~SiteName)


ggsave(l_rate_plot,filename = "fig/non_lin_loc_rate.pdf", width = 10, height = 6)

# Lin Local Component
lin_gia_pred_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$g_h_z_x_pred
l_gia_post_df_full <- data.frame(pred = colMeans(lin_gia_pred_post_full),
                                 lwr = apply(lin_gia_pred_post_full,2,quantile, probs = 0.025),
                                 upr = apply(lin_gia_pred_post_full,2,quantile, probs = 0.975),
                                 lwr_50 = apply(lin_gia_pred_post_full,2,quantile, probs = 0.25),
                                 upr_50 = apply(lin_gia_pred_post_full,2,quantile, probs = 0.75),
                                 res_ni_gam_decomp$data_grid)

l_gia_post_df <- l_gia_post_df_full 
l_gia_plot <- ggplot()+
  geom_line(data = l_gia_post_df, aes(x = Age, y = pred,colour = "mean"))+
  geom_ribbon(data = l_gia_post_df,aes(x = Age,ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(data = l_gia_post_df,aes(x = Age,ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c(
                               "CI" = ggplot2::alpha("#5bac06", 0.2)
                             ),
                             labels = c(
                               CI = "50% & 95% Credible Interval")
  )+
  ggplot2::scale_colour_manual("",
                               values = c(
                                 "mean" = "#5bac06"
                               ),
                               labels = c("Posterior Fit")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(legend.box = "horizontal",
                 legend.position = "bottom",
                 axis.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10)) +
  ggplot2::labs(
    x = "Year (CE)",
    y = "Sea Level (m)",
    colour = ""
  )+
  facet_wrap(~SiteName)


ggsave(l_gia_plot,filename = "fig/lin_loc_gia_plot.pdf", width = 10, height = 6)



# Regional component
r_pred_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$r_pred
r_post_df <- data.frame(pred = colMeans(r_pred_post_full),
                        lwr = apply(r_pred_post_full,2,quantile, probs = 0.025),
                        upr = apply(r_pred_post_full,2,quantile, probs = 0.975),
                        lwr_50 = apply(r_pred_post_full,2,quantile, probs = 0.25),
                        upr_50 = apply(r_pred_post_full,2,quantile, probs = 0.75),
                        res_ni_gam_decomp$data_grid)

# # Correlation between reg components using the spline coefficient:
# b_r_south <-  colMeans(global_output$noisy_model_run_output$BUGSoutput$sims.list$b_t[,1,])
# b_r_north <-  colMeans(global_output$noisy_model_run_output$BUGSoutput$sims.list$b_t[,2,])
# cor_sn <-cor(b_r_south,b_r_north)
# #corrplot(cor_sn,method = 'number')

# Plot
reg_plot <- ggplot(data = r_post_df, aes(x = Age, y = pred))+
  geom_line(aes(colour = section))+
  geom_ribbon(aes(ymin = lwr,ymax=upr,fill = section),alpha = 0.2)+
  geom_ribbon(aes(ymin = lwr_50,ymax=upr_50,fill = section),alpha = 0.3)+
  labs(x = "Year (CE)", y = "Sea Level (m)",colour = "",fill = "")+
  theme_bw()+
  ggplot2::scale_colour_manual("",
                               values = c("deepskyblue2","blue4")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c("deepskyblue2","blue4")
  )
ggsave(reg_plot,filename = "fig/reg_eastwest.pdf", width = 10, height = 6)
# Compare with simulated regional component
compare_region_plot <- reg_plot + 
# Line for Region 1
geom_line(data = data.frame(Time = time, SeaLevel = region_sea_level1), 
          aes(x = Time, y = SeaLevel, colour = "Region 1"), 
          linetype = "dashed", size = 1.2) +
  
  # Line for Region 2
  geom_line(data = data.frame(Time = time, SeaLevel = region_sea_level2), 
            aes(x = Time, y = SeaLevel, colour = "Region 2"), 
            linetype = "dashed", size = 1.2) +
  # Custom legend
  scale_color_manual(
    name = "Components",
    values = c(
               "Region 1" = "blue", 
               "Region 2" = "purple")
  ) 
compare_region_plot
ggsave(compare_region_plot,filename = "fig/compare_reg_sim_eastwest.pdf", width = 10, height = 6)
# Rate of change Regional component
r_rate_pred_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$r_pred_deriv
r_rate_post_df <- data.frame(pred = colMeans(r_rate_pred_post_full),
                             lwr = apply(r_rate_pred_post_full,2,quantile, probs = 0.025),
                             upr = apply(r_rate_pred_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(r_rate_pred_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(r_rate_pred_post_full,2,quantile, probs = 0.75),
                             res_ni_gam_decomp$data_grid)
#global_output$data)
reg_rate_plot <- ggplot(data = r_rate_post_df, aes(x = Age, y = pred))+
  geom_line(aes(colour = section))+
  geom_ribbon(aes(ymin = lwr,ymax=upr,fill = section),alpha = 0.2)+
  geom_ribbon(aes(ymin = lwr_50,ymax=upr_50,fill = section),alpha = 0.3)+
  labs(x = "Year (CE)", y = "Rate of Change (mm/year)",colour = "",fill = "")+
  ggplot2::geom_hline(yintercept = 0) +
  theme_bw()+
  ggplot2::scale_colour_manual("",
                               values = c("deepskyblue2","blue4")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c("deepskyblue2","blue4")
  )

ggsave(reg_rate_plot,filename = "fig/rate_reg_eastwest.pdf", width = 10, height = 6)

# Difference Regional Rate component
diff_r_pred_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$diff_r_pred
diff_r_post_df <- data.frame(pred = colMeans(diff_r_pred_post_full),
                             lwr = apply(diff_r_pred_post_full,2,quantile, probs = 0.025),
                             upr = apply(diff_r_pred_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(diff_r_pred_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(diff_r_pred_post_full,2,quantile, probs = 0.75),
                             res_ni_gam_decomp$data_grid)
#global_output$data)
diff_reg_plot <- ggplot(data = diff_r_post_df, aes(x = Age, y = pred))+
  geom_line(aes(colour = "mean"))+
  geom_ribbon(aes(ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(aes(ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  labs(x = "Year (CE)", y = "South - North (m)",colour = "",fill = "")+
  theme_bw()+
  ggplot2::scale_colour_manual("",
                               values = c("mean"="darkmagenta"),
                               labels = "Posterior Fit"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white")),
    legend.box = "horizontal",
    legend.position = "bottom",
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 10)
  ) +
  ggplot2::scale_fill_manual("",
                             values= c("CI"="darkmagenta"),
                             labels = "50% & 95% Credible Interval"
  )

ggsave(diff_reg_plot,filename = "fig/diff_reg_eastwest.pdf", width = 10, height = 6)


# # Difference Regional Rate component
diff_r_rate_pred_post_full <- res_ni_gam_decomp$noisy_model_run_output$BUGSoutput$sims.list$diff_r_pred_deriv
diff_r_rate_post_df <- data.frame(pred = colMeans(diff_r_rate_pred_post_full),
                                  lwr = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.025),
                                  upr = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.975),
                                  lwr_50 = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.25),
                                  upr_50 = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.75),
                                  res_ni_gam_decomp$data_grid)
#global_output$data)
diff_reg_rate_plot <- ggplot(data = diff_r_rate_post_df, aes(x = Age, y = pred))+
  geom_line(aes(colour = "mean"))+
  ggplot2::geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(aes(ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  labs(x = "Year (CE)", y = "Rate of Change (mm/year)",colour = "",fill = "")+
  theme_bw()+
  ggplot2::scale_colour_manual("",
                               values = c("mean"="darkmagenta"),
                               labels = "Posterior Fit"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(
      alpha = c(0.3), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = c("white")),
    legend.box = "horizontal",
    legend.position = "bottom",
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 10)
  ) +
  ggplot2::scale_fill_manual("",
                             values= c("CI"="darkmagenta"),
                             labels = "50% & 95% Credible Interval"
  )

ggsave(diff_reg_rate_plot,filename = "fig/diff_reg_rate_eastwest.pdf", width = 10, height = 6)



