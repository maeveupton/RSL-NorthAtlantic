# Code to simulate data for simulation study
setwd("/Users/maeve.upton/Dropbox/PhD_Work_2019_2023/RSL_North_Atlantic_chp3/Codes_1Nov/simulation_study")
# Load required packages
if (!require(splines)) install.packages("splines")
if (!require(MASS)) install.packages("MASS")  # for mvrnorm function
if (!require(ggplot2)) install.packages("ggplot2")
library(splines)
library(MASS)
library(ggplot2)

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
region_specific_spline_knots <- seq(xl, xr, by = 30)  # Fewer knots for site-specific splines
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
region_index <- c(rep(1,num_sites/2),rep(2,num_sites/2)) # 1st 5 sites to region 1 and 2nd 5 to region 2
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
ggsave(sim_data,filename = "simulate_data.pdf", width = 10, height = 7)


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

ggsave(component_sim_plot,filename = "component_simulate_data.pdf", width = 10, height = 7)


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

common_res_plot <- plot(res_ni_gam_decomp, plot_type = "regional_plot") +
  geom_line(data = data.frame(Time = time, SeaLevel = common_sea_level), 
            aes(x = Time, y = SeaLevel), 
            color = "red", linetype = "dashed", size = 1.2)
ggsave(common_res_plot,filename = "common_res_plot.pdf", width = 10, height = 7)

