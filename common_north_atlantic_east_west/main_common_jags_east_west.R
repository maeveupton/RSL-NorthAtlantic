# New approach for RSL analysis of East and West Atlantic

# f(x,t) = c(t) + r_East(t)+ r_West(t) + gt + l(x,t) + Error
# Model 1: Fit the global/common trend in just time
# Model 2: Fit the regional component with fixed priors from previous run
# Model 3: Strong priors on both common, regional and estimate the non-linear local component


# Data from Andy's updated data
# Clear workspace
rm(list = ls())
#---------Set working directory--------------
setwd("/users/research/mupton/3. RSL_North_Atlantic/1.common_north_atlantic_analysis_east_west")
#----------Load packages--------------------
library(devtools)
#install_github("maeveupton/reslr", force = TRUE)
# install.packages("reslr", dependencies = TRUE, INSTALL_opts = '--no-lock')
# install.packages("reslr")
library(reslr)
library(tidyverse)
library(ggplot2)
# remotes::install_github("tidyverse/purrr")
library(geosphere) # distm
library("rnaturalearth") # n_states for map
library(ggtext) # for the map
library(ggrepel) # labels_repel
library(ggspatial) # annotation_scale
library(readr) # read csv file
library(R2jags)# jags.parallel
library(xtable) # latex table
# R functions-------
source("R/plot_data.R") # Plotting separate east & west data & map plot
source("R/internal_functions.R") # New basis functions for east & west data
source("R/reslr_mcmc_fun.R") # New basis functions for east & west data
source("R/reslr_load_fun.R")

# Read in Andy's updated data:
global_data <- read_csv("https://www.dropbox.com/s/kqv3o10dnnbx38w/CommonEra2023_new.csv?dl=1")
colnames(global_data) <- c(
  "Basin", "Region", "Site", "Reference", "Indicator",
  "Latitude", "Longitude", "RSL", "RSL_err_upr", "RSL_err_lwr", "Age", "Age_2_err_upr", "Age_2_err_lwr"
)
# Updating the column names:
global_df <- global_data %>%
  dplyr::mutate(
    RSL_err = (RSL_err_upr + RSL_err_lwr) / 2,
    Age_err = ((Age_2_err_upr / 2) + (Age_2_err_lwr / 2)) / 2
  ) %>%
  # Only looking at Common era
  dplyr::filter(Age >= 0)
SL_df <- global_df %>%
  dplyr::filter(Basin == "North Atlantic")

# Removing Spain as it is problem with linear rate
SL_df <- SL_df %>%
  dplyr::filter(!Region == "Norway") %>%
  dplyr::filter(!Site %in% c("Muskiz Estuary", "Plentzia Estuary")) %>%
  dplyr::select(!c("RSL_err_upr", "RSL_err_lwr", "Age_2_err_upr", "Age_2_err_lwr", "Indicator"))

# New Irish dataset:
ireland_data <- read.csv("https://www.dropbox.com/scl/fi/kqa8xyb49v7egd5270129/Ireland_data_2023.csv?rlkey=myt1k5o1kl7nsblm3i7c8h2t4&dl=1")
ireland_df <- ireland_data %>% 
  dplyr::filter(!Site %in% c("Bull Island", "Cromane"))
SL_df <- rbind(SL_df, ireland_df)
# Remove index points
SL_df <- SL_df %>%
  dplyr::group_by(Site) %>%
  dplyr::filter(dplyr::n() > 2) %>%
  # Combine Sand Point Russian & VC
  dplyr::mutate(Site = replace(Site, Site == "Sand Point Russian", "Sand Point")) %>%
  dplyr::mutate(
    Site = replace(Site, Site == "Sand Point VC", "Sand Point"),
    Site = replace(Site, Site == " Timoleague", "Timoleague")
  )  %>% 
  # Filter greenland as it has no tg:
  dplyr::filter(!Region == "Greenland") %>% 
  # Filter out timoleague to see if it is weird
  dplyr::filter(!Site %in% c("Timoleague"))



# We include tide gauges for each proxy site--------------------
north_at_reslr <- reslr_load_fun(
  data = SL_df,
  include_tide_gauge = TRUE,
  include_linear_rate = TRUE,
  #all_TG_1deg = TRUE,
  TG_minimum_dist_proxy = TRUE,
  # Including tg longer than 75 years in europe
  list_preferred_TGs = c(
    "DUBLIN",
    # Germany & Denmark extra 100 years data
    #"CUXHAVEN 2", "ESBJERG",
    # Amsterdam extra 100 & 75 years
    "VLISSINGEN", "MAASSLUIS", "HOEK VAN HOLLAND",
    "IJMUIDEN","DEN HELDER", "HARLINGEN","WEST-TERSCHELLING",
    # Belgium extra 100 & 75 years
    "NIEUWPOORT", "OOSTENDE","ZEEBRUGGE",
    # France extra 100 & 75 years --> not working 
    "DUNKERQUE", "BOULOGNE","CALAIS",
    "BREST", "LE HAVRE", "ST. NAZAIRE",
    "ST JEAN DE LUZ (SOCOA)","LA ROCHELLE-LA PALLICE",
    # UK extra 100 & 75 years--> site not working:"HOLYHEAD"
    #"LIVERPOOL (GLADSTONE DOCK)",
    "SHEERNESS", "TOWER PIER", "SOUTHEND",
    "NORTH SHIELDS", "NEWLYN",
    # Portugal & Spain extra 100 & 75 years --> sites not working: "LAGOS"
    "CASCAIS", "CADIZ II", "CEUTA","VIGO",#"LA CORUÑA I"#,--> problem site
    # # North American sites
    "CHARLESTON I",
    "FORT PULASKI",
    "DAYTONA BEACH",
    "WILMINGTON",
    "SPRINGMAID PIER",
    "NAPLES",
    "FORT MYERS",
    # 50 years
    #"DEAL",
    #"PORTSMOUTH",
    "DEVONPORT",
    "HEYSHAM",
    "MILLPORT",
    "LERWICK",
    "TORSHAVN",
    "DIEPPE",
    "ROSCOFF",
    "LE CONQUET",
    "CONCARNEAU",
    "LES SABLES D OLONNE",
    "PORT BLOC",
    "ARCACHON-EYRAC",
    "BOUCAU",
    "LA CORUÑA II",
    "LEIXOES",
    "GIBRALTAR"
  ),
  prediction_grid_res = 30
)
saveRDS(north_at_reslr, file = "reslr_inputs/full_dataset/north_at_reslr_TG_LR.rds")
north_at_reslr <- readRDS("reslr_inputs/full_dataset/north_at_reslr_TG_LR.rds")

# Setting up 2 sections for East & west coast of north atlantic coast
# West data: North america
west_data <- north_at_reslr$data %>%
  filter(Longitude < -40) %>%
  filter(Latitude < 53) %>%
  mutate(section = "West")
# East data: europe
east_data <- anti_join(north_at_reslr$data, west_data)
east_data <- east_data %>% mutate(section = "East")
data <- rbind(east_data, west_data)

# Plotting raw data
data_plot_west <- plot_data(west_data)
ggsave(data_plot_west, filename = "fig/full_dataset/data_plot_west.pdf", height = 6, width = 10)
data_plot_east <- plot_data(east_data)
ggsave(data_plot_east, filename = "fig/full_dataset/data_plot_east.pdf", height = 6, width = 10)
# Plot map of sites:
plot_map(
  data = north_at_reslr$data,
  save_name = "fig/full_dataset/NA_map.pdf"
)

# West data grid
west_data_grid <- north_at_reslr$data_grid %>%
  filter(Longitude < -40) %>%
  filter(Latitude < 53) %>%
  mutate(section = "West")
# East data grid
east_data_grid <- anti_join(north_at_reslr$data_grid, west_data_grid)
east_data_grid <- east_data_grid %>% mutate(section = "East")
data_grid <- rbind(east_data_grid, west_data_grid)

# # Latex table of the data sites
# tg_tb<- data %>% filter(data_type_id == "TideGaugeData")%>%
#   dplyr::select(SiteName,Longitude, Latitude) %>% unique() %>% arrange(SiteName)
# print(xtable(tg_tb),
#       include.rownames = FALSE, 
#       include.colnames = FALSE)

## Running the ni-gam using updated reslr code------------------
# global_output <-
#   reslr_mcmc_fun(
#     data = data,
#     data_grid = data_grid,
#     n_iterations = 8000,
#     n_burnin = 1000,
#     n_thin = 5,
#     n_chains = 2,
#     spline_nseg_t = 8,
#     spline_nseg_st = 6,
#     spline_nseg_c = 20
#   )
# saveRDS(global_output,file = "reslr_outputs/full_dataset/global_reslr_output.rds")
global_output <- readRDS("reslr_outputs/full_dataset/global_reslr_output.rds")

# Total model fit
data_grid <-  global_output$data_grid
data <- global_output$data
tot_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$mu_pred
tot_post_df_full <- data.frame(pred = colMeans(tot_pred_post_full),
                               lwr = apply(tot_pred_post_full,2,quantile, probs = 0.025),
                               upr = apply(tot_pred_post_full,2,quantile, probs = 0.975),
                               lwr_50 = apply(tot_pred_post_full,2,quantile, probs = 0.25),
                               upr_50 = apply(tot_pred_post_full,2,quantile, probs = 0.75),
                               data_grid)
tot_post_df <- tot_post_df_full %>% filter(data_type_id == "ProxyRecord")
data <- data %>% filter(data_type_id == "ProxyRecord")
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

ggsave(tot_plot,filename = "fig/full_dataset/tot_mod_fit.pdf", width = 10, height = 6)

# Rate of change for Common Component
tot_pred_rate_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$mu_pred_deriv
tot_rate_post_df_full <- data.frame(pred = colMeans(tot_pred_rate_post_full),
                                    lwr = apply(tot_pred_rate_post_full,2,quantile, probs = 0.025),
                                    upr = apply(tot_pred_rate_post_full,2,quantile, probs = 0.975),
                                    lwr_50 = apply(tot_pred_rate_post_full,2,quantile, probs = 0.25),
                                    upr_50 = apply(tot_pred_rate_post_full,2,quantile, probs = 0.75),
                                    data_grid)
tot_rate_post_df <- tot_rate_post_df_full %>% filter(data_type_id == "ProxyRecord")

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


ggsave(tot_rate_plot,filename = "fig/full_dataset/tot_mod_fit_rate.pdf", width = 10, height = 6)


# Common Component
c_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$c_pred
c_post_df <- data.frame(pred = colMeans(c_pred_post_full),
                        lwr = apply(c_pred_post_full,2,quantile, probs = 0.025),
                        upr = apply(c_pred_post_full,2,quantile, probs = 0.975),
                        lwr_50 = apply(c_pred_post_full,2,quantile, probs = 0.25),
                        upr_50 = apply(c_pred_post_full,2,quantile, probs = 0.75),
                        global_output$data_grid)

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
    y = "Relative Sea Level (m)",
    colour = ""
  )
ggsave(c_plot,filename = "fig/full_dataset/common.pdf", width = 10, height = 6)

# Rate of change for Common Component
c_pred_rate_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$c_pred_deriv
c_rate_post_df <- data.frame(pred = colMeans(c_pred_rate_post_full),
                             lwr = apply(c_pred_rate_post_full,2,quantile, probs = 0.025),
                             upr = apply(c_pred_rate_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(c_pred_rate_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(c_pred_rate_post_full,2,quantile, probs = 0.75),
                             global_output$data_grid)

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

ggsave(c_rate_plot,filename = "fig/full_dataset/common_rate.pdf", width = 10, height = 6)

# Non Lin Local Component
local_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$l_pred
l_post_df_full <- data.frame(pred = colMeans(local_pred_post_full),
                             lwr = apply(local_pred_post_full,2,quantile, probs = 0.025),
                             upr = apply(local_pred_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(local_pred_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(local_pred_post_full,2,quantile, probs = 0.75),
                             data_grid)
l_post_df <- l_post_df_full %>% filter(data_type_id == "ProxyRecord")
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


ggsave(l_plot,filename = "fig/full_dataset/non_lin_loc.pdf", width = 10, height = 6)

# Rate of change of Non Lin Local Component
local_pred_rate_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$l_pred_deriv
l_post_rate_df_full <- data.frame(pred = colMeans(local_pred_rate_post_full),
                                  lwr = apply(local_pred_rate_post_full,2,quantile, probs = 0.025),
                                  upr = apply(local_pred_rate_post_full,2,quantile, probs = 0.975),
                                  lwr_50 = apply(local_pred_rate_post_full,2,quantile, probs = 0.25),
                                  upr_50 = apply(local_pred_rate_post_full,2,quantile, probs = 0.75),
                                  data_grid)
l_post_rate_df <- l_post_rate_df_full %>% filter(data_type_id == "ProxyRecord")
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


ggsave(l_rate_plot,filename = "fig/full_dataset/non_lin_loc_rate.pdf", width = 10, height = 6)

# Lin Local Component
lin_gia_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$g_h_z_x_pred
l_gia_post_df_full <- data.frame(pred = colMeans(lin_gia_pred_post_full),
                                 lwr = apply(lin_gia_pred_post_full,2,quantile, probs = 0.025),
                                 upr = apply(lin_gia_pred_post_full,2,quantile, probs = 0.975),
                                 lwr_50 = apply(lin_gia_pred_post_full,2,quantile, probs = 0.25),
                                 upr_50 = apply(lin_gia_pred_post_full,2,quantile, probs = 0.75),
                                 data_grid)
l_gia_post_df <- l_gia_post_df_full %>% filter(data_type_id == "ProxyRecord")
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


ggsave(l_gia_plot,filename = "fig/full_dataset/lin_loc_gia_plot.pdf", width = 10, height = 6)



# Regional component
r_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$r_pred
r_post_df <- data.frame(pred = colMeans(r_pred_post_full),
                        lwr = apply(r_pred_post_full,2,quantile, probs = 0.025),
                        upr = apply(r_pred_post_full,2,quantile, probs = 0.975),
                        lwr_50 = apply(r_pred_post_full,2,quantile, probs = 0.25),
                        upr_50 = apply(r_pred_post_full,2,quantile, probs = 0.75),
                        global_output$data_grid)
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
      alpha = c(0.2), # , 0.4),
      size = 1
    )),
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1),
      shape = c(NA),
      size = 2
    ))
  ) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 7,face = "bold"),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c("deepskyblue2","blue4")
  )
ggsave(reg_plot,filename = "fig/full_dataset/reg_east_west_eastonly.pdf", width = 10, height = 6)


# Rate of change Regional component
r_rate_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$r_pred_deriv
r_rate_post_df <- data.frame(pred = colMeans(r_rate_pred_post_full),
                             lwr = apply(r_rate_pred_post_full,2,quantile, probs = 0.025),
                             upr = apply(r_rate_pred_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(r_rate_pred_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(r_rate_pred_post_full,2,quantile, probs = 0.75),
                             data_grid)
view(r_rate_post_df %>% filter(section == "East"))
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
    strip.text.x = ggplot2::element_text(size = 7,face = "bold"),
    strip.background = ggplot2::element_rect(fill = c("white"))
  ) +
  ggplot2::scale_fill_manual("",
                             values = c("deepskyblue2","blue4")
  )

ggsave(reg_rate_plot,filename = "fig/full_dataset/rate_reg_east_west.pdf", width = 10, height = 6)

# Difference Regional Rate component
diff_r_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$diff_r_pred
diff_r_post_df <- data.frame(pred = colMeans(diff_r_pred_post_full),
                             lwr = apply(diff_r_pred_post_full,2,quantile, probs = 0.025),
                             upr = apply(diff_r_pred_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(diff_r_pred_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(diff_r_pred_post_full,2,quantile, probs = 0.75),
                             data_grid)
diff_reg_plot <- ggplot(data = diff_r_post_df, aes(x = Age, y = pred))+
  geom_line(aes(colour = "mean"))+
  geom_ribbon(aes(ymin = lwr,ymax=upr,fill = "CI"),alpha = 0.2)+
  geom_ribbon(aes(ymin = lwr_50,ymax=upr_50,fill = "CI"),alpha = 0.3)+
  labs(x = "Year (CE)", y = "East - West (m)",colour = "",fill = "")+
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
  )+
  scale_x_continuous(minor_breaks = seq(0, 2023, 100))

ggsave(diff_reg_plot,filename = "fig/full_dataset/diff_reg_east_west.pdf", width = 10, height = 6)

write_csv(diff_r_post_df,"Upton_reg_diff_east_west.csv")

# # Difference Regional Rate component
diff_r_rate_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$diff_r_pred_deriv
diff_r_rate_post_df <- data.frame(pred = colMeans(diff_r_rate_pred_post_full),
                                  lwr = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.025),
                                  upr = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.975),
                                  lwr_50 = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.25),
                                  upr_50 = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.75),
                                  data_grid)
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
  )+
  scale_x_continuous(minor_breaks = seq(0, 2023, 100))

ggsave(diff_reg_rate_plot,filename = "fig/full_dataset/diff_reg_rate_east_west.pdf", width = 10, height = 6)
write_csv(diff_r_rate_post_df,"Upton_rate_for_reg_diff_east_west.csv")
