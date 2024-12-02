# New approach for RSL analysis of North South of Cape Hatteras

# f(x,t) = c(t) + r_North(t)+ r_South(t) + gt + l(x,t) + Error
# Model 1: Fit the global/common trend in just time
# Model 2: Fit the regional component with fixed priors from previous run
# Model 3: Strong priors on both common, regional and estimate the non-linear local component


# Data from Andy's updated data
# Clear workspace
rm(list = ls())
#---------Set working directory--------------
setwd("/Users/maeve.upton/Desktop/GitHub/RSL-NorthAtlantic/common_north_atlantic_south_north")
#----------Load packages--------------------
library(devtools)
#install_github("maeveupton/reslr", force = TRUE)
#install.packages("reslr", dependencies = TRUE, INSTALL_opts = '--no-lock')
#install.packages("reslr")
library(reslr)
library(tidyverse)
library(ggplot2)
#remotes::install_github("tidyverse/purrr")
library(geosphere) #distm
library("rnaturalearth")#n_states for map
library(ggtext)# for the map
library(ggrepel)# labels_repel
library(ggspatial)# annotation_scale
library(readr)# read csv file
# R functions-------
source("R/plot_data.R")# Plotting separate north & south data & map plot
source("R/internal_functions.R")# New basis functions for north & south data
source("R/reslr_mcmc_fun.R")# New basis functions for north & south data
source("R/reslr_load_fun.R")# New basis functions for north & south data


# Read in Andy's updated data:
global_data <- read_csv("https://www.dropbox.com/s/kqv3o10dnnbx38w/CommonEra2023_new.csv?dl=1")
colnames(global_data)<- c("Basin","Region","Site","Reference","Indicator",
                          "Latitude","Longitude","RSL","RSL_err_upr","RSL_err_lwr","Age","Age_2_err_upr","Age_2_err_lwr")
# Updating the column names:
global_df <- global_data %>% 
  mutate(RSL_err = (RSL_err_upr + RSL_err_lwr)/2,
         Age_err = ((Age_2_err_upr/2) + (Age_2_err_lwr/2))/2) %>% 
  # Only looking at Common era
  filter(Age >= 0) 
SL_df <- global_df %>%
  filter(Basin == "North Atlantic") 

# Removing Spain as it is problem with linear rate
SL_df <- SL_df %>%
  dplyr::filter(!Region == "Norway") %>% 
  dplyr::filter(!Site %in% c("Muskiz Estuary","Plentzia Estuary")) 



# Remove index points
SL_df <- SL_df %>%
  dplyr::group_by(Site)%>%
  dplyr::filter(dplyr::n()>2) %>%
  # Combine Sand Point Russian & VC
  mutate(Site = replace(Site, Site == "Sand Point Russian", "Sand Point")) %>%
  mutate(Site = replace(Site, Site == "Sand Point VC", "Sand Point"))

# North American sites only------------
SL_df <- SL_df %>% filter(Region %in% c("Connecticut",
                                        "Florida",
                                        "Massachusetts",
                                        "North Carolina","New Jersey",
                                        "Magdelen Islands","Newfoundland",
                                        "Quebec","Nova Scotia","Maine","New York","Rhode Island"))

# We include tide gauges for each proxy site--------------------
north_at_reslr <- reslr_load_fun(data = SL_df,
                                 include_tide_gauge = TRUE,
                                 include_linear_rate = TRUE,
                                 list_preferred_TGs =
                                   c(
                                     #"LAKE WORTH PIER",
                                     "CHARLESTON I",
                                     "FORT PULASKI",
                                     "DAYTONA BEACH",
                                     "WILMINGTON",
                                     "SPRINGMAID PIER",
                                     "NAPLES",
                                     "FORT MYERS"),
                                 TG_minimum_dist_proxy = TRUE,
                                 prediction_grid_res = 30)
saveRDS(north_at_reslr,file = "reslr_inputs/north_at_reslr_TG_LR.rds")
north_at_reslr <- readRDS("reslr_inputs/north_at_reslr_TG_LR.rds")

# Setting up 2 sections for north and south cape Hatteras:
data <- north_at_reslr$data %>%
  mutate(section = ifelse(Latitude > 35.3,"North","South"),
         SiteName = as.character(SiteName)) %>% 
  mutate(SiteName = factor(SiteName))

data_grid <- north_at_reslr$data_grid %>%
  mutate(section = ifelse(Latitude > 35.3,"North","South"),
         SiteName = as.character(SiteName)) %>% 
  mutate(SiteName = factor(SiteName))

# north data
north_data <- data %>% 
  filter(section == "North")
# south data 
south_data <- anti_join(data,north_data)

# Plotting raw data
data_plot_south <- plot_data(south_data)
ggsave(data_plot_south,filename = "fig/data_plot_south.pdf",height = 6, width= 10)
data_plot_north <- plot_data(north_data)
ggsave(data_plot_north,filename = "fig/data_plot_north.pdf",height = 6, width= 10)

# Plot Map of data
plot_map(data)

# north data grid
north_data_grid <- data_grid %>% 
  filter(section == "North")
# south data 
south_data_grid <- anti_join(data_grid,north_data_grid)

# Running the ni-gam using updated reslr code------------------
# global_output <-
#   reslr_mcmc_fun(
#     data = data,
#     data_grid = data_grid,
#     n_iterations = 8000,
#     n_burnin = 2000,
#     n_thin = 4,
#     n_chains = 2,
#     spline_nseg_t = 4,#6,#NULL,#5,#6,# 3 = NULL,#10,
#     spline_nseg_st = 6,
#     spline_nseg_c = 20
#   )
# saveRDS(global_output,file = "reslr_outputs/global_reslr_output.rds")

# Outputs
global_output <- readRDS("reslr_outputs/global_reslr_output.rds")


# Plot outputs--------
# Total model fit
tot_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$mu_pred

tot_post_df_full <- data.frame(pred = colMeans(tot_pred_post_full),
                               lwr = apply(tot_pred_post_full,2,quantile, probs = 0.025),
                               upr = apply(tot_pred_post_full,2,quantile, probs = 0.975),
                               lwr_50 = apply(tot_pred_post_full,2,quantile, probs = 0.25),
                               upr_50 = apply(tot_pred_post_full,2,quantile, probs = 0.75),
                               global_output$data_grid)

tot_post_df <- tot_post_df_full %>% 
  filter(data_type_id == "ProxyRecord") 
data <- global_output$data %>% 
  filter(data_type_id == "ProxyRecord")
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
                                    global_output$data_grid)

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
    y = "Sea Level (m)",
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
                             global_output$data_grid)
#global_output$data)
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
                                  global_output$data_grid)
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
                                 global_output$data_grid)

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

# Correlation between reg components using the spline coefficient:
b_r_south <-  colMeans(global_output$noisy_model_run_output$BUGSoutput$sims.list$b_t[,1,])
b_r_north <-  colMeans(global_output$noisy_model_run_output$BUGSoutput$sims.list$b_t[,2,])
cor_sn <-cor(b_r_south,b_r_north)
#corrplot(cor_sn,method = 'number')
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
ggsave(reg_plot,filename = "fig/full_dataset/reg_ns.pdf", width = 10, height = 6)


# Rate of change Regional component
r_rate_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$r_pred_deriv
r_rate_post_df <- data.frame(pred = colMeans(r_rate_pred_post_full),
                             lwr = apply(r_rate_pred_post_full,2,quantile, probs = 0.025),
                             upr = apply(r_rate_pred_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(r_rate_pred_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(r_rate_pred_post_full,2,quantile, probs = 0.75),
                             global_output$data_grid)
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

ggsave(reg_rate_plot,filename = "fig/full_dataset/rate_reg_ns.pdf", width = 10, height = 6)

# Difference Regional Rate component
diff_r_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$diff_r_pred
diff_r_post_df <- data.frame(pred = colMeans(diff_r_pred_post_full),
                             lwr = apply(diff_r_pred_post_full,2,quantile, probs = 0.025),
                             upr = apply(diff_r_pred_post_full,2,quantile, probs = 0.975),
                             lwr_50 = apply(diff_r_pred_post_full,2,quantile, probs = 0.25),
                             upr_50 = apply(diff_r_pred_post_full,2,quantile, probs = 0.75),
                             global_output$data_grid)
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

ggsave(diff_reg_plot,filename = "fig/full_dataset/diff_reg_ns.pdf", width = 10, height = 6)


# # Difference Regional Rate component
diff_r_rate_pred_post_full <- global_output$noisy_model_run_output$BUGSoutput$sims.list$diff_r_pred_deriv
diff_r_rate_post_df <- data.frame(pred = colMeans(diff_r_rate_pred_post_full),
                                  lwr = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.025),
                                  upr = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.975),
                                  lwr_50 = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.25),
                                  upr_50 = apply(diff_r_rate_pred_post_full,2,quantile, probs = 0.75),
                                  global_output$data_grid)
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

ggsave(diff_reg_rate_plot,filename = "fig/full_dataset/diff_reg_rate_ns.pdf", width = 10, height = 6)

