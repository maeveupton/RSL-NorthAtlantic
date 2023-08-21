# Cross validation for new approach for RSL analysis of East West of North Atlantic

# f(x,t) = c(t) + r_east(t)+ r_west(t) + gt + l(x,t) + Error
# Model 1: Fit the global/common trend in just time
# Model 2: Fit the regional component with fixed priors from previous run
# Model 3: Strong priors on both common, regional and estimate the non-linear local component


# Data from Andy's updated data
# Clear workspace
rm(list = ls())
#---------Set working directory--------------
setwd("/users/research/mupton/3. RSL_North_Atlantic/cross_val_common_north_atlantic_analysis_east_west")
#----------Load packages--------------------
library(devtools)
library(reslr)
library(tidyverse)
library(ggplot2)
remotes::install_github("tidyverse/purrr")
library(geosphere) # distm
library(readr) # read csv file
library(xtable)# Latex tables
# R functions-------
source("R/plot_data.R") # Plotting separate north & south data & map plot
source("R/internal_functions.R") # New basis functions for north & south data
source("R/reslr_mcmc_fun.R") # New basis functions for north & south data
source("R/cross_val_check_fun.R") # Cross validation
source("R/reslr_load_fun.R") # updated reslr_load


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
  ) %>%
  # Filter greenland as it has no tg:
  dplyr::filter(!Region == "Greenland") %>%
  # Filter out timoleague to see if it is weird
  dplyr::filter(!Site %in% c("Timoleague"))



# We include gia rates before cross val--------------------
SL_df_gia <- reslr_load_fun(
  data = SL_df,
  include_tide_gauge = TRUE,
  include_linear_rate = TRUE,
  # all_TG_1deg = TRUE,
  TG_minimum_dist_proxy = TRUE,
  # Including tg longer than 75 years in europe
  list_preferred_TGs = c(
    "DUBLIN",
    # Germany & Denmark extra 100 years data
    # "CUXHAVEN 2", "ESBJERG",
    # Amsterdam extra 100 & 75 years
    "VLISSINGEN", "MAASSLUIS", "HOEK VAN HOLLAND",
    "IJMUIDEN", "DEN HELDER", "HARLINGEN", "WEST-TERSCHELLING",
    # Belgium extra 100 & 75 years
    "NIEUWPOORT", "OOSTENDE", "ZEEBRUGGE",
    # France extra 100 & 75 years --> not working
    "DUNKERQUE", "BOULOGNE", "CALAIS",
    "BREST", "LE HAVRE", "ST. NAZAIRE",
    "ST JEAN DE LUZ (SOCOA)", "LA ROCHELLE-LA PALLICE",
    # UK extra 100 & 75 years--> site not working:"HOLYHEAD"
    # "LIVERPOOL (GLADSTONE DOCK)",
    "SHEERNESS", "TOWER PIER", "SOUTHEND",
    "NORTH SHIELDS", "NEWLYN",
    # Portugal & Spain extra 100 & 75 years --> sites not working: "LAGOS"
    "CASCAIS", "CADIZ II", "CEUTA", "VIGO", # "LA CORUÑA I"#,--> problem site
    # # North American sites
    "CHARLESTON I",
    "FORT PULASKI",
    "DAYTONA BEACH",
    "WILMINGTON",
    "SPRINGMAID PIER",
    "NAPLES",
    "FORT MYERS",
    # 50 years
    # "DEAL",
    # "PORTSMOUTH",
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

# Filling in the NA columns for the dataset:
input_gia_data <- 
  SL_df_gia$data %>%
  mutate(Basin = ifelse(is.na(Basin) == TRUE,"North Atlantic",Basin),
         Region = ifelse(is.na(Region) == TRUE, "",Region),
         Reference = ifelse(is.na(Reference) == TRUE, "PSMSL database", Reference),
         SiteName = as.character(SiteName),
        Site = ifelse(is.na(Site) == TRUE,SiteName,Site),
        SiteName = as.factor(SiteName))


# # Running cross validation for just proxy records------------------
# east_west_cross_val_list <- cross_val_check_fun(
#   data = input_gia_data,
#   data_grid = SL_df_gia$data_grid,
#   n_fold = 10,
#   spline_nseg_c = 20,
#   spline_nseg_t = 8,
#   spline_nseg_st = 6,
#   n_iterations = 8000,
#   n_burnin = 1000,
#   n_thin = 5,
#   n_chains = 2,
#   seed = 38233,
#   prediction_grid_res = 30
# )

#saveRDS(east_west_cross_val_list, file = "reslr_outputs/cross_val_ew_without_extra_sites.rds")

cross_val_ew <- readRDS("reslr_outputs/cross_val_ew_without_extra_sites.rds")

cross_val_ew$ME_MAE_RSME_fold_site
cross_val_ew$ME_MAE_RSME_overall
true_pred_plt <- cross_val_ew$true_pred_plot_proxy
true_pred_plt <- cross_val_ew$true_pred_plot_tg
ggsave(true_pred_plt,filename = "fig/true_pred_ew_proxy.pdf", width = 10, height = 6)
cross_val_ew$total_coverage
CV_model_df_proxy <- cross_val_ew$CV_model_df %>% 
  filter(!data_type_id == "TideGaugeData") 

CV_model_df_tg <- cross_val_ew$CV_model_df %>% 
  filter(data_type_id == "TideGaugeData") %>% 
  dplyr::mutate(SiteName = gsub(",\n ","",SiteName))

CV_model_df<- rbind(CV_model_df_proxy,CV_model_df_tg)

true_pred_plot_full <- ggplot2::ggplot(data = CV_model_df, ggplot2::aes(
  x = true_RSL,
  y = y_post_pred,
  colour = "PI"
)) +
  ggplot2::geom_errorbar(
    data = CV_model_df,
    ggplot2::aes(
      x = true_RSL,
      ymin = lwr_PI,
      ymax = upr_PI
    ),
    colour = "red3",
    width = 0, alpha = 0.5
  ) +
  ggplot2::geom_point() +
  ggplot2::geom_abline(
    data = CV_model_df,
    ggplot2::aes(intercept = 0, slope = 1, colour = "True = Predicted")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 9, face = "bold"),
    axis.text = ggplot2::element_text(size = 9),
    strip.background = ggplot2::element_rect(fill = c("white")),
    strip.text = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7),
    legend.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 8)
  ) +
  ggplot2::theme(legend.box = "horizontal", 
                 legend.position = "bottom") +
  ggplot2::labs(
    x = "True Relative Sea Level (m)",
    y = "Predicted Relative Sea Level (m)"
  ) +
  ggplot2::scale_colour_manual("",
                               values = c(
                                 c(
                                   "PI" = "red3",
                                   # "True = Predicted" = "black")
                                   "True = Predicted" = "black"
                                 )
                               ),
                               labels = c(
                                 "PI" = paste0(unique(CV_model_df$CI), " Prediction Interval"),
                                 "True = Predicted" = "True = Predicted"
                               )
  ) +
  ggplot2::facet_wrap(~SiteName, scales = "free") +
  ggplot2::guides(
    colour = ggplot2::guide_legend(override.aes = list(
      linetype = c(1, 1),
      shape = c(NA, NA),
      size = 2
    ))
  )+
  theme(
    panel.spacing = unit(0,'lines')
  )+
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.title=element_blank(),
              legend.margin=margin(c(1,5,5,5)))
true_pred_plot_full
ggsave(true_pred_plot_full,filename = "fig/true_pred_ew.pdf", width = 10, height = 7)


# Prediction Interval size
prediction_interval_size <- CV_model_df %>%
  dplyr::group_by(SiteName) %>%
  dplyr::reframe(PI_width = abs(unique(mean(upr_PI - lwr_PI))),
                 coverage_by_site =
                   unique(length(which(obs_in_PI == "TRUE")) / dplyr::n()),
                 RSME = unique(sqrt((sum(true_RSL - pred_RSL)^2) / dplyr::n())))
print(xtable(prediction_interval_size,type = "latex",digits=4),
      file = "reslr_outputs/coverage_by_site_RSME.tex",
      include.rownames=FALSE)
