# Other proxy data tests

# Clear workspace
rm(list = ls())
#---------Set working directory--------------
setwd("/users/research/mupton/3. RSL_North_Atlantic/Caesar2021paper_plots")
library(tidyverse)
library(ggplot2)
library(zoo)
library(cowplot)
library(reticulate) # calling in Python script
library(akima) # cubic interpolation
library(ggrepel)# labelling different lines in final plot


# Ramhsoft 2015 data:
rahm_data <- read.csv("https://www.dropbox.com/s/t9gbld343tf2fcx/rahmsoft2015.csv?dl=1")
# y is an amoc index
names(rahm_data) <- c("year", "Amoc", "2sigma", "lower", "upper", "ID")
# Dataframe form LC Python code
rahm_df_update <- read.csv("LC_outputs/lowess_outputs/rahm2015.csv")
colnames(rahm_df_update) <- c("index", "year", "mean","lower","upper")

rahm_df <- rahm_df_update %>%
  dplyr::select(!index) %>%
  mutate(
    ID = "Rahmstorf et al. 2015"
    #upper = rahm_data$upper,#mean + unique(rahm_data$`2sigma`),
    #lower = rahm_data$lower#mean - unique(rahm_data$`2sigma`),
  )


# Smoothing the data by 50 years using lowess filter --> 50/1075
# smooth_rahm <- lowess(rahm_df$year, rahm_df$Amoc, f = 50 / 1075)

# rahm_smooth_df <- data.frame(
#   year = smooth_rahm$x,
#   Amoc = smooth_rahm$y,
#   upper_smooth = smooth_rahm$y + rahm_df$`2sigma`,
#   lower_smooth = smooth_rahm$y - rahm_df$`2sigma`,
#   ID = rahm_df$ID
# )


rahm_plot <-
  ggplot() +
  geom_line(data = rahm_df, aes(x = year, y = mean, colour = ID)) +
  geom_ribbon(
    data = rahm_df,
    aes(x = year, ymax = lower, ymin = upper, fill = ID), alpha = 0.4
  ) +
  labs(x = "Year (CE)", y = "Temperature anomaly (K)", colour = "", fill = "") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2005)) +
  scale_fill_manual(values = c("cyan4")) +
  scale_colour_manual(values = c("cyan4")) +
  theme(text = element_text(size = 7))
rahm_plot
ggsave(rahm_plot, filename = "fig/rahmstorf2015.pdf", width = 6, height = 10)

# Spooner 2020--------------
spooner_df <- read.csv("https://www.dropbox.com/scl/fi/9looo0teo061ayyn6wz95/Spooner2020.csv?rlkey=uwhrs1250cifp6byhvsxe28qb&dl=1")
colnames(spooner_df) <- c("year", "mean", "lower", "upper")
# Dataframe form LC Python code
spooner_df_update <- read.csv("LC_outputs/lowess_outputs/Spooner.csv")
colnames(spooner_df_update) <- c("index", "year", "mean", "lower", "upper")
spooner_df_new <- spooner_df_update %>%
  dplyr::select(!index) %>%
  mutate(ID = "Spooner et al. 2020")

# spooner_df <- spooner_df %>%
#   arrange(year) %>%
#   mutate(
#     ID = "Spooner et al. 2020",
#     year = round(year),
#     # Not sure about this?
#     `2sigma` = mean - lower,
#     sigma_test = 2 * (sd(mean) / sqrt(n()))
#   )

# # Sequence of full years in dataframe:
# year <- data.frame(year = seq(392,2013,by = 1))
# #spooner_long_df <- merge(year,spooner_df, by = "year",all.x = TRUE) %>%
# #  mutate(ID = "Spooner et al. 2020")
#
#
# # Smoothing the data by 50 years using lowess filter--> range : #50/1620.94
# smooth_spooner <- lowess(spooner_df$year, spooner_df$mean,
#   f = 0.05#50/1621
# )
# spooner_smooth_df <- data.frame(
#   year = smooth_spooner$x,
#   mean = smooth_spooner$y,
#   upper_smooth = smooth_spooner$y + spooner_df$sigma_test,#`2sigma`, # sigma_test,#spooner_df$`2sigma`,
#   lower_smooth = smooth_spooner$y - spooner_df$sigma_test,#`2sigma`, # sigma_test,
#   ID = spooner_df$ID
# )

spooner_plot <-
  ggplot() +
  geom_line(data = spooner_df_new, aes(x = year, y = mean, colour = ID)) +
  geom_ribbon(
    data = spooner_df_new,
    aes(x = year, ymax = lower, ymin = upper, fill = ID), alpha = 0.4
  ) +
  labs(x = "Year (CE)", y = "Abundance of \n high productivity \n species (%)", colour = "", fill = "") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2005)) +
  scale_fill_manual(values = c("orange")) +
  scale_colour_manual(values = c("orange")) +
  theme(text = element_text(size = 7))
spooner_plot
ggsave(spooner_plot, filename = "fig/spooner2020.pdf", width = 6, height = 10)

# Thornally 2018 Silt 1----------
thornally_silt1_data <- read.csv("https://www.dropbox.com/scl/fi/gtndau5ebns1gujjnw0ld/Thornally2018_s_silt1.csv?rlkey=p2qmaplp0kl2iute1jkxxvs9n&dl=1")
colnames(thornally_silt1_data) <- c("year", "mean", "smooth_mean", "2sigma", "lower", "upper")
# Dataframe form LC Python code
thornally_silt1_df_update <- read.csv("LC_outputs/lowess_outputs/thorn56.csv")
colnames(thornally_silt1_df_update) <- c("index", "year", "mean")

thornally_silt1_df <- thornally_silt1_df_update %>%
  dplyr::select(!index) %>%
  mutate(
    ID = "Thornally et al. 2018 56JPC",
    upper = mean + 0.8, # unique(thornally_silt1_data$`2sigma`),
    lower = mean - 0.8 # unique(thornally_silt1_data$`2sigma`),
  )

# Sequence of full years in dataframe:
# year <- data.frame(year = seq(1474,2004,by = 1))
# test_thorn_slit1_df <- merge(year,thornally_silt1_df, by = "year",all.x = TRUE) %>%  mutate(ID = "Thornally et al. 2018 56JPC")
# Smoothing the data by 50 years using lowess filter
# 50/528.7
# smooth_thornally_silt1 <- lowess(thornally_silt1_df$year,
#                                 thornally_silt1_df$mean, f = 0.0945)#0.05)# f??
# smooth_thornally_silt1 <- lowess(test_thorn_slit1_df$year,test_thorn_slit1_df$mean,f = 0.0945)#0.04)#0.0945)#0.05)# f??
# thornally_silt1_smooth_df <- data.frame(
#   year = smooth_thornally_silt1$x,
#   mean = smooth_thornally_silt1$y,
#   upper_smooth = smooth_thornally_silt1$y + thornally_silt1_df$`2sigma`,
#   lower_smooth = smooth_thornally_silt1$y - thornally_silt1_df$`2sigma`,
#   ID = thornally_silt1_df$ID
# )

thornally_plot_silt1 <-
  ggplot() +
  geom_line(data = thornally_silt1_df, aes(x = year, y = mean, colour = ID)) +
  geom_ribbon(
    data = thornally_silt1_df,
    aes(x = year, ymax = lower, ymin = upper, fill = ID), alpha = 0.4
  ) +
  # geom_line(data = thornally_silt1_df, aes(x = year, y = mean, colour = ID)) +
  # geom_ribbon(
  #  data = thornally_silt1_df,
  #  aes(x = year, ymax = lower, ymin = upper, fill = ID), alpha = 0.4
  # ) +
  labs(x = "Year (CE)", y = "56JPC s\u0305s\u0305 (\U003BCm)", colour = "", fill = "") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2005)) +
  scale_y_continuous(limits = c(29, 34.5)) +
  scale_fill_manual(values = "darkseagreen4") +
  scale_colour_manual(values = "darkseagreen4") +
  theme(text = element_text(size = 7))
thornally_plot_silt1
ggsave(thornally_plot_silt1, filename = "fig/thornally2018_silt1.pdf", width = 6, height = 10)

# Thornally 2018 silt 2------
thornally_silt2_df <- read.csv("https://www.dropbox.com/scl/fi/alj227ll6ziionz5gnkif/Thornally2018_s_silt2.csv?rlkey=vcci4tluwij73vubfkp75dnpi&dl=1")
colnames(thornally_silt2_df) <- c("year", "mean", "smooth_mean", "2sigma", "lower", "upper")
thornally_silt2_df <- thornally_silt2_df %>%
  mutate(ID = "Thornally et al. 2018 48JPC") %>%
  # Remove row of NA
  drop_na()

# Dataframe form LC Python code
thornally_silt2_df_update <- read.csv("LC_outputs/lowess_outputs/thorn48.csv")
colnames(thornally_silt2_df_update) <- c("index", "year", "mean")

thornally_silt2_df <- thornally_silt2_df_update %>%
  dplyr::select(!index) %>%
  mutate(
    ID = "Thornally et al. 2018 48JPC",
    upper = mean + unique(thornally_silt2_df$`2sigma`),
    lower = mean - unique(thornally_silt2_df$`2sigma`),
  )

## Smoothing the data by 50 years using lowess filter
# smooth_thornally_silt2 <- lowess(thornally_silt2_df$year, thornally_silt2_df$mean, f = 0.09) # 0.01)#0.09)# 0.030977)#0.05)# f??

# thornally_silt2_smooth_df <- data.frame(
#  year = smooth_thornally_silt2$x,
#  mean = smooth_thornally_silt2$y,
#  upper_smooth = smooth_thornally_silt2$y + thornally_silt2_df$`2sigma`,
#  lower_smooth = smooth_thornally_silt2$y - thornally_silt2_df$`2sigma`,
#  ID = thornally_silt2_df$ID
# )


thornally_plot_silt2 <-
  ggplot() +
  geom_line(data = thornally_silt2_df, aes(x = year, y = mean, colour = ID)) +
  geom_ribbon(
    data = thornally_silt2_df,
    aes(x = year, ymax = lower, ymin = upper, fill = ID), alpha = 0.4
  ) +
  # geom_line(data = thornally_silt2_df, aes(x = year, y = `smooth.mean.SS..µm.`,colour = ID))+
  # geom_ribbon(data = thornally_silt2_df,
  #            aes(x = year,ymax = `lower.bound.95.`,ymin = `upper.bound.95.`,fill = ID),alpha = 0.4)+
  labs(x = "Year (CE)", y = "48JPC s\u0305s\u0305 (\U003BCm)", colour = "", fill = "") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2005)) +
  scale_y_continuous(limits = c(27.5, 32.5)) +
  scale_fill_manual(values = "darkgoldenrod4") +
  scale_colour_manual(values = "darkgoldenrod4") +
  theme(text = element_text(size = 7))
thornally_plot_silt2
ggsave(thornally_plot_silt2, filename = "fig/thornally2018_silt2.pdf", width = 6, height = 10)

# Thibodeau 2018---------
thibodeau_data <- read.csv("https://www.dropbox.com/scl/fi/evpqashn818i7s9rs6yht/Thibodeau2018_md99_220.csv?rlkey=gqh1e5l2ri0z8xcutd2w0cl1q&dl=1")
colnames(thibodeau_data) <- c("year", "mean", "2sigma", "lower", "upper")

# Dataframe form LC Python code
# thibodeau_df_update <- read.csv("LC_outputs/lowess_outputs/thibcr2_2018.csv")
thibodeau_df_update <- read.csv("LC_outputs/lowess_outputs/thibmd99_2018.csv")
colnames(thibodeau_df_update) <- c("index", "year", "mean")

thibodeau_df <- thibodeau_df_update %>%
  dplyr::select(!index) %>%
  mutate(
    ID = "Thibodeau et al. 2018",
    upper = mean + unique(thibodeau_data$`2sigma`),
    lower = mean - unique(thibodeau_data$`2sigma`),
  )

# # Smoothing the data by 50 years using lowess filter
#
# smooth_thibodeau <- lowess(thibodeau_df$year, thibodeau_df$mean, f = 50 / range(thibodeau_df$year)) # 0.05)# f??
# thibodeau_smooth_df <- data.frame(
#   year = smooth_thibodeau$x,
#   mean = smooth_thibodeau$y,
#   upper_smooth = smooth_thibodeau$y + thibodeau_df$`2sigma`,
#   lower_smooth = smooth_thibodeau$y - thibodeau_df$`2sigma`,
#   ID = thibodeau_df$ID
# )


thibodeau_plot <-
  ggplot() +
  geom_line(data = thibodeau_df, aes(x = year, y = mean, colour = ID)) +
  geom_ribbon(
    data = thibodeau_df,
    aes(
      x = year, ymax = lower,
      ymin = upper, fill = ID
    ), alpha = 0.4
  ) +
  # geom_ribbon(data = thibodeau_df,
  #            aes(x = year,ymax = `lower.bound.95.`,ymin = `upper.bound.95.`,fill = ID),alpha = 0.4)+
  labs(x = "Year (CE)", y = paste(expression("\u03B4^18"), "O", "\n", "MD99-2220 (%o VPDB)"), colour = "", fill = "") +
  theme_bw() +
  # scale_x_continuous(limits = c(0,2005))+
  scale_x_continuous(limits = c(500, 2000)) +
  scale_y_continuous(limits = c(2.75, 3.3)) +
  scale_fill_manual(values = "coral4") +
  scale_colour_manual(values = "coral4") +
  theme(text = element_text(size = 7))
thibodeau_plot
ggsave(thibodeau_plot, filename = "fig/thibodeau2018.pdf", width = 6, height = 10)

# Sherwood 2011
sherwood_data <- read.csv("https://www.dropbox.com/scl/fi/11hnpl10jnua414cv5ftj/Sherwood2011.csv?rlkey=co85c633ckbsa2qxx4vnst8rg&dl=1")
colnames(sherwood_data) <- c("year", "mean", "lower", "upper")

# Dataframe form LC Python code
sherwood_df_update <- read.csv("LC_outputs/lowess_outputs/Sherwood2011.csv")
colnames(sherwood_df_update) <- c("index", "year", "mean", "lower", "upper")

sherwood_df <- sherwood_df_update %>%
  dplyr::select(!index) %>%
  mutate(
    ID = "Sherwood et al. 2011"
  )

# # Smoothing the data by 50 years using lowess filter
# #
# smooth_sherwood <- lowess(sherwood_df$year, sherwood_df$mean, f = 0.05) # f??
# sherwood_smooth_df <- data.frame(
#   year = smooth_sherwood$x,
#   mean = smooth_sherwood$y,
#   upper_smooth = smooth_sherwood$y + sherwood_df$`2sigma`,
#   lower_smooth = smooth_sherwood$y - sherwood_df$`2sigma`,
#   ID = sherwood_df$ID
# )

sherwood_point_df <- read.csv("https://www.dropbox.com/scl/fi/ram6zvii8q6ee5tokpi22/sherwood2011_points.csv?rlkey=s08nnzwx9bez3i5odhljp58xe&dl=1")
colnames(sherwood_point_df) <- c("year_min", "year_max", "year", "mean", "2sigma", "lower", "upper")
sherwood_point_df <- sherwood_point_df %>%
  mutate(ID = "Sherwood et al. 2011")

sherwood_both_plot <-
  ggplot() +
  geom_line(data = sherwood_df, aes(x = year, y = mean, colour = ID)) +
  geom_ribbon(
    data = sherwood_df,
    aes(x = year, ymax = lower, ymin = upper, fill = ID), alpha = 0.4
  ) +
  # Points
  geom_point(data = sherwood_point_df, aes(x = year, y = mean, colour = ID)) +
  geom_errorbar(data = sherwood_point_df, aes(
    x = year, y = mean,
    ymin = lower, ymax = upper, colour = ID
  )) +
  geom_errorbarh(data = sherwood_point_df, aes(
    y = mean,
    xmin = year_min, xmax = year_max, colour = ID
  )) +
  labs(x = "Year (CE)", y = paste("Bulk", expression("\u03B4^15"), "N", "(%o)"), colour = "", fill = "") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2005)) +
  # scale_y_continuous(limits = c(25,35))+
  scale_fill_manual(values = c("darkolivegreen3", "darkolivegreen3")) +
  scale_colour_manual(values = c("darkolivegreen3", "darkolivegreen3")) +
  theme(text = element_text(size = 7))
sherwood_both_plot
ggsave(sherwood_both_plot, filename = "fig/sherwood2011both.pdf", width = 6, height = 10)


# Upton 2023
upton_df <- read.csv("/users/research/mupton/3. RSL_North_Atlantic/1.common_north_atlantic_analysis_east_west/Upton_reg_diff_east_west.csv")
upton_df <- upton_df %>% mutate(ID = "Upton et al. 2023")
upton_plot <-
  ggplot() +
  geom_line(data = upton_df, aes(x = Age, y = pred, colour = ID)) +
  geom_ribbon(
    data = upton_df,
    aes(x = Age, ymax = upr, ymin = lwr, fill = ID), alpha = 0.4
  ) +
  labs(x = "Year (CE)", y = "East - West Sea Level (m)", colour = "", fill = "") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2005)) +
  # scale_y_continuous(limits = c(25,35))+
  scale_fill_manual(values = ("darkmagenta")) +
  scale_colour_manual(values = ("darkmagenta")) +
  theme(text = element_text(size = 7))
upton_plot
ggsave(upton_plot, filename = "fig/upton2023.pdf", width = 6, height = 10)


# Joining the plots together vertically:
compare_plot <- plot_grid(upton_plot,
  rahm_plot,
  thibodeau_plot,
  sherwood_both_plot,
  spooner_plot,
  thornally_plot_silt2,
  thornally_plot_silt1,
  ncol = 1, align = "hv"
)
ggsave(compare_plot, filename = "fig/compare_proxy_data.pdf", height = 15, width = 10)


# Put all the dataframes together with same x axis
upton_df_update <- upton_df %>%
  dplyr::select(pred, lwr, upr, ID, Age) %>%
  rename(mean = pred, lower = lwr, upper = upr, year = Age)

# Adding column with y label for each df
upton_df_update <- upton_df_update %>%
  mutate(y_lab_name = "East - West Sea Level (m)")
rahm_df <- rahm_df %>%
  mutate(y_lab_name = "Temperature anomaly (K)")
spooner_df_new <- spooner_df_new %>%
  mutate(y_lab_name = "Abundance of \n high productivity \n species (%)")
thibodeau_df <- thibodeau_df %>%
  mutate(y_lab_name = paste(expression("\u03B4^18"), "O", "\n", "MD99-2220 (%o VPDB)"))
thornally_silt1_df <- thornally_silt1_df %>%
  mutate(y_lab_name = "48JPC s\u0305s\u0305 (\U003BCm)")
thornally_silt2_df <- thornally_silt2_df %>%
  mutate(y_lab_name = "56JPC s\u0305s\u0305 (\U003BCm)")
sherwood_df <- sherwood_df %>%
  mutate(y_lab_name = paste("Bulk", expression("\u03B4^15"), "N", "(%o)"))
sherwood_point_df <- sherwood_point_df %>%
  mutate(y_lab_name = paste("Bulk", expression("\u03B4^15"), "N", "(%o)"))

full_df <- rbind(
  upton_df_update,
  rahm_df,
  spooner_df_new,
  thibodeau_df,
  thornally_silt1_df,
  sherwood_df,
  thornally_silt2_df
)
# sherwood_point_df
# )

my_plot_labs <- as_labeller(
  c(
    "Upton et al. 2023" = "East - West Sea Level (m)",
    "Rahmstorf et al. 2015" = "Temperature anomaly (K)",
    "Spooner et al. 2020" = "Abundance of \n high productivity \n species (%)",
    "Thibodeau et al. 2018" = paste(expression("\u03B4^{18}"), "O", "\n", "MD99-2220 (%o VPDB)"),
    "Thornally et al. 2018 48JPC" = "48JPC s\u0305s\u0305 (\U003BCm)",
    "Thornally et al. 2018 56JPC" = "56JPC s\u0305s\u0305 (\U003BCm)",
    "Sherwood et al. 2011" = paste("Bulk", expression("\u03B4^{15}"), "N", "(%o)")
  )
)
# Subset for labels
subset_labs <- full_df %>% 
  group_by(ID) %>%
  slice_min(year,with_ties = FALSE)

# Plotting all data together:
all_plot <-
  ggplot() +
  geom_point(data = sherwood_point_df, aes(x = year, y = mean, colour = ID)) +
  geom_errorbar(data = sherwood_point_df, aes(
    x = year, y = mean,
    ymin = lower, ymax = upper, colour = ID
   )) +
  geom_errorbarh(data = sherwood_point_df, aes(
    y = mean,
    xmin = year_min, xmax = year_max, colour = ID
  )) +
  geom_line(data = full_df, aes(x = year, y = mean, colour = ID)) +
  geom_ribbon(
    data = full_df,
    aes(x = year, ymax = upper, ymin = lower, fill = ID), alpha = 0.4
  ) +
  # labs(x = "Year (CE)", y = "East - West Sea Level (m)", colour = "", fill = "") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2005)) +
  # scale_y_continuous(limits = c(25,35))+
  # scale_fill_manual(values = ("darkmagenta")) +
  # scale_colour_manual(values = ("darkmagenta")) +
  theme(text = element_text(size = 7)) +
  facet_wrap(ID ~ .,
    scales = "free_y",
    ncol = 1, # strip.position="right",
    strip.position = "left",labeller = my_plot_labs
  ) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside", 
    axis.title.x = element_text(size=8, face="bold", color = "black"),
    strip.text.y = element_text(size=8, face="bold", color = "black")
  ) +
  # geom_label_repel(data = full_df,aes(x = year, y = mean,label = ID),
  #                  nudge_x = 1,
  #                  na.rm = TRUE)+
  #geom_text(data = subset_labs,aes(x = 0, y = mean,label = ID,colour = ID))+
  ylab(NULL)+
  xlab("Year (CE)")+
  labs(colour = NULL,fill = NULL)
all_plot
ggsave(all_plot, filename = "fig/all_plot.pdf", 
       height = 15, width = 10,dpi = 1200,
       device =cairo_pdf)
