#' Plot raw data with measurement uncertainty.
plot_data <- function(input_data, title = "",
                      xlab = "Year (CE)",
                      ylab = "Relative Sea Level (m)",
                      plot_tide_gauges = FALSE,
                      plot_proxy_records = TRUE,
                      plot_caption = TRUE) {

  # Input data-------
  data <- input_data
  # data <- input_data$data
  # data_grid <- input_data$data_grid
  n_sites <- length(data$SiteName %>% unique())
  n_proxy <- data %>%
    dplyr::filter(data_type_id == "ProxyRecord") %>%
    dplyr::select(SiteName, data_type_id) %>%
    unique() %>%
    nrow()
  # Plotting only Proxy Record
  if (plot_proxy_records == TRUE & plot_tide_gauges == FALSE) {
    data <- data %>%
      dplyr::filter(data_type_id == "ProxyRecord")
  }
  # Plotting tide gauge only
  if (plot_proxy_records == FALSE & plot_tide_gauges == TRUE) {
    data <- data %>%
      dplyr::filter(data_type_id == "TideGaugeData")
  }
  # Raw data plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = data, ggplot2::aes(
      xmin = Age - Age_err, xmax = Age + Age_err,
      ymin = RSL - RSL_err, ymax = RSL + RSL_err,
      fill = "gray",
    ), alpha = 0.7) +
    ggplot2::geom_point(
      data = data,
      ggplot2::aes(y = RSL, x = Age, colour = "black"), size = 0.3
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_bw() +
    ggplot2::labs(colour = "") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 7),
      strip.background = ggplot2::element_rect(fill = c("white"))
    ) +
    ggplot2::scale_fill_manual("",
      values = "grey",
      labels = expression(paste("1-sigma Error")),
      guide = ggplot2::guide_legend(override.aes = list(alpha = 0.7))
    ) +
    ggplot2::scale_colour_manual(
      values = c("black"),
      labels = c("Data")
    ) +
    ggplot2::facet_wrap(~SiteName) +
    ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_blank()
    )

   # Plotting both TG & proxy
  if (plot_proxy_records == TRUE & plot_tide_gauges == TRUE) {
    p <- p + ggplot2::facet_wrap(~SiteName, scales = "free")
  }

  # Informed caption
  if (plot_caption == TRUE) {
    p <- p + ggplot2::labs(caption = paste0(
      "No. proxy sites:", n_proxy,
      "\n No. tide gauge sites:", n_sites - n_proxy
    ))
  } else {
    p <- p
  }


  
  return(p)
}

# Map of the data sites
plot_map <- function(data){
  # Plot map of data:
  world <- ne_countries(scale = "medium", returnclass = "sf")
  SL_df_unique <- data %>%  
    dplyr::select(SiteName,Longitude,Latitude,data_type_id) %>% 
    unique() %>% 
    mutate(SiteName = as.character(SiteName))
  map_locations <- ggplot(data = world) +
    geom_sf(color="darkgrey",fill="darkseagreen3")+
    geom_point(data=SL_df_unique,aes(x = Longitude,y = Latitude,
                                     colour = data_type_id,shape = data_type_id))+
    coord_sf(xlim = c(-98, 20), ylim = c(20, 80), expand = FALSE)+
    labs(x="Longitude",y="Latitude")+
    annotation_scale(location = "br", width_hint = 0.5)+
    #geom_text(aes(x = -40,y=35),label = "North Atlantic Ocean",
    #          fontface="italic",size = 6)+
    annotation_north_arrow(location = "br", which_north = "true",
                           pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                           style = north_arrow_fancy_orienteering) +
    theme(plot.title = element_text(size=18,face="bold"),
          #panel.grid.major = element_line(color = "lightblue1", linetype = "dashed", size = 0.05),
          panel.background = element_rect(fill = "aliceblue"),
          axis.title=element_text(size=10,face="bold"),
          legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    scale_colour_manual(values = c("black","saddlebrown"),
                        labels=c("Proxy Record Sites","Tide Gauge Sites"),
                        guide = TRUE)+
    scale_shape_manual(values = c(16,17),
                       labels=c("Proxy Record Sites","Tide Gauge Sites"),
                       guide = TRUE)+
    theme(legend.position=c(0,1),legend.justification = c(0,1), 
          legend.box.background = element_rect(color="black",fill="white", linewidth=0.5),
          #legend.box = "horizontal",legend.margin = margin(6, 6, 6, 6),
          legend.title=element_blank(),
          legend.text = element_text(size = 10,face="bold"))+
    guides(color = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3)))
  
  
  map_locations
  ggsave(map_locations,filename = "fig/map_NA_no_labels.pdf", width = 10, height = 6)
  
}