#' Loading in data
#'
#' @param data The input data as a dataframe.
#' @param input_age_type The inputted age in years "CE" or year "BP". The default is "CE" and is the preferred structure of the package. The package has the ability to use Before Present ("BP") observations.
#' @param prediction_grid_res Resolution of grid. Predictions over every 50 years(default) can vary based on user preference, as larger values will reduce computational run time.
#' @param include_linear_rate User decides to include linear_rate and linear_rate_err associated. This relates to linear rate which corresponds to an important physical process that impacts sea level observations which is glacial isostatic adjustment (GIA). For the linear_rate and its associated linear_rate_err, the user can provide these values as additional columns in the input dataframe. If they prefer, the package will calculate the linear_rate and the linear_rate_err using the data.
#' @param include_tide_gauge Including tide gauge data from PSMSL website. The tide gauge data undergo a cleaning process in this function where flagged stations are removed as recommended by the online database. Next, the tide gauge data is averaged over period defined by sediment_average_TG which default is 10 years corresponding to accumulation rates of proxy records. Then, the user selects their preferred tide gauge based on three criteria: 1.nearest tide gauge to the proxy site; 2. User supplies a list of names of preferred tide gauges; 3. all tide gauges within 1 degree are chosen.
#' @param sediment_average_TG Average the tide gauge data to make it comparable to accumulation rates of proxy records. The default averaging period for tide gauges is 10 years and the user can alter this.
#' @param list_preferred_TGs The user can supply the name or names of the preferred tide gauges from the PSMSL database.
#' @param TG_minimum_dist_proxy The package finds the tide gauge closest to the proxy site
#' @param all_TG_1deg The package finds all tide gauges within 1 degree of the proxy site
#' @param detrend_data Detrend the data using the linear rate provided to remove this component.
#' @param core_col_year The year the sediment core was collected in order to the data to be detrended.
#' @param cross_val For the spline in time, spline in space time and the NIGAM the user can undertake cross validation to examine the
#' @param test_set The test set dataframe for cross validation
#'
#' @return A list containing data frame of data and prediction grid. The output of this function is two data frames, one with the data and one with the data_grid which represent a grid with evenly spaced time points. If tide gauge data is used, an ID column is included in the two output dataframes displaying the data source, "ProxyRecord" or "TideGaugeData".
reslr_load_fun <- function(data,
                       prediction_grid_res = 30,
                       include_tide_gauge = FALSE,
                       include_linear_rate = FALSE,
                       list_preferred_TGs = NULL,
                       TG_minimum_dist_proxy = FALSE,
                       all_TG_1deg = FALSE,
                       sediment_average_TG = 10,
                       detrend_data = FALSE,
                       core_col_year = NULL,
                       cross_val = FALSE,
                       test_set) {
  Age <- Age_BP <- RSL <- Age_err <- RSL_err <- SiteName <- max_Age <- min_Age<- bounds <- Longitude <- Latitude <- Site <- Region <- data_type_id <- ICE5_GIA_slope <- linear_rate_err <- linear_rate <- n <- obs_index <- x_4_upr <-x_lwr_box<- x_upr_box<- y_1_lwr<-y_lwr<- y_upr<- NULL
  # Dividing Age & Age_err by 1000 for easier calculations-----
  data <- data %>%
    dplyr::mutate(Age = Age / 1000) %>%
    dplyr::mutate(Age_err = Age_err / 1000)
  # If minus signs for Long and Lat maybe an issue
  
  # Tidy Original data-------------------------------
  if (!("SiteName" %in% colnames(data))) {
    data <- data %>%
      dplyr::mutate(SiteName = as.factor(paste0(Site, ",", "\n", " ", Region)))
  } else {
    message("Error: User must provide a column with site name(s) and a column with region name(s). \n")
    stop()
  }

  
  # Including no TG or no linear rates
  if (include_tide_gauge == FALSE & include_linear_rate == FALSE) {
    data <- data %>% dplyr::mutate(data_type_id = "ProxyRecord")
  }
  
  # Including TG & no linear rates but forget to include TG method
  if (include_tide_gauge == TRUE &
      include_linear_rate == FALSE &
      is.null(list_preferred_TGs) == TRUE &
      TG_minimum_dist_proxy == FALSE &
      all_TG_1deg == FALSE) {
    # # If the user wants to use their own tide gauge data
    # if("data_type_id" %in% colnames(data) & "TideGaugeData" %in% data$data_type_id){
    #   data <- data
    # }
    #else{
    stop("Error: No tide gauge selection method chosen, please provide criteria for choosing preferred tide gauge. \n")
  }
  
  
  # Including TGs and no linear rates
  if (include_tide_gauge == TRUE & include_linear_rate == FALSE) {
    # # If the user wants to use their own tide gauge data
    # if("data_type_id" %in% colnames(data) & "TideGaugeData" %in% data$data_type_id){
    #   data <- data
    # }
    # else{
    data <- data %>% dplyr::mutate(data_type_id = "ProxyRecord")
    data <- clean_tidal_gauge_data(
      data = data,
      list_preferred_TGs = list_preferred_TGs,
      TG_minimum_dist_proxy = TG_minimum_dist_proxy,
      all_TG_1deg = all_TG_1deg,
      sediment_average_TG = sediment_average_TG
    )
  }
  
  
  
  # Including linear rates & no TG
  if (include_linear_rate == TRUE & include_tide_gauge == FALSE) {
    # Checking if user provided GIA rates----------
    if (!("linear_rate" %in% colnames(data) & "linear_rate_err" %in% colnames(data))) {
      lm_data_rates <- linear_reg_rates(data)
      data <- dplyr::left_join(data, lm_data_rates, by = "SiteName")
      data <- data %>% dplyr::mutate(data_type_id = "ProxyRecord")
    } else {
      data <- data
      data <- data %>% dplyr::mutate(data_type_id = "ProxyRecord")
    }
  }
  
  # Including TG & linear rates but forget to include TG method
  if (include_tide_gauge == TRUE &
      include_linear_rate == TRUE &
      is.null(list_preferred_TGs) == TRUE &
      TG_minimum_dist_proxy == FALSE &
      all_TG_1deg == FALSE) {
    stop("Error: No tide gauge selection method chosen. Select criteria to chose your prefered tide gauge")
  }
  
  # Including linear rates & TG data
  if (include_linear_rate == TRUE & include_tide_gauge == TRUE) {
    data <- data %>% dplyr::mutate(data_type_id = "ProxyRecord")
    # Checking if user provided GIA rates----------
    if (!("linear_rate" %in% colnames(data) & "linear_rate_err" %in% colnames(data))) {
      lm_data_rates <- linear_reg_rates(data)
      data <- dplyr::left_join(data, lm_data_rates, by = "SiteName")
    } else {
      data <- data
    }
    data <- clean_tidal_gauge_data(
      data = data,
      list_preferred_TGs = list_preferred_TGs,
      TG_minimum_dist_proxy = TG_minimum_dist_proxy,
      all_TG_1deg = all_TG_1deg,
      sediment_average_TG = sediment_average_TG
    )
    #---Adding linear rates from ICE5G for TG-----
    data <- add_linear_rate(data = data)
    data <- data %>%
      dplyr::mutate(
        linear_rate = ifelse(data_type_id == "TideGaugeData", ICE5_GIA_slope, linear_rate),
        linear_rate_err = ifelse(data_type_id == "TideGaugeData", 0.3, linear_rate_err)
      )
  }
  
  # Prediction dataframe-------------------------------------
  sites <- data %>%
    dplyr::select(
      Longitude,
      Latitude,
      SiteName,
      dplyr::contains("linear_rate"),
      dplyr::contains("linear_rate_err"),
      dplyr::contains("data_type_id")
    ) %>%
    unique()
    
  times <- rep(seq(min(data$Age),
                   max(data$Age),
                   by = prediction_grid_res / 1000
  ), nrow(sites))
  sites <- sites[rep(seq_len(nrow(sites)),
                     each = length(times %>% unique())
  ), ]
  data_grid_full <- dplyr::tibble(
    sites,
    Age = times
  ) 

  data_age_boundary_max <- data %>%
    dplyr::group_by(SiteName) %>%
    dplyr::slice_max(Age, with_ties = FALSE) %>%
    #dplyr::
    dplyr::reframe(
      #Round here
      max_Age = round(Age + Age_err,digits = 3) 
    ) 
  data_age_boundary_min <- data %>%
    dplyr::group_by(SiteName) %>%
    dplyr::slice_min(Age, with_ties = FALSE) %>%
    dplyr::reframe(
      # Testing--> Fix here
      min_Age = round(Age,digits = 3)
      #min_Age = ifelse(Age>0, Age - Age_err,Age - Age_err)
    )
  # Max and min Ages to replace the prediction grid max and mins
  data_age_boundary <- dplyr::left_join(data_age_boundary_max,
                                        data_age_boundary_min,
                                        by = "SiteName")

  # # Filtering prediction grids to just cover the data
  # data_grid <- data_grid_full %>%
  #   dplyr::left_join(data_age_boundary, by = "SiteName") %>%
  #   dplyr::group_by(SiteName) %>%
  #   dplyr::filter(Age >= (min_Age) & Age <= (max_Age)) %>%
  #   dplyr::tibble() %>%
  #   dplyr::group_by(SiteName) %>%
  #   dplyr::mutate(Age = replace(Age, Age == min(Age), unique(min_Age))) %>%
  #   dplyr::mutate(Age = replace(Age, Age == max(Age), unique(max_Age)))
  
  # Filtering prediction grids to just cover the data
  data_grid <- data_grid_full %>%
    dplyr::left_join(data_age_boundary, by = "SiteName") %>%
    dplyr::group_by(SiteName) %>%
    dplyr::mutate(Age = replace(Age, Age == min(Age), unique(min_Age))) %>%
    dplyr::mutate(Age = replace(Age, Age == max(Age), unique(max_Age))) %>% 
    dplyr::filter(Age >= (min_Age) & Age <= (max_Age))

  
  # # Filtering prediction grids to just cover the data
  # # data_grid_filter <- data_grid_full %>%
  # #   dplyr::left_join(data_age_boundary, by = "SiteName") %>%
  # #   dplyr::group_by(SiteName) %>%
  # #   dplyr::mutate(Age_min_cover = ifelse(min(Age) >= (min_Age), TRUE, FALSE)) %>% 
  # #   dplyr::mutate(Age_max_cover = ifelse(max(Age) <= (max_Age), TRUE, FALSE)) %>% 
  # #   dplyr::tibble()
  #   
  # data_grid<- data_grid_full %>% 
  #   dplyr::left_join(data_age_boundary, by = "SiteName") %>%
  #   dplyr::group_by(SiteName) %>%
  #     # dplyr::mutate(Age = case_when(min(Age) == min(Age) ~ min_Age,
  #     #                               max(Age) == max(Age) ~ max_Age),
  #     #               )
  #   dplyr::mutate(Age = replace(Age, Age == min(Age), unique(min_Age))) %>%
  #   dplyr::mutate(Age = replace(Age, Age == max(Age), unique(max_Age)))
  
  # Ensuring SiteName is a factor
  data <- data %>%
    dplyr::mutate(
      SiteName = as.factor(SiteName),
      data_type_id = as.factor(data_type_id)
    )%>%
    # Arrange by age
    dplyr::arrange(Age)

  data_grid <- data_grid %>%
    dplyr::mutate(
      SiteName = as.factor(SiteName),
      data_type_id = as.factor(data_type_id)
    ) %>%
    dplyr::select(!c(max_Age, min_Age)) %>%
    dplyr::arrange(Age)
  
  # Multiply by 1000 just keep it in right units
  data <- data %>%
    dplyr::mutate(Age = Age*1000, Age_err = Age_err*1000)
  data_grid <- data_grid %>%
    dplyr::mutate(Age = Age*1000)

 
    input_data <- base::list(
      data = data,
      data_grid = data_grid,
      prediction_grid_res = prediction_grid_res
    )
    class(input_data) <- "reslr_input"
  
  return(input_data)
}