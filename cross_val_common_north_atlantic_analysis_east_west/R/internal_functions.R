#' Function to create the dataframes for plotting
#'
#' @param noisy_model_run_output The JAGS output
#' @param rate_grid If rate of change is included in the dataframe
#' @param decomposition Is the full model decomposition included in dataframe
#' @param CI Size of the credible intervals. The default in 0.95.
#' @noRd
create_output_df <- function(noisy_model_run_output,
                             data_grid,
                             CI = 0.95,
                             CI) {
  
  mu_post_pred <- noisy_model_run_output$BUGSoutput$sims.list$mu_pred
  upr <- apply(mu_post_pred, 2, stats::quantile, probs = (1 - CI) / 2)
  lwr <- apply(mu_post_pred, 2, stats::quantile, probs = 1 - ((1 - CI) / 2))
  # Cross Validation
  y_post_pred <- noisy_model_run_output$BUGSoutput$sims.list$y_pred
  upr_PI <- apply(y_post_pred, 2, stats::quantile, probs = (1 - CI) / 2)
  lwr_PI <- apply(y_post_pred, 2, stats::quantile, probs = 1 - ((1 - CI) / 2))
  
  output_dataframes <- data.frame(
    data_grid,
    pred = apply(mu_post_pred, 2, mean),
    upr = upr,
    lwr = lwr,
    y_post_pred = apply(y_post_pred, 2, mean),
    upr_PI = upr_PI,
    lwr_PI = lwr_PI,
    CI = paste0(CI * 100, "%"))
  
  return(output_dataframes)
}
#' A function to create folds for cross validation, recreated from the dismo package to avoid dependency issues.
#'
#' @param x A vector containing the input data
#' @param k The number of folds
#' @param by A vector or factor used to identify sub-groups in the data.
#'
#' @return A vector of indices corresponding to the k-folds.
#' @noRd
kfold_fun <- function(x, k=5, by=NULL) {
  
  singlefold <- function(obs, k) {
    if (k==1) {
      return(rep(1, obs))
    } else {
      i <- obs / k
      if (i < 1) {
        stop('insufficient records:', obs, ', with k=', k)
      }
      i <- round(c(0, i * 1:(k-1), obs))
      times = i[-1] - i[-length(i)]
      
      group <- c()
      for (j in 1:(length(times))) {
        group <- c( group, rep(j, times=times[j]) )
      }
      
      r <- order(stats::runif(obs))
      return(group[r])
    }
  }
  
  if (is.vector(x)) {
    if (length(x) == 1) {
      if (x > 1) {
        x <- 1:x
      }
    }
    obs <- length(x)
  }
  # else if (inherits(x, 'Spatial')) {
  #   if (inherits(x, 'SpatialPoints')) {
  #     obs <- nrow(sf:coordinates(x))
  #   } else {
  #     obs <- nrow(x@data)
  #   }
  # } else {
  #   obs <- nrow(x)
  # }
  if (is.null(by)) {
    return(singlefold(obs, k))
  }
  
  by = as.vector(as.matrix(by))
  if (length(by) != obs) {
    stop('by should be a vector with the same number of records as x')
  }
  un <- unique(by)
  group <- vector(length=obs)
  for ( u in un ) {
    i = which(by==u)
    kk = min(length(i), k)
    if (kk < k) warning('lowered k for by group: ', u  ,'  because the number of observations was  ',  length(i))
    group[i] <- singlefold(length(i), kk)
  }
  return(group)
}

#' If the user decides to include tide gauge data, this function adds the linear rate and the associated linear rate error for those sites.
#' The linear rate comes from a physical model known as an Earth-ice models which use a representation of the physical Earth structure (such as lithospheric thickness and properties such as mantle viscosity) topredict changes in GIA that occur through loading and unloading of ice, and provide estimates of GIA rates.
#' One such example of an Earth-ice physical model is the ICE5G VM2-90 created by Peltier, 2004 and is used as the source of the linear rates for the tide gauges.
#' Engelhart et al 2009 demonstrated the associated uncertainty for the linear rate for tide gauges to be 0.3mm per year.
#' Hence, this function calculates these rates and uncertainty values for tide gauge data.
#'
#' @param data Input data
#' @noRd
add_linear_rate <- function(data) {
  
  # GIA DATA from Peltier Website ICE5G----------------
  # Set up the URL for downloading the data
  #url <- "https://www.atmosp.physics.utoronto.ca/~peltier/datasets/GRID/dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc"
  
  # Create a temporary file
  #temp_file <- tempfile()
  
  # # Download the file and save it to the temporary file
  #  utils::download.file(url,
  #   destfile = temp_file,
  #   method = "libcurl",
  #   mode = "wb",
  #  quiet = TRUE
  # )
  
  # File is stored in the package:
  ice5g_data <- system.file("extdata", "dsea.1grid_O512.nc", package = "reslr")
  # Opening the files
  nc_data <- ncdf4::nc_open(ice5g_data) # ICE5G: better fit for data
  
  
  # Rounding to 1 decimal point to reduce number of spatial options--
  dat_lon <- round(data$Longitude, 1)
  dat_lat <- round(data$Latitude, 1)
  
  # Get lon and lat from GIA model output
  lon <- round(ncdf4::ncvar_get(nc_data, "Lon"), 1)
  lat <- round(ncdf4::ncvar_get(nc_data, "Lat"), 1)
  
  # Needs to be sorted for the match.closest function() below
  # Note need the index for the unsorted coordinate for pulling the correct SL rate later
  gia_lat <- dplyr::tibble(index = 1:length(lat), lat) %>% dplyr::arrange(lat)
  gia_lon <- dplyr::tibble(index = 1:length(lon), lon) %>% dplyr::arrange(lon)
  
  # Matching closest long & lat values
  lat_index <- gia_lat$index[match.closest(dat_lat, gia_lat$lat)]
  lon_index <- gia_lon$index[match.closest(360 + dat_lon, gia_lon$lon)] # change data lon to degrees east
  
  # GIA rates
  SL <- ncdf4::ncvar_get(nc_data, "Dsea_250")
  # Replicating to match dim of data
  linear_slope <- rep(NA, nrow(data))
  for (i in 1:nrow(data)) {
    linear_slope[i] <- (SL[lon_index[i], lat_index[i]])
  }
  
  # Combining linear with other dataset
  data <- cbind(data, ICE5_GIA_slope = linear_slope) # mm/yr
  
  # Remove the temporary file and directory
  #unlink(temp_file)
  #unlink(temp_dir, recursive = TRUE)
  
  return(data)
}

#' Match.closest: Will be used when matching similar long & lat values from both datasets
#'
#' @param x Input data
#' @param table Table of options
#' @param tolerance Identifying the nearby sites
#' @param nomatch No match
#' @noRd
match.closest <- function(x, table, tolerance = Inf, nomatch = NA_integer_) {
  #------Match Closest function doesn't exist on version R 3.6.3----
  #--Will be used when matching similar long & lat values from both datasets--
  lIdx <- findInterval(x, table, rightmost.closed = FALSE, all.inside = TRUE)
  rIdx <- lIdx + 1L
  lIdx[lIdx == 0L] <- 1L
  lDiff <- abs(table[lIdx] - x)
  rDiff <- abs(table[rIdx] - x)
  d <- which(lDiff >= rDiff)
  lIdx[d] <- rIdx[d]
  if (any(is.finite(tolerance))) {
    if (any(tolerance < 0L)) {
      warning(sQuote("tolerance"), " < 0 is meaningless. Set to zero.")
      tolerance[tolerance < 0L] <- 0L
    }
    if (length(nomatch) != 1L) {
      stop("Length of ", sQuote("nomatch"), " has to be one.")
    }
    tolerance <- rep_len(tolerance, length(table))
    lDiff[d] <- rDiff[d]
    lIdx[lDiff > tolerance[lIdx]] <- nomatch
  }
  lIdx
}
#' In this function, tide gauge data from the Permanent Service for Mean Sea Level online database is accessed in a temporary path.
#' The tide gauge data undergo a cleaning process in this function where flagged stations are removed as recommended by the online database.
#' Next, the data is averaged using a rolling window over a decade to ensure it is comparable with proxy data and the tide gauge data is given an RSL uncertainty with is the standard deviation of the data over the decade and an Age error of 5 years corresponding to half a decade.
#' Then, the user selects their preferred tide gauge based on three criteria: 1.nearest tide gauge to the proxy site; 2. User supplies a list of names of preferred tide gauges; 3. all tide gauges within 1 degree are chosen.
#' The tide gauge dataframe is joined with the proxy dataframe with an ID column for data source, "ProxyRecord" or "TideGaugeData"
#'
#' @param data Input data
#' @param list_preferred_TGs user can supply the name or names of the preferred tide gauges
#' @param TG_minimum_dist_proxy The user wants the tide gauge closest to the proxy site
#' @param all_TG_1deg The user wants all tide gauges within 1 degree of the proxy site
#' @param sediment_average_TG Average tide gauge data to make it comparable to accumulation rates of proxy records. The default averaging period for tide gauges is 10 years and the user can alter this.
#' @noRd
clean_tidal_gauge_data <- function(data,
                                   list_preferred_TGs = NULL,
                                   TG_minimum_dist_proxy = FALSE,
                                   all_TG_1deg = FALSE,
                                   sediment_average_TG) {
  Age_epoch_id <- LongLat <- rolling_avg <- median <- nearest_proxy_site <- RSL_annual <- TG_min_dist1 <- minimum_dist <- nearest_TG <- rows_site <- site <- min_dist1 <- stationflag <- name <- sd <- sd_TG <- n_obs_by_site <- RSL_offset <- data_type_id <- decade <- decade_meanRSL <- Age <- RSL <- Age_err <- RSL_err <- linear_rate <- linear_rate_err <- SiteName <- Longitude <- Latitude <- id <- NULL
  # Using data from PSMSL website for annual tide gauge data----------------------------------
  # Set up the URL for downloading the data
  url <- "https://psmsl.org/data/obtaining/rlr.annual.data/rlr_annual.zip"
  
  # Create a temporary file
  temp_file <- tempfile()
  
  # Download the file and save it to the temporary file
  utils::download.file(url,
                       destfile = temp_file,
                       quiet = TRUE)
  
  # Back up if an error
  # # Create a temporary file
  # temp_file <- paste0(tempdir(),'/tg',format(Sys.Date(), "%Y%m%d"))
  # # Download the file and save it to the temporary file
  # if(file.exists(temp_file) | override_download) {
  #   utils::download.file(url,
  #                        destfile = temp_file,
  #                        quiet = TRUE)
  # }
  
  # Unzip the data file to a temporary directory
  temp_dir <- tempfile()
  utils::unzip(temp_file, exdir = temp_dir)
  
  
  
  ### ------------Loop to open all RSL & Age data files------------
  read_plus <- function(flnm) {
    # fread quicker way to read in & allows for ; to be used
    data.table::fread(flnm, sep = ";") %>%
      # allows you to include the file name as id
      dplyr::mutate(filename = flnm)
  }
  # Warnings: there are some files without data
  suppressWarnings(
    temp_SL <-
      list.files(
        path = file.path(temp_dir, "rlr_annual", "data"),
        pattern = "*.rlrdata",
        full.names = T
      ) %>%
      purrr::map_df(~ read_plus(.)) %>%
      dplyr::tibble()
  )
  
  colnames(temp_SL) <- c("Age", "RSL", "flag_attention_1", "flag_attention_2", "id")
  temp_SL$id <- stringr::str_extract(basename(temp_SL$id), "[0-9]+")
  
  # Access the individual data files within the 'rlr_annual' folder
  file_path <- file.path(temp_dir, "rlr_annual", "filelist.txt")
  file_list <- utils::read.csv(file_path, stringsAsFactors = FALSE, header = F, sep = ";")
  colnames(file_list) <- c(
    "id", "Latitude", "Longitude", "name",
    "coastline", "stationcode", "stationflag"
  )
  
  # Removing white space in the name of each site
  file_list$name <- stringr::str_trim(file_list$name, side = "both")
  file_list$stationflag <- stringr::str_trim(file_list$stationflag, side = "both")
  
  # Data from the PSMSL website
  data_TG <- temp_SL %>%
    # Pulling out the file number from string so that it matches the name from other files
    dplyr::mutate(id = stringr::str_extract(basename(temp_SL$id), "[0-9]+")) %>%
    # Cases where bad data was collected
    dplyr::filter(!RSL == -99999) # %>%
  # dplyr::group_by(id) %>%
  # 2000-2018 used as the tidal epoch
  # dplyr::mutate(Age_epoch_id = ifelse(dplyr::between(Age, 2000, 2018), TRUE, FALSE))
  
  # Removing offset based on the location---
  # Offset value is the mean of RSL over the tidal epoch
  # Setting 2000-2018 as the tidal epoch
  # Age_epoch_ref <- data_TG %>%
  #  dplyr::select(RSL, Age_epoch_id) %>%
  #  dplyr::filter(Age_epoch_id == TRUE) %>%
  #  dplyr::summarise(RSL_offset = unique(mean(RSL)))
  # data_TG <- merge(data_TG, Age_epoch_ref, by = "id", all = TRUE)
  # Cases where no data between 2000-2018 set the offset to 7000
  # data_TG$RSL_offset[is.na(data_TG$RSL_offset)] <- 7000
  
  # Updating the RSL to the shifted RSL value using PSMSL instructions
  data_TG$RSL <- data_TG$RSL - 7000 # data_TG$RSL_offset #
  
  #--Joining SL data with location names--
  annual_SL_tide_df <- base::merge(data_TG, file_list, by = "id", all = TRUE)
  #-- Removing sites which have a station flag raised as they are poor sites---
  annual_SL_tide_df <- annual_SL_tide_df %>%
    dplyr::filter(!stationflag == "Y") %>%
    tidyr::drop_na()
  
  # Remove the temporary file and directory
  unlink(temp_file)
  unlink(temp_dir, recursive = TRUE)
  
  # Annual Tidal Gauge data----
  annual_tidal_gauge_data_df <- annual_SL_tide_df %>%
    dplyr::select(
      Age, RSL, Latitude, name,
      # RSL_offset, Age_epoch_id,
      Longitude
    ) %>%
    dplyr::rename(SiteName = name) %>%
    # from mm --> m
    dplyr::mutate(RSL = RSL / 1000) %>%
    # Reordering by group
    dplyr::group_by(SiteName) %>%
    dplyr::arrange(SiteName, .by_group = TRUE) %>%
    dplyr::arrange(Age) %>%
    dplyr::mutate(data_type_id = "TideGaugeData")
  
  
  # Set the window size for the moving average (in this case, 10 years)
  # Rate of sedimentation for proxies when using continuous cores
  window_size <- sediment_average_TG
  
  # # Version A: Create a new column for the decade based on the midpoint of the rolling window
  # annual_tidal_gauge_data_df$rolling_avg <- zoo::rollapply(annual_tidal_gauge_data_df$RSL,
  #   width = window_size,
  #   FUN = mean,
  #   align = "right",
  #   fill = NA
  # )
  # # Calculate the decadal averages based on the rolling average
  # decadal_averages_TG <- annual_tidal_gauge_data_df %>% tidyr::drop_na()
  
  # Version B: Decadal Averages using simple method------
  decadal_averages_TG <-
    annual_tidal_gauge_data_df %>%
    dplyr::mutate(decade = (Age - 1) %/% window_size) %>%
    dplyr::group_by(decade, SiteName) %>%
    dplyr::summarise(
      # decade_meanRSL = mean(RSL)#,
      rolling_avg = mean(RSL),
      Age = max(Age) # ,
      # rows_site = dplyr::n()
    ) # Age=min(Age)
  
  # Using standard deviation of RSL over the decade as uncertainty----
  decadal_averages_TG <- decadal_averages_TG %>%
    dplyr::group_by(SiteName) %>%
    # dplyr::mutate(sd_TG = sd(decade_meanRSL))
    dplyr::mutate(sd_TG = sd(rolling_avg))
  
  #----- New df with decadal averages for tide gauges-----
  tidal_gauge_average_10_df <- base::merge(
    decadal_averages_TG,
    annual_tidal_gauge_data_df
  )
  
  #---Rsl & Age error for tidal gauge data----
  tidal_gauge_full_df <- tidal_gauge_average_10_df %>%
    dplyr::mutate(
      Age_err = 1, # years --> half a year/half a decade
    ) %>%
    # dplyr::mutate(sd_TG = ifelse(is.na(sd_TG), 0.001, sd_TG)) %>%
    dplyr::group_by(SiteName) %>%
    dplyr::mutate(
      RSL_err = 0.001 # sd_TG
      # RSL_err = sd_TG
    )
  
  tidal_gauge_full_df <- tidal_gauge_full_df %>%
    dplyr::mutate(Age = Age / 1000) %>%
    dplyr::mutate(Age_err = Age_err / 1000) %>%
    dplyr::mutate(RSL_annual = RSL) %>%
    # Change for Version A
    dplyr::mutate(RSL = rolling_avg) %>%
    # Change for Version B
    # dplyr::mutate(RSL = decade_meanRSL) %>%
    dplyr::select(!c(
      decade, # Age_epoch_id,RSL_offset,
      rolling_avg,
      RSL_annual
    ))
  
  # No user option here -> this is a must: Removing sites with only 2 points (20 years of data)-----
  decadal_TG_df <-
    tidal_gauge_full_df %>%
    dplyr::group_by(SiteName) %>%
    dplyr::filter(dplyr::n() > 2)
  
  #-----Uniting original dataset and model run to give a site index to model_result data set-----
  SL_site_df <- data %>%
    dplyr::mutate(Longitude = round(Longitude, 1)) %>%
    dplyr::mutate(Latitude = round(Latitude, 1)) %>%
    # Uniting 2 columns
    tidyr::unite("LongLat", Latitude:Longitude, remove = FALSE) %>%
    dplyr::mutate(site = sprintf("%02d", as.integer(as.factor(LongLat)))) %>%
    dplyr::mutate(data_type_id = "ProxyRecord") %>%
    dplyr::group_by(SiteName) %>%
    dplyr::mutate(
      Longitude = dplyr::first(Longitude),
      Latitude = dplyr::first(Latitude)
    ) %>%
    dplyr::group_by(SiteName) %>%
    dplyr::mutate(n_obs_by_site = dplyr::n()) %>%
    dplyr::ungroup()
  
  SL_tide_site_df <- decadal_TG_df %>%
    dplyr::mutate(Longitude = round(Longitude, 1)) %>%
    dplyr::mutate(Latitude = round(Latitude, 1)) %>%
    # Uniting 2 columns
    tidyr::unite("LongLat", Latitude:Longitude, remove = FALSE) %>%
    dplyr::mutate(site = sprintf("%02d", as.integer(as.factor(LongLat)))) %>%
    dplyr::mutate(data_type_id = "TideGaugeData") %>%
    dplyr::group_by(SiteName) %>%
    dplyr::mutate(n_obs_by_site = dplyr::n()) %>%
    dplyr::ungroup()
  
  #------Joining proxy sites to gauges----
  SL_proxy_unique <- SL_site_df %>%
    dplyr::select(SiteName, Longitude, 
                  Latitude, data_type_id, n_obs_by_site, LongLat) %>%
    unique() %>%
    as.data.frame()
  SL_tide_unique <- SL_tide_site_df %>%
    dplyr::select(SiteName, Longitude, 
                  Latitude, data_type_id, n_obs_by_site,LongLat) %>%
    unique() %>%
    as.data.frame()
  
  #---Distance Matrix for each site to each other---
  mat.distance <- geosphere::distm(SL_proxy_unique[, 2:3], SL_tide_unique[, 2:3])
  mat.distance_m <- as.matrix(mat.distance)
  #--finding row mins & corresponding tidal gauge--
  rownames(mat.distance) <- SL_proxy_unique$SiteName
  colnames(mat.distance) <- SL_tide_unique$LongLat#SL_tide_unique$SiteName
  #--finding row mins & corresponding tidal gauge--
  dist_TG_proxy <- t(sapply(seq(nrow(mat.distance)), function(z) {
    js <- order((mat.distance[z, ]))[1:5]
    c(
      rownames(mat.distance)[z], colnames(mat.distance)[js[1]], mat.distance[z, js[1]],
      colnames(mat.distance)[js[2]], mat.distance[z, js[2]],
      colnames(mat.distance)[js[3]], mat.distance[z, js[3]],
      colnames(mat.distance)[js[4]], mat.distance[z, js[4]],
      colnames(mat.distance)[js[5]], mat.distance[z, js[5]]
    )
  }))
  
  dist_TG_proxy <- as.data.frame(dist_TG_proxy)
  colnames(dist_TG_proxy) <- c(
    "nearest_proxy_site",
    "TG_site_1", "TG_min_dist1",
    "TG_site_2", "TG_min_dist2",
    "TG_site_3", "TG_min_dist3",
    "TG_site_4", "TG_min_dist4",
    "TG_site_5", "TG_min_dist5"
  )
  # Sorting the minimum distances from lowest to highest
  dist_TG_proxy <- dist_TG_proxy %>%
    dplyr::arrange(dplyr::desc(TG_min_dist1))
  
  dist_TG_proxy_long_1 <- dist_TG_proxy %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with(c("TG_min_dist")),
      values_to = c("minimum_distance")
    )
  dist_TG_proxy_long_2 <- dist_TG_proxy %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with(c("TG_site")),
      values_to = c("nearest_TG")
    )
  # obs_sites <- SL_tide_unique %>%
  #   dplyr::filter(SiteName %in% dist_TG_proxy_long_2$nearest_TG) %>%
  #   dplyr::select(n_obs_by_site)
  
  dist_TG_proxy_df_no_sitename <- data.frame(
    nearest_proxy_site = dist_TG_proxy_long_1$nearest_proxy_site,
    nearest_TG = dist_TG_proxy_long_2$nearest_TG,
    #nearest_TG_long = SL_tide_unique$LongLat,
    #nearest_proxy_long = SL_proxy_unique$LongLat,
    minimum_dist = as.numeric(dist_TG_proxy_long_1$minimum_distance)
  ) %>% rename(LongLat = nearest_TG)
  
  # Matching the long and lat with name of TG
  long_lat_name_match <- SL_tide_unique %>% 
    dplyr::select(SiteName,LongLat) %>% 
    filter(LongLat %in% dist_TG_proxy_df_no_sitename$LongLat)
  
  dist_TG_proxy_df_new <- left_join(dist_TG_proxy_df_no_sitename,
                                    long_lat_name_match, by = "LongLat",
                                    relationship = "many-to-many") %>% 
    rename(nearest_TG = SiteName)
  
  
  # Criteria 1: User provides a list of TGs------------------------
  # If list true: combining all dataframes & removing duplicates----
  if (TG_minimum_dist_proxy == FALSE & is.null(list_preferred_TGs) == FALSE & all_TG_1deg == FALSE) {
    # Check if TG exists in the list
    check_TG <- all(list_preferred_TGs %in% unique(decadal_TG_df$SiteName))
    if (check_TG == FALSE) {
      message("Warning: Tide Gauge provided does not exist or may contain a misprint in the name.\n")
      stop()
    }
    
    decadal_TG_df_filter <- base::subset(decadal_TG_df, SiteName %in% list_preferred_TGs)
    # There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_list <- plyr::rbind.fill(
      SL_site_df,
      decadal_TG_df_filter
    )
    data_tide_proxy <-data_tide_proxy_list
  }
  # Criteria 2: Minimum distance to proxy site----------------------
  if (TG_minimum_dist_proxy == TRUE & is.null(list_preferred_TGs) == TRUE & all_TG_1deg == FALSE) {
    
    # Finding the closest TG
    all_nearest_TG_closest <- dist_TG_proxy_df_new %>%
      dplyr::group_by(nearest_proxy_site) %>%
      dplyr::filter(minimum_dist == min(minimum_dist)) %>%
      # Removing any duplicate tide gauge sites.
      dplyr::distinct(nearest_TG, .keep_all = TRUE)
    
    # Joining the selected TG sites back with the original data
    join_new_index_tide_df <- SL_tide_site_df %>%
      #dplyr::filter(SiteName %in% all_nearest_TG_closest$nearest_TG)
      dplyr::filter(LongLat %in% all_nearest_TG_closest$LongLat)
    
    #--There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_min_TG <- plyr::rbind.fill(
      SL_site_df,
      # stacking rows
      join_new_index_tide_df
    )
    data_tide_proxy <- data_tide_proxy_min_TG
  }
  
  
  # Criteria 3: All tide gauges within 1 degree away from proxy site-------------
  if (TG_minimum_dist_proxy == FALSE & is.null(list_preferred_TGs) == TRUE & all_TG_1deg == TRUE) {
    # 1 degree away from proxy site is 111.1km
    all_nearest_TG_closest <- dist_TG_proxy_df_new %>%
      dplyr::filter(minimum_dist <= 111100) %>%
      # Removing any duplicate tide gauge sites.
      dplyr::distinct(nearest_TG, .keep_all = TRUE)
    
    # Joining the selected TG sites back with the original data
    join_new_index_tide_df <- SL_tide_site_df %>%
      dplyr::filter(LongLat %in% all_nearest_TG_closest$LongLat)
    #dplyr::filter(SiteName %in% all_nearest_TG_closest$nearest_TG)
    
    #--There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_all_TG <- plyr::rbind.fill(
      SL_site_df,
      join_new_index_tide_df
    ) # stacking rows
    
    data_tide_proxy <- data_tide_proxy_all_TG
  }
  
  # If all options true: combining all dataframes & removing duplicates----
  if (TG_minimum_dist_proxy == TRUE & is.null(list_preferred_TGs) == FALSE & all_TG_1deg == TRUE) {
    
    # Criteria 1: List of names - Check if TG exists in the list
    check_TG <- all(list_preferred_TGs %in% unique(decadal_TG_df$SiteName))
    if (check_TG == FALSE) {
      message("Warning: Tide Gauge provided does not exist or may contain a misprint in the name.\n")
      stop()
    }
    
    decadal_TG_df_filter <- base::subset(decadal_TG_df, SiteName %in% list_preferred_TGs)
    # There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_list <- plyr::rbind.fill(
      SL_site_df,
      decadal_TG_df_filter
    )
    
    # Criteria 2: Minimum distance to proxy site - Finding the closest TG
    all_nearest_TG_closest <- dist_TG_proxy_df_new %>%
      dplyr::group_by(nearest_proxy_site) %>%
      dplyr::filter(minimum_dist == min(minimum_dist)) %>%
      # Removing any duplicate tide gauge sites.
      dplyr::distinct(nearest_TG, .keep_all = TRUE)
    
    # Joining the selected TG sites back with the original data
    join_new_index_tide_df <- SL_tide_site_df %>%
      #dplyr::filter(SiteName %in% all_nearest_TG_closest$nearest_TG)
      dplyr::filter(LongLat %in% all_nearest_TG_closest$LongLat)
    
    #--There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_min_TG <- plyr::rbind.fill(
      SL_site_df,
      # stacking rows
      join_new_index_tide_df
    )
    # Criteria 3:  All TG - 1 degree away from proxy site is 111.1km
    all_nearest_TG_closest <- dist_TG_proxy_df_new %>%
      dplyr::filter(minimum_dist <= 111100) %>%
      # Removing any duplicate tide gauge sites.
      dplyr::distinct(nearest_TG, .keep_all = TRUE)
    
    # Joining the selected TG sites back with the original data
    join_new_index_tide_df <- SL_tide_site_df %>%
      #dplyr::filter(SiteName %in% all_nearest_TG_closest$nearest_TG)
      dplyr::filter(LongLat %in% all_nearest_TG_closest$LongLat)
    
    #--There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_all_TG <- plyr::rbind.fill(
      SL_site_df,
      # stacking rows
      join_new_index_tide_df
    ) 
    
    
    # If you have multiple options together will need to merge them 
    data_tide_proxy <- base::merge(data_tide_proxy_all_TG,
                                   data_tide_proxy_min_TG,all = TRUE)
    # Merge only works with 2 entries
    data_tide_proxy <- base::merge(data_tide_proxy, data_tide_proxy_list,all = TRUE)
  }
  
  # If TG_min & list true: combining all dataframes & removing duplicates----
  if (TG_minimum_dist_proxy == TRUE & is.null(list_preferred_TGs) == FALSE & all_TG_1deg == FALSE) {
    
    # Criteria 1: List of names - Check if TG exists in the list
    check_TG <- all(list_preferred_TGs %in% unique(decadal_TG_df$SiteName))
    if (check_TG == FALSE) {
      message("Warning: Tide Gauge provided does not exist or may contain a misprint in the name.\n")
      stop()
    }
    
    decadal_TG_df_filter <- base::subset(decadal_TG_df, SiteName %in% list_preferred_TGs)
    # There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_list <- plyr::rbind.fill(
      SL_site_df,
      decadal_TG_df_filter
    )
    
    # Criteria 2: Minimum distance to proxy site - Finding the closest TG
    all_nearest_TG_closest <- dist_TG_proxy_df_new %>%
      dplyr::group_by(nearest_proxy_site) %>%
      dplyr::filter(minimum_dist == min(minimum_dist)) %>%
      # Removing any duplicate tide gauge sites.
      dplyr::distinct(nearest_TG, .keep_all = TRUE)
    
    # Joining the selected TG sites back with the original data
    join_new_index_tide_df <- SL_tide_site_df %>%
      #dplyr::filter(SiteName %in% all_nearest_TG_closest$nearest_TG)
      dplyr::filter(LongLat %in% all_nearest_TG_closest$LongLat)
    
    #--There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_min_TG <- plyr::rbind.fill(
      SL_site_df,
      # stacking rows
      join_new_index_tide_df
    )
    # Multiple join together
    data_tide_proxy <- base::merge(
      data_tide_proxy_min_TG,
      data_tide_proxy_list, all = TRUE)
  }
  
  # If TG_min & all tg true: combining all dataframes & removing duplicates----
  if (TG_minimum_dist_proxy == TRUE & is.null(list_preferred_TGs) == TRUE & all_TG_1deg == TRUE) {
    # Criteria 2: Minimum distance to proxy site -  Finding the closest TG
    all_nearest_TG_closest <- dist_TG_proxy_df_new %>%
      dplyr::group_by(nearest_proxy_site) %>%
      dplyr::filter(minimum_dist == min(minimum_dist)) %>%
      # Removing any duplicate tide gauge sites.
      dplyr::distinct(nearest_TG, .keep_all = TRUE)
    
    # Joining the selected TG sites back with the original data
    join_new_index_tide_df <- SL_tide_site_df %>%
      #dplyr::filter(SiteName %in% all_nearest_TG_closest$nearest_TG)
      dplyr::filter(LongLat %in% all_nearest_TG_closest$LongLat)
    
    #--There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_min_TG <- plyr::rbind.fill(
      SL_site_df,
      # stacking rows
      join_new_index_tide_df
    )
    # Criteria 3:  All TG - 1 degree away from proxy site is 111.1km
    all_nearest_TG_closest <- dist_TG_proxy_df_new %>%
      dplyr::filter(minimum_dist <= 111100) %>%
      # Removing any duplicate tide gauge sites.
      dplyr::distinct(nearest_TG, .keep_all = TRUE)
    
    # Joining the selected TG sites back with the original data
    join_new_index_tide_df <- SL_tide_site_df %>%
      #dplyr::filter(SiteName %in% all_nearest_TG_closest$nearest_TG)
      dplyr::filter(LongLat %in% all_nearest_TG_closest$LongLat)
    
    #--There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_all_TG <- plyr::rbind.fill(
      SL_site_df,
      join_new_index_tide_df
    ) # stacking rows
    
    # Joining all together
    data_tide_proxy <- base::merge(data_tide_proxy_all_TG,
                                   data_tide_proxy_min_TG,all = TRUE)
  }
  
  
  # If list & all tg true: combining all dataframes & removing duplicates----
  if (TG_minimum_dist_proxy == FALSE & is.null(list_preferred_TGs) == FALSE & all_TG_1deg == TRUE) {
    
    # Criteria 3: All tide gauges within 1 degree away from proxy site
    # 1 degree away from proxy site is 111.1km
    all_nearest_TG_closest <- dist_TG_proxy_df_new %>%
      dplyr::filter(minimum_dist <= 111100) %>%
      # Removing any duplicate tide gauge sites.
      dplyr::distinct(nearest_TG, .keep_all = TRUE)
    
    # Joining the selected TG sites back with the original data
    join_new_index_tide_df <- SL_tide_site_df %>%
      #dplyr::filter(SiteName %in% all_nearest_TG_closest$nearest_TG)
      dplyr::filter(LongLat %in% all_nearest_TG_closest$LongLat)
    
    #--There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_all_TG <- plyr::rbind.fill(
      SL_site_df,
      join_new_index_tide_df
    ) # stacking rows
    
    # Criteria 1: List of names -  Check if TG exists in the list
    check_TG <- all(list_preferred_TGs %in% unique(decadal_TG_df$SiteName))
    if (check_TG == FALSE) {
      message("Warning: Tide Gauge provided does not exist or may contain a misprint in the name.\n")
      stop()
    }
    
    decadal_TG_df_filter <- base::subset(decadal_TG_df, SiteName %in% list_preferred_TGs)
    # There will be NAs were the proxy data doesn't have a corresponding index--
    data_tide_proxy_list <- plyr::rbind.fill(
      SL_site_df,
      decadal_TG_df_filter
    )
    
    # If you have multiple options together will need to merge them 
    data_tide_proxy <- base::merge(data_tide_proxy_all_TG,
                                   data_tide_proxy_list,all = TRUE)
  }
  
  
  # If all true: combining all dataframes & removing duplicates----
  # if (TG_minimum_dist_proxy == FALSE & is.null(list_preferred_TGs) == FALSE & all_TG_1deg == FALSE) {
  #   stop("Error: No tide gauge selection method chosen. Select criteria to chose your prefered tide gauge")
  # }
  
  # Ensuring the SiteName is a factor
  data <- data_tide_proxy %>%
    dplyr::select(!c(
      n_obs_by_site, site, sd_TG
    )) %>%
    dplyr::mutate(
      SiteName = as.factor(SiteName),
      data_type_id = as.factor(data_type_id)
    )
  return(data)
}
#' Linear rate estimated using the data
#'
#' @param data Input data
#' @noRd
linear_reg_rates <- function(data) {
  Age <- RSL <- Age_err <- RSL_err <- linear_rate <- linear_rate_err <- SiteName <- Longitude <- Latitude <- NULL
  data_filter <- data %>%
    # Ignoring recent human influences to SL rise
    dplyr::filter(!Age > 1.800)
  # # Remove index points
  # data_filter <- data %>%
  #   dplyr::group_by(Site)%>%
  #   dplyr::filter(dplyr::n()>1)
  
  # Doing linear regression on rest of data
  data_lm <- data_filter %>%
    dplyr::group_by(SiteName) %>%
    dplyr::mutate(
      # linear_rate = round(stats::lm(RSL ~ Age)$coefficients[[2]], 2),
      linear_rate = stats::lm(RSL ~ Age)$coefficients[[2]],
      linear_rate_err = base::summary(stats::lm(RSL ~ Age))$coefficients[2, 2]
    )
  
  
  # Table of GIA rate vs lm rate from proxy data
  lm_slopes <- data_lm %>%
    dplyr::select(SiteName, linear_rate, linear_rate_err) %>%
    unique()
  return(lm_slopes)
}
#' Adding Noisy Input to the dataframe
#'
#' @param model_run JAGS output
#' @param data Input data
#' @param data_grid Input data grid
#' @param spline_nseg_t Number of segments used to create basis functions for NIGAM temporal component
#' @param spline_nseg_st Number of segments used to create basis functions for NIGMA spatial temporal component
add_noisy_input <- function(model_run,
                            jags_data,
                            data,
                            data_grid,
                            spline_nseg_t,
                            spline_nseg_st,
                            spline_nseg_c) {
    #-----Get posterior samples for SL-----
    intercept_post <- model_run$BUGSoutput$sims.list$intercept
    b_c_post <- model_run$BUGSoutput$sims.list$b_c
    b_g_post <- model_run$BUGSoutput$sims.list$b_g
    
    pred_mean_calc <- function(t_new) {
      # Create the regional basis functions
      B_c <- bs_bbase_c(t_new,
                        xl = min(data$Age),
                        xr = max(data$Age),
                        spline_nseg_c = spline_nseg_c,
                        data = data
      )
      #----Deriv----
      return(intercept_post[data$SiteName] + B_c %*% colMeans(b_c_post)+ b_g_post[data$SiteName] * (t_new))
  
    }
    #-------Now create derivatives---
    h <- 0.0000001#0.01
    t <- data$Age
    deriv <- (pred_mean_calc(t + h) - pred_mean_calc(t - h)) / (2 * h)
    # Predicted data
    pred_mean_calc_grid <- function(t_new) {
      # Create the regional basis functions
      B_c <- bs_bbase_c(t_new,
                        xl = min(data$Age),
                        xr = max(data$Age),
                        spline_nseg_c = spline_nseg_c,
                        data = data
      )
      #----Deriv----
      return(intercept_post[data_grid$SiteName] + B_c %*% colMeans(b_c_post) + b_g_post[data_grid$SiteName] * (t_new))
    }
    t_grid <- data_grid$Age
    deriv_grid <- (pred_mean_calc_grid(t_grid + h) - pred_mean_calc_grid(t_grid - h)) / (2 * h)
  
  # Add this new term in - this is the extra standard deviation on each term----
  data$NI_var_term <- sqrt(deriv^2 %*% data$Age_err^2)[, 1]
  data_grid$NI_var_grid_term <- sqrt(deriv_grid^2 %*% data$Age_err^2)[, 1]
  
  # Writing new dataframe with noisy extra column------
  data <- data.frame(data)
  data_grid <- data.frame(data_grid)
  update_input_df <- list(data = data,
                          data_grid = data_grid)
  return(update_input_df)
}
  

#' Creating basis function for splines
#'
#' @param data Input data
#' @param data_grid Prediction data
#' @param model_type Type of model
#' @param spline_nseg_t Number of segments for the creation of the basis functions
#' @param spline_nseg_t Number of segments for the creation of the basis functions for NIGAM temporal component
#' @param spline_nseg_st Number of segments for the creation of the basis functions for NIGAM spatial temporal component
#' @param xl Minimum value for the basis function for splines
#' @param xr Maximum value for the basis function for splines
#' @noRd


spline_basis_fun <- function(data,
                             data_grid,
                             spline_nseg,
                             spline_nseg_c,
                             spline_nseg_st,
                             spline_nseg_t,
                             xl,
                             xr) {
  
  # Basis functions in time for common component-----------------------
  B_c <- bs_bbase_c(data$Age,
                    xl = min(data$Age),
                    xr = max(data$Age),
                    spline_nseg_c = spline_nseg_c,
                    data = data
  )
  # Finding derivative  of basis functions using first principals-----------
  first_deriv_calc <- function(t_new) {
    # Create the regional basis functions
    B_c <- bs_bbase_c(t_new,
                      xl = min(data$Age),
                      xr = max(data$Age),
                      spline_nseg_c = spline_nseg_c,
                      data = data
    )
    return(B_c)
  }
  # Now create derivatives----------------------
  h <- 0.0000001#0.0000001#0.00001 # h <- 0.001
  t <- data$Age
  first_deriv_step1 <- first_deriv_calc(t + h)
  first_deriv_step2 <- first_deriv_calc(t - h)
  B_c_deriv <- (first_deriv_step1 - first_deriv_step2) / (2 * h)
  
  # Basis functions in time using prediction data frame-----------------------
  B_c_pred <- bs_bbase_c(data_grid$Age,
                         xl = min(data$Age),
                         xr = max(data$Age),
                         spline_nseg_c = spline_nseg_c,
                         data = data
  )
  # Now create derivatives----------------------
  h <- 0.0000001#0.00001 # h <- 0.001
  t_pred <- data_grid$Age
  first_deriv_step1 <- first_deriv_calc(t_pred + h)
  first_deriv_step2 <- first_deriv_calc(t_pred - h)
  B_c_pred_deriv <- (first_deriv_step1 - first_deriv_step2) / (2 * h)
  
  # Basis functions in time for regional component-----------------------
    B_t <- bs_bbase_t(data$Age,
                      xl = min(data$Age),
                      xr = max(data$Age),
                      spline_nseg_t = spline_nseg_t,
                      data = data
    )
    # Finding derivative  of basis functions using first principals-----------
    first_deriv_calc <- function(t_new) {
      # Create the regional basis functions
      B_t <- bs_bbase_t(t_new,
                        xl = min(data$Age),
                        xr = max(data$Age),
                        spline_nseg_t = spline_nseg_t,
                        data = data
      )
      return(B_t)
    }
    # Now create derivatives----------------------
    h <- 0.0000001#0.00001 # h <- 0.001
    t <- data$Age
    first_deriv_step1 <- first_deriv_calc(t + h)
    first_deriv_step2 <- first_deriv_calc(t - h)
    B_t_deriv <- (first_deriv_step1 - first_deriv_step2) / (2 * h)
    
    # Basis functions in time using prediction data frame-----------------------
    B_t_pred <- bs_bbase_t(data_grid$Age,
                           xl = min(data$Age),
                           xr = max(data$Age),
                           spline_nseg_t = spline_nseg_t,
                           data = data
    )
    # Now create derivatives----------------------
    h <- 0.0000001#0.00001 # h <- 0.001
    t_pred <- data_grid$Age
    first_deriv_step1 <- first_deriv_calc(t_pred + h)
    first_deriv_step2 <- first_deriv_calc(t_pred - h)
    B_t_pred_deriv <- (first_deriv_step1 - first_deriv_step2) / (2 * h)
    
    # Basis functions in space time for data-----------------------
    B_time <- bs_bbase_st(data$Age,
                          xl = min(data$Age),
                          xr = max(data$Age),
                          spline_nseg_st = spline_nseg_st,
                          data = data
    )
    B_space_1 <- bs_bbase_st(data$Latitude,
                             xl = min(data$Latitude),
                             xr = max(data$Latitude),
                             spline_nseg_st = spline_nseg_st,
                             data = data
    )
    B_space_2 <- bs_bbase_st(data$Longitude,
                             xl = min(data$Longitude),
                             xr = max(data$Longitude),
                             spline_nseg_st = spline_nseg_st,
                             data = data
    )
    
    B_st_full <- matrix(NA,
                        ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1),
                        nrow = nrow(data)
    )
    regional_knots_loc <- rep(NA,
                              ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1)
    )
    count <- 1
    for (i in 1:ncol(B_time)) {
      for (j in 1:ncol(B_space_1)) {
        for (k in 1:ncol(B_space_2)) {
          regional_knots_loc[count] <- i
          B_st_full[, count] <- B_time[, i] * B_space_1[, j] * B_space_2[, k]
          count <- count + 1
        }
      }
    }
    
    # Get rid of all the columns which are just zero
    B_st <- B_st_full[, -which(colSums(B_st_full) < 0.1)]
    
    # Find the index here that you remove then use this in the derivative
    remove_col_index <- which(colSums(B_st_full) < 0.1)
    
    first_deriv_calc <- function(t_new) {
      # Now the local basis functions
      B_time <- bs_bbase_st(t_new,
                            xl = min(data$Age),
                            xr = max(data$Age),
                            spline_nseg_st = spline_nseg_st,
                            data = data
      )
      B_space_1 <- bs_bbase_st(data$Latitude,
                               xl = min(data$Latitude),
                               xr = max(data$Latitude),
                               spline_nseg_st = spline_nseg_st,
                               data = data
      )
      B_space_2 <- bs_bbase_st(data$Longitude,
                               xl = min(data$Longitude),
                               xr = max(data$Longitude),
                               spline_nseg_st = spline_nseg_st,
                               data = data
      )
      
      B_st_full <- matrix(NA,
                          ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1),
                          nrow = nrow(data)
      )
      regional_knots_loc <- rep(NA, ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1))
      count <- 1
      for (i in 1:ncol(B_time)) {
        for (j in 1:ncol(B_space_1)) {
          for (k in 1:ncol(B_space_2)) {
            regional_knots_loc[count] <- i
            B_st_full[, count] <- B_time[, i] * B_space_1[, j] * B_space_2[, k]
            count <- count + 1
          }
        }
      }
      
      # Get rid of all the columns which are just zero
      B_st <- B_st_full[, -which(colSums(B_st_full) < 0.1)]
      return(B_st)
    }
    #-------Now create derivatives----
    h <- 0.0000001#0.00001
    t <- data$Age
    
    first_deriv_step1 <- first_deriv_calc(t + h)
    first_deriv_step2 <- first_deriv_calc(t - h)
    B_st_deriv <- (first_deriv_step1 - first_deriv_step2) / (2 * h)
    
    
    # Basis functions in space time using prediction data frame-----------------------
    B_pred_time <- bs_bbase_st(data_grid$Age,
                               xl = min(data$Age),
                               xr = max(data$Age),
                               spline_nseg_st = spline_nseg_st,
                               data = data
    )
    B_space_1 <- bs_bbase_st(data_grid$Latitude,
                             xl = min(data$Latitude),
                             xr = max(data$Latitude),
                             spline_nseg_st = spline_nseg_st,
                             data = data
    )
    B_space_2 <- bs_bbase_st(data_grid$Longitude,
                             xl = min(data$Longitude),
                             xr = max(data$Longitude),
                             spline_nseg_st = spline_nseg_st,
                             data = data
    )
    
    suppressWarnings({
      B_st_pred_full <- matrix(NA,
                               ncol = ncol(B_pred_time) * ncol(B_space_1) * ncol(B_space_1),
                               nrow = nrow(data_grid)
      )
      regional_knots_loc <- rep(NA,
                                ncol = ncol(B_pred_time) * ncol(B_space_1) * ncol(B_space_1)
      )
      count <- 1
      for (i in 1:ncol(B_pred_time)) {
        for (j in 1:ncol(B_space_1)) {
          for (k in 1:ncol(B_space_2)) {
            regional_knots_loc[count] <- i
            B_st_pred_full[, count] <- B_pred_time[, i] * B_space_1[, j] * B_space_2[, k]
            count <- count + 1
          }
        }
      }
      
      # Get rid of all the columns which are just zero corresponding to the previous basis functions
      B_st_pred <- B_st_pred_full[, -remove_col_index]
    })
    #-------Now create derivatives for prediciton----
    first_deriv_calc <- function(t_new) {
      # Now the local basis functions
      B_time <- bs_bbase_st(t_new,
                            xl = min(data$Age),
                            xr = max(data$Age),
                            spline_nseg_st = spline_nseg_st,
                            data = data
      )
      B_space_1 <- bs_bbase_st(data_grid$Latitude,
                               xl = min(data$Latitude),
                               xr = max(data$Latitude),
                               spline_nseg_st = spline_nseg_st,
                               data = data
      )
      B_space_2 <- bs_bbase_st(data_grid$Longitude,
                               xl = min(data$Longitude),
                               xr = max(data$Longitude),
                               spline_nseg_st = spline_nseg_st,
                               data = data
      )
      
      B_st_full <- matrix(NA,
                          ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1),
                          nrow = nrow(data_grid)
      )
      regional_knots_loc <- rep(NA, ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1))
      count <- 1
      for (i in 1:ncol(B_time)) {
        for (j in 1:ncol(B_space_1)) {
          for (k in 1:ncol(B_space_2)) {
            regional_knots_loc[count] <- i
            B_st_full[, count] <- B_time[, i] * B_space_1[, j] * B_space_2[, k]
            count <- count + 1
          }
        }
      }
      
      # Get rid of all the columns which are just zero
      B_st <- B_st_full[, -remove_col_index]
      return(B_st)
    }
    h <- 0.0000001#0.00001
    t_pred <- data_grid$Age
    
    first_deriv_step1 <- first_deriv_calc(t_pred + h)
    first_deriv_step2 <- first_deriv_calc(t_pred - h)
    B_st_deriv_pred <- (first_deriv_step1 - first_deriv_step2) / (2 * h)
    
    # Derivative of Basis function for the total model fit-----------------
    first_deriv_calc <- function(t_new) {
      # Create the regional basis functions
      B_t <- bs_bbase_t(t_new,
                        xl = min(data$Age),
                        xr = max(data$Age),
                        spline_nseg_t = spline_nseg_t,
                        data = data
      )
      colnames(B_t) <- c(paste("B_t", 1:ncol(B_t), sep = ""))
      
      # Now the local basis functions
      B_time <- bs_bbase_st(t_new,
                            xl = min(data$Age),
                            xr = max(data$Age),
                            spline_nseg_st = spline_nseg_st,
                            data = data
      )
      B_space_1 <- bs_bbase_st(data$Latitude,
                               xl = min(data$Latitude),
                               xr = max(data$Latitude),
                               spline_nseg_st = spline_nseg_st,
                               data = data
      )
      B_space_2 <- bs_bbase_st(data$Longitude,
                               xl = min(data$Longitude),
                               xr = max(data$Longitude),
                               spline_nseg_st = spline_nseg_st,
                               data = data
      )
      
      B_st_full <- matrix(NA,
                          ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1),
                          nrow = nrow(data)
      )
      regional_knots_loc <- rep(NA, ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1))
      count <- 1
      for (i in 1:ncol(B_time)) {
        for (j in 1:ncol(B_space_1)) {
          for (k in 1:ncol(B_space_2)) {
            regional_knots_loc[count] <- i
            B_st_full[, count] <- B_time[, i] * B_space_1[, j] * B_space_2[, k]
            count <- count + 1
          }
        }
      }
      
      # Get rid of all the columns which are just zero
      B_st <- B_st_full[, -which(colSums(B_st_full) < 0.1)]
      colnames(B_st) <- c(paste("B_st", 1:ncol(B_st), sep = ""))
      # Dummy matrix for intercept & GIA
      B_h <- fastDummies::dummy_cols(data.frame(data$SiteName))
      B_h <- B_h[, -1]
      colnames(B_h) <- c(paste("B_h", 1:ncol(B_h), sep = ""))
      B_g <- B_h * t_new
      colnames(B_g) <- c(paste("B_g", 1:ncol(B_g), sep = ""))
      # Basis function matrix with B_local & B_regional
      output_B_tot <- cbind(
        B_h, B_g,
        B_t, B_st
      )
      
      
      return(output_B_tot)
    }
    #-------Now create derivatives----
    h <- 0.0000001#0.00001
    t <- data$Age
    
    first_deriv_step1 <- first_deriv_calc(t + h)
    first_deriv_step2 <- first_deriv_calc(t - h)
    B_tot_deriv <- (first_deriv_step1 - first_deriv_step2) / (2 * h)
    
    # New basis function for site specific vertical offset----------
    B_h_deriv <- as.matrix(B_tot_deriv[, grepl("B_h", names(B_tot_deriv))])
    # New Basis Functions for Random linear local component-Site Specific slope-------
    B_deriv_g_re <- as.matrix(B_tot_deriv[, grepl("B_g", names(B_tot_deriv))])
    # New Basis Functions for spline in time-------------------------
    B_deriv_t <- as.matrix(B_tot_deriv[, grepl("B_t", names(B_tot_deriv))])
    # New Basis Functions in Space time-------------------------------
    B_deriv_st <- as.matrix(B_tot_deriv[, grepl("B_st", names(B_tot_deriv))])
    
    # All Basis Functions
    spline_basis_fun_list <- list(
      B_c = B_c,
      B_c_pred=B_c_pred,
      B_c_deriv = B_c_deriv,
      B_c_pred_deriv = B_c_pred_deriv,
      remove_col_index = remove_col_index,
      B_st = B_st,
      B_st_deriv = B_st_deriv,
      B_st_pred = B_st_pred,
      B_st_deriv_pred = B_st_deriv_pred,
      B_t = B_t,
      B_t_deriv = B_t_deriv,
      B_t_pred = B_t_pred,
      B_t_pred_deriv = B_t_pred_deriv,
      B_h_deriv = B_h_deriv,
      B_deriv_g_re = B_deriv_g_re,
      B_tot_deriv = B_tot_deriv,
      B_deriv_t = B_deriv_t,
      B_deriv_st = B_deriv_st
    )
  
  
  
  return(spline_basis_fun_list)
}

#' Creating spline basis functions for NIGAM for the common component
# Basis function approach
bs_bbase_c <- function(x,
                       xl =  min(x),
                       xr = max(x),
                       deg = 3,
                       spline_nseg_c = NULL,
                       data = data) {
  # Create basis functions------------------------------------------------------
  if (is.null(spline_nseg_c)) {
    spline_nseg_c <- round(deg / (1 + deg / length(data$Age)))
  }
  
  # Compute the length of the partitions
  dx <- (xr - xl) / spline_nseg_c
  # Create equally spaced knots
  knots <- seq(xl - deg * dx,
               xr + deg * dx,
               by = dx
  )
  
  # Use bs() function to generate the B-spline basis
  bs_matrix <- matrix(
    splines::bs(x,
                degree = deg,
                knots = knots,
                Boundary.knots = c(knots[1], knots[length(knots)])
    ),
    nrow = length(x)
  )
  return(bs_matrix)
}


#' Creating spline basis functions for NIGAM for the temporal component
# Basis function approach
bs_bbase_t <- function(x,
                       xl =  min(x),
                       xr = max(x),
                       deg = 3,
                       spline_nseg_t = NULL,
                       data = data) {
  # Create basis functions------------------------------------------------------
  if (is.null(spline_nseg_t)) {
    spline_nseg_t <- round(deg / (1 + deg / length(data$Age)))
  }

  # Compute the length of the partitions
  dx <- (xr - xl) / spline_nseg_t
  # Create equally spaced knots
  knots <- seq(xl - deg * dx,
               xr + deg * dx,
               by = dx
  )
  
  # Use bs() function to generate the B-spline basis
  bs_matrix <- matrix(
    splines::bs(x,
                degree = deg,
                knots = knots,
                Boundary.knots = c(knots[1], knots[length(knots)])
    ),
    nrow = length(x)
  )
  return(bs_matrix)
}

#' Creating spline basis functions for NIGAM spatial temporal component
# Basis function approach
bs_bbase_st <- function(x,
                        xl = min(x),
                        xr = max(x),
                        deg = 3,
                        spline_nseg_st = NULL,
                        data = data) {
  # Create basis functions------------------------------------------------------
  if (is.null(spline_nseg_st)) {
    spline_nseg_st <- round(deg / (1 + deg / length(data$Age)))
  }
  
  # Compute the length of the partitions
  dx <- (xr - xl) / spline_nseg_st
  # Create equally spaced knots
  knots <- seq(xl - deg * dx,
               xr + deg * dx,
               by = dx
  )
  
  # Use bs() function to generate the B-spline basis
  bs_matrix <- matrix(
    splines::bs(x,
                degree = deg,
                knots = knots,
                Boundary.knots = c(knots[1], knots[length(knots)])
    ),
    nrow = length(x)
  )
  return(bs_matrix)
}