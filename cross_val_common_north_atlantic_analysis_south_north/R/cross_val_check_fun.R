#' Cross validation check
cross_val_check_fun <- function(data, data_grid,
                                prediction_grid_res = 30,
                                spline_nseg_c = 8,
                                spline_nseg_t = 20,
                                spline_nseg_st = 6,
                                n_iterations = 10,
                                n_burnin = 10,
                                n_thin = 1,
                                n_chains = 3,
                                n_fold = 5,
                                seed = NULL) {

  # Random selection
  base::set.seed(seed)

  # Checking the number of observations > n_folds and give a warning
  # if(any(nrow_by_site$nrow_site < n_fold) == TRUE){
  # stop("Not enough observations in each site for the number of folds required.")
  # }

  # Checking the number of observations > n_folds & separating them for dataset
  nrow_by_site <- data %>%
    group_by(SiteName) %>%
    reframe(nrow_site = n()) %>%
    filter(nrow_site < n_fold + 2)
  # Put sites with less than 12 observations into the training set
  more_input_data <- data %>%
    filter(SiteName %in% nrow_by_site$SiteName) %>%
    mutate(CV_fold = "None")

  # Checking the number of observations > n_folds & separating them in data grid
  nrow_by_site_grid <- data_grid %>%
    group_by(SiteName) %>%
    reframe(nrow_site = n()) %>%
    filter(nrow_site < n_fold + 2)
  # Put sites with less than 12 observations into the training set
  more_input_data_grid <- data_grid %>%
    filter(SiteName %in% nrow_by_site$SiteName)

  # Filtering the sites with less than the number of folds but will put them back into training
  data_filt <- anti_join(data, more_input_data)
  data_filt_grid <- anti_join(data_grid, more_input_data_grid)

  # Setting up my fold index
  df_split_index <- kfold_fun(data_filt$Age,
    k = n_fold,
    by = data_filt$SiteName
  )
  # Is this in the wrong place?
  data_filt$CV_fold <- df_split_index

  # Empty list for model runs
  model_run_list <- list()
  for (i in 1:n_fold) {
    # Segment your data by fold using the which() function
    CV_fold <- base::which(df_split_index == i, arr.ind = TRUE)

    # Remove SiteName column as it causes errors in next step
    data_red <- data_filt %>% dplyr::select(!SiteName)
    # Test set--------
    test_set_full <- data_red[CV_fold, ]
    test_set<- test_set_full
    # Including the tide gauge sites aswell
    #more_input_data_clean <- more_input_data %>% dplyr::select(!SiteName)
    #test_set <- rbind(test_set_full, more_input_data_clean)

    # Training-----
    training_set_full <- data_red[-CV_fold, ]
    training_set <- training_set_full
    # Include extra tg data here
    #more_input_data_clean <- more_input_data %>% dplyr::select(!SiteName)
    #training_set <- rbind(training_set_full, more_input_data_clean)

    # reslr_load
    input_train <- reslr_load_fun(
      data = training_set,
      prediction_grid_res = prediction_grid_res,
      cross_val = TRUE,
      test_set = test_set,
      include_linear_rate = TRUE,
      include_tide_gauge = FALSE
    )
    # Joining in the rest of the dataset into the training set:
    # input_train_data <- cbind(input_train$data,more_input_data)
    
    # Setting up 2 sections for South & North Atlantic coast of North America
    data <-   input_train$data %>% 
      mutate(section = ifelse(Latitude > 35.3,"North","South"),
                    SiteName = as.character(SiteName)) %>% 
      mutate(SiteName = factor(SiteName))
    # Data grid
    data_grid <-   input_train$data_grid %>% 
      mutate(section = ifelse(Latitude > 35.3,"North","South"),
             SiteName = as.character(SiteName)) %>% 
      mutate(SiteName = factor(SiteName))

    # reslr_mcmc
    train_output <- reslr_mcmc_fun(data, data_grid,
      spline_nseg_t = spline_nseg_t,
      spline_nseg_st = spline_nseg_st,
      spline_nseg_c = spline_nseg_c,
      n_iterations = n_iterations,
      n_burnin = n_burnin,
      n_thin = n_thin,
      n_chains = n_chains,
      CI = 0.95
    )

    output_df <- train_output$output_dataframes
    # Column to identify the fold number in the loop
    output_df$CV_fold_number <- as.character(i)
    # Append this df into a list to combine to do the tests
    model_run_list[i] <- list(output_df)
    saveRDS(model_run_list,paste0("reslr_outputs/model_run_lists/mod_run_fold_",i,".rds"))
  }
  saveRDS(model_run_list,"reslr_outputs/model_run_lists/mod_run_v1.rds")

  # Combining all the dataframes
  CV_model_run_df <- suppressWarnings(
    dplyr::bind_rows(model_run_list)
  )
  # Removing rows without the test set:
  CV_model_df <- CV_model_run_df %>%
    dplyr::filter(is.na(CV_fold) == FALSE)

  # Mean Error & Mean Absolute Error & Root mean square error for each fold:
  ME_MAE_RSME_fold <- CV_model_df %>%
    dplyr::mutate(CV_fold_number = as.factor(CV_fold_number)) %>%
    dplyr::group_by(CV_fold_number) %>%
    dplyr::reframe(
      RSME = unique(sqrt((sum(RSL - pred)^2) / dplyr::n())),
      MAE = unique(sum(abs(RSL - pred)) / dplyr::n()),
      ME = unique(mean(RSL - pred))
    )
  # Mean Error & Mean Absolute Error & Root mean square error overall
  ME_MAE_RSME_overall <- CV_model_df %>%
    dplyr::reframe(
      RSME = unique(sqrt((sum(RSL - pred)^2) / dplyr::n())),
      MAE = unique(sum(abs(RSL - pred)) / dplyr::n()),
      ME = unique(mean(RSL - pred))
    )

  # Mean Error & Mean Absolute Error & Root mean square error for each fold & site:
  ME_MAE_RSME_fold_site <- CV_model_df %>%
    dplyr::mutate(CV_fold_number = as.factor(CV_fold_number)) %>%
    dplyr::group_by(SiteName, CV_fold_number) %>%
    dplyr::reframe(
      RSME = unique(sqrt((sum(RSL - pred)^2) / dplyr::n())),
      MAE = unique(sum(abs(RSL - pred)) / dplyr::n()),
      ME = unique(mean(RSL - pred))
    )

  # Mean Error & Mean Absolute Error & Root mean square error for each site:
  ME_MAE_RSME_site <- CV_model_df %>%
    dplyr::group_by(SiteName) %>%
    dplyr::reframe(
      RSME = unique(sqrt((sum(RSL - pred)^2) / dplyr::n())),
      MAE = unique(sum(abs(RSL - pred)) / dplyr::n()),
      ME = unique(mean(RSL - pred))
    )

  # Model dataframe CV
  CV_model_df <- CV_model_df %>%
    dplyr::rename(
      true_RSL = RSL,
      pred_RSL = pred
    )
  # Creating the prediction intervals outside JAGS to include Age Error
  # Overall Empirical Coverage
  CV_model_df <- CV_model_df %>%
    dplyr::mutate(
      obs_in_PI =
        ifelse(dplyr::between(
          true_RSL, upr_PI,
          lwr_PI
        ),
        TRUE, FALSE
        )
    )
  # Total coverage is trues/ number of rows with the prediction interval
  total_coverage <-
    length(which(CV_model_df$obs_in_PI == "TRUE")) / nrow(CV_model_df)
  # Coverage is trues/ number of rows with the prediction interval by site
  coverage_by_site <- CV_model_df %>%
    dplyr::group_by(SiteName) %>%
    dplyr::reframe(
      coverage_by_site =
        unique(length(which(obs_in_PI == "TRUE")) / dplyr::n())
    )

  # Prediction Interval size
  prediction_interval_size <- CV_model_df %>%
    dplyr::group_by(SiteName) %>%
    dplyr::reframe(PI_width = unique(mean(upr_PI - lwr_PI)))


  # True vs Predicted plot for tide gauges
  if ("TideGaugeData" %in% CV_model_df$data_type_id == TRUE) {
    CV_model_df_tg <- CV_model_df %>% filter(data_type_id == "TideGaugeData")
    true_pred_plot_tg <- ggplot2::ggplot(data = CV_model_df_tg, ggplot2::aes(
      x = true_RSL,
      y = y_post_pred,
      colour = "PI"
    )) +
      ggplot2::geom_errorbar(
        data = CV_model_df_tg,
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
        data = CV_model_df_tg,
        ggplot2::aes(intercept = 0, slope = 1, colour = "True = Predicted")
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 9, face = "bold"),
        axis.text = ggplot2::element_text(size = 9),
        strip.background = ggplot2::element_rect(fill = c("white")),
        strip.text = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 7),
        legend.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 8),
        axis.text.y = ggplot2::element_text(size = 8)
      ) +
      ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
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
      )
  }
  else{
    true_pred_plot_tg <- "No plot"
  }
  # True vs Predicted plot for proxy records
  CV_model_df_proxy <- CV_model_df %>% filter(data_type_id == "ProxyRecord")
  true_pred_plot_proxy <- ggplot2::ggplot(data = CV_model_df_proxy, ggplot2::aes(
    x = true_RSL,
    y = y_post_pred,
    colour = "PI"
  )) +
    ggplot2::geom_errorbar(
      data = CV_model_df_proxy,
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
      data = CV_model_df_proxy,
      ggplot2::aes(intercept = 0, slope = 1, colour = "True = Predicted")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 9, face = "bold"),
      axis.text = ggplot2::element_text(size = 9),
      strip.background = ggplot2::element_rect(fill = c("white")),
      strip.text = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 8)
    ) +
    ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
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
    )



  # Return a list of CV tests
  cross_validation_tests <- list(
    ME_MAE_RSME_fold_site = ME_MAE_RSME_fold_site,
    ME_MAE_RSME_site = ME_MAE_RSME_site,
    ME_MAE_RSME_overall = ME_MAE_RSME_overall,
    ME_MAE_RSME_fold = ME_MAE_RSME_fold,
    true_pred_plot_proxy = true_pred_plot_proxy,
    true_pred_plot_tg = true_pred_plot_tg,
    CV_model_df = CV_model_df,
    total_coverage = total_coverage,
    prediction_interval_size = prediction_interval_size,
    coverage_by_site = coverage_by_site
  )

  return(cross_validation_tests)
}
