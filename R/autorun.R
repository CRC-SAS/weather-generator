
#' @title Autorun from yaml
#' @description This function allow a autorun from a yaml named autorun_parameters.yml.
#' @export
autorun_from_yaml <- function(params_yaml_filename = NULL) {

    # -----------------------------------------------------------------------------#
    # --- STEP 1. Read YML parameter file ----

    # If the YML file wasn't provided as parameter, it'll be loaded the YML file provided with the gamwgen package
    if (is.null(params_yaml_filename)) {
        use_package_provided_files <- T
        params_yaml_filename <- system.file("autorun", "autorun_parameters.yml", package = "gamwgen")
    }

    # Read the YML file and load the parameters for the autorun_from_yaml function
    if (! file.exists(params_yaml_filename)) {
        stop("The parameters file named autorun_parameters.yml don't exist.\n")
    } else {
        cat(glue::glue("Reading parameters file autorun_parameters.yml...\n"))
        params <- yaml::yaml.load_file(params_yaml_filename)
    }

    # Set correct path to files provided with gamwgen package
    if (use_package_provided_files) {
        params$fit$climate  <- system.file(glue::glue("autorun/{params$model}"), "climate.csv",  package = "gamwgen")
        params$fit$stations <- system.file(glue::glue("autorun/{params$model}"), "stations.csv", package = "gamwgen")
        params$fit$seasonal_covariates  <- system.file(glue::glue("autorun/{params$model}"), "seasonal_covariates_fit.csv", package = "gamwgen")
        params$sim$simulation_locations <- system.file(glue::glue("autorun/{params$model}"), "simulation_locations.csv",    package = "gamwgen")
        params$sim$seasonal_covariates  <- system.file(glue::glue("autorun/{params$model}"), "seasonal_covariates_sim.csv", package = "gamwgen")
    }

    # ------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------#
    # --- STEP 2. Select model ----

    if (params$model == 'local') {
        fit_ctrl_fn <- gamwgen::local_fit_control
        fit_fn      <- gamwgen::local_calibrate
        sim_ctrl_fn <- gamwgen::local_simulation_control
        sim_fn      <- gamwgen::local_simulation
    }
    if (params$model == 'spatial') {
        fit_ctrl_fn <- gamwgen::spatial_fit_control
        fit_fn      <- gamwgen::spatial_calibrate
        sim_ctrl_fn <- gamwgen::spatial_simulation_control
        sim_fn      <- gamwgen::spatial_simulation
    }

    #-------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------#
    # --- STEP 3. Run fit ----

    if (params$fit$run) {

        fit_ctrl   <- fit_ctrl_fn(prcp_occurrence_threshold = params$fit$control$prcp_occurrence_threshold,
                                  avbl_cores = params$fit$control$avbl_cores,
                                  planar_crs_in_metric_coords = params$fit$control$planar_crs_in_metric_coords)
        fitted_obj <- fit_fn(climate = readr::read_csv(params$fit$climate, col_types = "Diddd") %>% tibble::as_tibble(),
                             stations = readr::read_csv(params$fit$stations, col_types = list(station_id = "i")) %>%
                                 sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
                                 sf::st_transform(crs = params$fit$control$planar_crs_in_metric_coords),
                             seasonal_covariates = readr::read_csv(params$fit$seasonal_covariates, col_types = "iiiddd") %>% tibble::as_tibble(),
                             control = fit_ctrl, verbose = params$fit$verbose)

        if (!is.null(params$fit$fitted_object))
            base::saveRDS(fitted_obj, params$fit$fitted_object)

    }

    #-------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------#
    # --- STEP 4. Run sim ----

    if (params$sim$run) {

        sim_ctrl   <- do.call(sim_ctrl_fn,
                              c(list(nsim = params$sim$control$nsim,
                                     seed = params$sim$control$seed,
                                     avbl_cores = params$sim$control$avbl_cores),
                                list(bbox_offset = params$sim$control$bbox_offset,
                                     sim_loc_as_grid = params$sim$control$sim_loc_as_grid)[base::eval(params$model == 'spatial')],
                                list(use_spatially_correlated_noise = params$sim$control$use_spatially_correlated_noise,
                                     use_temporary_files_to_save_ram = params$sim$control$use_temporary_files_to_save_ram,
                                     remove_temp_files_used_to_save_ram = params$sim$control$remove_temp_files_used_to_save_ram,
                                     manage_parallelization_externally = F)))
        simted_obj <- sim_fn(model = if(base::exists("fitted_obj")) fitted_obj else base::readRDS(params$sim$fitted_model),
                             simulation_locations = readr::read_csv(params$sim$simulation_locations, col_types = list(station_id = "i")) %>%
                                 sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
                                 sf::st_transform(crs = params$fit$control$planar_crs_in_metric_coords),
                             start_date = as.Date(params$sim$start_date),
                             end_date = as.Date(params$sim$end_date),
                             control = sim_ctrl,
                             output_folder = if(params$sim$output_folder == "") base::getwd() else params$sim$output_folder,
                             output_filename = params$sim$output_filename,
                             seasonal_covariates = readr::read_csv(params$sim$seasonal_covariates, col_types = "iiddd") %>% tibble::as_tibble(),
                             verbose = params$sim$verbose)

        if (!is.null(params$sim$simulated_object))
            base::saveRDS(simted_obj, params$sim$simulated_object)

    }

    #-------------------------------------------------------------------------------

}
