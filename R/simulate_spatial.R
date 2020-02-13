
## OBS: modo de ejecución (son tres)
## 1- locations a simular coinciden con los puntos ajustados (fit$stations), no necesitan estar todos los puntos ajustados pero todos los puntos a simular debieron haberse usado en el ajuste
## 1.1- ruido Santi
## 1.2- ruido nuevo (espacialmente correlacionado)
## 2- locations a simular no coinciden pero tampoco son una grilla regular (listo)
## 3- locations a simular es una grilla regular (listo)


#' @title Simulations control configuration
#' @description Provides fine control of different parameters that will be used to create new weather series.
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).#'
#'          Either NULL or an integer that will be used in a call to set.seed before simulating the response vectors.
#'          If set, the value is saved as the "seed" attribute of the returned value. The default, NULL will not change the random generator state.
#' @param avbl_cores ...
#' @export
spatial_simulation_control <- function(nsim = 1,
                                       seed = NULL,
                                       avbl_cores = 2,
                                       bbox_offset = 100000,
                                       sim_loc_as_grid = T,
                                       use_spatially_correlated_noise = T,
                                       use_temporary_files_to_save_ram = T,
                                       remove_temp_files_used_to_save_ram = T,
                                       manage_parallelization_externally = F) {

    prcp_noise_generating_function = gamwgen:::random_field_noise_prcp
    temperature_noise_generating_function = gamwgen:::random_field_noise_temperature

    if(!use_spatially_correlated_noise) {
        prcp_noise_generating_function = gamwgen:::not_spatially_correlated_random_field_noise_prcp
        temperature_noise_generating_function = gamwgen:::not_spatially_correlated_random_field_noise_temperature
    }

    return(list(nsim = nsim, seed = seed, avbl_cores = avbl_cores,
                bbox_offset = bbox_offset, sim_loc_as_grid = sim_loc_as_grid,
                use_spatially_correlated_noise = use_spatially_correlated_noise,
                use_temporary_files_to_save_ram = use_temporary_files_to_save_ram,
                remove_temp_files_used_to_save_ram = remove_temp_files_used_to_save_ram,
                manage_parallelization_externally = manage_parallelization_externally,
                prcp_noise_generating_function = prcp_noise_generating_function,
                temperature_noise_generating_function = temperature_noise_generating_function))
}



#' @title Simulates new weather trajectories in stations
#' @description Simulates new weather trajectories.
#' @param model A gamwgen model.
#' @param simulation_locations a sf object with the points at which weather should be simulated.
#'          If not set, the locations used to fit the model will be used.
#' @param start_date a start date in text format (will be converted using as.Date) or a date object.
#' @param end_date an end date in text format (will be converted using as.Date) or a date object.
#' @param control a gamwgen simulation control list.
#' @import dplyr
#' @import foreach
#' @export
spatial_simulation <- function(model, simulation_locations, start_date, end_date,
                               control = gamwgen:::spatial_simulation_control(),
                               output_folder = getwd(), output_filename = "sim_results.nc",
                               seasonal_covariates = NULL, verbose = F) {

    ## Object to be returned
    gen_climate <- list()

    ###############################################################

    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("foreach"))

    ###############################################################

    # Check that the fitted object is of the right class
    if(class(model) != 'gamwgen')
        stop(glue::glue('Received a model of class {class(model)} and a model of class "gamwgen" was expected.'))

    # Check that the locations to be simulated exists
    if (is.null(simulation_locations))
        stop("The parameter simulation_locations can't be null!")

    # Check that the input objects are of the right class
    gamwgen:::check.simulation.input.data(simulation_locations, seasonal_covariates)

    ###############################################################

    # Check that the locations to be simulated are projected in the same coordinate
    # system as the fitted data. Otherwise, convert it
    if(sf::st_crs(simulation_locations) != sf::st_crs(model$crs_used_to_fit)) {
        simulation_locations <- simulation_locations %>%
            sf::st_transform(sf::st_crs(model$crs_used_to_fit))
        warning('The crs used to fit and the crs of simulation_locations are not equals. ',
                'Se transforma simulation_locations al crs {del ajuste}')
    }

    ## Create a bounding box to check whether the simulation location are englufed by it
    sl_bbox <- sf::st_bbox(model$stations)

    sl_bbox_offset <- sf::st_bbox(
        obj = c(xmin = sl_bbox[['xmin']] - control$bbox_offset,
                ymin = sl_bbox[['ymin']] - control$bbox_offset,
                xmax = sl_bbox[['xmax']] + control$bbox_offset,
                ymax = sl_bbox[['ymax']] + control$bbox_offset),
        crs = sf::st_crs(model$crs_used_to_fit))

    polygon_offset <- sf::st_as_sfc(sl_bbox_offset)
    sf::st_crs(polygon_offset) <- sf::st_crs(model$crs_used_to_fit)

    ## Simulation locations should be inside the bounding box
    if(!all(sf::st_contains(polygon_offset, simulation_locations, sparse = F)))
        stop('Alguno de los puntos ajustados se encuentra fuera del bounding box ',
             'creado a partir de los puntos a simular y con un offset de 10 km.')

    ## Non spatially correlated noise can only be used if the weather stations
    ## to be simulated where used in the fitting process. The residuals of those stations are
    ## needed to model the multivariate distributions.
    if(!control$use_spatially_correlated_noise)
        if(any(lapply(sf::st_equals(simulation_locations, model$stations), length) != 1))
            stop('Los ruidos NO correlacionados espacialmente solo pueden ser usados cuando ',
                 'los puntos a ser simulados son los mismos que fueron utilizados en el ajuste!')


    ###############################################################

    # Check that the input location to be simulated is a valid sf object
    if(!all(sf::st_is_valid(simulation_locations)))
        stop('simulation_locations is not a valid sf object.')

    ###############################################################

    # Check consistency between start and end date of the simulation period
    if(end_date <= start_date)
        stop('End date should be greater than start date')

    ###############################################################

    # Check that the number of realization is equal to or larger than
    if(control$nsim < 1)
        stop('Number of simulations should be one or greater than one')

    ###############################################################

    # Check that the number of cores to be used is valid
    if(is.na(control$avbl_cores) || is.null(control$avbl_cores))
        stop('The control parameter avbl_cores must be at least 1.')

    ###############################################################

    if(!foreach::getDoParRegistered() && control$manage_parallelization_externally)
        stop('The control parameter manage_parallelization_externally was set as True, but',
             'neither sequential nor parallel backend was registered for the foreach package.')

    ###############################################################

    ## Control de uso correcto de covariables

    # Esquema de uso de covariables!!
    # cov ajuste    |     cov simulación
    #  interna      |      interna
    #  interna      |      externo
    #  externa      |      externa  ->  en ambos casos seasonal_covariates es el mismo tibble
    #  externa      |      externa  ->  diferentes (considerar siempre diferentes y listo)
    #  externa      |      interna  (no corresponde)
    # generalmente el ajuste es con covariables internas

    ## If seasonal_covariates is not NULL, model fit should have been done using covariates.
    if(!is.null(seasonal_covariates) & is.null(model$seasonal_covariates))
        stop('El ajuste fue hecho sin covariables, por lo tanto, la simulación ',
             'también debe hacerse sin covariables y no es valido utilizar el ',
             'parámetro seasonal_covariates!!')

    ## If model fit was done using an external set of covariates, the simulation
    ## should be done with an external set of variates as well.
    if(is.null(seasonal_covariates) & !is.null(model$seasonal_covariates))
        stop('El ajuste se hizo utilizando un archivo de covariables (parametro ',
             'seasonal_covariates), por lo tanto, la simulación también debe hacerse ',
             'con un archivo de covariables (parametro seasonal_covariates).')

    ## OBS:
    # Ocurre lo siguiente: si ajuste se hace sin covariables, la simulación también debe
    # hacerse sin covariables, y si el ajuste se hace con covariables, la simulación
    # también debe hacerse con convariables!!

    ###############################################################

    ## Check the presence of seasonal covariables time series as long as the simulation period.
    ## This control is only performed if the simulation will use covariables otherwise, it will be skipped.
    if(!is.null(seasonal_covariates)) {
        if(gamwgen:::repeat_seasonal_covariates(seasonal_covariates)) {
            sim_dates_control <- tidyr::crossing(year = base::seq.int(lubridate::year(start_date), lubridate::year(end_date)),
                                                 season = as.integer(c(1, 2, 3, 4)))
            seasonal_cov_ctrl <- seasonal_covariates %>% dplyr::select(year, season) %>% dplyr::distinct()
        } else {
            sim_dates_control <- tidyr::crossing(seasonal_covariates %>% dplyr::distinct(station_id),
                                                 year = base::seq.int(lubridate::year(start_date), lubridate::year(end_date)),
                                                 season = as.integer(c(1, 2, 3, 4)))
            seasonal_cov_ctrl <- seasonal_covariates %>% dplyr::select(station_id, year, season) %>% dplyr::distinct()
        }
        if (!all(do.call(paste0, sim_dates_control) %in% do.call(paste0, seasonal_cov_ctrl)))
            stop("Simulation years aren't in seasonal_covariates!")
    }

    ###############################################################



    ############################
    ## INITIAL CONFIGURATIONS ##
    ############################

    # Configuration of the RandomFields package in order to produce the expected results
    RandomFields::RFoptions(spConform=FALSE)

    # Para ...
    if(!is.null(control$seed))
        set.seed(control$seed)

    # Create simulation seeds to replicate results in different simulations
    realizations_seeds <- NULL
    if(!is.null(control$seed)) {
        realizations_seeds <- list()
        cant_dias_sim <- as.numeric(end_date - start_date) + 1
        for (r in seq(1, control$nsim, 1)) {
            realizations_seeds[[r]] <- list(general = ceiling(runif(min = 1, max = 10000000, n = control$nsim)), # uno por realizacion
                                            prcp_occ = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim)), # uno por día por realizacion
                                            prcp_amt = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim)), # uno por día por realizacion
                                            temp_dry = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim)), # uno por día por realizacion
                                            temp_wet = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim)), # uno por día por realizacion
                                            retries = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim))) # uno por día por realizacion
        }
    }

    ####################################
    ## Parallelization initialization ##
    ####################################


    ## Variable that indicate if it's necessary to remove
    ## the parallelization configuration
    remove_parallelization_conf <- F

    ## Register a sequential or a parallel backend for the foreach package
    if(!control$manage_parallelization_externally) {
        if (control$avbl_cores == 1) {
            foreach::registerDoSEQ()
        } else {
            remove_parallelization_conf <- T
            # Create cluster
            cluster  <- snow::makeCluster(type = "SOCK",
                                          spec = rep('localhost', length.out = control$avbl_cores),
                                          outfile = ifelse(verbose, "", snow::defaultClusterOptions$outfile))
            # Register cluster as backend for the %dopar% function
            doSNOW::registerDoSNOW(cluster)
        }
    }

    ## Register the number of workers to decide
    ## how manage progress bar in child processes
    nworkers <- foreach::getDoParWorkers()



    ##################################
    ## PREPARACIÓN DE LA SIMULACIÓN ##
    ##################################


    ############################################
    ## Process output_folder and output_filename
    output_folder <- sub('/$', '', output_folder)
    output_filename <- sub('\\.([^.]*)$', '', output_filename)


    ####################################
    ## Generate simulation dates
    simulation_dates <-
        tibble::tibble(date = seq.Date(from = as.Date(start_date),
                                       to = as.Date(end_date),
                                       by = "days")) %>%
        dplyr::mutate(year   = lubridate::year(date),
                      month  = lubridate::month(date),
                      season = lubridate::quarter(date, fiscal_start = 12))
    ## Numbers of days to be simulated
    ndates <- nrow(simulation_dates)


    ##################################
    ## Matrix with the locations to be simulated
    simulation_points <- simulation_locations %>%
        sf::st_transform(sf::st_crs(model$crs_used_to_fit)) %>%
        dplyr::mutate(point_id = dplyr::row_number(),
                      longitude = sf::st_coordinates(geometry)[,'X'],
                      latitude  = sf::st_coordinates(geometry)[,'Y'])


    ##################################
    ## Raster with the locations to be simulated.
    if (control$sim_loc_as_grid) {
        pnts_dist_matrix  <- gamwgen:::make_distance_matrix(simulation_locations)
        min_dist_btw_pnts <- floor(min(pnts_dist_matrix[upper.tri(pnts_dist_matrix)]))
        raster_resolution <- c(min_dist_btw_pnts, min_dist_btw_pnts)
        simulation_raster <- raster::rasterFromXYZ(
            xyz = sf::st_coordinates(simulation_points),
            res = raster_resolution,
            crs = sf::st_crs(simulation_points))
    }
    if(!control$sim_loc_as_grid) {
        stns_dist_matrix  <- gamwgen:::make_distance_matrix(model$stations)
        min_dist_btw_stns <- floor(min(stns_dist_matrix[upper.tri(stns_dist_matrix)]))
        raster_resolution <- c(min_dist_btw_stns, min_dist_btw_stns)
        simulation_raster <- raster::raster(
            xmn = sl_bbox_offset[['xmin']],
            xmx = sl_bbox_offset[['xmax']],
            ymn = sl_bbox_offset[['ymin']],
            ymx = sl_bbox_offset[['ymax']],
            resolution = raster_resolution,
            crs = sf::st_crs(simulation_points))
    }


    ################################################################
    ## Obtaining values for the day before the simulation start day
    start_date_prev_day_climatology <-
        gamwgen:::get_start_climatology(model, simulation_points, start_date, control)


    ######################################################################################
    ## Si no se recibe un seasonal_covariates externo, se utiliza el generado en el ajuste
    # Si por ahí fallan los controles, con esto se garantiza que: si se usaron covariables
    # al momento de realizar el ajuste, también se utilize covariables en la simulación, y
    # lo mismo para el caso inverso, es decir, cuando no se usó covariables ene l ajuste!!
    if(is.null(seasonal_covariates))
        seasonal_covariates <- model$seasonal_covariates


    #################################################################
    ## Se pueden excluir los registros de años que no serán simulados
    if (!is.null(seasonal_covariates))
        seasonal_covariates <- seasonal_covariates %>%
            dplyr::filter(year %in% unique(simulation_dates$year))


    #############################################################
    ## Obtención de covariables, si van a ser utilizadas, sino no
    if (!is.null(seasonal_covariates))
        seasonal_covariates <-
            gamwgen:::get_covariates(model, simulation_points, seasonal_covariates, simulation_dates, control)


    #########################################
    ## Global paramteres for noise generators
    if(control$use_spatially_correlated_noise)
        gen_noise_params <- gamwgen:::generate_month_params(
            residuals = model$models_residuals,
            observed_climate = model$models_data,
            stations = model$stations)
    if(!control$use_spatially_correlated_noise)
        gen_noise_params <- gamwgen:::generate_residuals_statistics(
            models_residuals = model$models_residuals)


    ############################################
    ## Matriz de clasificacion de días lluviosos
    clasification_matrix <- matrix(c(-Inf, 0, 0, 0, Inf, 1),
                                   ncol = 3,
                                   byrow = TRUE)


    #######################################################################################
    ## Se crea una matriz de simulación, esta va a contener todos los datos necesarios para
    ## la simulación de cada día a simular, si se utilizan covariables, se las incluye aquí
    simulation_matrix <- simulation_points %>%
        dplyr::select(point_id, longitude, latitude)


    ################################################################################################
    ## Add seasonal covariates if they were used when model fitting (only years in simulation_dates)
    if (!is.null(seasonal_covariates))
        simulation_matrix <- simulation_matrix %>% sf::st_join(seasonal_covariates) %>%
            dplyr::mutate(ST1 = dplyr::if_else(season == 1, seasonal_prcp, 0),
                          ST2 = dplyr::if_else(season == 2, seasonal_prcp, 0),
                          ST3 = dplyr::if_else(season == 3, seasonal_prcp, 0),
                          ST4 = dplyr::if_else(season == 4, seasonal_prcp, 0),
                          SN1 = dplyr::if_else(season == 1, seasonal_tmin, 0),
                          SN2 = dplyr::if_else(season == 2, seasonal_tmin, 0),
                          SN3 = dplyr::if_else(season == 3, seasonal_tmin, 0),
                          SN4 = dplyr::if_else(season == 4, seasonal_tmin, 0),
                          SX1 = dplyr::if_else(season == 1, seasonal_tmax, 0),
                          SX2 = dplyr::if_else(season == 2, seasonal_tmax, 0),
                          SX3 = dplyr::if_else(season == 3, seasonal_tmax, 0),
                          SX4 = dplyr::if_else(season == 4, seasonal_tmax, 0))


    #########################################
    ## Thresholds for performing retries. i.e.: if simulated values are above/below the
    ## max/min range, the simulation for that day will be repeated.
    temperature_range_thresholds <-
        gamwgen:::get_temperature_thresholds(model$stations, simulation_points,
                                             model$estadisticos_umbrales, control)


    ##########################################
    ## AQUI EMPIEZA REALMENTE LA SIMULACIÓN ##
    ##########################################


    #################################
    ## Control de tiempo de ejecución
    t.sim <- proc.time()


    #######################################################
    ## Set the progress bar to know how long this will take
    pb <- progress::progress_bar$new(
        format = paste0(ifelse(nworkers == 1, " realization:", " finished realizations:"),
                        " :r / ", control$nsim, ifelse(nworkers == 1, paste0(" | day: :d / ", ndates), ""),
                        " | progress: :bar :percent (in :elapsed) | eta: :eta"),
        total = ifelse(nworkers == 1, control$nsim*ndates, control$nsim),
        clear = FALSE, width= 90, show_after = 0)

    ## For print something until first realization finish
    if(nworkers > 1 && !verbose)
        pb$tick(0, tokens = list(r = 0))

    ## For manage the progress bar in parallel executions
    progress_pb <- function(r) {
        pb$tick(1, tokens = list(r = r))
    }


    ###########################################
    ## Crear objeto para guardar los resultados
    gamwgen:::CrearNetCDF(nc_file = glue::glue("{output_folder}/{output_filename}.nc"),
                          num_realizations = control$nsim,
                          sim_dates = simulation_dates$date,
                          coordinates = raster::rasterToPoints(simulation_raster),
                          coord_ref_system = sf::st_crs(simulation_points))


    ######
    ##
    nsim_gen_clim <- foreach::foreach(r = 1:control$nsim,
                                      .combine = dplyr::bind_rows,
                                      .export = c('output_folder', 'output_filename'),
                                      .packages = c('dplyr'),
                                      .options.snow = list(progress = progress_pb),
                                      .verbose=verbose) %dopar% {

        ######################################################################
        ## Para que las funciones de RandomFields devuelvan lo esperado!! ----
        RandomFields::RFoptions(spConform=FALSE)


        ##################################################
        ## Para cuando necesitamos repetir resultados ----
        set.seed(realizations_seeds[[r]]$general[[r]])


        ################################################################################
        ## Cuando se ejecuta el código en paralelo, simulation_matrix no es un sf válido
        if(nworkers > 1)
            simulation_matrix <- simulation_matrix %>%
                sf::st_as_sf(coords = c('longitude', 'latitude'), crs = sf::st_crs(simulation_points))


        #################################################################################
        ## Creacion de los puntos de simulacion para el dia i (eq. daily covariates) ----
        #################################################################################
        simulation_matrix.d <- simulation_matrix %>%
            # Si se usan covariables, simulation_matrix tiene las covariables
            # para cada season de cada year en simulation_dates! Por lo tanto,
            # se debe filtrar por season y year para acelerar el sf::st_join!
            {if (is.null(seasonal_covariates)) dplyr::filter(.)
             else dplyr::filter(., year == simulation_dates$year[1], season == simulation_dates$season[1])} %>%
            # Luego se agrega la climatología inicial, es decir, la del día previo al primer día
            # a simular. Las columnas incorporadas por esta acción son: prcp_occ, tmax y tmin,
            # por lo tanto, deben ser renombradas a: prcp_occ_prev, tmax_prev y tmin_prev
            sf::st_join(start_date_prev_day_climatology) %>%
            dplyr::rename(prcp_occ_prev = prcp_occ, tmax_prev = tmax, tmin_prev = tmin) %>%
            # a esto se debe agregar tipo_dia_prev y prcp_amt_prev
            dplyr::mutate(tipo_dia_prev = factor(prcp_occ_prev, levels = c(1, 0), labels = c('Lluvioso', 'Seco')),
                          prcp_amt_prev = NA_real_) %>%  # no se tiene la amplitud de prcp para el día previo al 1er día a simular!
            # Se agregan date, time y doy del primer día a simular
            dplyr::mutate(date = simulation_dates$date[1],
                          time = as.numeric(date)/1000,
                          doy = lubridate::yday(date),
                          month = lubridate::month(date)) %>%
            # Se crean columnas para almacenar los resultados de la simulación
            dplyr::mutate(prcp_occ = NA_integer_,
                          tmax = NA_real_,
                          tmin = NA_real_,
                          tipo_dia = NA_character_,
                          prcp_amt = NA_real_) %>%
            # para control de paralelización
            dplyr::mutate(nsim = r)


        #######################################
        ## Tiempos a tomar por cada realización
        tiempos <- tibble::tibble(tiempo.gen_clim = list(),
                                  tiempo.save_rds = list())
        tiempos <- tiempos %>% dplyr::add_row()


        #################################
        ## Control de tiempo de ejecución
        t.daily_gen_clim <- proc.time()


        #####################################################################################
        ## Antes se usaba un foreach para paralelizar esto, pero no se puede ser paralelizado
        ## porque simulation_matrix.d no toma los valores correctos al paralelizar!!
        ## Ver version anterior para más detalles (commit: 1898e5a)
        daily_gen_clim <- purrr::map_dfr(1:ndates, function(d) {

            ##############################################################
            ## Temporal index for each month of the simulation/realization
            current_month <- simulation_dates$month[d]



            ###########################################
            ## Precipitation occurrence (prcp_occ) ----
            ###########################################

            # Simulation of precipitation occurrence
            SIMocc <- mgcv::predict.bam(model$fitted_models$prcp_occ_fit,
                                        newdata = simulation_matrix.d,
                                        #cluster = cluster,  # empeora el tiempo para grillas grandes
                                        newdata.guaranteed = TRUE) # una optimizacion

            # Raster con los valores "climáticos"
            SIMocc_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMocc = !!SIMocc) %>%
                gamwgen:::sf2raster('SIMocc', simulation_raster)

            # Raster con los valores de "ruido"
            SIMocc_points_noise.d <- control$prcp_noise_generating_function(
                simulation_points = simulation_points,
                gen_noise_params = gen_noise_params,
                month_number = current_month,
                selector = 'prcp',
                seed = realizations_seeds[[r]]$prcp_occ[[d]]) %>%
                gamwgen:::sf2raster('prcp_residuals', simulation_raster)

            # Raster con los valores simulados
            SIMocc_points.d <- SIMocc_points_climate.d + SIMocc_points_noise.d

            # Raster con los valores simulados reclasificados a 0 y/o 1
            SIMocc_points.d <- raster::reclassify(SIMocc_points.d, clasification_matrix)

            # Agregar valores de ocurrencia a la grilla de simulacion
            simulation_matrix.d <- simulation_matrix.d %>%
                dplyr::mutate(prcp_occ = raster::extract(SIMocc_points.d, simulation_points),
                              tipo_dia = factor(prcp_occ, levels = c(1, 0),
                                                labels = c('Lluvioso', 'Seco')))



            #########################################
            ## Temperature (both, tmax and tmin) ----
            #########################################

            #  Raster con los valores de "ruido" para días secos
            temperature_random_fields_dry <-
                control$temperature_noise_generating_function(
                    simulation_points = simulation_points,
                    gen_noise_params = gen_noise_params,
                    month_number = current_month,
                    selector = c('tmax_dry', 'tmin_dry'),
                    seed = realizations_seeds[[r]]$temp_dry[[d]])

            # Procesamiento de residuos para dias secos
            rasters_secos.d <- purrr::map(
                .x = c("tmin", "tmax"),
                .f = function(variable, objeto_sf) {
                    return (gamwgen:::sf2raster(objeto_sf, paste0(variable, '_residuals'), simulation_raster))
                }, objeto_sf = temperature_random_fields_dry
            )
            names(rasters_secos.d) <- c("tmin", "tmax")

            # Raster con los valores de "ruido" para días lluviosos
            temperature_random_fields_wet <-
                control$temperature_noise_generating_function(
                    simulation_points = simulation_points,
                    gen_noise_params = gen_noise_params,
                    month_number = current_month,
                    selector = c('tmax_wet', 'tmin_wet'),
                    seed = realizations_seeds[[r]]$temp_wet[[d]])

            # Procesamiento de residuos para dias humedos
            rasters_humedos.d <- purrr::map(
                .x = c("tmin", "tmax"),
                .f = function(variable, objeto_sf) {
                    return (gamwgen:::sf2raster(objeto_sf, paste0(variable, '_residuals'), simulation_raster))
                }, objeto_sf = temperature_random_fields_wet
            )
            names(rasters_humedos.d) <- c("tmin", "tmax")

            #Ahora vamos a generar 2 rasters: uno para tmax y otro para tmin
            #Cada raster tiene los residuos de las variables correspondientes
            #considerando la ocurrencia de dia seco o humedo
            SIMmax_points_noise.d <-
                gamwgen:::ensamblar_raster_residuos(rasters_humedos.d$tmax, rasters_secos.d$tmax, SIMocc_points.d)
            SIMmin_points_noise.d <-
                gamwgen:::ensamblar_raster_residuos(rasters_humedos.d$tmin, rasters_secos.d$tmin, SIMocc_points.d)



            ##################################
            ## Maximum temperature (tmax) ----
            ##################################

            # Simulation of maximum temperature (SIMmax)
            SIMmax <- mgcv::predict.bam(model$fitted_models$tmax_fit,
                                        newdata = simulation_matrix.d,
                                        #cluster = cluster,  # no mejora mucho el tiempo
                                        newdata.guaranteed = TRUE) # una optimizacion

            # Raster con los valores "climáticos"
            SIMmax_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMmax = !!SIMmax) %>%
                gamwgen:::sf2raster('SIMmax', simulation_raster)

            # Raster con los valores simulados
            SIMmax_points.d <- SIMmax_points_climate.d + SIMmax_points_noise.d

            # Agregar valores de temperatura mínima a los puntos de simulación
            simulation_matrix.d <- simulation_matrix.d %>%
                dplyr::mutate(tmax = raster::extract(SIMmax_points.d, simulation_points))



            ##################################
            ## Minimum temperature (tmin) ----
            ##################################

            # Simulation of minimum temperature (SIMmin)
            SIMmin <- mgcv::predict.bam(model$fitted_models$tmin_fit,
                                        newdata = simulation_matrix.d,
                                        #cluster = cluster,  # no mejora mucho el tiempo
                                        newdata.guaranteed = TRUE) # una optimizacion

            # Raster con los valores "climáticos"
            SIMmin_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMmin = !!SIMmin) %>%
                gamwgen:::sf2raster('SIMmin', simulation_raster)

            # Raster con los valores simulados
            SIMmin_points.d <- SIMmin_points_climate.d + SIMmin_points_noise.d

            # Agregar valores de temperatura mínima a los puntos de simulación
            simulation_matrix.d <- simulation_matrix.d %>%
                dplyr::mutate(tmin = raster::extract(SIMmin_points.d, simulation_points))



            #################################################
            ## Check Temperatures (both, tmax and tmin) ----
            #################################################

            # Simulated maximum and minium temperature values will be valid the daily
            # temperature range (tmax - tmin) falls between the thresholds estimated
            # from the original data. If the simulated daily temperature is outside the
            # thresholds, i.e.: above the maximum range or beneath the minimum range,
            # the daily value will be resimulated until the condition is satisfied.


            #################
            ## Dry thresholds
            # Maximum range
            dry_max_range <- temperature_range_thresholds %>%
                dplyr::filter(month == current_month & prcp_occ == 0) %>%
                sf::st_as_sf(., coords = c('longitude', 'latitude')) %>%
                gamwgen:::sf2raster('max.range', simulation_raster)
            # Minimum range
            dry_min_range <- temperature_range_thresholds %>%
                dplyr::filter(month == current_month & prcp_occ == 0) %>%
                sf::st_as_sf(., coords = c('longitude', 'latitude')) %>%
                gamwgen:::sf2raster('min.range', simulation_raster)


            ################
            # Wet thresholds
            # Maximum range
            wet_max_range <- temperature_range_thresholds %>%
                dplyr::filter(month == current_month & prcp_occ == 1) %>%
                sf::st_as_sf(., coords = c('longitude', 'latitude')) %>%
                gamwgen:::sf2raster('max.range', simulation_raster)
            # Minimum range
            wet_min_range <- temperature_range_thresholds %>%
                dplyr::filter(month == current_month & prcp_occ == 0) %>%
                sf::st_as_sf(., coords = c('longitude', 'latitude')) %>%
                gamwgen:::sf2raster('min.range', simulation_raster)


            ################
            ## Merge rasters
            # Maximum range
            maximum_daily_range_raster <-
                gamwgen:::ensamblar_raster_residuos(dry_max_range, wet_max_range, SIMocc_points.d)
            # Minimum range
            minimum_daily_range_raster <-
                gamwgen:::ensamblar_raster_residuos(dry_min_range, wet_min_range, SIMocc_points.d)
            # Daily range
            daily_range_points.d <- base::abs(SIMmax_points.d - SIMmin_points.d)


            # Perform the test. If the simulate temperature range is above the minimum observed range and
            # below the maximum observed range, the simulations are valid. Otherwise, re-simulate
            daily_retries <- 0
            while ( daily_retries < 100 && (any(raster::getValues(SIMmax_points.d < SIMmin_points.d), na.rm = T) ||
                                            any(raster::getValues(daily_range_points.d > maximum_daily_range_raster), na.rm = T) ||
                                            any(raster::getValues(daily_range_points.d < minimum_daily_range_raster), na.rm = T)) ) {
                daily_retries <- daily_retries  + 1

                # Raster con los valores de "ruido" para días secos
                temperature_random_fields_dry <-
                    control$temperature_noise_generating_function(
                        simulation_points = simulation_points,
                        gen_noise_params = gen_noise_params,
                        month_number = current_month,
                        selector = c('tmax_dry', 'tmin_dry'),
                        seed = if (is.null(control$seed)) NULL else realizations_seeds[[r]]$retries[[d]] + daily_retries)
                #}, times = 10) # 180 milisegundos

                # Procesamiento de residuos para dias secos
                rasters_secos.d <- purrr::map(
                    .x = c("tmin", "tmax"),
                    .f = function(variable, objeto_sf) {
                        return (gamwgen:::sf2raster(objeto_sf, paste0(variable, '_residuals'), simulation_raster))
                    }, objeto_sf = temperature_random_fields_dry
                )
                names(rasters_secos.d) <- c("tmin", "tmax")

                # Raster con los valores de "ruido" para días lluviosos
                temperature_random_fields_wet <-
                    control$temperature_noise_generating_function(
                        simulation_points = simulation_points,
                        gen_noise_params = gen_noise_params,
                        month_number = current_month,
                        selector = c('tmax_wet', 'tmin_wet'),
                        seed = if (is.null(control$seed)) NULL else realizations_seeds[[r]]$retries[[d]] + daily_retries)

                # Procesamiento de residuos para dias humedos
                rasters_humedos.d <- purrr::map(
                    .x = c("tmin", "tmax"),
                    .f = function(variable, objeto_sf) {
                        return (gamwgen:::sf2raster(objeto_sf, paste0(variable, '_residuals'), simulation_raster))
                    }, objeto_sf = temperature_random_fields_wet
                )
                names(rasters_humedos.d) <- c("tmin", "tmax")

                #Ahora vamos a generar 2 rasters: uno para tmax y otro para tmin
                #Cada raster tiene los residuos de las variables correspondientes
                #considerando la ocurrencia de dia seco o humedo
                SIMmax_points_noise.d <-
                    gamwgen:::ensamblar_raster_residuos(rasters_humedos.d$tmax, rasters_secos.d$tmax, SIMocc_points.d)
                SIMmin_points_noise.d <-
                    gamwgen:::ensamblar_raster_residuos(rasters_humedos.d$tmin, rasters_secos.d$tmin, SIMocc_points.d)

                # Generar nuevas temperaturas
                # Maximum temperature
                new_tmax <- SIMmax_points_climate.d + SIMmax_points_noise.d
                # Minimum temperature
                new_tmin <- SIMmin_points_climate.d + SIMmin_points_noise.d

                # Actualizar temperaturas simuladas
                SIMmax_points.d <- new_tmax
                SIMmin_points.d <- new_tmin

                # Daily range
                daily_range_points.d <- SIMmax_points.d - SIMmin_points.d

                #################################################
                ## Progress Bar (for non parallel execution) ----
                if(nworkers == 1) # for report retries!!
                    pb$tick(0, tokens = list(r = r, d = d, t = daily_retries))
            }


            ##########################
            ## Report retries problems
            if(daily_retries >= 100)
                warning("Failed to simulate random noise that doesn't violate the constraint of max. temp. > min. temp.")


            ###############################################################
            ## Update simulation_matrix.d with new values for tmax and tmin
            if (daily_retries > 0) {
                # Population of the simulation matrix with the simulated maximum temperature data
                simulation_matrix.d <- simulation_matrix.d %>%
                    dplyr::mutate(tmax = raster::extract(SIMmax_points.d, simulation_points))
                # Population of the simulation matrix with the simulated minimum temperature data
                simulation_matrix.d <- simulation_matrix.d %>%
                    dplyr::mutate(tmin = raster::extract(SIMmin_points.d, simulation_points))
            }



            ########################################
            ## Precipitation amounts (prcp_amt) ----
            ########################################

            # Filtrar el modelo a usar por el mes en curso
            prcp_amt_fit <- model$fitted_models$prcp_amt_fit[[current_month]]

            # Estimación del parametro de forma
            alphaamt <- MASS::gamma.shape(prcp_amt_fit)$alpha
            # Estimación de los parametros de escala
            betaamt <- base::exp(mgcv::predict.bam(prcp_amt_fit,
                                                   newdata = simulation_matrix.d,
                                                   #cluster = cluster,  # no mejora mucho el tiempo
                                                   newdata.guaranteed = TRUE))/alphaamt

            # Raster con los valores de "ruido"
            SIMamt_points_noise.d <- control$prcp_noise_generating_function(
                simulation_points = simulation_points,
                gen_noise_params = gen_noise_params,
                month_number = current_month,
                selector = 'prcp',
                seed = realizations_seeds[[r]]$prcp_amt[[d]]) %>%
                gamwgen:::sf2raster('prcp_residuals', simulation_raster)

            # Simulacion de montos
            SIMamt <- stats::qgamma(stats::pnorm(raster::extract(SIMamt_points_noise.d, simulation_points)),
                                    shape = rep(alphaamt, length(betaamt)), scale = betaamt)

            # Raster con los valores "climáticos"
            SIMamt_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMamt = !!SIMamt) %>%
                gamwgen:::sf2raster('SIMamt', simulation_raster)

            # Enmascarar pixeles sin ocurrencia de lluvia
            SIMamt_points.d <- SIMamt_points_climate.d * SIMocc_points.d

            # Agregar valores de los montos de prcp a los puntos de simulación
            simulation_matrix.d <- simulation_matrix.d %>%
                dplyr::mutate(prcp_amt = raster::extract(SIMamt_points.d, simulation_points))



            #########################################################################
            ## Preparar simulation_matrix.d para la simulación del siguiente día ----
            #########################################################################

            current_sim_matrix  <- simulation_matrix.d %>%
                sf::st_drop_geometry() %>% tibble::as_tibble()
            current_sim_results <- current_sim_matrix %>%
                dplyr::select(dplyr::one_of("station_id", "point_id"), prcp_occ, tmax, tmin, prcp_amt, tipo_dia)
            # OJO: se usa el operador <<- para utilizar los resultados el siguiente día
            simulation_matrix.d <<- simulation_matrix %>%
                # Si se usan covariables, simulation_matrix tiene las covariables
                # para cada season de cada year en simulation_dates! Por lo tanto,
                # se debe filtrar por season y year para acelerar el sf::st_join!
                {if (is.null(seasonal_covariates)) dplyr::filter(.)
                else dplyr::filter(., year == simulation_dates$year[d+1], season == simulation_dates$season[d+1])} %>%
                # Ya no se agrega la climatología inicial (la del día previo al primer día a simular),
                # sino que, como climatología previa se usan los resultados del día en curso
                dplyr::inner_join(current_sim_results, by = c("point_id")) %>%
                # Finalmente, se hacen las actualizaciones necesarias para que simulation_matrix.d
                # pueda ser utilizada en la siguiente iteración, es decir para el siguiente día a simular
                dplyr::mutate(prcp_occ_prev = prcp_occ,
                              tmax_prev = tmax,
                              tmin_prev = tmin,
                              tipo_dia_prev = tipo_dia,
                              prcp_amt_prev = prcp_amt,
                              date = simulation_dates$date[d+1],
                              time = as.numeric(date)/1000,
                              doy = lubridate::yday(date),
                              month = lubridate::month(date),
                              prcp_occ = NA_integer_,
                              tmax = NA_real_,
                              tmin = NA_real_,
                              tipo_dia = NA_character_,
                              prcp_amt = NA_real_,
                              nsim = r)



            #################################################
            ## Progress Bar (for non parallel execution) ----
            if(nworkers == 1)
                pb$tick(1, tokens = list(r = r, d = d))



            ###########################
            ## Devolver resultados ----
            return (
                tibble::tibble(
                    date = simulation_dates$date[d],
                    raster_prcp = list(SIMamt_points.d),
                    raster_tmax = list(SIMmax_points.d),
                    raster_tmin = list(SIMmin_points.d)
                ))
        })



        ##############################################
        ## Tomar tiempo de generación del clima diario
        tiempos <- dplyr::mutate(tiempos, tiempo.gen_clim = list(proc.time() - t.daily_gen_clim))



        ###################################################################
        ## Se guarda en disco el tibble con los rasters de la realización
        rds_path <- glue::glue("{output_folder}/{output_filename}_realization_{r}.rds")
        if(control$use_temporary_files_to_save_ram) {
            t.saveRDS <- proc.time()
            base::saveRDS(daily_gen_clim, rds_path)
            tiempos <- dplyr::mutate(tiempos, tiempo.save_rds = list(proc.time() - t.saveRDS))
        }



        ######################
        ## Liberar memoria RAM
        if(control$use_temporary_files_to_save_ram) {
            rm(daily_gen_clim)
            invisible(gc())
        }



        #################
        ## Retorno final
        return (tibble::tibble(nsim = r,
                               nsim_gen_climate = ifelse(control$use_temporary_files_to_save_ram, list(rds_path), list(daily_gen_clim))
                               ) %>% dplyr::bind_cols(tiempos))
    }



    ##############################################
    ## Guardar realizacion en archivo de salida ##
    ##############################################


    ####################################################
    ## Tomar tiempos de generación del archivo de salida
    ctrl_output_file <- purrr::map2_dfr(
        nsim_gen_clim %>% dplyr::pull(nsim),
        nsim_gen_clim %>% dplyr::pull(nsim_gen_climate),
        function(r, daily_gen_clim) {

            #######################################
            ## Tiempos a tomar por cada realización
            tiempos <- tibble::tibble(tiempo.read_rds = list(),
                                      tiempo.gen_file = list())
            tiempos <- tiempos %>% dplyr::add_row()

            ######################################
            ## Leer archivos con rasters generados
            if(control$use_temporary_files_to_save_ram) {
                t.read_rds <- proc.time()
                rds_path <- daily_gen_clim
                daily_gen_clim <- base::readRDS(rds_path)
                tiempos <- dplyr::mutate(tiempos, tiempo.read_rds = list(proc.time() - t.read_rds))
            }

            #############################
            ## Generar archivos de salida
            t.gen_file <- proc.time()
            gamwgen:::GuardarRealizacionNetCDFfromRasters(nc_file = glue::glue("{output_folder}/{output_filename}.nc"),
                                                          numero_realizacion = r,
                                                          sim_dates = simulation_dates$date,
                                                          raster_tmax = daily_gen_clim$raster_tmax,
                                                          raster_tmin = daily_gen_clim$raster_tmin,
                                                          raster_prcp = daily_gen_clim$raster_prcp)
            tiempos <- dplyr::mutate(tiempos, tiempo.gen_file = list(proc.time() - t.gen_file))

            ######################
            ## Liberar memoria RAM
            if(control$use_temporary_files_to_save_ram) {
                rm(daily_gen_clim)
                invisible(gc())
            }

            #################
            ## Retorno final
            return (tibble::tibble(nsim = r) %>% dplyr::bind_cols(tiempos))
        })



    #################################
    ## Control de tiempo de ejecución
    tiempo.sim <- proc.time() - t.sim



    ##############################
    ## Preparar datos de salida ##
    ##############################


    gen_climate[['nsim']] <- control$nsim
    gen_climate[['seed']] <- control$seed
    gen_climate[['realizations_seeds']] <- realizations_seeds
    gen_climate[['simulation_points']] <- simulation_points
    gen_climate[['output_file_with_results']] <- glue::glue("{output_folder}/{output_filename}.nc")
    gen_climate[['output_file_fomart']] <- "NetCDF4"

    fitted_stations <- model$stations;  climate <- model$climate
    fsc_filename <- glue::glue("{output_folder}/fitted_stations_and_climate.RData")
    save(fitted_stations, climate, file = fsc_filename)
    gen_climate[['rdata_file_with_fitted_stations_and_climate']] <- fsc_filename
    rm(fsc_filename); invisible(gc())

    names(nsim_gen_clim$tiempo.gen_clim) <- paste0("sim_", nsim_gen_clim$nsim)
    gen_climate[['exec_times']][["gen_clim_time"]] <- nsim_gen_clim$tiempo.gen_clim

    if(control$use_temporary_files_to_save_ram) {
        names(nsim_gen_clim$tiempo.save_rds) <- paste0("sim_", nsim_gen_clim$nsim)
        gen_climate[['exec_times']][["rds_save_time"]] <- nsim_gen_clim$tiempo.save_rds
    }
    if(control$use_temporary_files_to_save_ram) {
        names(ctrl_output_file$tiempo.read_rds) <- paste0("sim_", ctrl_output_file$nsim)
        gen_climate[['exec_times']][["rds_read_time"]] <- ctrl_output_file$tiempo.read_rds
    }

    names(ctrl_output_file$tiempo.gen_file) <- paste0("sim_", ctrl_output_file$nsim)
    gen_climate[['exec_times']][["gen_output_time"]] <- ctrl_output_file$tiempo.gen_file
    gen_climate[['exec_times']][["exec_total_time"]] <- tiempo.sim

    class(gen_climate) <- c(class(gen_climate), 'gamwgen.climate')



    #########################
    ## FINALIZAR EJECUCIÓN ##
    #########################


    ## Cerrar progress bar
    pb$terminate()


    ## Remove parallelization conf, if necessary
    if(remove_parallelization_conf) {
        foreach::registerDoSEQ()
        snow::stopCluster(cluster)
    }


    ## Remove temporary files
    if(control$use_temporary_files_to_save_ram && control$remove_temp_files_used_to_save_ram)
        purrr::walk( nsim_gen_clim %>% dplyr::pull(nsim_gen_climate),
                     function(filename) { file.remove(filename); invisible(gc()) } )


    ## Return result
    gen_climate
}
