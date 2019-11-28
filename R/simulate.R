
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
#' @param coord_ref_system ...
#' @param avbl_cores ...
#' @export
glmwgen_simulation_control <- function(nsim = 1, seed = NULL, avbl_cores = 2,
                                       bbox_offset = 100000,
                                       use_spatially_correlated_noise = T) {

    prcp_noise_genrating_function = glmwgen:::random_field_noise_prcp
    temperature_noise_genrating_function = glmwgen:::random_field_noise_temperature

    # if(!use_spatially_correlated_noise) {
    #     prcp_noise_genrating_function = glmwgen:::random_field_noise_prcp
    #     temperature_noise_genrating_function = glmwgen:::not_spatially_correlated_random_field_noise_temperature
    # }

    return(list(nsim = nsim, seed = seed, avbl_cores = avbl_cores, bbox_offset = bbox_offset,
                use_spatially_correlated_noise = use_spatially_correlated_noise,
                prcp_noise_genrating_function = prcp_noise_genrating_function,
                temperature_noise_genrating_function = temperature_noise_genrating_function))
}



#' @title Simulates new weather trajectories in stations
#' @description Simulates new weather trajectories.
#' @param model A glmwgen model.
#' @param simulation_locations a SpatialPoints object with the points at which weather should be simulated.
#'          If not set, the locations used to fit the model will be used.
#' @param start_date a start date in text format (will be converted using as.Date) or a date object.
#' @param end_date an end date in text format (will be converted using as.Date) or a date object.
#' @param control a glmwgen simulation control list.
#' @import dplyr
#' @import foreach
#' @export
sim.glmwgen <- function(model, simulation_locations, start_date, end_date,
                        control = glmwgen:::glmwgen_simulation_control()) {

    ## Objeto a ser retornado
    gen_climate <- list()

    ###############################################################

    if(class(model) != 'glmwgen')
        stop(paste('Received a model of class', class(model), 'and a model of class "glmwgen" was expected.'))

    glmwgen:::check.simulation.input.data(simulation_locations)

    ###############################################################

    if(sf::st_crs(simulation_locations) != sf::st_crs(model$crs_used_to_fit)) {
        simulation_locations %>% simulation_locations %>%
            sf::st_transform(sf::st_crs(model$crs_used_to_fit))
        warning('The crs used to fit and the crs of simulation_locations are not equals. Se transformar simulation_locations al crs {del ajuste}')
    }

    ## Crear bounding_box con un buffer de 10km y verificar que los puntos ajustados estén dentro!!
    sl_bbox <- sf::st_bbox(model$stations)

    sl_bbox_offset <- sf::st_bbox(
        obj = c(xmin = sl_bbox[['xmin']] - control$bbox_offset,
                ymin = sl_bbox[['ymin']] - control$bbox_offset,
                xmax = sl_bbox[['xmax']] + control$bbox_offset,
                ymax = sl_bbox[['ymax']] + control$bbox_offset),
        crs = sf::st_crs(model$crs_used_to_fit))

    polygon_offset <- sf::st_as_sfc(sl_bbox_offset)
    sf::st_crs(polygon_offset) <- sf::st_crs(model$crs_used_to_fit)

    # los puntos a simular tienen que estar dentro del boundign box del ajuste
    if(!all(sf::st_contains(polygon_offset, simulation_locations, sparse = F)))
        stop('Alguno de los puntos ajustados se encuentra fuera del bounding box ',
             'creado a partir de los puntos a simular y con un offset de 10 km.')

    ## Solo pueden usarse ruidos NO correlacionados espacialmente
    ## para simular en los puntos utilizados en el ajuste
    if(!control$use_spatially_correlated_noise)
        if(any(!diag(sf::st_intersects(simulation_locations, model$stations, sparse = F))))
            stop('Los ruidos NO correlacionados espacialmente solo pueden ser usados cuando ',
                 'los puntos a ser simulados son los mismos que fueron utilizados en el ajuste!')


    ###############################################################

    if(!all(sf::st_is_valid(simulation_locations)))
        stop('simulation_locations is not a  valid sf object.')

    ###############################################################

    if(end_date <= start_date)
        stop('End date should be greater than start date')

    ###############################################################

    years_in_sim_dates <- base::seq.int(lubridate::year(start_date), lubridate::year(end_date))
    years_in_senal_cov <- dplyr::distinct(model$seasonal_data, year) %>% dplyr::pull()
    if (!all(is.element(years_in_sim_dates, years_in_senal_cov)))
        stop("Simulation years aren't in model$seasonal_data!")

    ###############################################################

    if(control$nsim < 1)
        stop('Number of simulations should be one or greater than one')

    ###############################################################



    ############################
    ## INITIAL CONFIGURATIONS ##
    ############################

    # Para que las funciones de RandomFields devuelvan lo esperado!!
    RandomFields::RFoptions(spConform=FALSE)

    # Para ...
    if(!is.null(control$seed))
        set.seed(control$seed)

    # Para ...
    realizations_seeds <- NULL
    if(!is.null(control$seed))
        realizations_seeds <- ceiling(runif(min = 1, max = 10000000, n = control$nsim))



    ####################################
    ## Parallelization initialization ##
    ####################################


    ## Variable that indicate if it's necessary to remove
    ## the parallelization configuration
    remove_parallelization_conf <- F

    ## Register a sequential backend if the user didn't register a parallel
    ## in order to avoid a warning message if we use %dopar%.
    if(!foreach::getDoParRegistered()) {
        if(is.na(control$avbl_cores) | is.null(control$avbl_cores)) {
            foreach::registerDoSEQ()
        } else if (control$avbl_cores <= 1) {
            foreach::registerDoSEQ()
        } else if (control$avbl_cores > 1) {
            remove_parallelization_conf <- T
            #doMC::registerDoMC(control$avbl_cores)
            cluster <- parallel::makeCluster(control$avbl_cores)
            doSNOW::registerDoSNOW(cluster)
            #doParallel::registerDoParallel(cl)
        }
    }



    ##################################
    ## PREPARACIÓN DE LA SIMULACIÓN ##
    ##################################


    ####################################
    ## Generate simulation dates
    simulation_dates <-
        tibble::tibble(date = seq.Date(from = as.Date(start_date),
                                       to = as.Date(end_date),
                                       by = "days"))
    ## Numbers of days to be simulated
    ndates <- nrow(simulation_dates)


    ####################################
    ## Matriz con los puntos a simular
    simulation_points <- simulation_locations %>%
        sf::st_transform(model$crs_used_to_fit) %>%
        dplyr::mutate(point_id = dplyr::row_number(),
                      longitude = sf::st_coordinates(geometry)[,'X'],
                      latitude  = sf::st_coordinates(geometry)[,'Y'])
    ## Tibble to be used as grid for many functions
    tibble_of_points <- tibble::as_tibble(simulation_points) %>%
        dplyr::select(longitude, latitude)


    #####################################
    ## Interpolación de valores iniciales
    if(control$use_spatially_correlated_noise)
        interpolated_data_start_date_prev <-
            glmwgen:::interpolate_month_day(model, simulation_points,
                                            seed = control$seed,
                                            month = lubridate::month(start_date-1),
                                            day = lubridate::day(start_date-1))


    ###########################
    ## Global months paramteres
    if(control$use_spatially_correlated_noise)
        month_params <- glmwgen:::generate_month_params(
            residuals = model$models_residuals,
            observed_climate = model$models_data,
            stations = model$stations)


    ##############################
    ## Global residuals statistics
    if(!control$use_spatially_correlated_noise)
        residuals_statistics <- glmwgen:::generate_residuals_statistics(
            residuals = model$models_residuals,
            observed_climate = model$models_data
        )

    ############################################
    ## Matriz de clasificacion de días lluviosos
    clasification_matrix <- matrix(c(-Inf, 0, 0, 0, Inf, 1),
                                   ncol = 3,
                                   byrow = TRUE)


    ##########################################
    ## AQUI EMPIEZA REALMENTE LA SIMULACIÓN ##
    ##########################################


    #################################
    ## Control de tiempo de ejecución
    t.sim <- proc.time()


    #######################################################
    ## Set the progress bar to know how long this will take
    pb <- progress::progress_bar$new(
        format = paste(" realization: :r /", control$nsim ,
                       " | progress: :percent (in :elapsed) | eta: :eta"),
        total = control$nsim, clear = TRUE, width= 80, show_after = 0)

    invisible(pb$tick(0, tokens = list(r = 1)))

    progress_pb <- function(r) {
        pb$tick(1, tokens = list(r = r))
    }


    ###########################################
    ## Crear objeto para guardar los resultados
    glmwgen:::CrearNetCDF('prueba.nc',
                num_realizations = control$nsim,
                simulation_points = simulation_points,
                sim_dates = simulation_dates$date)


    ######
    ##
    #microbenchmark::microbenchmark({
    ctrl_sim <- foreach::foreach(r = 1:control$nsim, .combine = dplyr::bind_rows, .packages = c('dplyr'),
                                 .options.snow = list(progress = progress_pb)) %dopar% {

        ######################################################################
        ## Para que las funciones de RandomFields devuelvan lo esperado!! ----
        RandomFields::RFoptions(spConform=FALSE)


        #################################################################################
        ## Creacion de los puntos de simulacion para el dia i (eq. daily covariates) ----
        simulation_points.d <- simulation_points %>%
            sf::st_join(interpolated_data_start_date_prev) %>%
            dplyr::rename(prcp_occ_prev = prcp_occ,
                          tmax_prev = tmax,
                          tmin_prev = tmin) %>%
            # Columnas para poder ejecutar los predict
            dplyr::mutate(date = simulation_dates$date[1],
                          time = as.numeric(date)/1000,
                          doy = 1) %>%
            # 1- colnames(model$fitted_models$prcp_occ_fit$model)
            dplyr::mutate(tipo_dia_prev = factor(prcp_occ_prev, levels = c(1, 0),
                                                 labels = c('Lluvioso', 'Seco')),
                          prcp_occ = NA_integer_) %>%
            # 2- colnames(model$fitted_models$prcp_amt_fit$Jan$model)
            dplyr::mutate(prcp_amt = NA_real_) %>%
            # 3- colnames(model$fitted_models$tmax_fit$model)
            dplyr::mutate(tmax = NA_real_) %>%
            # 4- colnames(model$fitted_models$tmin_fit$model)
            dplyr::mutate(tmin = NA_real_) %>%
            # para control de paralelización
            dplyr::mutate(nsim = r)


        #################################
        ## Control de tiempo de ejecución
        t.rasters <- proc.time()


        #####################################################################################
        ## Antes se usaba un foreach para paralelizar esto, pero no se puede ser paralelizado
        ## porque simulation_points.d no toma los valores correctos al paralelizar!!
        ## Ver version anterior para más detalles (commit: 1898e5a)
        #microbenchmark::microbenchmark({
        rasters <- purrr::map_dfr(1:ndates, function(d) {
            ##################################################
            ## Control de tiempo de ejecución para un día ----
            #microbenchmark::microbenchmark({


            #######################
            ## Indice de meses ----
            #microbenchmark::microbenchmark({
            current_month <- lubridate::month(simulation_dates$date[d])
            #}, times = 10) # 40 milisegundos


            #######################################
            ## Ocurrencia de lluvia (prcp_occ) ----

            # Simulacion de ocurrencia de lluvia
            #microbenchmark::microbenchmark({
            SIMocc <- mgcv::predict.bam(model$fitted_models$prcp_occ_fit,
                                        newdata = simulation_points.d,
                                        #cluster = cluster,  # empeora el tiempo para grillas grandes
                                        newdata.guaranteed = TRUE) # una optimizacion
            #}, times = 10) # 480 milisegundos

            # Raster con los valores "climáticos"
            #microbenchmark::microbenchmark({
            SIMocc_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMocc = !!SIMocc) %>%
                glmwgen:::sf2raster('SIMocc')
            #}, times = 10) # 8 milisegundos

            # Raster con los valores de "ruido"
            #microbenchmark::microbenchmark({
            SIMocc_points_noise.d <- glmwgen:::random_field_noise_prcp(
                month_parameters = month_params,
                grid = tibble_of_points,
                month_number = current_month,
                var_name = 'prcp',
                coord_ref_system = model$crs_used_to_fit,
                seed = realizations_seeds[d]) %>%
                glmwgen:::sf2raster('prcp_residuals')
            #}, times = 10) # 45 milisegundos

            # Raster con los valores simulados
            #microbenchmark::microbenchmark({
            SIMocc_points.d <- SIMocc_points_climate.d + SIMocc_points_noise.d
            #}, times = 10) # 4 milisegundos

            # Raster con los valores simulados reclasificados a 0 y/o 1
            #microbenchmark::microbenchmark({
            SIMocc_points.d <- raster::reclassify(SIMocc_points.d, clasification_matrix)
            #}, times = 10) # 4 milisegundos

            # Agregar valores de ocurrencia a la grilla de simulacion
            #microbenchmark::microbenchmark({
            simulation_points.d <- simulation_points.d %>%
                dplyr::mutate(prcp_occ = raster::extract(SIMocc_points.d, simulation_points),
                              tipo_dia = factor(prcp_occ, levels = c(1, 0),
                                                labels = c('Lluvioso', 'Seco')))
            #}, times = 10) # 15 milisegundos



            ########################################
            ## Temperatura (ambas, tmax y tmin) ----

            #  Raster con los valores de "ruido" para días secos
            #microbenchmark::microbenchmark({
            temperature_random_fields_dry <-
                glmwgen:::random_field_noise_temperature(
                    month_parameters = month_params,
                    month_number = current_month,
                    grid = tibble_of_points,
                    var_name = c('tmax_dry', 'tmin_dry'),
                    coord_ref_system = sf::st_crs(model$crs_used_to_fit),
                    seed = realizations_seeds[d])
            #}, times = 10) # 180 milisegundos

            # Procesamiento de residuos para dias secos
            #microbenchmark::microbenchmark({
            rasters_secos.d <- purrr::map(
                .x = c("tmin", "tmax"),
                .f = function(variable, objeto_sf) {
                    return (glmwgen:::sf2raster(objeto_sf, paste0(variable, '_residuals')))
                }, objeto_sf = temperature_random_fields_dry
            )
            #}, times = 10) # 14 milisegundos
            #microbenchmark::microbenchmark({
            names(rasters_secos.d) <- c("tmin", "tmax")
            #}, times = 10) # 2 milisegundos

            # Raster con los valores de "ruido" para días lluviosos
            #microbenchmark::microbenchmark({
            temperature_random_fields_wet <-
                glmwgen:::random_field_noise_temperature(
                    month_parameters = month_params,
                    month_number = current_month,
                    grid = tibble_of_points,
                    var_name = c('tmax_wet', 'tmin_wet'),
                    coord_ref_system = sf::st_crs(model$crs_used_to_fit),
                    seed = realizations_seeds[d])
            #}, times = 10) # 180 milisegundos

            # Procesamiento de residuos para dias humedos
            #microbenchmark::microbenchmark({
            rasters_humedos.d <- purrr::map(
                .x = c("tmin", "tmax"),
                .f = function(variable, objeto_sf) {
                    return (glmwgen:::sf2raster(objeto_sf, paste0(variable, '_residuals')))
                }, objeto_sf = temperature_random_fields_wet
            )
            #}, times = 10) # 14 milisegundos
            #microbenchmark::microbenchmark({
            names(rasters_humedos.d) <- c("tmin", "tmax")
            #}, times = 10) # 2 milisegundos

            #Ahora vamos a generar 2 rasters: uno para tmax y otro para tmin
            #Cada raster tiene los residuos de las variables correspondientes
            #considerando la ocurrencia de dia seco o humedo
            #microbenchmark::microbenchmark({
            SIMmax_points_noise.d <-
                glmwgen:::ensamblar_raster_residuos(rasters_humedos.d$tmax, rasters_secos.d$tmax, SIMocc_points.d)
            #}, times = 10) # 15 milisegundos
            #microbenchmark::microbenchmark({
            SIMmin_points_noise.d <-
                glmwgen:::ensamblar_raster_residuos(rasters_humedos.d$tmin, rasters_secos.d$tmin, SIMocc_points.d)
            #}, times = 10) # 15 milisegundos



            #################################
            ## Temperatura Máxima (tmax) ----

            # Simulacion de temperatura máxima
            #microbenchmark::microbenchmark({
            SIMmax <- mgcv::predict.bam(model$fitted_models$tmax_fit,
                                        newdata = simulation_points.d,
                                        #cluster = cluster,  # no mejora mucho el tiempo
                                        newdata.guaranteed = TRUE) # una optimizacion
            #}, times = 100) # 540 milisegundos

            # Raster con los valores "climáticos"
            #microbenchmark::microbenchmark({
            SIMmax_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMmax = !!SIMmax) %>%
                glmwgen:::sf2raster('SIMmax')
            #}, times = 10) # 8 milisegundos

            # Raster con los valores simulados
            #microbenchmark::microbenchmark({
            SIMmax_points.d <- SIMmax_points_climate.d + SIMmax_points_noise.d
            #}, times = 10) # 4 milisegundos

            # Agregar valores de temperatura mínima a los puntos de simulación
            #microbenchmark::microbenchmark({
            simulation_points.d <- simulation_points.d %>%
                dplyr::mutate(tmax = raster::extract(SIMmax_points.d, simulation_points))
            #}, times = 10) # 15 milisegundos



            #################################
            ## Temperatura Mínima (tmin) ----

            # Simulacion de temperatura mínima
            #microbenchmark::microbenchmark({
            SIMmin <- mgcv::predict.bam(model$fitted_models$tmin_fit,
                                        newdata = simulation_points.d,
                                        #cluster = cluster,  # no mejora mucho el tiempo
                                        newdata.guaranteed = TRUE) # una optimizacion
            #}, times = 100) # 580 milisegundos

            # Raster con los valores "climáticos"
            #microbenchmark::microbenchmark({
            SIMmin_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMmin = !!SIMmin) %>%
                glmwgen:::sf2raster('SIMmin')
            #}, times = 10) # 8 milisegundos

            # Raster con los valores simulados
            #microbenchmark::microbenchmark({
            SIMmin_points.d <- SIMmin_points_climate.d + SIMmin_points_noise.d
            #}, times = 10) # 4 milisegundos

            # Agregar valores de temperatura mínima a los puntos de simulación
            #microbenchmark::microbenchmark({
            simulation_points.d <- simulation_points.d %>%
                dplyr::mutate(tmin = raster::extract(SIMmin_points.d, simulation_points))
            #}, times = 10) # 15 milisegundos



            ########################################
            ## Montos de lluvia (prcp_amt) ----

            # Filtrar el modelo a usar por el mes en curso
            #microbenchmark::microbenchmark({
            prcp_amt_fit <- model$fitted_models$prcp_amt_fit[[current_month]]
            #}, times = 10) # 5 milisegundos

            # Estimación del parametro de forma
            #microbenchmark::microbenchmark({
            alphaamt  <- MASS::gamma.shape(prcp_amt_fit)$alpha
            #}, times = 10) # 20 milisegundos
            # Estimación de los parametros de escala
            #microbenchmark::microbenchmark({
            betaamt <- base::exp(mgcv::predict.bam(prcp_amt_fit,
                                   newdata = simulation_points.d,
                                   #cluster = cluster,  # no mejora mucho el tiempo
                                   newdata.guaranteed = TRUE))/alphaamt
            #}, times = 100) # 360 milisegundos

            # Simulacion de montos
            #microbenchmark::microbenchmark({
            SIMamt <- stats::qgamma(stats::pnorm(raster::extract(SIMocc_points_noise.d, simulation_points)),
                             shape = rep(alphaamt, length(betaamt)), scale = betaamt)
            #}, times = 10) # 15 milisegundos

            # Raster con los valores "climáticos"
            #microbenchmark::microbenchmark({
            SIMamt_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMamt = !!SIMamt) %>%
                glmwgen:::sf2raster('SIMamt')
            #}, times = 10) # 8 milisegundos

            # Enmascarar pixeles sin ocurrencia de lluvia
            #microbenchmark::microbenchmark({
            SIMamt_points.d <- SIMamt_points_climate.d * SIMocc_points.d
            #}, times = 10) # 4 milisegundos

            # Agregar valores de los montos de prcp a los puntos de simulación
            #microbenchmark::microbenchmark({
            simulation_points.d <- simulation_points.d %>%
                dplyr::mutate(prcp_amt = raster::extract(SIMamt_points.d, simulation_points))
            #}, times = 10) # 15 milisegundos



            #########################################################################
            ## Preparar simulation_points.d para la simulación del siguiente día ----
            current_sim_points  <- simulation_points.d
            # OJO: se usa el operador <<- para utilizar los resultados el siguiente día
            simulation_points.d <<- simulation_points.d %>%
                dplyr::mutate(date = simulation_dates$date[d+1],
                              doy = lubridate::yday(date),
                              time = as.numeric(date)/1000,
                              prcp_occ = NA_integer_,
                              tipo_dia = NA_character_,
                              prcp_amt = NA_real_,
                              tmax = NA_real_,
                              tmin = NA_real_,
                              prcp_occ_prev = raster::extract(SIMocc_points.d, simulation_points),
                              tipo_dia_prev = factor(prcp_occ_prev, levels = c('0', '1'),
                                                     labels = c("Seco", "Lluvioso")),
                              tmax_prev = raster::extract(SIMmax_points.d, simulation_points),
                              tmin_prev = raster::extract(SIMmin_points.d, simulation_points),
                              nsim = r)



            #################################################
            ## Control de tiempo de ejecución para un día ----
            #}, times = 10) # 2.5 segundos



            ###########################
            ## Devolver resultados ----
            return (
                tibble::tibble(
                    date = simulation_dates$date[d],
                    raster_prcp = list(SIMocc_points.d),
                    raster_tmax = list(SIMmax_points.d),
                    raster_tmin = list(SIMmin_points.d),
                    csim_points = list(current_sim_points)
                ))
        })
        #}, times = 1) # con ndates=365, 926 segundos



        #################################
        ## Control de tiempo de ejecución
        tiempo.rasters <- proc.time() - t.rasters



        ##############
        # Progress Bar
        #pb$tick(1, tokens = list(r = r, d = d))


        return (tibble::tibble(nsim = r, rasters = list(rasters),
                               tiempo = list(tiempo.rasters)))
    }
    #}, times = 2)


    ########################################
    ## Guardar realizacion en archivo NetCDF
    purrr::walk2(
        ctrl_sim %>% dplyr::pull(nsim),
        ctrl_sim %>% dplyr::pull(rasters),
        function(r, rasters) {
            glmwgen:::GuardarRealizacionNetCDF('prueba.nc',
                                               numero_realizacion = r,
                                               sim_dates = simulation_dates$date,
                                               raster_tmax = rasters$raster_tmax,
                                               raster_tmin = rasters$raster_tmin,
                                               raster_prcp = rasters$raster_prcp)
        })


    #################################
    ## Control de tiempo de ejecución
    tiempo.sim <- proc.time() - t.sim


    ###########################
    ## Preparar datos de salida

    gen_climate[['nsim']] <- control$nsim
    gen_climate[['seed']] <- control$seed
    gen_climate[['realizations_seeds']] <- realizations_seeds
    gen_climate[['simulation_coordinates']] <- simulation_points

    names(ctrl_sim$tiempo) <- paste0("sim_", ctrl_sim$nsim)
    gen_climate[['exec_times']][["gen_rast_time"]] <- ctrl_sim$tiempo
    gen_climate[['exec_times']][["exec_tot_time"]] <- tiempo.sim

    class(gen_climate) <- c(class(gen_climate), 'glmwgen.climate')



    #########################
    ## FINALIZAR EJECUCIÓN ##
    #########################


    ## Remove parallelization conf, if necessary
    if(remove_parallelization_conf) {
        foreach::registerDoSEQ()
        parallel::stopCluster(cluster)
    }


    ## Return result
    gen_climate
}
