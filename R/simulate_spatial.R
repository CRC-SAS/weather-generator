
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
spatial_simulation_control <- function(nsim = 1, seed = NULL, avbl_cores = 2,
                                       bbox_offset = 100000, sim_loc_as_grid = T,
                                       use_spatially_correlated_noise = T,
                                       use_temporary_files_to_save_ram = T,
                                       remove_temp_files_used_to_save_ram = T,
                                       manage_parallelization_externally = F) {

    prcp_noise_generating_function = glmwgen:::random_field_noise_prcp
    temperature_noise_generating_function = glmwgen:::random_field_noise_temperature

    if(!use_spatially_correlated_noise) {
        prcp_noise_generating_function = glmwgen:::not_spatially_correlated_random_field_noise_prcp
        temperature_noise_generating_function = glmwgen:::not_spatially_correlated_random_field_noise_temperature
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
#' @param model A glmwgen model.
#' @param simulation_locations a sf object with the points at which weather should be simulated.
#'          If not set, the locations used to fit the model will be used.
#' @param start_date a start date in text format (will be converted using as.Date) or a date object.
#' @param end_date an end date in text format (will be converted using as.Date) or a date object.
#' @param control a glmwgen simulation control list.
#' @import dplyr
#' @import foreach
#' @export
spatial_simulation <- function(model, simulation_locations, start_date, end_date,
                               control = glmwgen:::glmwgen_simulation_control(),
                               output_filename = "sim_results.nc",
                               seasonal_climate = NULL, verbose = F) {

    # TODO:
    # 1- quitar rasters de ejecución día a día () --> opcional
    # 2- paralizar prcp_occ y prcp_amt con tmax y tmin (LISTO) --> no se puede porque la cantidad de
    #    procesos crece muy rápido, cada proceso hijo del for each crear dos hijos más, no es factible
    #    porque no hay control sobre la cantidad de procesos a ser creados y por ende sobre la cantidad
    #    de procesadores requeridos (además, crear procesos hijos solo haría todo más lento).
    # 3- agregar covariables (LISTO)
    # 4- cambiar ruidos de monto de precipitacion para no tener montos de lluvia tan altos (LISTO)
    # 5- Controles sobre archivo recibido con covariables, copiar de fit (LISTO).
    # 6- verificar que con el mismo seed se generen los mismos resultados para corridas diferentes (LISTO)
    # 7- Agregar capa raster con los station_id en la celdas que corresponden a las estaciones simuladas!!
    # 8- Preg a Alessio si los ruidos (ahora calculados con un sed por realización), no deben usar seed por día?
    # *- Tratar de cerrar esto lo más pronto posible

    ## Objeto a ser retornado
    gen_climate <- list()

    ###############################################################

    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("foreach"))

    ###############################################################

    if(class(model) != 'glmwgen')
        stop(glue::glue('Received a model of class {class(model)} and a model of class "glmwgen" was expected.'))

    # If
    if (is.null(simulation_locations))
        stop("The parameter simulation_locations can't be null!")

    # Se controlan que los datos recibidos tengan el formato correcto
    glmwgen:::check.simulation.input.data(simulation_locations, seasonal_climate)

    ###############################################################

    if(sf::st_crs(simulation_locations) != sf::st_crs(model$crs_used_to_fit)) {
        simulation_locations %>% simulation_locations %>%
            sf::st_transform(sf::st_crs(model$crs_used_to_fit))
        warning('The crs used to fit and the crs of simulation_locations are not equals. ',
                'Se transforma simulation_locations al crs {del ajuste}')
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

    ## los puntos a simular tienen que estar dentro del boundign box del ajuste
    if(!all(sf::st_contains(polygon_offset, simulation_locations, sparse = F)))
        stop('Alguno de los puntos ajustados se encuentra fuera del bounding box ',
             'creado a partir de los puntos a simular y con un offset de 10 km.')

    ## Solo pueden usarse ruidos NO correlacionados espacialmente
    ## para simular en los puntos utilizados en el ajuste
    if(!control$use_spatially_correlated_noise)
        if(any(lapply(sf::st_equals(simulation_locations, model$stations), length) != 1))
            stop('Los ruidos NO correlacionados espacialmente solo pueden ser usados cuando ',
                 'los puntos a ser simulados son los mismos que fueron utilizados en el ajuste!')


    ###############################################################

    if(!all(sf::st_is_valid(simulation_locations)))
        stop('simulation_locations is not a valid sf object.')

    ###############################################################

    if(end_date <= start_date)
        stop('End date should be greater than start date')

    ###############################################################

    if(control$nsim < 1)
        stop('Number of simulations should be one or greater than one')

    ###############################################################

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
    #  interna      |      interna  -->  me va a pasar Alessio el cod, hay un pequeño cambio
    #  interna      |      externo
    #  externa      |      externa  ->  iguales
    #  externa      |      externa  ->  diferentes (considerar siempre diferentes y listo)
    #  externa      |      interna  (no corresponde)
    # generalmente el ajuste es con covariables internas

    ## Si seasonal_climate != NULL se debió hacer el ajuste usando covariables
    if (!is.null(seasonal_climate) & !model$control$use_covariates)
        stop('El ajuste fue hecho sin covariables, por lo tanto, la simulación ',
             'también debe hacerse sin covariables y no es valido utilizar el ',
             'parámetro seasonal_climate!!')

    ## Si el ajuste se hizo utilizando un archivo de covariables externo,
    ## entonces la simulación también debe hacerse con un archivo externo
    if (model$control$use_external_seasonal_climate & is.null(seasonal_climate))
        stop('El ajuste se hizo utilizando un archivo de covariables externo (parametro ',
             'seasonal_climate), por lo tanto, la simulación también debe hacerse con un ',
             'archivo de covariables externo (parametro seasonal_climate).')

    ## OBS:
    # A pesar de que cuando el ajuste se hace sin covariables, la simulación también debe
    # hacerse sin covariables, y de que cuando el ajuste se hace con covariables, la simulación
    # también debe hacerse con convariables; esto no verifica explícitamente porque se toma
    # directamente el parametro de control use_covariates del ajuste para determinar si la
    # simulación se hace con covariables o se hace sin covariables!!

    ###############################################################

    ## Se verifica que hayan covariables suficientes para cubrir todas las fechas a simular,
    ## pero el control solo se hace si se van a utilizar covariables en la simulación!!
    sim_dates_control <- tidyr::crossing(model$seasonal_data %>% dplyr::distinct(station_id),
                                         year = base::seq.int(lubridate::year(start_date), lubridate::year(end_date)),
                                         season = as.integer(c(1, 2, 3, 4)))
    seasonal_cov_ctrl <- model$seasonal_data %>% dplyr::select(station_id, year, season) %>% dplyr::distinct()
    if (!all(do.call(paste0, sim_dates_control) %in% do.call(paste0, seasonal_cov_ctrl)) & model$control$use_covariates)
        stop("Simulation years aren't in model$seasonal_data!")

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
    if(!is.null(control$seed)) {
        realizations_seeds <- list()
        cant_dias_sim <- as.numeric(end_date - start_date) + 1
        for (r in seq(1, control$nsim, 1)) {
            realizations_seeds[[r]] <- list(general = ceiling(runif(min = 1, max = 10000000, n = control$nsim)), # uno por realizacion
                                            prcp_occ = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim)), # uno por día por realizacion
                                            prcp_amt = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim)), # uno por día por realizacion
                                            temp_dry = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim)), # uno por día por realizacion
                                            temp_wet = ceiling(runif(min = 1, max = 10000000, n = cant_dias_sim))) # uno por día por realizacion
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


    ###################################################################################
    ## Si no se recibe un seasonal_climate externo, se utiliza el generado en el ajuste
    if(is.null(seasonal_climate))
        seasonal_climate <- model$seasonal_data


    ############################################################
    ## Define name of the output netcdf4 file (as absolute path)
    netcdf_filename <- glue::glue("{getwd()}/{output_filename}")


    ####################################
    ## Generate simulation dates
    simulation_dates <-
        tibble::tibble(date = seq.Date(from = as.Date(start_date),
                                       to = as.Date(end_date),
                                       by = "days")) %>%
        dplyr::mutate(year = lubridate::year(date),
                      month = lubridate::month(date),
                      season = lubridate::quarter(date, fiscal_start = 12))
    ## Numbers of days to be simulated
    ndates <- nrow(simulation_dates)


    ##################################
    ## Matriz con los puntos a simular
    simulation_points <- simulation_locations %>%
        sf::st_transform(sf::st_crs(model$crs_used_to_fit)) %>%
        dplyr::mutate(point_id = dplyr::row_number(),
                      longitude = sf::st_coordinates(geometry)[,'X'],
                      latitude  = sf::st_coordinates(geometry)[,'Y'])


    ##################################
    ## Raster con los puntos a simular
    if (control$sim_loc_as_grid) {
        pnts_dist_matrix  <- glmwgen:::make_distance_matrix(simulation_locations)
        min_dist_btw_pnts <- floor(min(pnts_dist_matrix[upper.tri(pnts_dist_matrix)]))
        raster_resolution <- c(min_dist_btw_pnts, min_dist_btw_pnts)
        simulation_raster <- raster::rasterFromXYZ(
            xyz = sf::st_coordinates(simulation_points),
            res = raster_resolution,
            crs = sf::st_crs(simulation_points))
    }
    if(!control$sim_loc_as_grid) {
        stns_dist_matrix  <- glmwgen:::make_distance_matrix(model$stations)
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


    ############################################################################
    ## Obtención de valores para el día previo al día de inicio de la simulación
    start_date_prev_day_climatology <-
        glmwgen:::get_start_climatology(model, simulation_points, start_date, control)


    #############################################################
    ## Obtención de covariables, si van a ser utilizadas, sino no
    if (model$control$use_covariates)
        covariates <- glmwgen:::get_covariates(model, simulation_points, seasonal_climate,
                                               simulation_dates, control)


    #########################################
    ## Global paramteres for noise generators
    if(control$use_spatially_correlated_noise)
        gen_noise_params <- glmwgen:::generate_month_params(
            residuals = model$models_residuals,
            observed_climate = model$models_data,
            stations = model$stations)
    if(!control$use_spatially_correlated_noise)
        gen_noise_params <- glmwgen:::generate_residuals_statistics(
            models_residuals = model$models_residuals)


    ############################################
    ## Matriz de clasificacion de días lluviosos
    clasification_matrix <- matrix(c(-Inf, 0, 0, 0, Inf, 1),
                                   ncol = 3,
                                   byrow = TRUE)


    #######################################################################################
    ## Se crea una matriz de simulación, esta va a contener todos los datos necesarios para
    ## la simulación de cada día a simular, si se utilizan covariables, se las incluye aquí
    if (!model$control$use_covariates) {
        simulation_matrix <- simulation_points
    } else {
        simulation_matrix <- simulation_points %>%
            sf::st_join(covariates %>% {if (control$sim_loc_as_grid) dplyr::select(., dplyr::everything())
                                        else dplyr::select(. , -dplyr::one_of("station_id"))}) %>%
            dplyr::mutate(ST1 = if_else(season == 1, sum_prcp, 0),
                          ST2 = if_else(season == 2, sum_prcp, 0),
                          ST3 = if_else(season == 3, sum_prcp, 0),
                          ST4 = if_else(season == 4, sum_prcp, 0),
                          SN1 = if_else(season == 1, mean_tmin, 0),
                          SN2 = if_else(season == 2, mean_tmin, 0),
                          SN3 = if_else(season == 3, mean_tmin, 0),
                          SN4 = if_else(season == 4, mean_tmin, 0),
                          SX1 = if_else(season == 1, mean_tmax, 0),
                          SX2 = if_else(season == 2, mean_tmax, 0),
                          SX3 = if_else(season == 3, mean_tmax, 0),
                          SX4 = if_else(season == 4, mean_tmax, 0))
    }


    ##########################################
    ## AQUI EMPIEZA REALMENTE LA SIMULACIÓN ##
    ##########################################


    #################################
    ## Control de tiempo de ejecución
    t.sim <- proc.time()


    #######################################################
    ## Set the progress bar to know how long this will take
    pb <- progress::progress_bar$new(
        format = paste0(ifelse(nworkers == 1, " realization:", "finished realizations:"),
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
    glmwgen:::CrearNetCDF(netcdf_filename,
                num_realizations = control$nsim,
                sim_dates = simulation_dates$date,
                simulation_raster = simulation_raster,
                coord_ref_system = sf::st_crs(simulation_points))


    ######
    ##
    #microbenchmark::microbenchmark({
    ctrl_sim <- foreach::foreach(r = 1:control$nsim, .combine = dplyr::bind_rows, .packages = c('dplyr'),
                                 .options.snow = list(progress = progress_pb), .verbose=verbose) %dopar% {

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
        #microbenchmark::microbenchmark({
        simulation_matrix.d <- simulation_matrix %>%
            {if (!model$control$use_covariates) dplyr::filter(.)
             else dplyr::filter(., year == simulation_dates$year[1],
                                   season == simulation_dates$season[1])} %>%
            sf::st_join(start_date_prev_day_climatology) %>%
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
            # 5- tipo día actual, no previo
            dplyr::mutate(tipo_dia = NA_character_) %>%
            # para control de paralelización
            dplyr::mutate(nsim = r)
        #}, times = 100) # 17 milisegundos


        #######################################
        ## Tiempos a tomar por cada realización
        tiempos <- tibble::tibble(tiempo.gen_rast = list(),
                                  tiempo.save_rds = list())
        tiempos <- tiempos %>% dplyr::add_row()


        #################################
        ## Control de tiempo de ejecución
        t.rasters <- proc.time()


        #####################################################################################
        ## Antes se usaba un foreach para paralelizar esto, pero no se puede ser paralelizado
        ## porque simulation_matrix.d no toma los valores correctos al paralelizar!!
        ## Ver version anterior para más detalles (commit: 1898e5a)
        #microbenchmark::microbenchmark({
        rasters <- purrr::map_dfr(1:ndates, function(d) {

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
                                        newdata = simulation_matrix.d,
                                        #cluster = cluster,  # empeora el tiempo para grillas grandes
                                        newdata.guaranteed = TRUE) # una optimizacion
            #}, times = 10) # 480 milisegundos

            # Raster con los valores "climáticos"
            #microbenchmark::microbenchmark({
            SIMocc_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMocc = !!SIMocc) %>%
                glmwgen:::sf2raster('SIMocc', simulation_raster)
            #}, times = 10) # 8 milisegundos

            # Raster con los valores de "ruido"
            #microbenchmark::microbenchmark({
            SIMocc_points_noise.d <- control$prcp_noise_generating_function(
                simulation_points = simulation_points,
                gen_noise_params = gen_noise_params,
                month_number = current_month,
                selector = 'prcp',
                seed = realizations_seeds[[r]]$prcp_occ[[d]]) %>%
                glmwgen:::sf2raster('prcp_residuals', simulation_raster)
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
            simulation_matrix.d <- simulation_matrix.d %>%
                dplyr::mutate(prcp_occ = raster::extract(SIMocc_points.d, simulation_points),
                              tipo_dia = factor(prcp_occ, levels = c(1, 0),
                                                labels = c('Lluvioso', 'Seco')))
            #}, times = 10) # 15 milisegundos



            ########################################
            ## Temperatura (ambas, tmax y tmin) ----

            #  Raster con los valores de "ruido" para días secos
            #microbenchmark::microbenchmark({
            temperature_random_fields_dry <-
                control$temperature_noise_generating_function(
                    simulation_points = simulation_points,
                    gen_noise_params = gen_noise_params,
                    month_number = current_month,
                    selector = c('tmax_dry', 'tmin_dry'),
                    seed = realizations_seeds[[r]]$temp_dry[[d]])
            #}, times = 10) # 180 milisegundos

            # Procesamiento de residuos para dias secos
            #microbenchmark::microbenchmark({
            rasters_secos.d <- purrr::map(
                .x = c("tmin", "tmax"),
                .f = function(variable, objeto_sf) {
                    return (glmwgen:::sf2raster(objeto_sf, paste0(variable, '_residuals'), simulation_raster))
                }, objeto_sf = temperature_random_fields_dry
            )
            #}, times = 10) # 14 milisegundos
            #microbenchmark::microbenchmark({
            names(rasters_secos.d) <- c("tmin", "tmax")
            #}, times = 10) # 2 milisegundos

            # Raster con los valores de "ruido" para días lluviosos
            #microbenchmark::microbenchmark({
            temperature_random_fields_wet <-
                control$temperature_noise_generating_function(
                    simulation_points = simulation_points,
                    gen_noise_params = gen_noise_params,
                    month_number = current_month,
                    selector = c('tmax_wet', 'tmin_wet'),
                    seed = realizations_seeds[[r]]$temp_wet[[d]])
            #}, times = 10) # 180 milisegundos

            # Procesamiento de residuos para dias humedos
            #microbenchmark::microbenchmark({
            rasters_humedos.d <- purrr::map(
                .x = c("tmin", "tmax"),
                .f = function(variable, objeto_sf) {
                    return (glmwgen:::sf2raster(objeto_sf, paste0(variable, '_residuals'), simulation_raster))
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
                                        newdata = simulation_matrix.d,
                                        #cluster = cluster,  # no mejora mucho el tiempo
                                        newdata.guaranteed = TRUE) # una optimizacion
            #}, times = 100) # 540 milisegundos

            # Raster con los valores "climáticos"
            #microbenchmark::microbenchmark({
            SIMmax_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMmax = !!SIMmax) %>%
                glmwgen:::sf2raster('SIMmax', simulation_raster)
            #}, times = 10) # 8 milisegundos

            # Raster con los valores simulados
            #microbenchmark::microbenchmark({
            SIMmax_points.d <- SIMmax_points_climate.d + SIMmax_points_noise.d
            #}, times = 10) # 4 milisegundos

            # Agregar valores de temperatura mínima a los puntos de simulación
            #microbenchmark::microbenchmark({
            simulation_matrix.d <- simulation_matrix.d %>%
                dplyr::mutate(tmax = raster::extract(SIMmax_points.d, simulation_points))
            #}, times = 10) # 15 milisegundos



            #################################
            ## Temperatura Mínima (tmin) ----

            # Simulacion de temperatura mínima
            #microbenchmark::microbenchmark({
            SIMmin <- mgcv::predict.bam(model$fitted_models$tmin_fit,
                                        newdata = simulation_matrix.d,
                                        #cluster = cluster,  # no mejora mucho el tiempo
                                        newdata.guaranteed = TRUE) # una optimizacion
            #}, times = 100) # 580 milisegundos

            # Raster con los valores "climáticos"
            #microbenchmark::microbenchmark({
            SIMmin_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMmin = !!SIMmin) %>%
                glmwgen:::sf2raster('SIMmin', simulation_raster)
            #}, times = 10) # 8 milisegundos

            # Raster con los valores simulados
            #microbenchmark::microbenchmark({
            SIMmin_points.d <- SIMmin_points_climate.d + SIMmin_points_noise.d
            #}, times = 10) # 4 milisegundos

            # Agregar valores de temperatura mínima a los puntos de simulación
            #microbenchmark::microbenchmark({
            simulation_matrix.d <- simulation_matrix.d %>%
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
                                   newdata = simulation_matrix.d,
                                   #cluster = cluster,  # no mejora mucho el tiempo
                                   newdata.guaranteed = TRUE))/alphaamt
            #}, times = 100) # 360 milisegundos

            # Raster con los valores de "ruido"
            #microbenchmark::microbenchmark({
            SIMamt_points_noise.d <- control$prcp_noise_generating_function(
                simulation_points = simulation_points,
                gen_noise_params = gen_noise_params,
                month_number = current_month,
                selector = 'prcp',
                seed = realizations_seeds[[r]]$prcp_amt[[d]]) %>%
                glmwgen:::sf2raster('prcp_residuals', simulation_raster)
            #}, times = 10) # 45 milisegundos

            # Simulacion de montos
            #microbenchmark::microbenchmark({
            SIMamt <- stats::qgamma(stats::pnorm(raster::extract(SIMamt_points_noise.d, simulation_points)),
                             shape = rep(alphaamt, length(betaamt)), scale = betaamt)
            #}, times = 10) # 15 milisegundos

            # Raster con los valores "climáticos"
            #microbenchmark::microbenchmark({
            SIMamt_points_climate.d <- simulation_points %>%
                dplyr::mutate(SIMamt = !!SIMamt) %>%
                glmwgen:::sf2raster('SIMamt', simulation_raster)
            #}, times = 10) # 8 milisegundos

            # Enmascarar pixeles sin ocurrencia de lluvia
            #microbenchmark::microbenchmark({
            SIMamt_points.d <- SIMamt_points_climate.d * SIMocc_points.d
            #}, times = 10) # 4 milisegundos

            # Agregar valores de los montos de prcp a los puntos de simulación
            #microbenchmark::microbenchmark({
            simulation_matrix.d <- simulation_matrix.d %>%
                dplyr::mutate(prcp_amt = raster::extract(SIMamt_points.d, simulation_points))
            #}, times = 10) # 15 milisegundos



            #########################################################################
            ## Preparar simulation_matrix.d para la simulación del siguiente día ----
            current_sim_points  <- simulation_matrix.d %>% sf::st_drop_geometry()
            # OJO: se usa el operador <<- para utilizar los resultados el siguiente día
            simulation_matrix.d <<- simulation_matrix %>%
                {if (!model$control$use_covariates) dplyr::filter(.)
                 else dplyr::filter(., year == simulation_dates$year[1],
                                       season == simulation_dates$season[1])} %>%
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
                    raster_tmin = list(SIMmin_points.d),
                    csim_points = list(current_sim_points)
                ))
        })
        #}, times = 1) # con ndates=365, 926 segundos



        ############################################
        ## Tomar tiempo de generación de los rasters
        tiempos <- dplyr::mutate(tiempos, tiempo.gen_rast = list(proc.time() - t.rasters))



        ########################################
        ## Guardar realizacion en archivo NetCDF
        ## Falla cuando un sub-proceso intenta abrir el archivo,
        ## pero esté ya fue abierto por otro sub-proceso!!!
        # glmwgen:::GuardarRealizacionNetCDF(netcdf_filename,
        #                                    numero_realizacion = r,
        #                                    sim_dates = simulation_dates$date,
        #                                    raster_tmax = rasters$raster_tmax,
        #                                    raster_tmin = rasters$raster_tmin,
        #                                    raster_prcp = rasters$raster_prcp)


        ###################################################################
        ## Se guarda en memoria el tibble con los rasters de la realización
        rasters_path <- paste0(getwd(), "/realizacion_", r ,".rds")
        if(control$use_temporary_files_to_save_ram) {
            t.saveRDS <- proc.time()
            base::saveRDS(rasters, rasters_path)
            tiempos <- dplyr::mutate(tiempos, tiempo.save_rds = list(proc.time() - t.saveRDS))
        }



        ######################
        ## Liberar memoria RAM
        if(control$use_temporary_files_to_save_ram)
            rm(rasters)



        #################
        ## Retorno final
        return (tibble::tibble(nsim = r,
                               rasters = ifelse(control$use_temporary_files_to_save_ram,
                                                list(rasters_path), list(rasters))
                               ) %>% dplyr::bind_cols(tiempos))
    }
    #}, times = 2)


    ########################################
    ## Guardar realizacion en archivo NetCDF
    ctrl_gen <- purrr::map2_dfr(
        ctrl_sim %>% dplyr::pull(nsim),
        ctrl_sim %>% dplyr::pull(rasters),
        function(r, rasters) {

            #######################################
            ## Tiempos a tomar por cada realización
            tiempos <- tibble::tibble(tiempo.read_rds = list(),
                                      tiempo.gen_ncdf = list())
            tiempos <- tiempos %>% dplyr::add_row()

            ######################################
            ## Leer archivos con rasters generados
            if(control$use_temporary_files_to_save_ram) {
                t.read_rds <- proc.time()
                rasters <- base::readRDS(rasters)
                tiempos <- dplyr::mutate(tiempos, tiempo.read_rds = list(proc.time() - t.read_rds))
            }

            ###########################
            ## Generar archivos netcdf4
            t.gen_ncdf <- proc.time()
            glmwgen:::GuardarRealizacionNetCDF(netcdf_filename,
                                               numero_realizacion = r,
                                               sim_dates = simulation_dates$date,
                                               raster_tmax = rasters$raster_tmax,
                                               raster_tmin = rasters$raster_tmin,
                                               raster_prcp = rasters$raster_prcp)
            tiempos <- dplyr::mutate(tiempos, tiempo.gen_ncdf = list(proc.time() - t.gen_ncdf))

            ######################
            ## Liberar memoria RAM
            if(control$use_temporary_files_to_save_ram)
                rm(rasters)

            #################
            ## Retorno final
            return (tibble::tibble(nsim = r) %>% dplyr::bind_cols(tiempos))
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
    gen_climate[['netcdf4_file_with_results']] <- netcdf_filename

    stations <- model$stations;  climate <- model$climate
    save(stations, climate, file = "station_climate.RData")
    gen_climate[['rdata_file_with_station_climate']] <- "station_climate.RData"

    names(ctrl_sim$tiempo.gen_rast) <- paste0("sim_", ctrl_sim$nsim)
    gen_climate[['exec_times']][["gen_rast_time"]] <- ctrl_sim$tiempo.gen_rast

    if(control$use_temporary_files_to_save_ram) {
        names(ctrl_sim$tiempo.save_rds) <- paste0("sim_", ctrl_sim$nsim)
        gen_climate[['exec_times']][["rds_save_time"]] <- ctrl_sim$tiempo.save_rds
    }
    if(control$use_temporary_files_to_save_ram) {
        names(ctrl_gen$tiempo.read_rds) <- paste0("sim_", ctrl_gen$nsim)
        gen_climate[['exec_times']][["rds_read_time"]] <- ctrl_gen$tiempo.read_rds
    }

    names(ctrl_gen$tiempo.gen_ncdf) <- paste0("sim_", ctrl_gen$nsim)
    gen_climate[['exec_times']][["gen_ncdf_time"]] <- ctrl_gen$tiempo.gen_ncdf
    gen_climate[['exec_times']][["exec_tot_time"]] <- tiempo.sim

    class(gen_climate) <- c(class(gen_climate), 'glmwgen.climate')


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
        purrr::walk( ctrl_sim %>% dplyr::pull(rasters),
                     function(filename) { file.remove(filename) } )


    ## Return result
    gen_climate
}
