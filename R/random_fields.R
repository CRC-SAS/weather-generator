

# Definicion de funcion para la generación de campos gaussianos para ocurrencia de temperatura
random_field_noise_temperature <- function(simulation_points, gen_noise_params, month_number, selector, seed) {

    if (length(selector) == 1 && selector == "Wet")
        selector = c('tmax_wet', 'tmin_wet')

    if (length(selector) == 1 && selector == "Dry")
        selector = c('tmax_dry', 'tmin_dry')

    ## OBS:
    # Estamos trabajando con coordenadas métricas, sin embargo,
    # el paquete RandomFields trabaja mejor en metros, por la tanto,
    # pasamos grid a kilometros, dividiendo cada columna por 1000!!
    grid <- simulation_points %>%
        dplyr::select(longitude, latitude) %>%
        sf::st_drop_geometry() %>% tibble::as_tibble() %>%
        dplyr::mutate(longitude = longitude / 1000,
                      latitude = latitude / 1000)

    # Extraer parametros correspondientes a la variable y mes determinado
    month_params_vario_tmax <- gen_noise_params[[month_number]]$variogram_parameters[[selector[1]]]
    month_params_vario_tmin <- gen_noise_params[[month_number]]$variogram_parameters[[selector[2]]]

    # Crear modelo con los parametros
    model <- RandomFields::RMbiwm(nudiag = c(0.5, 0.5),
                                  nured12 = 1,
                                  cdiag = c(month_params_vario_tmax[2], month_params_vario_tmin[2]),
                                  rhored = gen_noise_params[[month_number]]$correlation,
                                  s = c(month_params_vario_tmax[3],
                                        month_params_vario_tmin[3],
                                        0.5*sum(month_params_vario_tmax[3], month_params_vario_tmin[3])))

    # Simular sobre una grilla regular
    campos_simulados <- RandomFields::RFsimulate(model, x = grid, grid = F, coord_sys="cartesian", coordunits="km", seed=seed)

    # Extraer resultados
    # campos.simulados.df <- RandomFields::RFspDataFrame2conventional(campos_simulados, data.frame = T)

    # Crear objeto sf (se vuelve, longitude y latitude, a metros)
    campo <- simulation_points %>%
        dplyr::mutate(id = dplyr::row_number(),
                      tmax_residuals = campos_simulados[, 1],
                      tmin_residuals = campos_simulados[, 2]) %>%
        dplyr::select(tmax_residuals, tmin_residuals)

    # Devolver resultados
    return(campo)

}


# Definicion de funcion para la generación de campos gaussianos para ocurrencia de precipitacion
random_field_noise_prcp <- function(simulation_points, gen_noise_params, month_number, selector, seed) {

    ## OBS:
    # Estamos trabajando con coordenadas métricas, sin embargo,
    # el paquete RandomFields trabaja mejor en metros, por la tanto,
    # pasamos grid a kilometros, dividiendo cada columna por 1000!!
    grid <- simulation_points %>%
        dplyr::select(longitude, latitude) %>%
        sf::st_drop_geometry() %>% tibble::as_tibble() %>%
        dplyr::mutate(longitude = longitude / 1000,
                      latitude = latitude / 1000)

    # Extraer parametros correspondientes a la variable y mes determinado
    month_params_vario_prcp <- gen_noise_params[[month_number]]$variogram_parameters[[selector]]

    # Crear modelo con los parametros
    model <-RandomFields::RMexp(var = month_params_vario_prcp[[2]],
                                scale = month_params_vario_prcp[[3]])

    # Simular sobre una grilla regular
    campos_simulados <- RandomFields::RFsimulate(model, x = grid, grid = F, coord_sys="cartesian", coordunits="km", seed=seed)

    # Extraer resultados
    # campos.simulados.df <- RandomFields::RFspDataFrame2conventional(campos_simulados, data.frame = T)

    # Crear objeto sf
    campo <- simulation_points %>%
        dplyr::mutate(id = dplyr::row_number(),
                      prcp_residuals = campos_simulados) %>%
        dplyr::select(prcp_residuals)

    # Devolver resultados
    return(campo)

}


# ...
not_spatially_correlated_random_field_noise_temperature <- function(simulation_points, gen_noise_params, month_number, selector, seed) {

    # para repetir resultados
    set.seed(seed)

    if (all(endsWith(selector, "_wet")))
        selector = "Wet"

    if (all(endsWith(selector, "_dry")))
        selector = "Dry"

    campos_simulados <- gen_noise_params %>%
        dplyr::filter(month == month_number, type_day == selector) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(ruidos =
            list(MASS::mvrnorm(mu = c(mean.tmax_residuals, mean.tmin_residuals),
                               Sigma = matrix(c(var.tmax_residuals, cov.residuals, cov.residuals, var.tmin_residuals), 2, 2),
                               empirical = FALSE) %>% #Alessio, falla con TRUE, podemos usar FALSE??
                     magrittr::set_names(c("tmax_noise", "tmin_noise"))
                 )) %>%
        dplyr::mutate(tmax_noise = ruidos["tmax_noise"],
                      tmin_noise = ruidos["tmin_noise"])

    # Crear objeto sf
    campo <- simulation_points %>%
        dplyr::left_join(campos_simulados, by = "station_id") %>%
        dplyr::rename(tmax_residuals = tmax_noise, tmin_residuals = tmin_noise) %>%
        dplyr::select(tmax_residuals, tmin_residuals)

    return (campo)
}


not_spatially_correlated_random_field_noise_prcp <- function(simulation_points, gen_noise_params, month_number, selector, seed) {

    # para repetir resultados
    set.seed(seed)

    campos_simulados <- simulation_points %>% dplyr::rowwise() %>%
        dplyr::mutate(prcp_residuals = stats::rnorm(n = 1, mean = 0, sd = 1)) %>%
        dplyr::select(station_id, prcp_residuals)

    # Crear objeto sf
    campo <- simulation_points %>%
        sf::st_join(.) %>%
        #dplyr::left_join(campos_simulados, by = "station_id") %>%
        dplyr::select(prcp_residuals)

    return (campo)
}

##############
## OLD METHODS

noise_method <- proto::proto(
    new = function(self, model, simulation_locations, control, noise_function = no_noise) {
        cache <- NULL
        cache_next_index <- NULL
        month_cache_max_size <- NULL

        if('cache_size' %in% names(control) && control$cache_size > 1) {
            month_cache_max_size <- control$cache_size

            cache <- array(data = NA,
                           dim = c(3, 12, month_cache_max_size, nrow(simulation_locations)))
            dimnames(cache)[[1]] <- c('tx', 'tn', 'prcp')
            dimnames(cache)[[2]] <- base::month.abb
            dimnames(cache)[[4]] <- rownames(simulation_locations)

            cache_next_index <- array(data = control$cache_size + 1,
                                      dim = c(3, 12))
            dimnames(cache_next_index)[[1]] <- c('tx', 'tn', 'prcp')
            dimnames(cache_next_index)[[2]] <- base::month.abb
        }

        wrapped_noise_f <- function(self, ...) noise_function(...)


        proto::proto(model = model, grid = simulation_locations,
                     control = control, noise_function = noise_function,
                     internal_noise_function = wrapped_noise_f,
                     cache = cache, cache_next_index = cache_next_index,
                     month_cache_max_size = month_cache_max_size)
    },

    generate_noise = function(self, month_number, var_name) {
        if(is.null(self$cache)) return(drop(self$internal_noise_function(self$model, self$grid, month_number, var_name)))

        self$cache_next_index[var_name, month_number] <- self$cache_next_index[var_name, month_number] + 1

        if(self$cache_next_index[var_name, month_number] > self$month_cache_max_size) {
            # Repopulate the cache for this variable and month.
            self$cache[var_name, month_number, , ] <- self$internal_noise_function(self$model, self$grid, month_number, var_name, nsim = self$month_cache_max_size)
            self$cache_next_index[var_name, month_number] <- 1
        }
        return(self$cache[var_name, month_number, self$cache_next_index[var_name, month_number], ])
    }
)

random_field_noise <- function(model, simulation_locations, month_number, var_name, nsim = 1) {
    month_params <- model$month_params[[month_number]]$variogram_parameters[[var_name]]
    rf_model     <- RandomFields::RMexp(var=month_params[2], scale=month_params[3])
    rf_simulate  <- RandomFields::RFsimulate(rf_model,
                                             distances = model$simulation_dist_matrix,
                                             dim = nrow(simulation_locations),
                                             n = nsim,
                                             printlevel = 0)
    return (t(rf_simulate))
}


rnorm_noise <- function(model, simulation_locations, month_number, var_name, nsim = 1) {
    sapply(model$month_params[[month_number]]$residuals_sd[rownames(simulation_locations), var_name], rnorm, mean = 0, n = nsim)
}


mvrnorm_noise <- function(model, simulation_locations, month_number, var_name, nsim = 1) {
    MASS::mvrnorm(mu = rep(0, nrow(simulation_locations)),
                  Sigma = model$month_params[[month_number]]$cov_matrix[[var_name]][rownames(simulation_locations), rownames(simulation_locations)],
                  n = nsim)
}


no_noise <- function(model, simulation_locations, month_number, var_name, nsim = 1) {
    return(0)
}

