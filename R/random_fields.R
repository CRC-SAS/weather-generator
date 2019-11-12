


# Definicion de funcion para la generación de campos gaussianos para ocurrencia de temperatura
random_field_noise_temperature <- function(month_parameters, grid, month_number, var_name, crs_value) {

    # Extraer parametros correspondientes a la variable y mes determinado
    month_params_vario_tmax <- month_parameters[[month_number]]$variogram_parameters[[var_name[1]]]
    month_params_vario_tmin <- month_parameters[[month_number]]$variogram_parameters[[var_name[2]]]

    # Crear modelo con los parametros
    model <- RandomFields::RMbiwm(nudiag = c(0.5, 0.5),
                                  nured12 = 1,
                                  cdiag = c(month_params_vario_tmax[2], month_params_vario_tmin[2]),
                                  rhored = 1,
                                  s=c(month_params_vario_tmax[3],
                                      month_params_vario_tmin[3],
                                      0.5*sum(month_params_vario_tmax[3], month_params_vario_tmin[3])))

    # Simular sobre una grilla regular
    campos_simulados <- RandomFields::RFsimulate(model, x = grid, grid = F)

    # Extraer resultados
    # campos.simulados.df <- RandomFields::RFspDataFrame2conventional(campos_simulados, data.frame = T)

    # Crear objeto sf
    campo <- grid %>%
        dplyr::mutate(x_coord = x_coord*1000, y_coord = y_coord*1000) %>%
        sf::st_as_sf(coords = c('x_coord', 'y_coord')) %>%
        sf::st_set_crs(value = crs_value) %>%
        dplyr::mutate(id = 1:n(),
                      tmax_residuals = campos_simulados[, 1],
                      tmin_residuals = campos_simulados[, 2]) %>%
        dplyr::select(tmax_residuals, tmin_residuals)

    # Devolver resultados
    return(campo)

}

# Definicion de funcion para la generación de campos gaussianos para ocurrencia de precipitacion
random_field_noise_prcp <- function(month_parameters, grid, month_number, var_name, crs_value) {

    # Extraer parametros correspondientes a la variable y mes determinado
    month_params_vario_prcp <- month_parameters[[month_number]]$variogram_parameters[[var_name]]

    # Crear modelo con los parametros
    model <-RandomFields::RMexp(var = month_params_vario_prcp[[2]],
                                scale = month_params_vario_prcp[[3]])

    # Simular sobre una grilla regular
    campos_simulados <- RandomFields::RFsimulate(model, x = grid, grid = F)

    # Extraer resultados
    # campos.simulados.df <- RandomFields::RFspDataFrame2conventional(campos_simulados, data.frame = T)

    # Crear objeto sf
    campo <- grid %>%
        dplyr::mutate(x_coord = x_coord*1000, y_coord = y_coord*1000) %>%
        sf::st_as_sf(coords = c('x_coord', 'y_coord')) %>%
        sf::st_set_crs(value = crs_value) %>%
        dply::mutate(id = 1:n(),
                     prcp_residuals = campos_simulados) %>%
        dplyr::select(prcp_residuals)

    # Devolver resultados
    return(campo)

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
        # ret_values <- self$cache[var_name, month_number, self$cache_next_index[var_name, month_number], ]
        # names(ret_values) <- dimnames(self$cache)[[4]]
        # return(ret_values)
        return(self$cache[var_name, month_number, self$cache_next_index[var_name, month_number], ])
    }
)
#
# make_noise_proto <- function(self, model, simulation_locations, control, noise_function = no_noise) {
#     cache <- NULL
#     cache_next_index <- NULL
#     month_cache_max_size <- NULL
#
#     if('cache_size' %in% names(control) && control$cache_size > 1) {
#         month_cache_max_size <- control$cache_size
#
#         cache <- array(data = NA,
#                        dim = c(3, 12, month_cache_max_size, nrow(simulation_locations)))
#         dimnames(cache)[[1]] <- c('tx', 'tn', 'prcp')
#         dimnames(cache)[[2]] <- base::month.abb
#         dimnames(cache)[[4]] <- rownames(simulation_locations)
#
#         cache_next_index <- array(data = control$cache_size + 1,
#                                   dim = c(3, 12))
#         dimnames(cache_next_index)[[1]] <- c('tx', 'tn', 'prcp')
#         dimnames(cache_next_index)[[2]] <- base::month.abb
#     }
#
#     wrapped_noise_f <- function(self, ...) noise_function(...)
#
#
#     proto::proto(model = model, grid = simulation_locations,
#                  control = control, noise_function = noise_function,
#                  internal_noise_function = wrapped_noise_f,
#                  cache = cache, cache_next_index = cache_next_index,
#                  month_cache_max_size = month_cache_max_size,
#                  generate_noise = function(self, month_number, var_name) {
#                      if(is.null(self$cache)) return(drop(self$internal_noise_function(self$model, self$grid, month_number, var_name)))
#
#
#                      self$cache_next_index[var_name, month_number] <- self$cache_next_index[var_name, month_number] + 1
#
#                      if(self$cache_next_index[var_name, month_number] > self$month_cache_max_size) {
#                          # Repopulate the cache for this variable and month.
#                          self$cache[var_name, month_number, , ] <- self$internal_noise_function(self$model, self$grid, month_number, var_name, nsim = self$month_cache_max_size)
#                          self$cache_next_index[var_name, month_number] <- 1
#                      }
#
#                      ret_values <- self$cache[var_name, month_number, self$cache_next_index[var_name, month_number], ]
#                      names(ret_values) <- dimnames(self$cache)[[4]]
#                      return(ret_values)
#                      # return(self$cache[var_name, month_number, self$cache_next_index[var_name, month_number], ])
#                  })
# }


#
# noise_generator <- noise_method$new(glmwgen_fit, glmwgen_fit$distance_matrix, glmwgen_simulation_control(cache_size = 1000), noise_function = gaussian_random_field)
# noise_generator$cache_next_index['tx', ]
# #
# noise_generator$generate_noise(2, 'tx')
#
# require(microbenchmark)
# microbenchmark(noise_generator$generate_noise(3, 'prcp'), times = 1500)

# gaussian_random_field <- function(model, simulation_locations, month_number, var_name, nsim = 1) {
#     month_params <- model$month_params[[month_number]]$variogram_parameters[[var_name]]
#     generated_noise <- t(suppressWarnings(geoR::grf(nrow(simulation_locations),
#                                           grid = simulation_locations,
#                                           cov.model = "exponential",
#                                           cov.pars = c(month_params[2], month_params[3]),
#                                           nugget = month_params[1],
#                                           mean = 0,
#                                           # method = control$grf_method,
#                                           messages = F,
#                                           nsim = nsim)$data))
#     colnames(generated_noise) <- rownames(simulation_locations)
#     # generated_noise - apply(generated_noise, 1, mean)
#     generated_noise
# }

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

# random_field_noise <- function(model, simulation_locations, month_number, var_name, nsim = 1) {
#     month_params <- model$month_params[[month_number]]$variogram_parameters[[var_name]]
#     t(RandomFields::RFsimulate(model = RandomFields::RMexp(var = month_params[2], scale = month_params[3]),
#                              x = simulation_locations,
#                              n = nsim,
#                              printlevel = 0)@data)
# }
#
# random_field_noise_2 <- function(model, simulation_locations, month_number, var_name, nsim = 1) {
#     month_params <- model$month_params[[month_number]]$variogram_parameters[[var_name]]
#     replicate(nsim, RandomFields::RFsimulate(model = RandomFields::RMexp(var = month_params[2], scale = month_params[3]),
#                                x = coordinates(simulation_locations),
#                                printlevel = 0)@data)
# }
#
# microbenchmark(
#     nsim = random_field_noise(model, simulation_locations, 1, 'tx', 10),
#     rep = random_field_noise_2(model, simulation_locations, 1, 'tx', 10),
#     times = 10)
#
# RandomFields::RFoptions(seed = 12454565, spConform = F)
#
# dist_noise <- RandomFields::RFsimulate(model = RandomFields::RMexp(var = 8, scale = 6.179372e+05),
#                          distances = d_object,
#                          dim = nrow(simulation_coordinates),
#                          n = 10)
#
# coords_noise <- RandomFields::RFsimulate(model = RandomFields::RMexp(var = 8, scale = 6.179372e+05),
#                                          x = simulation_locations,
#                                          n = 10)
# cholesky_random_field <- function(model, simulation_locations, month_number, var_name, nsim = 1) {
#     model$cSigma[[month_number]][[var_name]] %*% rnorm(nrow(simulation_locations)) * sqrt(model$month_params[[month_number]]$variogram_parameters[[var_name]][2])
# }

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

# # file MASS/R/mvrnorm.R
# # copyright (C) 1994-2015 W. N. Venables and B. D. Ripley
# #
# #  This program is free software; you can redistribute it and/or modify
# #  it under the terms of the GNU General Public License as published by
# #  the Free Software Foundation; either version 2 or 3 of the License
# #  (at your option).
# #
# #  This program is distributed in the hope that it will be useful,
# #  but WITHOUT ANY WARRANTY; without even the implied warranty of
# #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# #  GNU General Public License for more details.
# #
# #  A copy of the GNU General Public License is available at
# #  http://www.r-project.org/Licenses/
# #
# mvrnorm <-
#     function(n = 1, mu, Sigma, tol=1e-6, empirical = FALSE, EISPACK = FALSE)
# {
#     p <- length(mu)
#     if(!all(dim(Sigma) == c(p,p))) stop("incompatible arguments")
#     if(EISPACK) stop("'EISPACK' is no longer supported by R", domain = NA)
#     eS <- eigen(Sigma, symmetric = TRUE)
#     ev <- eS$values
#     if(!all(ev >= -tol*abs(ev[1L]))) stop("'Sigma' is not positive definite")
#     X <- matrix(rnorm(p * n), n)
#     if(empirical) {
#         X <- scale(X, TRUE, FALSE) # remove means
#         X <- X %*% svd(X, nu = 0)$v # rotate to PCs
#         X <- scale(X, FALSE, TRUE) # rescale PCs to unit variance
#     }
#     X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
#     nm <- names(mu)
#     if(is.null(nm) && !is.null(dn <- dimnames(Sigma))) nm <- dn[[1L]]
#     dimnames(X) <- list(nm, NULL)
#     if(n == 1) drop(X) else t(X)
# }
