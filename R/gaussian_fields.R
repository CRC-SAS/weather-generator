
# Definicion de funcion para la generacion de campos gaussianos para alterar los valores iniciales
setting_variograms_initial_values <- function(simulation_matrix, distance_matrix,
                                              day_of_year = 365, grid = NULL, seed = seed) {

    # Matriz de valores iniciales para inicializar la simulacion

    prcp_intial_values_interpolation <- simulation_matrix %>%
        dplyr::filter(., doy == day_of_year) %>%
        dplyr::select(., station_id, date, prcp_occ_prev) %>%
        tidyr::spread(., key = station_id, value = prcp_occ_prev) %>%
        dplyr::select(., colnames(distance_matrix))

    tmax_intial_values_interpolation <- simulation_matrix %>%
        dplyr::filter(., doy == day_of_year) %>%
        dplyr::select(., station_id, date, tmax_prev) %>%
        tidyr::spread(., key = station_id, value = tmax_prev) %>%
        dplyr::select(., colnames(distance_matrix))

    tmin_intial_values_interpolation <- simulation_matrix %>%
        dplyr::filter(., doy == day_of_year) %>%
        dplyr::select(., station_id, date, tmin_prev) %>%
        tidyr::spread(., key = station_id, value = tmin_prev) %>%
        dplyr::select(., colnames(distance_matrix))


    prcp_cor   <- stats::cor(prcp_intial_values_interpolation, use = "pairwise.complete")
    prcp_vario <- stats::var(prcp_intial_values_interpolation, use = "pairwise.complete") * (1 - prcp_cor)
    tmax_vario <- stats::cov(tmax_intial_values_interpolation, use = "pairwise.complete.obs")
    tmin_vario <- stats::cov(tmin_intial_values_interpolation, use = "pairwise.complete.obs")

    prcp_params <- stats::optim(par = c(0.01, 1, max(distance_matrix)), f = glmwgen:::partially_apply_LS(prcp_vario, distance_matrix))$par
    prcp_params[prcp_params < 0] <- 0
    prcp_params <- c(0, prcp_params[2], prcp_params[3])

    sill_initial_value_tmax <- mean(var(tmax_intial_values_interpolation, na.rm = T))
    tmax_params <- stats::optim(par = c(sill_initial_value_tmax, max(distance_matrix)), fn = glmwgen:::partially_apply_LS(tmax_vario, distance_matrix, base_p = c(0)))$par
    tmax_params <- c(0, tmax_params)

    sill_initial_value_tmin <- mean(var(tmin_intial_values_interpolation, na.rm = T))
    tmin_params <- stats::optim(par = c(sill_initial_value_tmin, max(distance_matrix)), fn = glmwgen:::partially_apply_LS(tmin_vario, distance_matrix, base_p = c(0)))$par
    tmin_params <- c(0, tmin_params)

    variogram_parameters <- list(prcp = prcp_params,
                                 tmax = tmax_params,
                                 tmin = tmin_params)

    models_parameters <- list(prcp = RandomFields::RMexp(var = variogram_parameters$prcp[2],
                                                         scale = variogram_parameters$prcp[3]),
                              tmax = RandomFields::RMexp(var = variogram_parameters$tmax[2],
                                                         scale = variogram_parameters$tmax[3]),
                              tmin = RandomFields::RMexp(var = variogram_parameters$tmin[2],
                                                         scale = variogram_parameters$tmin[3]))

    random_fields <- list(prcp = RandomFields::RFsimulate(models$prcp, x = grid, grid = F, seed = seed),
                          tmax = RandomFields::RFsimulate(models$tmax, x = grid, grid = F, seed = seed),
                          tmin = RandomFields::RFsimulate(models$tmin, x = grid, grid = F, seed = seed))

    return(list(variogram = variogram_parameters,
                models = models_parameters,
                random_fields = random_fields))

}
