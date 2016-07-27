gaussian_random_field <- function(model, simulation_locations, control, month_number, var_name) {
    month_params <- model$month_params[[month_number]][[var_name]]
    suppressWarnings(geoR::grf(nrow(simulation_locations),
                               grid = simulation_locations,
                               cov.model = "exponential",
                               cov.pars = c(month_params[2], month_params[3]),
                               nugget = month_params[1],
                               mean = 0,
                               # method = control$grf_method,
                               messages = FALSE)$data)
}

cholesky_random_field <- function(model, simulation_locations, control, month_number, var_name) {
    model$cSigma[[month_number]][[var_name]] %*% rnorm(nrow(simulation_locations)) * sqrt(model$month_params[[month_number]][[var_name]][2])
}
