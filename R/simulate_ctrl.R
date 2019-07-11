partial <- function(f, ...) {
    l <- list(...)
    function(...) {
        do.call(f, c(list(...), l))
    }
}

# get_seasonal_covariate <- function(year, season_number, season_values, quantile_probs = c(0, 0.33, 0.66, 1), levels_probabilities = c(1/3, 1/3, 1/3)) {
#     splitted <- cut(season_values, breaks = quantile(season_values, probs = quantile_probs), include.lowest = T)
#     values_level <- season_values[splitted == sample(levels(splitted), 1, prob = levels_probabilities)]
#     if(length(values_level) == 1) return(values_level)
#     return(base::sample(values_level, 1))
# }

get_seasonal_covariate <- function(years, season_number, seasonal_values, variable) {
    seasonal_values <- seasonal_values %>% filter(season == season_number) %>% select(year, value = ends_with(variable))
    if(all(years %in% seasonal_values$year)) {
        return(seasonal_values[match(years, seasonal_values$year), 'value'])
    }
    warning('Seasonal values not found for some years, starting series at the first year')
    return(seasonal_values$value[1:length(years)])
}

get_rainfall_seasonal_covariate <- function(years, season_number, seasonal_values) {
    get_seasonal_covariate(years, season_number, seasonal_values, 'prcp')
}

# get_temperatures_seasonal_covariate <- function(year, season_number, season_tx_values, season_tn_values, quantile_probs = c(0, 0.33, 0.66, 1),
#                                                levels_probabilities = c(0.25, 0.35, 0.4)) {
#     splitted_tx <- cut(season_tx_values, breaks = quantile(season_tx_values, probs = quantile_probs), include.lowest = T)
#     splitted_tn <- cut(season_tn_values, breaks = quantile(season_tn_values, probs = quantile_probs), include.lowest = T)
#
#     values_level_index <- sample(seq_along(levels(splitted_tx)), 1, prob = levels_probabilities)
#     tx_sampled_level <- season_tx_values[splitted_tx == levels(splitted_tx)[values_level_index]]
#     tn_sampled_level <- season_tn_values[splitted_tn == levels(splitted_tn)[values_level_index]]
#     return(list(tx = ifelse(length(tx_sampled_level) > 1, sample(tx_sampled_level, 1), tx_sampled_level),
#                 tn = ifelse(length(tn_sampled_level) > 1, sample(tn_sampled_level, 1), tn_sampled_level)))
# }
get_temperatures_seasonal_covariate <- function(years, season_number, seasonal_values) {
    return(list(tx = get_seasonal_covariate(years, season_number, seasonal_values, 'tx'),
                tn = get_seasonal_covariate(years, season_number, seasonal_values, 'tn')))
}

#' @title Simulations control configuration
#' @description Provides fine control of different parameters that will be used to create new weather series.
#' @export
glmwgen_simulation_control <- function(seasonal_temps_covariates_getter = get_temperatures_seasonal_covariate,
                                       seasonal_prcp_covariates_getter = get_rainfall_seasonal_covariate,
                                       random_fields_method = 'rf',
                                       minimum_temperatures_difference_threshold = 0.1,
                                       Rt = NULL,
                                       # grf_method = NULL,
                                       multicombine = F,
                                       always_krig_coefficients = T,
                                       interpolation_method = c('idw', 'kriging'),
                                       cache_size = 1) {
    rf_function <- NULL
    if(!is.function(random_fields_method)) {
        # if(startsWith(random_fields_method, 'gauss')) rf_function <- gaussian_random_field
        if(startsWith(random_fields_method, 'rf')) rf_function <- glmwgen:::random_field_noise
        if(startsWith(random_fields_method, 'chol')) rf_function <- glmwgen:::cholesky_random_field
        if(startsWith(random_fields_method, 'rn')) rf_function <- glmwgen:::rnorm_noise
        if(startsWith(random_fields_method, 'mv')) rf_function <- glmwgen:::mvrnorm_noise
    } else {
        rf_function <- random_fields_method
    }

    interpolation_method <- match.arg(interpolation_method)

    interpolation_method <- c('kriging' = glmwgen:::krige_covariate_automap,
                              'idw' = glmwgen:::idw_covariate)[[interpolation_method]]

    return(list(seasonal_temps_covariates_getter = seasonal_temps_covariates_getter,
                seasonal_prcp_covariates_getter = seasonal_prcp_covariates_getter,
                random_fields_method = rf_function,
                Rt = Rt,
                # grf_method = grf_method,
                always_krig_coefficients = always_krig_coefficients,
                multicombine = multicombine,
                interpolation_method = interpolation_method,
                cache_size = cache_size,
                minimum_temperatures_difference_threshold = minimum_temperatures_difference_threshold))
}
