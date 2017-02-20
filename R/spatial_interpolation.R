# krige_covariate <- function(model, station_coordinates_grid, simulation_coordinates_grid, covariate) {
#     variogram <- fields::vgram(loc = station_coordinates_grid, y = covariate,
#                                breaks = seq(0, max(model$distance_matrix),
#                                             by = max(model$distance_matrix) / 30))
#     variogram <- data.frame(y = variogram$stats[2, ], x = variogram$centers)
#
#     krigged_coefs <- tryCatch({
#         # variogram_fit <- with(variogram,
#         #                       stats::nls(formula = y ~ rho * (1 - exp(-x / theta)),
#         #                                  start = list(rho = max(y, na.rm = TRUE), theta = max(dist_values)/2),
#         #                                  control = list(maxiter = 500, minFactor = 1/10000, warnOnly = F)))
#         variogram_fit <- with(variogram,
#                               stats::nls(formula = y ~ rho * (1 - exp(-x / theta)),
#                                          algorithm = "port",
#                                          start = list(rho = max(y, na.rm = TRUE), theta = max(model$distance_matrix)/2),
#                                          lower = list(rho = 0, theta = 0),
#                                          upper = list(rho = max(y, na.rm = TRUE) * 1.5, theta = max(model$distance_matrix) * 1.5)))
#         krig_model <- fields::Krig(station_coordinates_grid, covariate, give.warnings = F,
#                                    sigma2 = 0, rho = coefficients(variogram_fit)['rho'], theta = coefficients(variogram_fit)['theta'])
#         fields::predict.Krig(krig_model, simulation_coordinates_grid)
#
#     }, error = function(...) {
#         fields::predict.Krig(fields::Krig(station_coordinates_grid, covariate, give.warnings = F), simulation_coordinates_grid)
#     })
# }


krige_covariate_automap <- function(model, stations_locations, simulation_locations, covariate) {
    covariate <- SpatialPointsDataFrame(coordinates(stations_locations),
                                        data.frame(value = covariate),
                                        proj4string = stations_locations@proj4string)
    automap::autoKrige(value ~ 1, input_data = covariate, new_data = simulation_locations, debug.level = 0)$krige_output$var1.pred
}

idw_covariate <- function(model, stations_locations, simulation_locations, covariate) {
    drop(model$idw_weights %*% covariate)
}


## TODO: cleanup.

# estimate model coefficients on grid using ordinary kriging (OK) occurrence

# dist_function <- ifelse(is_projected, 'rdist', 'rdist.earth')
#
#
#         krig_covariate <- function(covariate) {
#             # krig_model <- fields::Krig(coordinates(model$stations), covariate, give.warnings = F, Distance = dist_function)
#             krig_model <- fields::Krig(station_coordinates_grid, covariate, give.warnings = F)
#             fields::predict.Krig(krig_model, simulation_coordinates_grid)
#         }
#
#         krig_gstat <- function(covariate) {
#             covariate <- data.frame(covariate, coordinates(model$stations))
#             colnames(covariate) <- c('value', 'x', 'y')
#             coordinates(covariate) <- ~ x + y
#             covariate@proj4string <- model$stations@proj4string
#
#             max_range <- max(model$distance_matrix)
#
#             vgm_model <- c('Exp', 'Gau', 'Log')
#
#             result <- NULL
#             model_index <- 0
#
#             krige_internal <- function(vgm_model) {
#                 vario <- gstat::variogram(value ~ 1, data = covariate, cutoff = max_range, width = max_range / 20)
#                 # vario <- gstat::variogram(value ~ 1, data = covariate)
#                 vario_fit <- gstat::fit.variogram(vario, gstat::vgm(model = vgm_model[model_index],
#                                                                     psill = mean(vario$gamma, na.rm = T),
#                                                                     range = max_range / 10))
#                 krige_fit <- gstat::krige(value ~ 1, covariate, krige_locations, model = vario_fit, debug.level = 0)
#                 krige_fit
#             }
#
#             while(is.null(result) && model_index < length(vgm_model)) {
#                 model_index <- model_index + 1
#                 result <- tryCatch(krige_internal(vgm_model), warning = function(w) {
#                     NULL
#                 })
#             }
#
#             if(is.null(result)) {
#                 warning('Failed to build variogram')
#                 return(NA)
#             } else {
#                 covariate_idx <- apply(spDists(covariate, result), 1, which.min)
#                 return(result[covariate_idx, ]$var1.pred)
#             }
#         }
#
#         # krig_automap <- function(covariate) {
#         #     covariate <- data.frame(covariate, coordinates(model$stations))
#         #     colnames(covariate) <- c('value', 'x', 'y')
#         #     coordinates(covariate) <- ~ x + y
#         #     covariate@proj4string <- model$stations@proj4string
#         #     automap::autoKrige(value ~ 1, covariate, krige_locations)
#         # }
#
#         spline_covariate <- function(covariate) {
#             spline_model <- fields::Tps(coordinates(model$stations), covariate, give.warnings = F, lon.lat = !is_projected, miles = F)
#             fields::predict.Krig(spline_model, simulation_coordinates)
#         }
#
#         dist_values <- do.call(c, sapply(1:nrow(model$distance_matrix), function(r) model$distance_matrix[r, (r:nrow(model$distance_matrix))]))
#         dist_values <- dist_values[dist_values > 0]
#         #
#         # covariate_sp <- station_sp_gk
#         # covariate_sp@data$value <- covariate
#         # automap::autofitVariogram(value ~ 1, covariate_sp)
#
#         krig_covariate_2 <- function(covariate) {
#             # variogram <- fields::vgram(loc = station_coordinates_grid, d = dist_values, y = covariate,
#             #                            breaks = seq(0, max(dist_values), by = max(dist_values) / 20))
#             variogram <- fields::vgram(loc = station_coordinates_grid, y = covariate,
#                                        breaks = seq(0, max(model$distance_matrix), by = max(model$distance_matrix) / 30))
#             variogram <- data.frame(y = variogram$stats[2, ], x = variogram$centers)
#
#             krigged_coefs <- tryCatch({
#                 # variogram_fit <- with(variogram,
#                 #                       stats::nls(formula = y ~ rho * (1 - exp(-x / theta)),
#                 #                                  start = list(rho = max(y, na.rm = TRUE), theta = max(dist_values)/2),
#                 #                                  control = list(maxiter = 500, minFactor = 1/10000, warnOnly = F)))
#                 variogram_fit <- with(variogram,
#                                       stats::nls(formula = y ~ rho * (1 - exp(-x / theta)),
#                                                  algorithm = "port",
#                                                  start = list(rho = max(y, na.rm = TRUE), theta = max(model$distance_matrix)/2),
#                                                  lower = list(rho = 0, theta = 0),
#                                                  upper = list(rho = max(y, na.rm = TRUE) * 1.5, theta = max(model$distance_matrix) * 1.5)))
#                 krig_model <- fields::Krig(station_coordinates_grid, covariate, give.warnings = F,
#                                            sigma2 = 0, rho = coefficients(variogram_fit)['rho'], theta = coefficients(variogram_fit)['theta'])
#                 fields::predict.Krig(krig_model, simulation_coordinates_grid)
#
#             }, error = function(...) {
#                 fields::predict.Krig(fields::Krig(station_coordinates_grid, covariate, give.warnings = F), simulation_coordinates_grid)
#             })
#
#         }
