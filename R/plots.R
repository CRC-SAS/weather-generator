
plot_means <- function(in_simulated_climate, observed_climate = NULL, stations = NULL, title = NULL) {
    if(is.null(stations)) stations <- attributes(in_simulated_climate)$model$stations

    if(is.null(observed_climate)) observed_climate <- attributes(in_simulated_climate)$model$start_climatology

    simulation_locations <- sp::SpatialPoints(attributes(in_simulated_climate)$simulation_coordinates,
                                              proj4string = attributes(in_simulated_climate)$model$stations_proj4string)

    stations_dist <- sp::spDists(stations, simulation_locations)

    nearest_simulation_location <- apply(stations_dist, 1, which.min)
    names(nearest_simulation_location) <- stations$id


    observed_summarised_climate <- observed_climate %>% group_by(station) %>%
        summarise(tx = mean(tx, na.rm = T), tn = mean(tn, na.rm = T), prcp = sum(prcp, na.rm = T))

    # If the amount of stations is lower than the amount of simulation points,
    # filter only the matching locations.
    if(length(nearest_simulation_location) < dim(in_simulated_climate)[3]) {
        in_simulated_climate <- in_simulated_climate[, , nearest_simulation_location, ]
    }
    dimnames(in_simulated_climate)[[3]] <- names(nearest_simulation_location)

    summarised_simulated_climate <- apply(in_simulated_climate, c(1, 3, 4), mean)
    simulated_stations <- dimnames(in_simulated_climate)[[3]]
    observed_summarised_climate <- data.frame(observed_summarised_climate %>% filter(station %in% simulated_stations))
    observed_summarised_climate$station <- as.character(observed_summarised_climate$station)

    realization_means <- do.call(rbind, lapply(dimnames(summarised_simulated_climate)[[3]], function(weather_variable) {
        values <- summarised_simulated_climate[, , weather_variable]
        do.call(rbind, lapply(colnames(values), function(column_name) {
            cbind(station = column_name, variable = weather_variable, value = values[, column_name])
        }))
    }))

    realization_means <- data.frame(realization_means, stringsAsFactors = F)
    realization_means$value <- as.numeric(realization_means$value)

    realization_means <- realization_means %>% arrange(station, variable)
    realization_means$x <- rep(1:dim(in_simulated_climate)[1], length(simulated_stations) * length(unique(realization_means$variable)))

    weather_variables_colours <- c('tx' = '#DD2222', 'tn' = '#22DD22')

    require(ggplot2)
    require(gridExtra)

    tx_plot <- ggplot(realization_means %>% filter(variable == 'tx')) +
        geom_boxplot(aes(x = station, y = value)) +
        geom_point(data = observed_summarised_climate, aes(x = station, y = tx), colour = weather_variables_colours[['tx']]) +
        scale_y_continuous(name = 'Mean Maximum Temperature') +
        scale_x_discrete(name = NULL) +
        ggtitle(label = title) +
        theme_light()

    tn_plot <- ggplot(realization_means %>% filter(variable == 'tn')) +
        geom_boxplot(aes(x = station, y = value)) +
        geom_point(data = observed_summarised_climate, aes(x = station, y = tn), colour = weather_variables_colours[['tn']]) +
        scale_y_continuous(name = 'Mean Minimum Temperature') +
        scale_x_discrete(name = 'Weather Station ID') +
        theme_light()

    gridExtra::grid.arrange(tx_plot, tn_plot, ncol = 1)
}

plot_monthly_means <- function(in_simulated_climate, observed_climate) {
    observed_summarised_climate <- observed_climate %>% group_by(station, lubridate::month(date)) %>%
        summarise(tx = mean(tx, na.rm = T), tn = mean(tn, na.rm = T), prcp = sum(prcp, na.rm = T))

    summarised_simulated_climate <- apply(in_simulated_climate, c(1, 3, 4), mean)
    simulated_stations <- dimnames(in_simulated_climate)[[3]]
    observed_summarised_climate <- data.frame(observed_summarised_climate %>% filter(station %in% simulated_stations))
    observed_summarised_climate$station <- as.character(observed_summarised_climate$station)

    realization_means <- do.call(rbind, lapply(dimnames(summarised_simulated_climate)[[3]], function(weather_variable) {
        values <- summarised_simulated_climate[, , weather_variable]
        do.call(rbind, lapply(colnames(values), function(column_name) {
            cbind(station = column_name, variable = weather_variable, value = values[, column_name])
        }))
    }))

    realization_means <- data.frame(realization_means, stringsAsFactors = F)
    realization_means$value <- as.numeric(realization_means$value)

    realization_means <- realization_means %>% arrange(station, variable)
    realization_means$x <- rep(1:dim(in_simulated_climate)[1], length(simulated_stations) * length(unique(realization_means$variable)))

    weather_variables_colours <- c('tx' = '#DD2222', 'tn' = '#22DD22')

    require(ggplot2)
    # plot(ggplot(observed_summarised_climate) +
    #     geom_hline(aes(yintercept = tx), colour = 'red', alpha = 0.5) +
    #     geom_hline(aes(yintercept = tx * 1.01), lty = 2, alpha = 0.5) +
    #     geom_hline(aes(yintercept = tx * 0.99), lty = 2, alpha = 0.5) +
    #     geom_point(data = realization_means %>% filter(variable == 'tx'),
    #                aes(x = x, y = value, colour = variable), alpha = 0.5) +
    #     scale_color_manual(values = weather_variables_colours, guide = F) +
    #     scale_x_discrete(name = 'Realization N°') +
    #     scale_y_continuous(name = 'Mean Maximum Temperature') +
    #     facet_wrap(~ station))
    #
    # plot(ggplot(observed_summarised_climate) +
    #     geom_hline(aes(yintercept = tn), colour = 'green', alpha = 0.5) +
    #     geom_hline(aes(yintercept = tn * 1.01), lty = 2, alpha = 0.5) +
    #     geom_hline(aes(yintercept = tn * 0.99), lty = 2, alpha = 0.5) +
    #     geom_point(data = realization_means %>% filter(variable == 'tn'),
    #                aes(x = x, y = value, colour = variable), alpha = 0.5) +
    #     scale_color_manual(values = weather_variables_colours, guide = F) +
    #     scale_x_discrete(name = 'Realization N°') +
    #     scale_y_continuous(name = 'Mean Minimum Temperature') +
    #     theme_light() +
    #     facet_wrap(~ station))

    require(gridExtra)

    tx_plot <- ggplot(realization_means %>% filter(variable == 'tx')) +
        geom_boxplot(aes(x = station, y = value)) +
        geom_point(data = observed_summarised_climate, aes(x = station, y = tx), colour = weather_variables_colours[['tx']]) +
        scale_y_continuous(name = 'Mean Maximum Temperature') +
        scale_x_discrete(name = NULL) +
        theme_light()

    tn_plot <- ggplot(realization_means %>% filter(variable == 'tn')) +
        geom_boxplot(aes(x = station, y = value)) +
        geom_point(data = observed_summarised_climate, aes(x = station, y = tn), colour = weather_variables_colours[['tn']]) +
        scale_y_continuous(name = 'Mean Minimum Temperature') +
        scale_x_discrete(name = 'Weather Station ID') +
        theme_light()

    gridExtra::grid.arrange(tx_plot, tn_plot, ncol = 1)
}



plot_seasonal_means <- function(in_simulated_climate, observed_climate, show_boxplot = F, weather_variable = 'tx', ...) {

    simulated_dates <- as.Date(dimnames(in_simulated_climate)[[2]])

    season_name <- c('JFM', 'AMJ', 'JAS', 'OND')

    ## Boxplot de realizaciones y clima observado.
    ## Los valores están promediados para todas las estaciones simuladas.
    clima_observado <- observed_climate[observed_climate$date %in% simulated_dates, ]

    clima_por_estacion <- as.data.frame(clima_observado %>%
                                            group_by(año = lubridate::year(date), estacion = ((lubridate::month(date) - 0.5) %/% 3)+1) %>%
                                            summarise(tx = mean(tx, na.rm = T), tn = mean(tn, na.rm = T), prcp = sum(prcp, na.rm = T)))
    clima_por_estacion$año_estacion <- factor(paste(clima_por_estacion$año, clima_por_estacion$estacion, sep = '-'))

    summary_f <- c('tx' = mean, 'tn' = mean, 'prcp' = sum)
    groupped_df <- data.frame()

    for(r_num in 1:(dim(in_simulated_climate)[1])) {
        groupped <- as.data.frame(data.frame(r_num = r_num, idx = 1:dim(in_simulated_climate)[2], date = clima_observado$date) %>%
                                      group_by(r_num, año = lubridate::year(date), estacion = ((lubridate::month(date) - 0.5) %/% 3)+1) %>%
                                      summarise(idx = list(unique(idx))))
        groupped$v <- NA
        for(r_idx in 1:nrow(groupped)) {
            groupped[r_idx, 'v'] <- summary_f[[weather_variable]](in_simulated_climate[r_num, unlist(groupped[r_idx, 'idx']), , weather_variable], na.rm = T)
        }

        groupped_df <- rbind(groupped_df, groupped)
    }

    groupped_df$año <- factor(groupped_df$año)
    groupped_df$estacion <- factor(groupped_df$estacion)
    groupped_df$año_estacion <- paste(groupped_df$año, groupped_df$estacion, sep = '-')


    if(show_boxplot) {
        ## Boxplot de realizaciones y clima observado
        boxplot(v ~ año_estacion, data = groupped_df, main = sprintf('Variable: %s', weather_variable))
        lines(clima_por_estacion[, weather_variable], lwd = 1.5, type = 'l', col = rgb(0.9, 0.2, 0.2),
              ylab = weather_variable)
        points(clima_por_estacion[, weather_variable], lwd = 2, type = 'p', col = rgb(1, 0, 0), pch = 20, cex = 0.9)
    }


    ## Gráfico de diferencias
    simulated_seasonal_mean <- groupped_df %>% group_by(año_estacion) %>% summarise(v = mean(v)) %>% arrange(año_estacion)


    observed_seasonal_mean <- clima_por_estacion %>% arrange(año_estacion)
    observed_seasonal_mean <- observed_seasonal_mean[, c('año_estacion', weather_variable)]

    seasonal_means_diff <- observed_seasonal_mean
    seasonal_means_diff[, weather_variable] <- simulated_seasonal_mean[, 'v'] - seasonal_means_diff[, weather_variable]

    plot(x = ordered(seasonal_means_diff$año_estacion), y = seasonal_means_diff[, weather_variable], type = 'p',
         main = sprintf('Variable: %s', weather_variable), ...)
    lines(x = seasonal_means_diff$año_estacion, y = seasonal_means_diff[, weather_variable])
    abline(b = 0, a = mean(seasonal_means_diff[, weather_variable]), lty = 2, col = 'blue')
    abline(b = 0, a = 0, col = 'red')
}


