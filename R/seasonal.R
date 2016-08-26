
estimate_seasonal_covariates <- function(climate) {
    ## compute regional average seasonal total precip to use as predictor
    d_montot <- climate %>% group_by(station, year = lubridate::year(date), month = lubridate::month(date)) %>%
        summarise(montot = sum(prcp), maxmean = mean(tx), minmean = mean(tn))

    d_avg_montot <- d_montot %>% group_by(year, month) %>%
        summarise(avgmontot = mean(na.omit(montot)),
                  avgmaxmean = mean(na.omit(maxmean)),
                  avgminmean = mean(na.omit(minmean)))

    d_seatot <- d_avg_montot %>% group_by(year, season = ceiling(month/3)) %>%
        summarise(seatot = sum(avgmontot),
                  seamax = mean(avgmaxmean),
                  seamin = mean(avgminmean)) %>%
        arrange(year, season)  # Sort in ascending order.

    # sort( sapply(objects(),function(x){ format(object.size(get(x)), units='MB') }))

    ## create seasonal total vector
    season_len = c(sum(c(31, 28, 31)), sum(c(30, 31, 30)), sum(c(31, 31, 30)), sum(c(31, 30, 31)))
    season_len_leap = c(sum(c(31, 29, 31)), sum(c(30, 31, 30)), sum(c(31, 31, 30)), sum(c(31, 30, 31)))

    season_rainfall_covariates <- list(c(), c(), c(), c())
    season_max_temp_covariates <- list(c(), c(), c(), c())
    season_min_temp_covariates <- list(c(), c(), c(), c())

    # ST1 = ST2 = ST3 = ST4 = c() SMX1 = SMX2 = SMX3 = SMX4 = c() SMN1 = SMN2 = SMN3 = SMN4 = c()

    for (yr in unique(d_seatot$year)) {
        for (season in unique(d_seatot$season)) {
            season_data <- d_seatot[d_seatot$year == yr & d_seatot$season == season, ]
            season_len_data <- if (lubridate::leap_year(yr)) season_len_leap else season_len

            values_before <- if (season > 1)
                rep(0, times = sum(season_len_data[1:(season - 1)])) else NULL
            values_after <- if (season < 4)
                rep(0, times = sum(season_len_data[(season + 1):4])) else NULL

            season_rainfall_covariates[[season]] <- c(season_rainfall_covariates[[season]], values_before, rep(season_data$seatot, times = season_len_data[season]), values_after)

            season_max_temp_covariates[[season]] <- c(season_max_temp_covariates[[season]], values_before, rep(season_data$seamax, times = season_len_data[season]), values_after)

            season_min_temp_covariates[[season]] <- c(season_min_temp_covariates[[season]], values_before, rep(season_data$seamin, times = season_len_data[season]), values_after)
        }
    }

    return(list(tx = season_max_temp_covariates, tn = season_min_temp_covariates, prcp = season_rainfall_covariates))
}


# rm(season, values_after, values_before, yr, season_data, season_len_data, d_montot, d_avg_montot)
