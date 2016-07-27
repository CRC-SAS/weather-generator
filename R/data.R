climate.data.frame <- function(date, station, tx, tn, prcp, ...) {
    as.climate.data.frame(data.frame(date, station, tx, tn, prcp, ...))
}

as.climate.data.frame <- function(x, variables_names = c(date = "date", station = "station", tx = "tx", tn = "tn", prcp = "prcp")) {
    for (climvar in names(variables_names)) {
        print(climvar)
    }
}
