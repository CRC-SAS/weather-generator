
# ----------------------------------------------------------------------------- #
# ---- Funciones mensuales ----
# ----------------------------------------------------------------------------- #

# Definicion de funcion para graficos mensuales de Temperatura maxima
PearsonsResidualsMonthlyMaxTemp <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    residals.pearson <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.pearson = resid(model, type = 'pearson')) %>%
        dplyr::mutate(., month = lubridate::month(date, label = T)) %>%
        dplyr::group_by(., month) %>%
        dplyr::summarise(., mean = mean(residuos.pearson, na.rm = T),
                         std = sd(residuos.pearson, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = month, y = mean,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::labs(x = "Month", y = "Mean Pearson's residuals", title = "Maximum Temperature monthly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::theme_bw()

    # Standard deviation
    # Desvios
    pp <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = month, y = std,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = residals.pearson) +
        ggplot2::labs(x = "Month", y = "Std.Dev. Pearson's residuals", title = "Maximum Temperature monthly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(residals.pearson$std))-1, ceiling(mean(residals.pearson$std))+1), breaks = seq(from = floor(mean(residals.pearson$std))-1, to = ceiling(mean(residals.pearson$std))+1, 0.5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/Tx"))) {
            dir.create(paste0(path, "/", stn_name, "/Tx"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/Tx/Mean.Pearsons.Residuals.monthly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/Tx/Std.Pearsons.Residuals.monthly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# Definicion de funcion para graficos mensuales de Temperatura minima
PearsonsResidualsMonthlyMinTemp <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    residals.pearson <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.pearson = resid(model, type = 'pearson')) %>%
        dplyr::mutate(., month = lubridate::month(date, label = T)) %>%
        dplyr::group_by(., month) %>%
        dplyr::summarise(., mean = mean(residuos.pearson, na.rm = T),
                         std = sd(residuos.pearson, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = month, y = mean,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::labs(x = "Month", y = "Mean Pearson's residuals", title = "Maximum Temperature monthly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::theme_bw()

    # Standard deviation
    pp <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = month, y = std,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = residals.pearson) +
        ggplot2::labs(x = "Month", y = "Std.Dev. Pearson's residuals", title = "Minimum Temperature monthly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(residals.pearson$std))-1, ceiling(mean(residals.pearson$std))+1), breaks = seq(from = floor(mean(residals.pearson$std))-1, to = ceiling(mean(residals.pearson$std))+1, 0.5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/Tn"))) {
            dir.create(paste0(path, "/", stn_name, "/Tn"), recursive = T)
        }

        # Save graph of the mean
        gg +
            ggplot2::ggsave(paste0(path, "/", stn_name, "/Tn/Mean.Pearsons.Residuals.monthly.png"),
                            dpi = 600, device = 'png')

        # Save graph of the sd
        pp +
            ggplot2::ggsave(paste0(path, "/", stn_name, "/Tn/Std.Pearsons.Residuals.monthly.png"),
                            dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# Definicion de funcion para graficos mensuales de Ocurrencia de precipitación
PearsonsResidualsMonthlyPrcpOcc <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    residals.pearson <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.pearson = resid(model, type = 'pearson')) %>%
        dplyr::mutate(., month = lubridate::month(date, label = T)) %>%
        dplyr::group_by(., month) %>%
        dplyr::summarise(., mean = mean(residuos.pearson, na.rm = T),
                         std = sd(residuos.pearson, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = month, y = mean,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::labs(x = "Month", y = "Mean Pearson's residuals", title = "Precipitacion Occurrence monthly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::theme_bw()

    # Standard deviation
    pp <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = month, y = std,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = residals.pearson) +
        ggplot2::labs(x = "Month", y = "Std.Dev. Pearson's residuals", title = "Precipitation Occurrence monthly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(residals.pearson$std))-1, ceiling(mean(residals.pearson$std))+1), breaks = seq(from = floor(mean(residals.pearson$std))-1, to = ceiling(mean(residals.pearson$std))+1, 0.5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/PrcpOcc"))) {
            dir.create(paste0(path, "/", stn_name, "/PrcpOcc"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/PrcpOcc/Mean.Pearsons.Residuals.monthly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/PrcpOcc/Std.Pearsons.Residuals.monthly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# Definicion de funcion para graficos mensuales de Montos de precipitación
PearsonsResidualsMonthlyPrcpAmt <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    residals.pearson <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.pearson = resid(model, type = 'pearson')) %>%
        dplyr::mutate(., month = lubridate::month(date, label = T)) %>%
        dplyr::group_by(., month) %>%
        dplyr::summarise(., mean = mean(residuos.pearson, na.rm = T),
                         std = sd(residuos.pearson, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = month, y = mean,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::labs(x = "Month", y = "Mean Pearson's residuals", title = "Precipitacion Amounts monthly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::theme_bw()

    # Standard deviation
    pp <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = month, y = std,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = residals.pearson) +
        ggplot2::labs(x = "Month", y = "Std.Dev. Pearson's residuals", title = "Precipitation Amnounts monthly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(residals.pearson$std))-1, ceiling(mean(residals.pearson$std))+1), breaks = seq(from = floor(mean(residals.pearson$std))-1, to = ceiling(mean(residals.pearson$std))+1, 0.5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/PrcpAmt"))) {
            dir.create(paste0(path, "/", stn_name, "/PrcpAmt"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/PrcpAmt/Mean.Pearsons.Residuals.monthly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/PrcpAmt/Std.Pearsons.Residuals.monthly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# Definicion de funcion para graficos mensuales de Montos de precipitación de Residuos de Anscombe
AnscombeResidualsMonthlyPrcpAmt <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    anscombe.residuals <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.anscombe = wle::residualsAnscombe(y = data[,'prcp'], mu = model$fitted.values, family = stats::Gamma(link = log))) %>%
        dplyr::mutate(., month = lubridate::month(date, label = T)) %>%
        dplyr::group_by(., month) %>%
        dplyr::summarise(., mean = mean(residuos.anscombe, na.rm = T),
                         std = sd(residuos.anscombe, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = anscombe.residuals, ggplot2::aes(x = month, y = mean,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = median(anscombe.residuals$mean)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = anscombe.residuals, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = anscombe.residuals, linetype = 'dashed') +
        ggplot2::labs(x = "Month", y = "Mean Anscombe's residuals", title = "Precipitacion Amounts monthly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::theme_bw()

    # Standard deviation
    pp <- ggplot2::ggplot(data = anscombe.residuals, ggplot2::aes(x = month, y = std,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = anscombe.residuals) +
        ggplot2::labs(x = "Month", y = "Std.Dev. Anscombe's residuals", title = "Precipitation Amnounts monthly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(anscombe.residuals$std))-1, ceiling(mean(anscombe.residuals$std))+1), breaks = seq(from = floor(mean(anscombe.residuals$std))-1, to = ceiling(mean(anscombe.residuals$std))+1, 0.5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/PrcpAmt"))) {
            dir.create(paste0(path, "/", stn_name, "/PrcpAmt"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/PrcpAmt/Mean.Anscombe.Residuals.monthly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/PrcpAmt/Std.Anscombe.Residuals.monthly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# ------------------------------------------------------------------------------


# ----------------------------------------------------------------------------- #
# ---- Funciones anuales ----
# ----------------------------------------------------------------------------- #

# Definicion de funcion para graficos anuales de Temperatura maxima
PearsonsResidualsYearlyMaxTemp <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    residals.pearson <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.pearson = resid(model, type = 'pearson')) %>%
        dplyr::mutate(., year = lubridate::year(date)) %>%
        dplyr::group_by(., year) %>%
        dplyr::summarise(., mean = mean(residuos.pearson, na.rm = T),
                         std = sd(residuos.pearson, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = year, y = mean,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::labs(x = "Years", y = "Mean Pearson's residuals", title = "Maximum Temperature monthly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::scale_x_continuous(limits = c(min(residals.pearson$year), max(residals.pearson$year)),
                                    breaks = seq(from = min(residals.pearson$year), to = max(residals.pearson$year), by = 5)) +
        ggplot2::theme_bw()

    # Standard deviation
    # Desvios
    pp <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = year, y = std,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = residals.pearson) +
        ggplot2::labs(x = "Years", y = "Std.Dev. Pearson's residuals", title = "Maximum temperature yearly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(residals.pearson$std))-1, ceiling(mean(residals.pearson$std))+1), breaks = seq(from = floor(mean(residals.pearson$std))-1, to = ceiling(mean(residals.pearson$std))+1, 0.5)) +
        ggplot2::scale_x_continuous(limits = c(min(residals.pearson$year), max(residals.pearson$year)),
                                    breaks = seq(from = min(residals.pearson$year), to = max(residals.pearson$year), by = 5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/Tx"))) {
            dir.create(paste0(path, "/", stn_name, "/Tx"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/Tx/Mean.Pearsons.Residuals.yearly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/Tx/Std.Pearsons.Residuals.yearly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# Definicion de funcion para graficos mensuales de Temperatura minima
PearsonsResidualsYearlyMinTemp <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    residals.pearson <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.pearson = resid(model, type = 'pearson')) %>%
        dplyr::mutate(., year = lubridate::year(date)) %>%
        dplyr::group_by(., year) %>%
        dplyr::summarise(., mean = mean(residuos.pearson, na.rm = T),
                         std = sd(residuos.pearson, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = year, y = mean,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::labs(x = "Years", y = "Mean Pearson's residuals", title = "Minimum Temperature monthly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::scale_x_continuous(limits = c(min(residals.pearson$year), max(residals.pearson$year)),
                                    breaks = seq(from = min(residals.pearson$year), to = max(residals.pearson$year), by = 5)) +
        ggplot2::theme_bw()

    # Standard deviation
    pp <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = year, y = std,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = residals.pearson) +
        ggplot2::labs(x = "Years", y = "Std.Dev. Pearson's residuals", title = "Minimum temperature yearly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(residals.pearson$std))-1, ceiling(mean(residals.pearson$std))+1), breaks = seq(from = floor(mean(residals.pearson$std))-1, to = ceiling(mean(residals.pearson$std))+1, 0.5)) +
        ggplot2::scale_x_continuous(limits = c(min(residals.pearson$year), max(residals.pearson$year)),
                                    breaks = seq(from = min(residals.pearson$year), to = max(residals.pearson$year), by = 5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/Tn"))) {
            dir.create(paste0(path, "/", stn_name, "/Tn"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/Tn/Mean.Pearsons.Residuals.yearly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/Tn/Std.Pearsons.Residuals.yearly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# Definicion de funcion para graficos mensuales de Ocurrencia de lluvia
PearsonsResidualsYearlyPrcpOcc <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    residals.pearson <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.pearson = resid(model, type = 'pearson')) %>%
        dplyr::mutate(., year = lubridate::year(date)) %>%
        dplyr::group_by(., year) %>%
        dplyr::summarise(., mean = mean(residuos.pearson, na.rm = T),
                         std = sd(residuos.pearson, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = year, y = mean)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::labs(x = "Years", y = "Mean Pearson's residuals", title = "Precipitation Occurence yearly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::scale_x_continuous(limits = c(min(residals.pearson$year), max(residals.pearson$year)),
                                    breaks = seq(from = min(residals.pearson$year), to = max(residals.pearson$year), by = 5)) +
        ggplot2::theme_bw()

    # Standard deviation
    pp <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = year, y = std)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = residals.pearson) +
        ggplot2::labs(x = "Years", y = "Std.Dev. Pearson's residuals", title = "Precipitation Occurence yearly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(residals.pearson$std))-1, ceiling(mean(residals.pearson$std))+1), breaks = seq(from = floor(mean(residals.pearson$std))-1, to = ceiling(mean(residals.pearson$std))+1, 0.5)) +
        ggplot2::scale_x_continuous(limits = c(min(residals.pearson$year), max(residals.pearson$year)),
                                    breaks = seq(from = min(residals.pearson$year), to = max(residals.pearson$year), by = 5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/PrcpOcc"))) {
            dir.create(paste0(path, "/", stn_name, "/PrcpOcc"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/PrcpOcc/Mean.Pearsons.Residuals.yearly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/PrcpOcc/Std.Pearsons.Residuals.yearly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# Definicion de funcion para graficos mensuales de Montos de lluvia
PearsonsResidualsYearlyPrcpAmt <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    residals.pearson <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.pearson = resid(model, type = 'pearson')) %>%
        dplyr::mutate(., year = lubridate::year(date)) %>%
        dplyr::group_by(., year) %>%
        dplyr::summarise(., mean = mean(residuos.pearson, na.rm = T),
                         std = sd(residuos.pearson, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = year, y = mean)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = residals.pearson, linetype = 'dashed') +
        ggplot2::labs(x = "Years", y = "Mean Pearson's residuals", title = "Precipitation Amounts yearly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::scale_x_continuous(limits = c(min(residals.pearson$year), max(residals.pearson$year)),
                                    breaks = seq(from = min(residals.pearson$year), to = max(residals.pearson$year), by = 5)) +
        ggplot2::theme_bw()

    # Standard deviation
    pp <- ggplot2::ggplot(data = residals.pearson, ggplot2::aes(x = year, y = std)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = residals.pearson) +
        ggplot2::labs(x = "Years", y = "Std.Dev. Pearson's residuals", title = "Precipitation Amounts yearly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(floor(mean(residals.pearson$std))-1, ceiling(mean(residals.pearson$std))+1), breaks = seq(from = floor(mean(residals.pearson$std))-1, to = ceiling(mean(residals.pearson$std))+1, 0.5)) +
        ggplot2::scale_x_continuous(limits = c(min(residals.pearson$year), max(residals.pearson$year)),
                                    breaks = seq(from = min(residals.pearson$year), to = max(residals.pearson$year), by = 5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/PrcpAmt"))) {
            dir.create(paste0(path, "/", stn_name, "/PrcpAmt"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/PrcpAmt/Mean.Pearsons.Residuals.yearly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/PrcpAmt/Std.Pearsons.Residuals.yearly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# Definicion de funcion para graficos mensuales de Montos de precipitación de Residuos de Anscombe
AnscombeResidualsYearlyPrcpAmt <- function(stn_name = NULL, model = NULL, path = NULL) {
    # The model must be of class lm, lmrob, glm, glmrob
    base::stopifnot(class(model) %in% c('lm', 'lmrob', 'glm', 'glmrob'))
    # Estimation of Pearson's residuals
    anscombe.residuals <- data.frame(date = model$dates) %>%
        dplyr::mutate(., residuos.anscombe = wle::residualsAnscombe(y = data[,'prcp'], mu = model$fitted.values, family = stats::Gamma(link = log))) %>%
        dplyr::mutate(., year = lubridate::year(date)) %>%
        dplyr::group_by(., year) %>%
        dplyr::summarise(., mean = mean(residuos.anscombe, na.rm = T),
                         std = sd(residuos.anscombe, na.rm = T),
                         n = n()) %>%
        dplyr::mutate(., error = qnorm(0.975)*std/sqrt(n),
                      upper.bound = mean + error,
                      lower.bound = mean - error)
    # Mean
    gg <- ggplot2::ggplot(data = anscombe.residuals, ggplot2::aes(x = year, y = mean,  group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = median(anscombe.residuals$mean)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(lower.bound)), data = anscombe.residuals, linetype = 'dashed') +
        ggplot2::geom_hline(ggplot2::aes(yintercept = median(upper.bound)), data = anscombe.residuals, linetype = 'dashed') +
        ggplot2::labs(x = "Month", y = "Mean Anscombe's residuals", title = "Precipitacion Amounts yearly mean of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        ggplot2::scale_x_continuous(limits = c(min(anscombe.residuals$year), max(anscombe.residuals$year)),
                                    breaks = seq(from = min(anscombe.residuals$year), to = max(anscombe.residuals$year), by = 5)) +
        ggplot2::theme_bw()

    # Standard deviation
    pp <- ggplot2::ggplot(data = anscombe.residuals, ggplot2::aes(x = year, y = std)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean(std)), data = anscombe.residuals) +
        ggplot2::labs(x = "Month", y = "Std.Dev. Anscombe's residuals", title = "Precipitation Amnounts yearly standard deviation of the residuals", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.25)) +
        ggplot2::scale_x_continuous(limits = c(min(anscombe.residuals$year), max(anscombe.residuals$year)),
                                    breaks = seq(from = min(anscombe.residuals$year), to = max(anscombe.residuals$year), by = 5)) +
        ggplot2::theme_bw()


    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/PrcpAmt"))) {
            dir.create(paste0(path, "/", stn_name, "/PrcpAmt"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/PrcpAmt/Mean.Anscombe.Residuals.yearly.png"),
                        dpi = 600, device = 'png')

        # Save graph of the sd
        ggplot2::ggsave(plot = pp,
                        filename = paste0(path, "/", stn_name, "/PrcpAmt/Std.Anscombe.Residuals.yearly.png"),
                        dpi = 600, device = 'png')
    } else {
        list(gg, pp)
    }

}

# ------------------------------------------------------------------------------


# ---------------------------------------------------------------------------- #
# ---- Statistical diagnostics ----
# ---------------------------------------------------------------------------- #

# Definición de función para realizar Q-Q Plots de resiclduos estudentizados
qqplotStudentResiduesNormalMaxTemp <- function(stn_name = NULL, model = NULL, conf = 0.95, path = NULL) {
    # Si se introduce un objeto lm, glm se extraern los residuos studentizados
    if (class(model) == c('lm')) {
        # Extraer los residuos studentizados del modelo lineal
        student.residuals <- stats::rstudent(model)
    } else if (class(model) == c('lmrob')) {
        student.residuals <- MASS::studres(model)
    } else {
        warning(paste0("El formato de entrada debe ser un objeto lm o glm o un vector numérico"))
    }

    # Eliminar valores NAs si los hubiese
    student.residuals <- na.omit(student.residuals)
    # Determinar orden de los residuos
    ord <- order(student.residuals)
    # Determinada longitud del vector de residuos
    n <- length(student.residuals)
    # Crear plotting positions
    P <- ppoints(length(student.residuals))
    # Crear data frame con los residuos ordenados y los cuantiles teóricos de una normal
    residuos.data.frame <- data.frame(ord.res = student.residuals[ord], cuantiles.normales = qnorm(P))

    # Crear intervalos de confianza
    # Cuantiles de los residuos
    Q.x <- quantile(residuos.data.frame$ord.res, c(0.25, 0.75))
    # Cuantiles teóricos
    Q.z <- qnorm(c(0.25, 0.75))
    # Rrelación entre ambos
    beta <- diff(Q.x)/diff(Q.z)
    # Coeficiente de regresión
    coef <- c(Q.x[1] - beta * Q.z[1], beta)
    # Cuantil del intervalo de confianza elegido
    zz <- qnorm(1 - (1 - conf)/2)
    # Desvio del intervalo de confianza
    SE <- (coef[2]/dnorm(residuos.data.frame$cuantiles.normales)) * sqrt(P * (1 - P)/n)
    # Linea de gresión para crear la banda de confianza
    fit.value <- coef[1] + coef[2] * residuos.data.frame$cuantiles.normales
    # Umbral superior del intervalo de confianza
    residuos.data.frame$upper <- fit.value + zz * SE
    # Umbral inferior del intervalo de confianza
    residuos.data.frame$lower <- fit.value - zz * SE
    # Crear categrotía de datos normales o no si están dentro del intervalo de confianza
    residuos.data.frame <- residuos.data.frame %>%
        dplyr::mutate(., normal = if_else(ord.res >= lower & ord.res <= upper, 1, 0))

    # Creación del gráfico
    gg <- ggplot2::ggplot(residuos.data.frame, ggplot2::aes(x=cuantiles.normales, y=ord.res)) +
        ggplot2::geom_point(data = residuos.data.frame, ggplot2::aes(color = as.factor(normal))) +
        ggplot2::theme(legend.position = 'bottom') +
        ggplot2::coord_equal(ratio=1) +
        ggplot2::geom_qq_line(data = residuos.data.frame, ggplot2::aes(sample = ord.res)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha=0.2) +
        ggplot2::labs(title = "Maximum Temperature Q-Q Plot", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_x_continuous(limits = c(min(residuos.data.frame - 0.1), max(residuos.data.frame + 0.1)),
                                    breaks = seq(floor(min(residuos.data.frame - 0.1)), ceiling(max(residuos.data.frame + 0.1)),1),
                                    name = "Theoretical quantiles") +
        ggplot2::scale_y_continuous(limits = c(min(residuos.data.frame - 0.1), max(residuos.data.frame + 0.1)),
                                    breaks = seq(floor(min(residuos.data.frame - 0.1)), ceiling(max(residuos.data.frame + 0.1)),1),
                                    name = "Studentized residuals") +
        ggplot2::scale_color_brewer(palette = 'Set1', labels = c("Not normal", "Normal"),
                                    guide = guide_legend(reverse=TRUE),
                                    name = element_blank()) +
        ggplot2::theme_bw()

    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/Tx"))) {
            dir.create(paste0(path, "/", stn_name, "/Tx"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/Tx/Studentized.residuals.normal.png"),
                        dpi = 600, device = 'png')


    } else {
        return(gg)
    }

}

# Definición de función para realizar Q-Q Plots de resiclduos estudentizados
qqplotStudentResiduesNormalMinTemp <- function(stn_name = NULL, model = NULL, conf = 0.95, path = NULL) {
    # Si se introduce un objeto lm, glm se extraern los residuos studentizados
    if (class(model) == c('lm')) {
        # Extraer los residuos studentizados del modelo lineal
        student.residuals <- stats::rstudent(model)
    } else if (class(model) == c('lmrob')) {
        student.residuals <- MASS::studres(model)
    } else {
        warning(paste0("El formato de entrada debe ser un objeto lm o glm o un vector numérico"))
    }

    # Eliminar valores NAs si los hubiese
    student.residuals <- na.omit(student.residuals)
    # Determinar orden de los residuos
    ord <- order(student.residuals)
    # Determinada longitud del vector de residuos
    n <- length(student.residuals)
    # Crear plotting positions
    P <- ppoints(length(student.residuals))
    # Crear data frame con los residuos ordenados y los cuantiles teóricos de una normal
    residuos.data.frame <- data.frame(ord.res = student.residuals[ord], cuantiles.normales = qnorm(P))

    # Crear intervalos de confianza
    # Cuantiles de los residuos
    Q.x <- quantile(residuos.data.frame$ord.res, c(0.25, 0.75))
    # Cuantiles teóricos
    Q.z <- qnorm(c(0.25, 0.75))
    # Rrelación entre ambos
    beta <- diff(Q.x)/diff(Q.z)
    # Coeficiente de regresión
    coef <- c(Q.x[1] - beta * Q.z[1], beta)
    # Cuantil del intervalo de confianza elegido
    zz <- qnorm(1 - (1 - conf)/2)
    # Desvio del intervalo de confianza
    SE <- (coef[2]/dnorm(residuos.data.frame$cuantiles.normales)) * sqrt(P * (1 - P)/n)
    # Linea de gresión para crear la banda de confianza
    fit.value <- coef[1] + coef[2] * residuos.data.frame$cuantiles.normales
    # Umbral superior del intervalo de confianza
    residuos.data.frame$upper <- fit.value + zz * SE
    # Umbral inferior del intervalo de confianza
    residuos.data.frame$lower <- fit.value - zz * SE
    # Crear categrotía de datos normales o no si están dentro del intervalo de confianza
    residuos.data.frame <- residuos.data.frame %>%
        dplyr::mutate(., normal = if_else(ord.res >= lower & ord.res <= upper, 1, 0))

    # Creación del gráfico
    gg <- ggplot2::ggplot(residuos.data.frame, ggplot2::aes(x=cuantiles.normales, y=ord.res)) +
        ggplot2::geom_point(data = residuos.data.frame, ggplot2::aes(color = as.factor(normal))) +
        ggplot2::theme(legend.position = 'bottom') +
        ggplot2::coord_equal(ratio=1) +
        ggplot2::geom_qq_line(data = residuos.data.frame, ggplot2::aes(sample = ord.res)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha=0.2) +
        ggplot2::labs(title = "Minimum Temperature Q-Q Plot", subtitle = paste0("Station: ", stn_name)) +
        ggplot2::scale_x_continuous(limits = c(min(residuos.data.frame - 0.1), max(residuos.data.frame + 0.1)),
                                    breaks = seq(floor(min(residuos.data.frame - 0.1)), ceiling(max(residuos.data.frame + 0.1)),1),
                                    name = "Theoretical quantiles") +
        ggplot2::scale_y_continuous(limits = c(min(residuos.data.frame - 0.1), max(residuos.data.frame + 0.1)),
                                    breaks = seq(floor(min(residuos.data.frame - 0.1)), ceiling(max(residuos.data.frame + 0.1)),1),
                                    name = "Studentized residuals") +
        ggplot2::scale_color_brewer(palette = 'Set1', labels = c("Not normal", "Normal"),
                                    guide = guide_legend(reverse=TRUE),
                                    name = element_blank()) +
        ggplot2::theme_bw()

    if (!is.null(path)) {
        if (!dir.exists(paste0(path, "/", stn_name, "/Tn"))) {
            dir.create(paste0(path, "/", stn_name, "/Tn"), recursive = T)
        }

        # Save graph of the mean
        ggplot2::ggsave(plot = gg,
                        filename = paste0(path, "/", stn_name, "/Tn/Studentized.residuals.normal.png"),
                        dpi = 600, device = 'png')


    } else {
        return(gg)
    }

}

# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# ---- General functions ----
# ---------------------------------------------------------------------------- #

#' @title ApplyMaxTempDiagnostics
#' @description Function that runs Pearson's residuals graphical diagnostics on Max Temperature model.
#' @export
ApplyMaxTempDiagnostics <- function(stn_name = NULL, model = NULL, path = NULL) {
    graphs <- purrr::map(
        .x = c("PearsonsResidualsMonthly", "PearsonsResidualsYearly"),
        .f = function(test.name) {
            func.name <- paste0(test.name, "MaxTemp")
            tryCatch({
                diagnostics <- do.call(what = func.name, args = list(stn_name = stn_name, model = model, path = path))
            }, error = function(e) {
                cat(e$message, "\n")
                return (NULL)
            })
        }
    )
    if (is.null(path)) {
        return(graphs)
    }
}

#' @title ApplyMinTempDiagnostics
#' @description Function that runs Pearson's residuals graphical diagnostics on Min Temperature model.
#' @export
ApplyMinTempDiagnostics <- function(stn_name = NULL, model = NULL, path = NULL) {
    graphs <- purrr::map(
        .x = c("PearsonsResidualsMonthly", "PearsonsResidualsYearly"),
        .f = function(test.name) {
            func.name <- paste0(test.name, "MinTemp")
            tryCatch({
                diagnostics <- do.call(what = func.name, args = list(stn_name = stn_name, model = model, path = path))
            }, error = function(e) {
                cat(e$message, "\n")
                return (NULL)
            })
        }
    )
    if (is.null(path)) {
        return(graphs)
    }
}

#' @title ApplyPrcpOccDiagnostics
#' @description Function that runs Pearson's residuals graphical diagnostics on Precipitation occurrence model.
#' @export
ApplyPrcpOccDiagnostics <- function(stn_name = NULL, model = NULL, path = NULL) {
    graphs <- purrr::map(
        .x = c("PearsonsResidualsMonthly", "PearsonsResidualsYearly"),
        .f = function(test.name) {
            func.name <- paste0(test.name, "PrcpOcc")
            tryCatch({
                diagnostics <- do.call(what = func.name, args = list(stn_name = stn_name, model = model, path = path))
            }, error = function(e) {
                cat(e$message, "\n")
                return (NULL)
            })
        }
    )
    if (is.null(path)) {
        return(graphs)
    }
}

#' @title ApplyPrcpAmtDiagnostics
#' @description Function that runs Pearson's and Anscombe's residuals graphical diagnostics on Precipitation occurrence model.
#' @export
ApplyPrcpAmtDiagnostics <- function(stn_name = NULL, model = NULL, path = NULL) {
    graphs <- purrr::map(
        .x = c("PearsonsResidualsMonthly", "PearsonsResidualsYearly", "AnscombeResidualsMonthly", "AnscombeResidualsYearly"),
        .f = function(test.name) {
            func.name <- paste0(test.name, "PrcpAmt")
            tryCatch({
                diagnostics <- do.call(what = func.name, args = list(stn_name = stn_name, model = model, path = path))
            }, error = function(e) {
                cat(e$message, "\n")
                return (NULL)
            })
        }
    )
    if (is.null(path)) {
        return(graphs)
    }
}

#' @title SummaryTablePrecipitation
#' @description Summary precipitation model.
#' @export
SummaryTablePrecipitation <- function(prcp_occ_fit = NULL, prcp_amt_fit = NULL) {
    summary.table.precipitation <- do.call(rbind, lapply(list(prcp_occ_fit,
                                                              prcp_amt_fit), broom::glance)) %>%
        dplyr::mutate(., model = c("occ", "amt")) %>%
        dplyr::select(., model, null.deviance, df.null, logLik, AIC, BIC, deviance, df.residual)
    return(summary.table.precipitation)
}

#' @title SummaryTableTemperature
#' @description Summary temperature model.
#' @export
SummaryTableTemperature <- function(tx_fit = NULL, tn_fit = NULL) {
    summary.table.temperature <- do.call(rbind, lapply(list(tx_fit,
                                                            tn_fit), broom::glance)) %>%
        dplyr::mutate(., model = c("tx", "tn")) %>%
        dplyr::select(., model, r.squared, adj.r.squared, sigma, statistic, p.value, df, logLik, AIC, BIC, deviance)
    return(summary.table.temperature)
}

# ---------------------------------------------------------------------------- #
