

#' @title Recover partial runs
#' @description This feature enables recovery from unexpected situations.
#' @param run_output_folder the output_folder which was set to the execution to be recovered
#' @param run_output_filename the output_filename which was set to the execution to be recovered
#' @export
recover_data_of_previous_runs <- function(run_output_folder = getwd(), run_output_filename = "sim_results.csv", avbl_cores = 1) {

    ## Español: Carpeta de destino y nombre del archivo
    ## English: Process run_output_folder and run_output_filename
    run_output_folder <- sub('/$', '', run_output_folder)
    run_output_filename <- sub('\\.([^.]*)$', '', run_output_filename)

    ## Español: Definir archivos a recuperar
    ## English: Define files to be recovered
    files_pattern <- glue::glue("{run_output_filename}_realization_[0-9]+\\.rds")
    files_to_recover <- list.files(run_output_folder, pattern = files_pattern, full.names = T)

    ## Español: Definir nombre del archivo con los datos recuperados
    ## English: Define the name of the file with the recovered data
    recovery_output_filename <- glue::glue("{run_output_folder}/{run_output_filename}.csv")

    ## Español: Recuperar archivos
    ## English: Recover files
    for (rds_path in sort(files_to_recover)) {
        message(glue::glue("Recovering \"{rds_path}\" ({which(rds_path == files_to_recover)} of {length(files_to_recover)})"))

        r <- as.integer(sub("^.*_realization_([0-9]+)\\.rds$", "\\1", rds_path))

        tibble_with_data <- base::readRDS(rds_path) %>%
            dplyr::select(nsim, tidyselect::any_of(c("station_id", "point_id")), longitude, latitude,
                          date, tmax, tmin, prcp_occ, prcp_amt, type_day)

        gamwgen:::GuardarRealizacionEnCSV(filename = recovery_output_filename,
                                          numero_realizacion = r,
                                          tibble_with_data = tibble_with_data,
                                          avbl_cores = avbl_cores)
    }

    return (glue::glue("Output file (file with the recovered data): {recovery_output_filename}"))
}

