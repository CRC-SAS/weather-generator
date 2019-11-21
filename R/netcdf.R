
# Definición de función para la creación de archivos NetCDF
CrearNetCDF <- function(nc_file, num_realizations, sim_dates, simulation_points) {
  # Obtener valores para dimensiones X a Y
  coordenadas <- sf::st_coordinates(simulation_points)
  x_vals      <- sort(unique(coordenadas[, 1]))
  y_vals      <- sort(unique(coordenadas[, 2]))
  t_vals      <- julian(sim_dates)
  r_vals      <- seq(from = 1, to = num_realizations, by = 1)

  # Creacion de dimensiones
  dimX <- ncdf4::ncdim_def(name = "projection_x_coordinate",
                           units = "meters",
                           vals = x_vals,
                           unlim = FALSE,
                           create_dimvar = TRUE,
                           longname = "x-coordinate")
  dimY <- ncdf4::ncdim_def(name = "projection_y_coordinate",
                           units = "meters",
                           vals = y_vals,
                           unlim = FALSE,
                           create_dimvar = TRUE,
                           longname = "y-coordinate")
  dimT <- ncdf4::ncdim_def(name = "time",
                           units = 'days since 1970-01-01',
                           vals = t_vals,
                           unlim = FALSE,
                           calendar = "standard",
                           create_dimvar = TRUE,
                           longname = 'Days since 1970-01-01')
  dimR <- ncdf4::ncdim_def(name = "realization",
                           units = "(none)",
                           vals = r_vals,
                           unlim = FALSE,
                           create_dimvar = TRUE,
                           longname = 'Realization number')

  # Creacion de variables de indices
  def_vars <- list(
    list(name = "tmax", units = 'degrees Celsius', description = 'Maximum daily air temperature'),
    list(name = "tmin", units = 'degrees Celsius', description = 'Minimum daily air temperature'),
    list(name = "prcp", units = 'millimeters', description = 'Daily precipitation')
  )

  variables <- list()
  for (i in seq(from = 1, to = length(def_vars))) {
    variable_data  <- def_vars[[i]]
    variable       <- ncdf4::ncvar_def(name = variable_data$name,
                                       units = variable_data$units,
                                       dim = list(dimX, dimY, dimT, dimR),
                                       longname = variable_data$description,
                                       prec = "float",
                                       compression = 9,
                                       verbose = FALSE)
    variables[[i]] <- variable
  }

  # Creacion de archivo
  nc <- ncdf4::nc_create(filename = nc_file,
                         vars = variables,
                         force_v4 = TRUE,
                         verbose = FALSE)

  # Agregado de atributos
  ncdf4::ncatt_put(nc, varid = 0,
                   attname = "CRS",
                   attval = sf::st_crs(simulation_points)$proj4string,
                   verbose = FALSE)
  ncdf4::ncatt_put(nc, varid = 0,
                   attname = "Conventions",
                   attval = "CF-1.7",
                   verbose = FALSE)

  # Cierre de archivo
  ncdf4::nc_close(nc)

  return (nc_file)
}

# Definición de función para la inserción de datos de una realización en un NetCDF
GuardarRealizacionNetCDF <- function(nc_file, numero_realizacion, sim_dates, raster_tmax, raster_tmin, raster_prcp) {
  # Abrir NetCDF para escritura
  nc <- ncdf4::nc_open(filename = nc_file, write = TRUE)

  # Agregado de variables
  variables <- list(
    "tmax" = raster_tmax,
    "tmin" = raster_tmin,
    "prcp" = raster_prcp
  )
  for (variable in names(variables)) {
    raster_variable <- variables[[variable]]
    if (length(raster_variable) > 0) {
      stack_variable  <- raster::flip(raster::stack(raster_variable), direction = "y")
      array_variable  <- base::aperm(as.array(stack_variable), perm = c(2, 1, 3))
      ncdf4::ncvar_put(nc = nc,
                       varid = variable,
                       vals  = array_variable,
                       start = c(1, 1, 1, numero_realizacion),
                       count = c(-1, -1, -1, 1))
    }
  }

  # Cierre de archivo
  ncdf4::nc_close(nc)
}
