
# Funcion para pasar de objeto SF a raster
sf2raster <- function(objeto_sf, variable) {
    # Genero data frame con las coordenadas del objeto SF
    raster_objeto <- data.frame(sf::st_coordinates(objeto_sf)) %>%
        # Agrego una variable Z que corresponda a la variable de interes
        dplyr::mutate(z = dplyr::pull(objeto_sf, !!variable)) %>%
        # Transformo data frame a raster
        raster::rasterFromXYZ(.)
    # Defino CRS para el raster tomando el mismo del objeto SF
    raster::crs(raster_objeto) <- sf::st_crs(objeto_sf)$proj4string
    # Devuelvo el raster
    return (raster_objeto)
}

# Funcion para ensamblar rasters de residuos considerando dias secos y humedos
ensamblar_raster_residuos <- function(raster_residuos_dias_humedos, raster_residuos_dias_secos, raster_ocurrencia) {
    # Genera el raster para dias humedos (eliminando los dias secos, cuyo valor de ocurrencia es 0)
    raster_humedos <- raster::mask(x = raster_residuos_dias_humedos, mask = raster_ocurrencia, maskvalue = 0, updatevalue = 0)
    # Genera el raster para dias secos (eliminando los dias humedos, cuyo valor de ocurrencia es 1)
    raster_secos   <- raster::mask(x = raster_residuos_dias_secos, mask = raster_ocurrencia, maskvalue = 1, updatevalue = 0)
    # Suma ambos rasters
    return (raster_humedos + raster_secos)
}
