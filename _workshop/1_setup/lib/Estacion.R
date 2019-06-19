require(R6)
require(dplyr)
require(dbplyr)
require(magrittr)

EstacionFacade <- R6Class("EstacionFacade",
	inherit = Facade,
	public = list(
		initialize = function(con) {
			private$con <- con
		},
  
		buscar = function(omm_id = NULL, institucion_id = NULL, institucion_siglas = NULL, pais_id = NULL,
		                  tipo = NULL, solo.activas = FALSE) {
		  instituciones <- dplyr::tbl(src = private$con, dbplyr::in_schema('public', 'institucion')) %>%
		    dplyr::rename(institucion_id = id) %>%
		    dplyr::select(institucion_id, pais_id, siglas) %>%
		    dplyr::collect() 
		  filtros       <- list(omm_id = as.integer(omm_id), institucion_id = as.integer(institucion_id), tipo = tipo)
		  estaciones    <- private$filtrar(dplyr::tbl(src = private$con, dbplyr::in_schema('public', 'estacion')), filtros) %>%
		    dplyr::collect() %>%
		    dplyr::inner_join(instituciones, by = c("institucion_id"))
		  if (solo.activas) {
		    estaciones %<>% dplyr::filter(is.na(fecha_fin))
		  }
		  if (! is.null(pais_id)) {
        pais.id    <- pais_id
        estaciones %<>% dplyr::filter(pais_id == pais.id)
      }
      if (! is.null(institucion_siglas)) {
        estaciones %<>% dplyr::filter(siglas == institucion_siglas)
      }
		  
		  return (estaciones)
		},
		
		buscarVecinasPrecalculadas = function(omm_id, max.distancia = NULL, max.dif.elevacion = NULL, max.vecinas = NULL) {
		  # a. Buscar estaciones vecinas
		  filtros <- list(omm_id = as.integer(omm_id))
		  vecinas <- private$filtrar(dplyr::tbl(src = private$con, dbplyr::in_schema('public', 'estacion_vecino')), filtros) %>%
		    dplyr::select(omm_vecino_id, distancia, diferencia_elevacion) %>%
		    dplyr::rename(omm_id = omm_vecino_id) %>%
		    dplyr::collect()
		  
		  if (! is.null(max.distancia)) {
		    vecinas <- vecinas %>%
		      dplyr::filter(distancia <= max.distancia)
		  }
		  
		  if (! is.null(max.dif.elevacion)) {
		    vecinas <- vecinas %>%
		      dplyr::filter(diferencia_elevacion <= max.dif.elevacion)
		  }
		  
		  if (! is.null(max.vecinas) && (nrow(vecinas) > max.vecinas)) {
		    vecinas <- vecinas %>%
		      dplyr::arrange(distancia, diferencia_elevacion)
		    vecinas <- vecinas[1:max.vecinas,]
		  }
		  
		  # b. Hacer join con estacion y devolver todos los datos
		  results <- self$buscar() %>%
		    dplyr::inner_join(vecinas, by = c("omm_id")) %>%
		    dplyr::arrange(distancia, diferencia_elevacion)
		  return (results)
		},
		
		buscarVecinasExhaustivo = function(omm_id, max.distancia = NULL, max.dif.elevacion = NULL, max.vecinas = NULL) {
		  query   <- "select e1.omm_id as omm_original_id, e2.omm_id,
                    st_distance(st_makepoint(e2.lon_dec, e2.lat_dec, 4326), st_makepoint(e1.lon_dec, e1.lat_dec, 4326), true)/1000 distancia,
                    abs(e2.elev - e1.elev) as diferencia_elevacion
                  from estacion e1, estacion e2
                  where e1.omm_id != e2.omm_id
                  order by distancia, diferencia_elevacion"
		  filtros <- list(omm_id = omm_id)
		  vecinas <- DBI::dbGetQuery(private$con, query) %>%       # no hay riesgo de sql injection, porque query no recibe paramteros!!
		    dplyr::filter(omm_original_id == filtros$omm_id) %>%   # si los recibiera, se debería utilizar dbSendQuery en su lugar!!
		    dplyr::select(omm_id, distancia, diferencia_elevacion) # ver: https://db.rstudio.com/best-practices/run-queries-safely/#sql-injection-attack
		  
		  # Aplicar filtros por distancia y diferencia de elevacion
		  if (! is.null(max.distancia)) {
		    vecinas %<>% dplyr::filter(distancia <= max.distancia)
		  }
		  if (! is.null(max.dif.elevacion)) {
		    vecinas %<>% dplyr::filter(diferencia_elevacion <= max.dif.elevacion)
		  }
		  
		  # Applicar otros filtros y devolver resultados
		  results <- self$buscar() %>%
		    dplyr::inner_join(vecinas, by = c("omm_id")) %>%
		    dplyr::arrange(distancia, diferencia_elevacion)
		  if (! is.null(max.vecinas) && (nrow(results) > max.vecinas)) {
		    results <- results[seq(from = 1, to = max.vecinas), ]
		  }
		  return (results)
		},
		
		buscarVecinas = function(omm_id, max.distancia = NULL, max.dif.elevacion = NULL, max.vecinas = NULL, exhaustivo = FALSE) {
		  if (exhaustivo) {
		    return (self$buscarVecinasExhaustivo(omm_id, max.distancia, max.dif.elevacion, max.vecinas))
		  } else {
		    return (self$buscarVecinasPrecalculadas(omm_id, max.distancia, max.dif.elevacion, max.vecinas))
		  }
		}
	)
)
