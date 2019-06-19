require(R6)
require(dplyr)
require(dbplyr)
require(magrittr)
require(tidyr)

RegistroDiarioFacade <- R6Class("RegistroDiarioFacade",
	inherit = Facade,
	public = list(
		initialize = function(con) {
			private$con <- con
		},

		buscar = function(omm_id = NULL, variable_id = NULL, estado.registro = NULL,
		                  fecha.desde = NULL, fecha.hasta = NULL) {
		  # a. Buscar registros por estacion y estado de registro
		  filtros.ancho   <- list(omm_id = omm_id, estado = estado.registro)
		  registros.ancho <- private$filtrar(dplyr::tbl(src = private$con, dbplyr::in_schema('public', 'estacion_registro_diario')), 
		                                     filtros.ancho)
		  
		  # b. Filtrar por fecha desde y hasta de ser necesario
		  if (! is.null(fecha.desde)) {
		    registros.ancho %<>% dplyr::filter(fecha >= fecha.desde)
		  }
		  if (! is.null(fecha.hasta)) {
		    registros.ancho %<>% dplyr::filter(fecha <= fecha.hasta)
		  }
		  registros.ancho %<>% dplyr::collect()
		  
		  # c. Pasar a tabla larga y filtrar por variable si es necesario
		  filtros.largo   <- list(variable_id = variable_id)
		  registros.largo <- tidyr::gather(data = registros.ancho, key = 'variable_id', value = valor, 
		                                   -omm_id, -fecha, -estado, -num_observaciones) %>%
		    private$filtrar(., filtros.largo)

		  # d. Ordenar por omm_id, fecha
		  registros.largo %<>% dplyr::arrange(omm_id, fecha)

		  return (dplyr::collect(registros.largo))
		},

		# buscarResultadosTests = function(red_id = NULL, estacion_id = NULL, variable_id = NULL,
		#                                  fecha.desde = NULL, fecha.hasta = NULL, resultado.test = NULL) {
		#   filtros    <- list(red_id = red_id, estacion_id = as.character(estacion_id), variable_id = variable_id,
		#                      fecha.desde = fecha.desde, fecha.hasta = fecha.hasta, resultado.test = resultado.test)
		#   resultados <- private$filtrar(dplyr::tbl(src = private$con, dbplyr::in_schema('public', 'registro_diario_test')), filtros)
		#   return (resultados)
		# },

		primeraFechaRegistro = function(omm_id = NULL, estado.registro = NULL) {
		  filtros   <- list(omm_id = omm_id, estado = estado.registro)
		  registros <- private$filtrar(dplyr::tbl(src = private$con, dbplyr::in_schema('public', 'estacion_registro_diario')), 
		                               filtros) %>%
		    dplyr::summarise(primera.fecha = min(fecha, na.rm = TRUE)) %>%
		    dplyr::collect()
		  return (registros$primera.fecha)
		},

		ultimaFechaRegistro = function(omm_id = NULL, estado.registro = NULL) {
		  filtros   <- list(omm_id = omm_id, estado = estado.registro)
		  registros <- private$filtrar(dplyr::tbl(src = private$con, dbplyr::in_schema('public', 'estacion_registro_diario')), 
		                               filtros) %>%
	      dplyr::summarise(ultima.fecha = max(fecha, na.rm = TRUE)) %>%
	      dplyr::collect()
		  return (registros$ultima.fecha)
		},

		# registros.nuevos = (omm_id, fecha, <variables>, estado, num_observaciones)
		agregarRegistrosNuevos = function(registros.nuevos) {
		  # i. Se inicia una transaccion y se realizan las operaciones en un bloque tryCatch
		  DBI::dbBegin(conn = private$con)
		  
		  tryCatch({
		    # ii. Insertar datos de estacion_registro_diario
		    DBI::dbWriteTable(conn = private$con, name = "estacion_registro_diario", 
		                      value = registros.nuevos, append = TRUE, row.names = FALSE)
		    
		    # iii. Hacer commit
		    DBI::dbCommit(conn = private$con)
		  }, error = function(e) {
		    # Error en la transaccion. Hacer rollback y luego mostrar el error.
		    DBI::dbRollback(conn = private$con)
		    stop(e)
		  })
		},
		
		# omm_id = ID OMM de la estacion cuyo datos ser van a modificar
		# registros.pendientes = (omm_id, fecha, variable_id, valor, pasa_test.*, dictamen)
		# fecha.control = fecha de ejecucion del script
 		actualizarEstacion = function(omm_id, registros.pendientes, fecha.control) {
 		  # i. Generar registros para estacion_variable_estado
 		  #    estacion_variable_estado(omm_id, fecha, variable, valor, estado)
 		  #    No se guardan los aprobados 
 		  estacion.variable.estado <- registros.pendientes %>%
 		    dplyr::filter(dictamen != 'A') %>%
 		    dplyr::rename(estado = dictamen, variable = variable_id) %>%
 		    dplyr::select(omm_id, fecha, variable, valor, estado)
 		  
 		  # ii. Generar registros para estacion_variable_test
 		  #     estacion_variable_test = (omm_id, fecha, variable, test, valor, fecha_control) 
 		  #     Solo se guardan los tests fallidos
 		  estacion.variable.test <- registros.pendientes %>%
 		    dplyr::select(-dictamen) %>%
 		    dplyr::rename(variable = variable_id) %>%
 		    tidyr::gather(key = pasa_test, value = resultado_test, -omm_id, -fecha, -variable, -valor) %>%
 		    dplyr::filter(! is.na(resultado_test) & ! resultado_test) %>%
 		    tidyr::separate(col = pasa_test, into = c('dummy', 'test'), sep = '\\.') %>%
 		    dplyr::select(omm_id, fecha, variable, test, valor) %>%
 		    dplyr::mutate(fecha_control = fecha.control)
 		  
 		  # iii. Generar registros para estacion_registro_diario
 		  #      Las variables que figuran con estado N deben guardarse como faltantes
 		  #      Se debe generar el dictamen a nivel de fecha.
 		  #      Si hay al menos un dictamen S => el estado es D
 		  #      Sino, si hay al menos un dictamen N => el estado es R
 		  #      Sino, el estado es V
 		  fechas    <- registros.pendientes %>%
 		    dplyr::distinct(fecha) %>%
 		    dplyr::arrange(fecha) %>%
 		    dplyr::pull(fecha)
 		  variables <- dplyr::tbl(src = private$con, dbplyr::in_schema('public', 'variable')) %>%
 		    dplyr::arrange(orden_columna) %>%
 		    dplyr::collect() %>%
 		    dplyr::pull(id)
 		  
 		  # Primero genero los dictamenes para cada fecha
 		  estacion.registro.diario.estado <- purrr::map_dfr(
 		    .x = fechas,
 		    .f = function(una.fecha) {
 		      estados.variable <- registros.pendientes %>%
 		        dplyr::filter(fecha == una.fecha) %>%
 		        dplyr::pull(dictamen)
 		      if (length(estados.variable) > 0) {
 		        if (any('S' == estados.variable, na.rm = TRUE)) {
 		          estado.registro <- 'D'
 		        } else if (any('N' == estados.variable, na.rm = TRUE)) {
 		          estado.registro <- 'R'
 		        } else {
 		          estado.registro <- 'V'
 		        }
 		      } else {
 		        estado.registro <- 'V'
 		      }
 		      
 		      return (data.frame(omm_id = omm_id, fecha = una.fecha, estado = estado.registro))
 		    }
 		  ) %>% as.data.frame()
 		  
 		  # Ahora busco las variables con estado N para reemplazar su valor por NA
 		  estacion.registro.diario.faltantes <- registros.pendientes %>%
 		    dplyr::filter(dictamen == 'N') %>%
 		    dplyr::select(omm_id, fecha, variable_id) %>%
        as.data.frame()
 		  
 		  # iv. Ya se tiene toda la informacion necesaria preparada para la actualizacion
 		  #     Se inicia una transaccion y se realizan las operaciones en un bloque tryCatch
		  DBI::dbBegin(conn = private$con)
		  
		  tryCatch({
		    # v. Insertar datos de estaciones_variable_estado
		    if (nrow(estacion.variable.estado) > 0) {
  		    DBI::dbWriteTable(conn = private$con, name = "estacion_variable_estado", 
  		                      value = estacion.variable.estado, append = TRUE, row.names = FALSE)
		    }
		    
		    # vi. Insertar datos de estaciones_variable_test
		    if (nrow(estacion.variable.test) > 0) {
		      DBI::dbWriteTable(conn = private$con, name = "estacion_variable_test", 
		                        value = estacion.variable.test, append = TRUE, row.names = FALSE)
		    }
		    
		    # vii. Actualizar estados de registros diarios
		    for (i in seq_len(nrow(estacion.registro.diario.estado))) {
		      omm_id <- estacion.registro.diario.estado[i, "omm_id"]
		      fecha  <- estacion.registro.diario.estado[i, "fecha"]
		      estado <- estacion.registro.diario.estado[i, "estado"]
		      self$executeStatement(query = "UPDATE estacion_registro_diario SET estado = $1 WHERE omm_id = $2 and fecha = $3",
		                            parameters = list(estado, omm_id, fecha))  
		    }
		    
		    # viii. Actualizar estados de registros diarios
		    for (i in seq_len(nrow(estacion.registro.diario.faltantes))) {
		      omm_id      <- estacion.registro.diario.faltantes[i, "omm_id"]
		      fecha       <- estacion.registro.diario.faltantes[i, "fecha"]
		      variable_id <- estacion.registro.diario.faltantes[i, "variable_id"]
		      self$executeStatement(query = paste0("UPDATE estacion_registro_diario SET ", variable_id, " = $1 where omm_id = $2 and fecha = $3"),
		                            parameters = list(NA, omm_id, fecha)) 
		    }

  		  # ix. Finalizar transaccion
  		  DBI::dbCommit(conn = private$con)
		  }, error = function(e) {
		    # Error en la transaccion. Hacer rollback y luego mostrar el error.
		    DBI::dbRollback(conn = private$con)
		    stop(e)
		  })
 		}
	)
)
