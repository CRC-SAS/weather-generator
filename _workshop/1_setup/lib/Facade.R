require(R6)
require(dplyr)
require(dbplyr)

Facade <- R6Class("Facade",
	private = list(
	  con     = NULL,
	  filtrar = function(registros, filtros) {
	    for (col.name in names(filtros)) {
	      if (! is.null(filtros[[col.name]]) && (length(filtros[[col.name]]) > 0)) {
	        field.name <- rlang::sym(col.name)
	        if (length(filtros[[col.name]]) > 1) {
	          field.values <- filtros[[col.name]]
	          registros    <- registros %>% dplyr::filter(UQ(field.name) %in% field.values)
	        } else {
	          field.value <- filtros[[col.name]]
	          registros   <- registros %>% dplyr::filter(UQ(field.name) == field.value)
	        }
	      }
	    }
	    return (registros)
	  }
	),
	public = list(
	  getConnection = function() {
	    return (private$con)
	  },
	  setConnection = function(new.con) {
	    private$con <- new.con
	  },
	  executeStatement = function(query, parameters) {
	    res.query <- DBI::dbSendStatement(private$con, query)
	    DBI::dbBind(res.query, parameters)
	    DBI::dbClearResult(res.query)
	  }
	)
)
