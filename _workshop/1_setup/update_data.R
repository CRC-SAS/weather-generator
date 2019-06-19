
## set the working directory -- this is wherever you put the folder from the zip file
setwd("/home/dbonhaure/RStudioProjects/glmwgen/_workshop/1_setup/")

# -----------------------------------------------------------------------------#
# --- PASO 1. Cargar librerias propias e iniciar script ----

# a) Carga de clases de uso general
source(paste0("/home/dbonhaure/RStudioProjects/glmwgen/_workshop/1_setup/lib/", "Facade.R"), echo = FALSE)
source(paste0("/home/dbonhaure/RStudioProjects/glmwgen/_workshop/1_setup/lib/", "Estacion.R"), echo = FALSE)
source(paste0("/home/dbonhaure/RStudioProjects/glmwgen/_workshop/1_setup/lib/", "RegistroDiario.R"), echo = FALSE)

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 2. Conectar a base de datos e instanciar facades ----

# a. Conectar a la base de datos e indicar que el huso horario a utilizar es UTC
con <- DBI::dbConnect(drv = RPostgres::Postgres(),
                      dbname = "crcssa",
                      user = "crcssa_user",
                      password = "asdf1234",
                      host = "localhost")
DBI::dbExecute(con, "SET TIME ZONE 'UTC'") # Definir zona horaria como UTC

# b. Crear facade para manipular entidades (estaciones)
estacion.facade <- EstacionFacade$new(con)

# c. Crear facade para obtener datos de las estaciones
registro.diario.facade <- RegistroDiarioFacade$new(con)

# ------------------------------------------------------------------------------

PP2 = read.table('data/Salado_prcp.dat.original',header=T)
MX2 = read.table('data/Salado_tmax.dat.original',header=T)
MN2 = read.table('data/Salado_tmin.dat.original',header=T)

stations.omm_id <- c(87544,87448,87453,87534,87468,87484,87537,87548,87532,87540,87550,87616,87623,87624,87637,87640,87679)

# load(file='data/test/climate.RData') # climate.control <- climate
variables.id <- c('prcp', 'tmax', 'tmin')
anhos.desde  <- 1961
anhos.hasta  <- 2013
#
climate <- registro.diario.facade$buscar(omm_id = stations.omm_id, variable_id = variables.id) %>%
    dplyr::select(fecha, omm_id, variable_id, valor) %>%
    tidyr::spread(key = variable_id, value = valor) %>%
    dplyr::rename(date = fecha, station = omm_id, prcp = prcp, tx = tmax, tn = tmin) %>%
    dplyr::filter(date >= as.Date(glue::glue('{anhos.desde}-01-01'))) %>%
    dplyr::filter(date <= as.Date(glue::glue('{anhos.hasta}-12-31'))) %>%
    tidyr::complete(date = base::seq.Date(min(date), max(date), by = "days"), station) %>%
    glmwgen::exclude_incomplete_years() # solo se pueden tomar anhos completos!!

pp2 <- climate %>% dplyr::select(date, station, prcp) %>%
    dplyr::mutate(station = glue::glue("pr.{station}")) %>%
    tidyr::spread(key = station, value = prcp) %>%
    dplyr::select(colnames(PP2))
write.table(pp2, file='data/Salado_prcp.dat', sep='\t', row.names = F)

mx2 <-climate %>% dplyr::select(date, station, tx) %>%
    dplyr::mutate(station = glue::glue("tx.{station}")) %>%
    tidyr::spread(key = station, value = tx) %>%
    dplyr::select(colnames(MX2))
write.table(mx2, file='data/Salado_tmax.dat', sep='\t', row.names = F)

mn2 <-climate %>% dplyr::select(date, station, tn) %>%
    dplyr::mutate(station = glue::glue("tn.{station}")) %>%
    tidyr::spread(key = station, value = tn) %>%
    dplyr::select(colnames(MN2))
write.table(mn2, file='data/Salado_tmin.dat', sep='\t', row.names = F)

# -----------------------------------------------------------------------------#
# --- PASO 6. Finalizar script cerrando conexion a base de datos ----

# a) Cerrar conexion a base de datos
DBI::dbDisconnect(con)

# ------------------------------------------------------------------------------
