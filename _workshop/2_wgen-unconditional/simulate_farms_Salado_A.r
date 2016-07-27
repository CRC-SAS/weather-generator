
# ------------------------------------------------------------------------------
# --- Script to generate synthetic farm boundaries and other data
# --- for the Salado A basin.
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Clean up all objects ----

remove(list = ls()); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Install needed R packages ----

# update.packages()

if (!require(lubridate)) {install.packages("lubridate"); library(lubridate)}
if (!require(rgdal)) {install.packages("rgdal"); library(rgdal)}
if (!require(sp)) {install.packages("sp"); library(sp)}
if (!require(PBSmapping)) {install.packages("PBSmapping"); library(PBSmapping)}
if (!require(spdep)) {install.packages("spdep"); library(spdep)}
if (!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if (!require(spatstat)) {install.packages("spatstat"); library(spatstat)}
if (!require(rgeos)) {install.packages("rgeos"); library(rgeos)}
if (!require(rgdal)) {install.packages("rgdal"); library(rgdal)}
if (!require(maptools)) {install.packages("maptools"); library(maptools)}
if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if (!require(Hmisc)) {install.packages("Hmisc"); library(Hmisc)}
if (!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if (!require(RANN)) {install.packages("RANN"); library(RANN)}

#  update.packages(ask=FALSE)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Define function to trim number of farms within a department ----

# --- The list of synthetic farms to generate was derived from the AG Census
# --- and involves multiple assumptions.
# --- It may happen sometimes that the number of farms to generate
# --- within the department/partido exceeds the available area.
# --- This function compares the number of available 'unit tiles' in a partido
# --- with the number of tiles required by all farms to be simulated.
# --- The list of farms is trimmed to a certain proportion of the available tiles.

FilterFarms <- function(synth.farms, farm.comparable.field,
  max.aggregated.value,
  initial.seed = 0, max.rate = 0.95) {
  
  # synth.farms <- farm.table.this.dept
  # farm.comparable.field <- 'adjusted.farm.size'
  # max.aggregated.value <- area
  # initial.seed <- 0
  # max.rate <- 0.50
  
  set.seed(initial.seed)
  filtered.synth.farms <- synth.farms
  
  aggregated.farm.value <- sum(filtered.synth.farms[, farm.comparable.field])
 
   while (aggregated.farm.value > (max.rate * max.aggregated.value)) {
    row.number            <- sample(x = seq(from = 1, 
      to = nrow(filtered.synth.farms)), size = 1)
    filtered.synth.farms  <- filtered.synth.farms[-c(row.number), ]
    aggregated.farm.value <-
      sum(filtered.synth.farms[, farm.comparable.field])
  }
  
  return (filtered.synth.farms)
  
}
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Decide whether results will be plotted ----

plot.results <- FALSE
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Define Coordinate Reference Systems for various projections ----

gk.crs.string <- '+proj=tmerc +lat_0=-90 +lon_0=-60 +k=1 +x_0=5500000 +y_0=0
  +ellps=intl +twogs84=-148,136,90,0,0,0,0 +units=m +no_defs'

utm.crs.string <- '+proj=utm +zone=20H +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

ll.crs.string <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Read outline of Salado Basin A  ----

proj.dir <- getwd()

basin.dir <- "D:/ENSO/Projects/CNH3/data/shapes_cuenca/"

tt0 <- file.info(basin.dir)   # Info about directory where basin shapefile is located
if (!tt0$isdir)
  stop("ERROR: Specified directory for basin shapefile does not exist... check name...\n")              

setwd(basin.dir)

# --- Read shapefile with Basin A outline

cuenca.A.ll <- readOGR(".", layer="Salado_A_latlon", verbose = TRUE)

setwd(proj.dir)

# --- Convert basin shape to Gauss-Kruger coordinates (in meters)

cuenca.A.GK <- sp::spTransform(cuenca.A.ll,
  CRS = CRS(gk.crs.string)) 

# --- Convert basin shape to UTM coordinates (in meters)

cuenca.A.UTM <- sp::spTransform(cuenca.A.ll,
  CRS = CRS(utm.crs.string))

if (plot.results) { 
  plot(cuenca.A.ll, axes = TRUE)
  plot(cuenca.A.GK, axes = TRUE)
  plot(cuenca.A.UTM, axes = TRUE)
}

rm(proj.dir, basin.dir, tt0); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Create individual shapefiles for partidos ----

# --- Read corrected shapefiles for Argentina's departments

orig.dept.dir <- "D:/ENSO-Data/Other Data/ARG_departments/"

tt0 <- dir.create(orig.dept.dir, showWarnings = TRUE)

tt0 <- file.info(orig.dept.dir)   # Informacion sobre el directorio especificado
if (tt0$isdir) {    				      # Es un directorio?
  setwd(orig.dept.dir)						# Cambiar al directorio especificado
} else {
  stop("ERROR: Specified directory does not exist... verify the name...\n")
}          			
rm(tt0)

dept.archivo <- paste(orig.dept.dir, "deptscorrected.shp", sep="")
if (!file.exists(dept.archivo)) {
  stop("Specified SHP file for departments", dept.archivo, "does not exist...\n")
}

oldPath <- getwd() 
setwd(orig.dept.dir)
depts <- rgdal::readOGR(".", layer="deptscorrected", verbose = TRUE) # Read departments
setwd(oldPath)

rm(orig.dept.dir, dept.archivo); gc()

# --- Aislar la parte de datos de departamentos

pp2 <- depts@data

# --- Identificar los partidos/departamentos que contienen
# --- a la cuenca A del Salado

depts.en.cuenca <- c("9 DE JULIO","ADOLFO ALSINA","BOLIVAR","BRAGADO",
 "CARLOS CASARES","CARLOS TEJEDOR","COLON",
 "DAIREAUX","FLORENTINO AMEGHINO","GENERAL ARENALES",
	"GENERAL PINTO","GENERAL VIAMONTE","GENERAL VILLEGAS","GUAMINI",
	"HIPOLITO YRIGOYEN","JUNIN","LEANDRO N. ALEM", "LINCOLN","PEHUAJO",
	"PELLEGRINI","RIVADAVIA","SALLIQUELO","TRENQUE LAUQUEN","TRES LOMAS",
  "GENERAL ROCA","PRESIDENTE ROQUE SAENZ PENA",
	"ATREUCO","CATRILO","CHAPALEUFU","CONHELO",
	"MARACO","QUEMU QUEMU","TRENEL",
	"GENERAL LOPEZ")

codigos.deptos.en.cuenca <- c(6588,6007,6105,6112,
  6147,6154,6175,
  6231,6277,6294,
  6351,6385,6392,6399,
  6406,6413,6462,6469,6609,
  6616,6679,6721,6826,6847,
  14035,14084,
  42007,42028,42056,42035,
  42105,42119,42147,
  82042)

pp3 <- pp2$code %in% codigos.deptos.en.cuenca

deptos.en.cuenca <- pp2[pp3, ]  # Seleccionar departamentos dentro de cuenca

# --- Iterate to create a shape file for each department inside the basin

for (i in 1:nrow(deptos.en.cuenca)) {
  
  code <- deptos.en.cuenca[i, "code"]
  dept <- as.character(deptos.en.cuenca[i, "dept"])
  
  indiv.shape <- depts[depts@data$code == code, ]
  
  oldPath <- getwd()
  setwd("D:/ENSO/Projects/CNH3/data/shapes_depts/")
  writeOGR(indiv.shape, dsn=".", layer = dept, driver = "ESRI Shapefile",
    check_exists = TRUE, overwrite_layer = TRUE, verbose = TRUE)
  setwd(oldPath)   
}

# --- Merge all shapefiles for departments in the Salado A basin.
# --- Build a single "outer boundary" for all polygons.

library(rgeos); library(maptools)

all.depts <- depts[pp3, ]

dept.IDs <- all.depts@data$code

dissolved.depts <- maptools::unionSpatialPolygons(all.depts, IDs = dept.IDs)

outer.shape.ll <- maptools::unionSpatialPolygons(all.depts,
  IDs = rep(1, times= length(dept.IDs)))

outer.shape.GK <- sp::spTransform(outer.shape.ll,
  CRS = CRS(gk.crs.string) )

outer.shape.UTM <- sp::spTransform(outer.shape.ll,
  CRS = CRS(utm.crs.string) )

dept.dir <- "D:/ENSO/Projects/CNH3/data/deptos_en_cuenca/"

if (plot.results) { 
  
  Cairo::CairoPDF(file = paste0(dept.dir,"Salado_depts.pdf"),
    width = 6, height = 6, onefile = TRUE, family = "Arial")
  
  plot(dissolved.depts, axes = TRUE, border = "grey80")
  plot(cuenca.A.ll, add = TRUE, border = "tomato")
  plot(outer.shape.ll, add = TRUE, border = "steelblue")
  
  dev.off()
}

rm(all.depts,pp2,pp3,code,dept,indiv.shape,codigos.deptos.en.cuenca,
   i,oldPath,dept.IDs,dissolved.depts); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Read table with departments inside Salado Basin ----

file1 <- "Salado_A_departments.txt"
infile1 <- paste(dept.dir, file1, sep="")

if (file.exists(infile1)) {
  basin.depts <- read.table(infile1, sep = "\t", header = TRUE, as.is = TRUE)
  basin.depts <- dplyr::arrange(basin.depts, dept.code, dept.name)
} else {
  cat("File with list of Salado departments not found...\n")
}

dept.names <- basin.depts$dept.name

# --- Convert all upper case characters to upper and lower case

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
{s <- substring(s,2); if(strict) tolower(s) else s},
sep = "", collapse = " " )
sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

dept.names.lc <- capwords(dept.names, strict=TRUE) # Upper and lower case

dept.names.noblanks <- stringr::str_replace_all(dept.names, " ", "_")

# --- Let us make sure we have individual shape files for
# --- each department inside the basin.

dept.dir <- "D:/ENSO/Projects/CNH3/data/shapes_depts/"

dept.files <- list.files(path = dept.dir, pattern = ".shp$",
  full.names = TRUE)

dept.files.short <- list.files(path = dept.dir, pattern = ".shp$",
  full.names = FALSE)

dept.files.short <- sort(stringr::str_replace(dept.files.short,
  pattern = ".shp$",
  replacement = ''))

# --- Check that the department names read from a table
# --- have corresponding shapefiles.

if (!all.equal(sort(dept.names), dept.files.short)) {
  stop("Check departments' table and shapefile...\n")
}

n.of.depts <- length(dept.files)  # Number of depts in Salado Basin

remove(capwords, file1, infile1, dept.files.short); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Generate a regular grid inside the outer shape of departments ----

# --- Each grid point will represent the center of a tile or unit area.
# --- The grid generated extends outside the polygon; it will be trimmed later.
# --- The tiles will be the "building blocks" to assemble synthetic farms.
# --- NOTE THAT THE POINTS MAY HAVE DIFFERENT DELTA-X and DELTA-Y VALUES.
# --- That is, tiles created from these points may NOT be square, but rectangular.

tile.area <-  20                  # area represented by a point (in hectares)
point.spacing.x <- 1000  					# point separation along x-axis (in meters)
point.spacing.y <- 200						# point separation along y-axis (in meters)

if (((point.spacing.x * point.spacing.y) / 10000) != tile.area)
  stop("Error in point spacing along x and y dimensions...\n")

point.offset.x <- point.spacing.x / 2
point.offset.y <- point.spacing.y / 2

# --- Generate regular grid of points inside merged polygon
# --- The number of points in the grid is roughly equal to the area
# --- of the target region divided by the area encompassed by one grid point.

gg1 <- bbox(outer.shape.GK)		# bounding box of outer polygon
range.x <- range(pretty(range(gg1["x",]), n = 40, high.u.bias = 100))	# range of lon (x) values
range.y <- range(pretty(range(gg1["y",]), n = 40, high.u.bias = 100))	# range of lat (y) values

# --- The x- and y-values generated exceed the boundaries of the outer polygon
# --- so it can be easier to generate unit-area tiles afterwards...

x.values <- seq(from=(floor(range.x[1]) - (2 * point.spacing.x) + point.offset.x),
  to=(ceiling(range.x[2] + (2 * point.spacing.x))),
  by=point.spacing.x)

y.values <- seq(from=(floor(range.y[1]) - (2 * point.spacing.y) + point.offset.y),
  to=ceiling(range.y[2] + (2 * point.spacing.y)),
  by=point.spacing.y)

gg2 <- as.data.frame(expand.grid(lon = x.values, lat = y.values))

# --- Reorder the points so they are ordered along rows of constant latitude.
# --- Latitude is sorted in reverse so the first line is the northernmost lat.

gg2 <- dplyr::arrange(gg2, desc(lat), lon)
row.names(gg2) <- 1:nrow(gg2)	# matrix of coordinates for points in grid

# --- Create a spatial points object with what will be the centers
# --- of the tiles to be created later.

points.info <- data.frame(point.ID = row.names(gg2),
  x.coord = gg2[, "lon"], y.coord = gg2[, "lat"])

gg3 <- SpatialPointsDataFrame(coords = gg2, data = points.info,
  proj = CRS(gk.crs.string), match.ID = TRUE)

# --- If plot.results == TRUE, plot grid of points

if (plot.results) {
  plot(gg3, axes=TRUE, pch="16", col="tomato", cex=0.1)
  plot(outer.shape.GK, lwd=2, border="steelblue4", add=TRUE)
  title("Departments in Salado A Basin")
}

remove(gg1, gg2, range.x, range.y, x.values, y.values)
remove(point.offset.x, point.offset.y); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Create SpatialPolygons object with unit-area tiles ----

polys <- vector(mode="list", length = length(coordinates(gg3)[,1]))

# --- From the points, first we generate a grid topology,
# --- and then polygons around each point.

hh1 <- points2grid(gg3, round=NULL)

hh2 <- GridTopology(cellcentre.offset=hh1@cellcentre.offset,
  cellsize=hh1@cellsize,
  cells.dim=hh1@cells.dim)

hh3 <- as.SpatialPolygons.GridTopology(hh2,
  proj4string = CRS(gk.crs.string) )
gc()

hh4 <- sp::spChFIDs(hh3, x=as.character(points.info$point.ID))	# assign unique polygon IDs
gc()

hh5 <- points.info$point.ID

polys[hh5] <- slot(hh4[hh5],"polygons"); gc()

tiles.SP <- SpatialPolygons(polys,
  proj4string = CRS(gk.crs.string))
gc()

# --- Plot tiles and tile IDs

if (plot.results) {
  plot(tiles.SP, axes=TRUE, asp=1, border="grey70", lwd=0.5)
  title("Salado A Basin")
  plot(cuenca.A.UTM, add=TRUE, lwd=2, border="steelblue3")
  
  # plot(gg3,add=TRUE, pch=16, col="tomato", cex=0.5)
  # text(points.info[,"x.coord"], points.info[,"y.coord"],
  #	labels=points.info[,"point.ID"], cex=0.5)
}

rm(hh1,hh2,hh3,hh4,hh5,polys); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Identify and select points and tiles inside the outer shape of depts ----

# --- Select polygons that are COMPLETELY inside outer department shape. 

ccc <- rgeos::gWithin(tiles.SP, outer.shape.GK, byid = TRUE)
points.in <- which(ccc)
cc1 <- tiles.SP[points.in]       # Tiles inside outer shape

tiles.SP2 <- sp::spChFIDs(cc1,
  x = as.character(1:length(cc1)))	# assign unique polygon IDs

# --- Plot tiles inside target region...

if (plot.results) {
  plot(tiles.SP2[1:5000], axes = TRUE, col = "grey70", lwd = 0.5)
  plot(outer.shape.GK, lwd = 2, add = TRUE, border = "steelblue1")
  title("Salado A Basin")
}

remove(tiles.SP, ccc, cc1, points.in); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Read or build a 'master.tile.list' data frame ----

# --- This master list has information about all tiles in the outer shape,
# --- of departments, given a tile size and shape.
# --- If a text file exists, the list will be read at the start of each run,
# --- to avoid selecting  tiles used in a previous run.
# --- If the text file does not exist, create a 'master.tile.list' data frame
# --- with the tile information.

# --- Define name of text file with the master tile list. 

master.tile.outdir <- "D:/ENSO/Projects/CNH3/data/campos_sinteticos/tablas_datos/"
master.list.outfile <- paste(master.tile.outdir,"master_tile_data.txt", sep = "")

# --- Check if file exists. 
# --- If it exists, then read it... else, create the data frame 'master.tile.list'

# --- WARNING: If a department/partido has already been processed and
# --- the master list is read from a file, only few available tiles may be left
# --- for that department.
# --- If a department has to be re-run, make sure that master list is not read.

if (file.exists(master.list.outfile)) {
  # Master file list exists, then READ it....
  cat('Reading existing master tile file...')
  master.tile.list <- read.table(master.list.outfile,
    header = TRUE, sep = '\t', strip.white = TRUE)
} else {
  # Master file list DOES NOT exist, CREATE it....
  cat('Creating master tile data frame')
  tile.IDs <- as.numeric(sapply(slot(tiles.SP2, "polygons"),
    function(x) {slot(x,"ID")}))
  tile.coords <- coordinates(tiles.SP2)
  
  # --- Generate a data frame with ALL tiles inside outer shape of departments.
  # --- The list will be used to identify tiles already assigned to a farm.
  
  master.tile.list <- data.frame(tile.ID = tile.IDs,
    x.coord = tile.coords[, 1],
    y.coord = tile.coords[, 2],
    available = rep(TRUE, length.out = length(tile.IDs)),
    stringsAsFactors = FALSE)
}

table(master.tile.list$available, useNA = 'always')
# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# ***** Start working with farms *****
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# --- Read table with sizes of farms to be simulated ----

indir1 <- "D:/ENSO/R stuff/SaladoBasin/"    # Input dir for Federico's table
file1 <- "Salado_A_to_GIS_2002.txt"         # Input file for Federico's table

# --- Either read the table with farm sizes provided by Federico or,
# --- alternatively,
# --- read a test file with a few synthetic farm sizes for testing...

test.farm.sizes <- FALSE		# If TRUE, use test farm size data set

if (!test.farm.sizes) {			# Do not use test farm sizes
  # --- Table with farm sizes from Federico
  infile1 <- paste(indir1, file1, sep="")
  if (!exists("infile1")) {
    stop("ERROR: Specified input file", file1, "does not exist...\n")
  }
  farm.table <- read.table(infile1, sep="\t", na.strings=c("",NA), header=TRUE)
  # --- Retain only a few variables:   
  farm.table <- farm.table[,c("region.id","province.id",
    "department.id","farm.id","owner.id","eap.id",
    "original.farm.size", "adjusted.farm.size")]
} else {
  # --- Table with TEST farm sizes
  file1 <- "TEST.farm.sizes.txt"
  infile1 <- paste(indir1, file1, sep="")
  if (!exists("infile1")) {
    stop("ERROR: Specified input file", file1, "does not exist...\n")
  }
  farm.table <- read.table(infile1, sep="\t", na.strings=c("",NA), header=TRUE)
  # --- Retain only a few variables:
  farm.table <- farm.table[,c("region.id","farm.id","owner.id",
    "original.farm.size")]
}

rm(indir1, file1, infile1); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Adjust synthetic farm sizes ----
# --- so they are integer multiples of the unit tile size.

# --- If a farm is too large ( > 5000 has or 250 times the area of the unit tile)
# --- then truncate adjusted size to about 5000 has (final size depends on tile area).
# --- If a farm is too small ( < 20 has) then increase adjusted area to minimum area
# --- (one unit tile).

#  tile.area <- 20
max.farm.size <-  250 * tile.area # 5000 has with a tile of 20 has
min.farm.size <- tile.area        # about 20 has

farm.table <- dplyr::tbl_df(farm.table) %>%
  dplyr::mutate(adjusted.farm.size = round((original.farm.size /
    tile.area), 0) * tile.area) %>%
  dplyr::mutate(adjusted.farm.size = ifelse(adjusted.farm.size >= max.farm.size,
    max.farm.size, adjusted.farm.size)) %>%
  dplyr::mutate(adjusted.farm.size = ifelse(adjusted.farm.size <= min.farm.size,
    min.farm.size, adjusted.farm.size)) %>%
  dplyr::rename(prov.code = province.id, dept.code = department.id)

if (any(farm.table$adjusted.farm.size %% tile.area != 0)) {
  stop('There are adjusted areas not integer multiples of tile area')
}

# --- Check the areas of all adjusted farm sizes against
# --- the total area of departments considered.

sum.farms <- sum(farm.table$adjusted.farm.size)  # Sum of adjusted sizes

sum.depts <- ceiling(sum(basin.depts$area) * 100)  # sum of department areas in km^2, expressed in has

if (sum.depts >= sum.farms) {
  cat("Sum of adjusted farm areas is lower than sum of total departments area...\n")
} else {
  stop("Sum of adjusted farm areas is HIGHER than sum of total departments area!\n")
}

rm(min.farm.size, max.farm.size, sum.farms, sum.depts); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Identify departments with problems in farm areas ----

# --- Identify those departments for which area of farms to simulate is
# --- larger or equal than the total department area.

ddd <- dplyr::tbl_df(farm.table) %>%
  dplyr::group_by(prov.code, dept.code) %>%
  dplyr::summarise(count.a = n(),
    sum.o = round(sum(original.farm.size), 0),
    sum.a = sum(adjusted.farm.size) )

dd2 <- dplyr::inner_join(ddd, dplyr::tbl_df(basin.depts)) %>%
  dplyr::mutate(dept.area = (area * 100),
    prop = round(sum.a / (area * 100), 3)) %>%
  dplyr::select(prov.code, prov.name, dept.code, dept.name,
    count.a, sum.o, sum.a, dept.area, prop) %>%
  dplyr::arrange(prov.code, count.a)

dd3 <- dplyr::filter(dd2, prop >= 0.97)

if (any(dd2$prop >= 1)) {
  warning('Farm areas > dept size for SOME depts...\n')
  print(dd3)
}

rm(ddd, dd2, dd3); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Make a dotchart with number of farms to be simulated by department ----
# --- Se usa un dotplot algo mas sofisticado, usando funcion dotchart2()
# --- en paquete Hmisc. Los departamentos se agrupan por provincia.

dept.dir <- "D:/ENSO/Projects/CNH3/data/shapes_depts/"

if (plot.results) { 
  dd1 <- dplyr::tbl_df(farm.table) %>%
    dplyr::group_by(dept.code) %>%
    dplyr::summarise(n.farms = n()) %>%
    dplyr::arrange(n.farms)
  
  dd3 <- dplyr::inner_join(dd1, dplyr::tbl_df(basin.depts)) %>%
    dplyr::select(dept.name, prov.name, n.farms) %>%
    dplyr::arrange(prov.name, desc(n.farms))
  
  dd3$dept.name <-stringr::str_to_title(dd3$dept.name)
  
  Cairo::CairoPNG(filename = paste0(dept0.dir,"synth_farms_dotplot.png"),
    width = 840, height = 940,
    pointsize = 12, bg = "white",  res = NA)
  
  Hmisc::dotchart2(dd3$n.farms,
    labels = dd3$dept.name,
    groups = dd3$prov.name,
    auxdata = dd3$n.farms,
    auxtitle="Synth\n farms",
    xlab="Number of farms",
    main="Number of synthetic farms - Salado A Basin",
    dotsize = 1.1,
    cex.labels = 0.6,
    cex.group.labels = 0.8,
    sort. = FALSE,
    pch = 16, col = "tomato")
  
  dev.off()
  rm(dd1, dd2, dd3); gc()
  
} # Fin de if que define si se grafican resultados
# ------------------------------------------------------------------------------


# --- Iterate over each department inside the Salado A basin

# for (iii in 1:n.of.depts) {

for (iii in 18:34) {

# iii <- 20

# -----------------------------------------------------------------------------#
# --- Read shape file for THIS department ----

dept.dir <- "D:/ENSO/Projects/CNH3/data/shapes_depts/"

oldPath <- getwd()
setwd(dept.dir)
cat("Reading shape file for department", dept.names[iii], "...\n")
dept.boundary <- rgdal::readOGR(".", layer = dept.names[iii], verbose = TRUE)
setwd(oldPath)

dept.data <- slot(dept.boundary, "data") # Data associated with THIS department

dept.name <- dept.names[iii]
dept.name.lc <- dept.names.lc[iii]

rm(oldPath); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#  
# --- Convert boundary for THIS department into GK and UTM coordinates (in meters) ----

dept.boundary.GK <- sp::spTransform(dept.boundary,
    CRS = CRS(gk.crs.string))

dept.boundary.UTM <- sp::spTransform(dept.boundary,
  CRS = CRS(utm.crs.string))

# --- Plot outer boundary of THIS department in GK coordinates
  
if (plot.results) {
  plot(dept.boundary.GK, axes=TRUE)
  plot(dept.boundary.GK, border="steelblue4",
    col=rgb(246/255, 232/255, 195/255),
    lwd=2, add=T)
  title(dept.name.lc)
}

# --- Plot outer boundary of THIS department in lat-lon coordinates

if (plot.results) {
  plot(dept.boundary, axes=TRUE, border="steelblue4",
    col=rgb(246/255, 232/255, 195/255))
  title(dept.name.lc)
}
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Calculate area of department in hectares ----
# --- Divide by 10,000 because 1 ha = 10,000 m^2
# --- Also, 1 km^2 = 100 ha

# area <- (sapply(slot(dept.boundary.GK,"polygons"),
#   function(x) {slot(x,"area")} )) / 10000.
  
area <- ceiling(rgeos::gArea(dept.boundary.GK) / 10000.)
  
cat(paste("Area of department", dept.name.lc, "is", area, "hectares...\n"))
cat(paste("Area of department", dept.name.lc, "is", area/100, "km^2...\n"))
# ------------------------------------------------------------------------------ 


# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# --- NOW LET US GENERATE PSEUDO_FARMS
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# --- Identify and select tiles inside THIS department ----

tiles.SP3 <- tiles.SP2[dept.boundary.GK, ]     # Tiles inside polygon

tile.IDs.this.dept <- sapply(slot(tiles.SP3, "polygons"),
  function(x) {slot(x,"ID")})

# --- Tiles NOT available within entire region.
# --- Some tiles inside THIS department may have been used
# --- while processing previous departments.

tiles.unavailable.all.region <- as.numeric(master.tile.list$tile.ID[
  master.tile.list$available == FALSE])

uu1 <- tile.IDs.this.dept %in% tiles.unavailable.all.region

if (any(uu1)) {
  # Must eliminate tiles from this dept used before...
  cat('SOME tiles in this department have been used previously...\n')
  cat('The number of tiles inside this department will be reduced accordingly.\n')
  tiles.SP3 <- tiles.SP3[!uu1]
} else {
  # No need to eliminate any tiles wiuthin this dept
  cat('ALL tiles inside department are available...')
}

n.tiles.in.dept <- length(tiles.SP3)  # Number of tiles inside polygon

# --- Plot tiles inside target region...

if (plot.results) {
  sp::plot(dept.boundary.GK, axes = TRUE, border = 'steelblue1')
  sp::plot(tiles.SP3, add = TRUE)
  title(dept.name.lc)
}
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Build table of farms (trimmed, if needed) for THIS department ----

# --- Select farms inside THIS department (the one we are working on)

tt1 <- match(dept.name, basin.depts$dept.name)
this.dept.code <- as.numeric(basin.depts[tt1, "dept.code"])

farm.table.this.dept <- dplyr::tbl_df(farm.table) %>%
  dplyr::filter(dept.code == this.dept.code) %>%
  dplyr::mutate(tiles.required = round((adjusted.farm.size / tile.area), 0))  %>%
  dplyr::arrange(desc(adjusted.farm.size))

n.sim.farms <- nrow(farm.table.this.dept)  # number of simulated farms for THIS dept

n.tiles.required.this.dept  <- sum(farm.table.this.dept$tiles.required)

# --- Check if number of tiles within THIS dept is sufficient
# --- to place all farms in the dept.
# --- The number of tiles needed for all farms is compared to 
# --- the value of 'n.tiles.in.dept', calculated in the previous step.
# --- This number of tiles lists those tiles AVAILABLE for this
# --- department, as some tiles may have been used when assembling farms
# --- for an adjacent department.

cat(paste('Number of tiles in THIS department:', n.tiles.in.dept,'\n'))
cat(paste('Number of tiles needed for farms:', n.tiles.required.this.dept,'\n'))

# Proportion of available tiles in dept that can assigned to farms
max.prop.tiles <- 0.97 

if (n.tiles.required.this.dept > n.tiles.in.dept) {
  # Are there sufficient tiles inside the department?
  warning("Number of tiles required is > tiles in THIS department")
  warning("List of farms MUST be adjusted to fit in department")
  
  cat('Trimming list of farms within THIS department...\n')

  farm.table.this.dept <- FilterFarms(farm.table.this.dept,
    farm.comparable.field = 'tiles.required',
    max.aggregated.value = n.tiles.in.dept,
    initial.seed = 0,
    max.rate = max.prop.tiles)
  
  n.sim.farms <- nrow(farm.table.this.dept)  # number of simulated farms for THIS dept
  cat(paste('Trimmed number of farms within THIS department:',
    n.sim.farms, '\n'))
  
} # End of if that checks if total number of tiles required is > tiles in department

if (plot.results) {
  hist(farm.table.this.dept$adjusted.farm.size,
    breaks = c(0,100,250,500, 1000, 2500, 5000.01),
    freq = FALSE, col = "thistle",
    xlab = "Simulated farm areas (ha)")
  box()
}

table(farm.table.this.dept$tiles.required)

rm(tt1, this.dept.code, max.prop.tiles)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Build a list to store tile IDs used for each simulated farm ----
# --- The list has length equal to the number of farms to be simulated.

tile.IDs.used <- vector(mode = "list", length = n.sim.farms)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Build a data frame with info on tiles in THIS department ----
# --- The data frame lists the IDs of tiles, their x- and y-coordinates.
# --- A fourth column indicates if tiles are available to be used
# --- to construct a farm. Initially all tiles are set to TRUE.

tile.IDs <- sapply(slot(tiles.SP3, "polygons"), function(x) {slot(x, "ID")})
tile.coords <- coordinates(tiles.SP3)

tiles.available.this.dept <- data.frame(tile.ID = tile.IDs,
  x.coord = tile.coords[, 1],
  y.coord = tile.coords[, 2],
  available = rep(TRUE, length.out = length(tile.IDs)),
  stringsAsFactors = FALSE)

tiles.available.this.dept.2 <- rep(TRUE, length.out = length(tile.IDs),
  stringsAsFactors = FALSE)
names(tiles.available.this.dept.2) <- tile.IDs
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Find the k nearest neighbors for each tile inside THIS department ----
# --- inside the region considered.
# --- The number of neighbors found (k) is equal to the highest number
# --- of tiles required (i.e., the number of tiles needed to simulate
# --- the largest synthetic farm) plus an extra (here, 10%).
# --- The calculation returns a matrix that has (a) a number of rows equal to
# --- the number of tiles inside the region and (b) k columns (for the k nearest neigbors).
# --- NOTE: the values of the elements of the matrix of nearest neighbors are the _indices_
# --- of tiles inthe input object.

cc3 <- spdep::knearneigh(tile.coords,
  k = ceiling(max(farm.table.this.dept[, "tiles.required"]) * 1.3),
  longlat = FALSE,
  RANN = TRUE); gc()

nn.indices <- cc3$nn    # Nearest neighbors that are available
row.names(nn.indices) <- tiles.available.this.dept[, "tile.ID"]

# --- We add one column at the beginning of the matrix of neighbors
# --- that includes the index of the tile for which neighbors
# --- are found.

nn.indices.with.base.tile <- cbind((1:nrow(nn.indices)), nn.indices)

rm(cc3);gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Iterate over the number of synthetic farms to be generated ----

set.seed(256)

for (i in 1:n.sim.farms) {
 
  ##  i <- 1

  cat(paste("Simulating farm", i, "of", n.sim.farms, "...\n"))
  sim.farm.area <- farm.table.this.dept[i, "adjusted.farm.size"]  # area of THIS simulated farm
  n.tiles <- as.numeric(farm.table.this.dept[i, "tiles.required"])  # how many unit tiles  required for THIS farm?
  cat("Tiles needed for this farm:", as.numeric(n.tiles),"\n")
  
  # Find IDs of tiles that are still available for allocation to a farm
  tile.IDs.still.available <- dplyr::filter(tiles.available.this.dept, available == TRUE)
  
  if (nrow(tile.IDs.still.available) < n.tiles) {
    # There are no tiles available for this farm
    stop(paste("No tiles available for farm", i, "of", nrow(farm.table.this.dept)))
  
  } else {
    
    uu0 <- tile.IDs.still.available[, "tile.ID"]
    cat("Tiles still available:", length(uu0), "...\n")
    
    if (n.tiles == 1) {
      # Only ONE tile needed for this farm...
            
      # --- ID for the first (and, here, only) tile in a simulated farm
      first.tile.ID <- sample(uu0, size = 1)
      
    } else {
      # More than one tile needed for this farm...
      
      n.neighbors.still.available <- apply(X = nn.indices.with.base.tile, MARGIN = 1,
        FUN = function(x) {length(x[!is.na(x)])})
      indices.tiles.with.enough.neighbors <- which(n.neighbors.still.available >= n.tiles)
      names(indices.tiles.with.enough.neighbors) <- NULL
      
      # --- ID for the first tile in a simulated farm
      
      first.tile.ID <- sample(as.numeric(
        tiles.available.this.dept[indices.tiles.with.enough.neighbors, "tile.ID"]),
        size = 1)
      
      # Select the first tile selected plus "n.tiles - 1" nearest neighbors.
      # Remember the first column in "nn.indices.with.base.tile" is the index of the
      # first tile selected.
      
    } # End of else for only one tile needed
 
    # --- Find out the row number for the tile ID selected
    first.tile.row <- which(tiles.available.this.dept[, "tile.ID"] == first.tile.ID)
    
    oo1 <- nn.indices.with.base.tile[first.tile.row, ]
    oo2 <- oo1[!is.na(oo1)]
    # Indices and positions of all tiles used for THIS farm
    all.tiles <- oo2[1:n.tiles]
    all.tiles.IDs <- as.numeric(tiles.available.this.dept[all.tiles, "tile.ID"])
     
  } # End of else for tiles still available

  # --- Now that we have finished selecting tiles for THIS synthetic farm
  # --- let's save the tiles used for it in a list.
  
  tile.IDs.used[[i]] <- all.tiles.IDs  # store points used for this farm
  
  # --- Mark those tiles used for THIS farm
  # --- as unavailable for use in other farms
  tiles.available.this.dept[all.tiles, "available"] <- FALSE
  
  # --- In the list of available neighbors, mark tiles used as unavailable (NA)
  
  for (k in 1:nrow(nn.indices.with.base.tile)) {
    uu1 <- nn.indices.with.base.tile[k, ] %in% all.tiles
    nn.indices.with.base.tile[k, which(uu1)] <- NA
  }
  
  gc()

}   # End of iterating over number of farms to be simulated

# --- Update the master tile list (for entire area) setting to FALSE
# --- column 'available' for those rows which have been used in this department.

uuu <- as.numeric(tiles.available.this.dept[tiles.available.this.dept$available == FALSE, 'tile.ID'])
uu2 <- match(uuu, as.numeric(master.tile.list$tile.ID))
master.tile.list[uu2, 'available'] <- FALSE

table(master.tile.list$available)

remove(oo1,oo2,first.tile.ID,first.tile.row,
       tile.IDs.still.available, uuu, uu2)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Check that there are no duplicates among selected tile IDs ----

oo5 <- sort(as.numeric(unlist(tile.IDs.used)))
if (any(duplicated(oo5))) {stop("ERROR: Duplicated tile IDs...\n")}

oo5[which(duplicated(oo5))]

# --- Check that areas of generated farms coincide with those specified

ww1 <- lapply(tile.IDs.used, FUN=length)
ww2 <- unlist(ww1) * tile.area

if (any(ww2 != farm.table.this.dept[,"adjusted.farm.size"])) {
  stop("ERROR: Specified and simulated farm areas do not coincide...\n")
  ww3 <- which(ww2 != farm.table.this.dept[, "adjusted.farm.size"])
  #ww2[ww3]
  #sim.farms[ww3,"area"]
  rm(ww3)
}

remove(ww1,ww2,oo5); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Plot tiles used for each simulated farm and their center points ----


if (plot.results) {
  
  uu1 <- unlist(tile.IDs.used)
  uu2 <- which(uu1 %in% tile.IDs)
  
  library(RColorBrewer)
  farm.colors <- sample(brewer.pal(8, "Pastel2"),
    size = n.sim.farms, replace = TRUE)

  plot(dept.boundary.GK, axes=TRUE, asp=1)
  plot(dept.boundary.GK, lwd=2, border="steelblue4", add=T)
  title(dept.name.lc)
  
  for (i in 1:n.sim.farms) {
    ppp <- sort(tile.IDs.used[[i]])
    pp2 <- match(ppp, tile.IDs)
    plot(tiles.SP3[pp2, ], add=TRUE,
      col = farm.colors[i], border = farm.colors[i], lwd=0.01)
  }	# End of iteration through simulated farms
  
  remove(ppp, pp2, farm.colors)
  
}   # End of "if plot.results"

rm(uu1, uu2); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Merge all unit polygons for each farm into a single polygon ----
# --- the unionSpatialPolygons() function produces an object of
# --- SpatialPolygon class, so we extract the polygon.

list.of.polygons <- vector("list", length = n.sim.farms)  # store farm polygons here

for (i in 1:n.sim.farms) {
  
  #  i <- 25
  sim.farm.ID <- as.numeric(farm.table.this.dept[i, "farm.id"])
  
  cat("Working on farm",i,"of", n.sim.farms,"... Farm ID:", sim.farm.ID,"...\n")
 
  pp1 <- sort(as.numeric(tile.IDs.used[[i]]))	# Tile IDs used for THIS farm
  pp2 <- match(pp1, as.numeric(tiles.available.this.dept[, "tile.ID"]))	# Row numbers in tiles object
  
  pp3 <- maptools::unionSpatialPolygons(tiles.SP3[pp2],
    IDs = as.numeric(rep(sim.farm.ID, times = length(pp2))),
    threshold = NULL)
  
  chk <- sapply(slot(pp3, "polygons"), function(x) slot(x, "ID"))
  
  if (length(chk) != 1)
    stop("ERROR: More than one polygon...")
  
  pp5 <- pp3@polygons   			# extract polygon
  list.of.polygons[i] <- pp5  # store polygon for THIS farm into list
  
}	# End of iterating through simulated farms to merge all unit polygons

remove(pp1, pp2, pp3, pp5, chk); gc()
remove(nn.indices, tile.coords, tiles.available.this.dept); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Assemble polygons for all simulated farms ----
# --- into a single SpatialPolygons object

farms.SP <- SpatialPolygons(list.of.polygons)
slot(farms.SP, "proj4string") <- CRS(gk.crs.string)

# --- In some cases, the same polygon may have two Polygons inside
# --- after performing the union of the individual SpatialPolygons.
# --- This may mean a farm formed by tiles that are geographically disjoint.
# --- WARNING: There may be issues with calculation of centroids in these cases.

nn1 <- sapply(farms.SP@polygons, function(x) {length(slot(x,"Polygons")) })
nn2 <- which(nn1 > 1)

# --- Check if more than one polygon with same ID

chk <- sapply(slot(farms.SP, "polygons"), function(x) slot(x, "ID"))
chk2 <- any(duplicated(chk))
if (chk2) stop("ERROR: More than one polygon in farms.SP...\n")

remove(nn1,nn2,chk,chk2,list.of.polygons); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Convert polygon coordinates from Gauss-Kruger to lat-lon and UTM ----

# Latitude and longitude

farms.SP2 <- sp::spTransform(farms.SP, CRS(ll.crs.string))

# UTM

farms.SP3 <- sp::spTransform(farms.SP, CRS(utm.crs.string))
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Find centroids for farms in disjoint polygons ----
# --- There are a few farms that involve disjointed polygons (see above).
# --- In these cases, the centroid does not coincide with the "labpt" slot.

# --- First do it in Gauss-Kruger coordinates

farms.centroids.SP <- rgeos::gCentroid(farms.SP, byid = TRUE)

# --- Now do it in lat-lon coordinates

farms.centroids.SP2 <- rgeos::gCentroid(farms.SP2, byid = TRUE)

# --- Now do it in UTM coordinates

farms.centroids.SP3 <- rgeos::gCentroid(farms.SP3, byid = TRUE)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Plot simulated farms ----
# --- using Spatial Polygons objects generated so far...

if (plot.results) {
  
  library(RColorBrewer)
  farm.colors <- sample(brewer.pal(8, "Pastel2"), size=n.sim.farms, replace=TRUE)
  
  # --- First plot in GK coordinates
  
  plot(dept.boundary.GK, axes = TRUE, border = "steelblue", lwd = 2)
  plot(farms.SP, lwd = 0.2, col = farm.colors, border = "darkgrey", add = TRUE)
  plot(farms.centroids.SP, add=TRUE, pch=16, cex=0.3, col="tomato")
  plot(dept.boundary.GK, add = TRUE, border = "steelblue", lwd = 2)
  title(paste("Synthetic Farms -", dept.name.lc))
  
  # --- Then plot in lat-lon coordinates
  
  plot(dept.boundary, axes = TRUE, border = "steelblue", lwd = 2)
  plot(farms.SP2, lwd = 0.2, col = farm.colors, border = "darkgrey", add = TRUE)
  plot(farms.centroids.SP2, add=TRUE, pch=16, cex=0.5, col="tomato")
  plot(dept.boundary, add = TRUE, border = "steelblue", lwd = 2)
  title(paste("Synthetic Farms -", dept.name.lc))
  
  # --- Then plot in UTM coordinates
  
  plot(dept.boundary.UTM, axes = TRUE, border = "steelblue", lwd = 2)
  plot(farms.SP3, lwd = 0.2, col = farm.colors, border = "darkgrey", add = TRUE)
  plot(farms.centroids.SP3, add=TRUE, pch=16, cex=0.5, col="tomato")
  plot(dept.boundary.UTM, add = TRUE, border = "steelblue", lwd = 2)
  title(paste("Synthetic Farms -", dept.name.lc))
  
  #farm.id <- sapply(slot(farms.SP2, "polygons"), function(x) {slot(x, "ID")})
  #pp1 <- t(sapply(slot(farms.SP2, "polygons"), function(x) {slot(x, "labpt")}))
  #text(pp1, labels=farm.id, cex=0.4)
  
}

gc()
# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# --- UP TO THIS POINT WE HAVE CREATED SPATIAL POLYGONS
# --- WITH SYNTHETIC FARM BORDERS.
# --- NOW WE NEED TO ASSEMBLE A DATA FRAME THAT WILL GO WITH THE SPATIAL DATA
# --- AND THAT WILL HAVE ADDITIONAL INFO ABOUT EACH SYNTHETIC FARM.
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# --- Create an empty data frame to store info about each simulated farm ----

farms.DF <- data.frame(farmId = rep(NA, n.sim.farms),
  eapId = rep(NA, n.sim.farms),
  regnId = rep(NA, n.sim.farms),
  provId = rep(NA, n.sim.farms),
  provnm = rep(NA, n.sim.farms),
  deptId = rep(NA, n.sim.farms),
  deptnm = rep(NA, n.sim.farms),
  ownerId = rep(NA, n.sim.farms),
  orgFrmSz = rep(NA, n.sim.farms),
  adjFrmSz = rep(NA, n.sim.farms),
  GKlon = rep(NA, n.sim.farms),
  GKlat = rep(NA, n.sim.farms),
  deglon = rep(NA, n.sim.farms),
  deglat = rep(NA, n.sim.farms),
  UTMlon = rep(NA, n.sim.farms),
  UTMlat = rep(NA, n.sim.farms),
  perim = rep(NA, n.sim.farms),
  soilUnit = rep(NA, n.sim.farms),
  unitType = rep(NA, n.sim.farms),
  prodIndx = rep(NA, n.sim.farms),
  pctgSl1 = rep(NA, n.sim.farms),
  sgrpSl1 = rep(NA, n.sim.farms),
  pctgSl2 = rep(NA, n.sim.farms),
  sgrpSl2 = rep(NA, n.sim.farms),
  pctgSl3 = rep(NA, n.sim.farms),
  sgrpSl3 = rep(NA, n.sim.farms),
  creaGpId = rep(NA, n.sim.farms),
  advsrId= rep(NA, n.sim.farms) )
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Start filling in data frame with info about each simulated farm ----

# --- First, complete information using the synthetic farms' Spatial Polygons objects.

farms.DF[, "farmId"] <- as.numeric(sapply(slot(farms.SP2, "polygons"),
  function(x) {slot(x, "ID")}))
farms.DF[, "adjFrmSz"] <- sapply(slot(farms.SP, "polygons"),
  function(x) {slot(x, "area")}) / 10000
# Extract Gauss-Kruger coordinates for farm centroids from farms.centroids.SP3 object
farms.DF[, "GKlon"] <- round(coordinates(farms.centroids.SP)[,1], 0)
farms.DF[, "GKlat"] <- round(coordinates(farms.centroids.SP)[,2], 0)
# Extract lat-lon coordinates for farm centroids from farms.centroids.SP2 object
farms.DF[, "deglon"] <- round(coordinates(farms.centroids.SP2)[,1], 3)
farms.DF[, "deglat"] <- round(coordinates(farms.centroids.SP2)[,2], 3)
# Extract UTM coordinates for farm centroids from farms.centroids.SP object
farms.DF[, "UTMlon"] <- coordinates(farms.centroids.SP3)[,1]
farms.DF[, "UTMlat"] <- coordinates(farms.centroids.SP3)[,2]

# --- Second, complete  information using "farm.table.this.dept"

tt11 <- farm.table.this.dept$farm.id
tt12 <- match(farms.DF[,"farmId"], tt11)

farms.DF[, "eapId"] <- farm.table.this.dept[tt12,"eap.id"]									# Fill in EAP ID
farms.DF[, "regnId"] <- farm.table.this.dept[tt12,"region.id"]							# Fill in region ID
farms.DF[, "ownerId"] <- farm.table.this.dept[tt12,"owner.id"]							# Fill in farm owner ID
farms.DF[, "orgFrmSz"] <- round(farm.table.this.dept[tt12,"original.farm.size"], 1)		# Fill in farm size in original table

remove(tt11, tt12); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Calculate perimeter of each farm and store in data frame ----
# --- To do this, we must convert data into polygons used in PBSmapping library.
# --- In some cases, a farm includes separate polygons (e.g., a separate farm lot), so we need
# --- to sum up perimeters for all polygons coresponding to a given farm ID.

farms.PS <- SpatialPolygons2PolySet(farms.SP)			# Convert spatialpolygons into PBS PolySet
farms.perimeter <- calcLength(farms.PS, rollup=1)	# calculate perimeter (in meters)
farms.perimeter2 <- tapply(farms.perimeter[,"length"], INDEX=farms.perimeter[,"PID"], FUN=sum)	# add up perimeter for isolated polygons

farms.DF[, "perim"] <- as.vector(round(farms.perimeter2, 0))
# table(farms.DF$adjFrmSz, farms.DF$perim)

remove(farms.PS, farms.perimeter, farms.perimeter2); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Fetch province corresponding to each farm centroid ----
# --- First, read CORRECTED shape files with Argentine provinces

prov.dir <- "D:/ENSO-Data/Other Data/ARG_provinces/"
prov.file <- paste(prov.dir, "provs.corrected.shp", sep="")
if (!exists("prov.file")) {
  stop("ERROR: Specified SHP file for provinces", prov.file, "does not exist...\n") }

oldPath <- getwd()
setwd(prov.dir)
pp1 <- readOGR(".", layer="provs.corrected", verbose = TRUE)
pp1 <- spTransform(pp1, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
setwd(oldPath)

# --- Isolate the data part of provinces file

pp2 <- pp1@data

# --- Use centroids of farms to extract province data

pp6 <- over(farms.centroids.SP2, pp1)

farms.DF[, "provId"] <- as.numeric(pp6[, "indec.code"])		# Fill in INDEC province code
farms.DF[, "provnm"] <- as.character(pp6[,"prov.name"])		# Fill in province name

remove(prov.dir, prov.file, oldPath, pp1, pp2, pp6); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Fetch department corresponding to each farm centroid ----
# --- First, read CORRECTED shape files with Argentine departments

dept.dir <- "D:/ENSO-Data/Other Data/ARG_departments/"
dept.file <- paste(dept.dir, "deptscorrected.shp", sep="")
if (!exists("dept.file")) {
  stop("ERROR: Specified SHP file for departments", dept.file, "does not exist...\n") }

oldPath <- getwd()
setwd(dept.dir)
pp1 <- readOGR(".", layer="deptscorrected", verbose = TRUE)
pp1 <- spTransform(pp1, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
setwd(oldPath)

# --- Isolate the data part of department file

pp2 <- pp1@data

# --- Use centroids of farms to extract department data

pp6 <- over(farms.centroids.SP2, pp1)

dept.code1 <- as.numeric(pp6[, "code"])
dept.name1 <- as.character(pp6[, "dept"])

# --- Convert all upper case characters to upper and lower case

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
  {s <- substring(s,2); if(strict) tolower(s) else s},
  sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

dept.name1 <- capwords(dept.name1, strict=TRUE)

farms.DF[, "deptId"] <- dept.code1		# Fill in department code
farms.DF[, "deptnm"] <- dept.name1		# Fill in department name

remove(dept.dir, dept.file, oldPath, pp1, pp2, pp6, dept.code1, dept.name1); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Add soil data for each synthetic farm from INTA soil atlas ----

# --- First, read data files downloaded from INTA.
# --- URL: http://geointa.inta.gov.ar/publico/INTA_SUELOS/suelos_500000_v8.zip.
# --- These files have soils by province, NOT by department.
# --- WARNING: "soils by department" file has different column names.

soils.dir <- "D:/ENSO-Data/Other Data/ARG_soils_by_province/"
oldPath <- getwd()
setwd(soils.dir)
soils <- readOGR(".", layer="suelos_500000_v9", verbose = TRUE)
soils <- spTransform(soils, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
setwd(oldPath)

# --- Get table with soil data

tt1 <- soils@data

# --- Extract soil data for farm centroids

pp7 <- over(farms.centroids.SP2, soils)

# --- Fill in dataframe with soil data

farms.DF[, "soilUnit"] <- as.character(pp7[, "SIMBC"])
farms.DF[, "unitType"] <- as.character(pp7[, "TIPO_UC"])
farms.DF[, "prodIndx"] <- as.numeric(pp7[, "IND_PROD"])
farms.DF[, "pctgSl1"] <- as.numeric(pp7[, "PORC_SUE1"])
farms.DF[, "pctgSl2"] <- as.numeric(pp7[, "PORC_SUE2"])
farms.DF[, "pctgSl3"] <- as.numeric(pp7[, "PORC_SUE3"])
farms.DF[, "sgrpSl1"] <- as.character(pp7[, "SGRUP_SUE1"])
farms.DF[, "sgrpSl2"] <- as.character(pp7[, "SGRUP_SUE2"])
farms.DF[, "sgrpSl3"] <- as.character(pp7[, "SGRUP_SUE3"])

remove(soils.dir,oldPath,soils,tt1,pp7); gc()
# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# --- WE ARE DONE WITH GENERATION OF THE ANCILLARY RECORDS FOR EACH FARM.
# --- NOW JOIN THE POLYGONS WITHIN THE DATA FRAME AND SAVE THE RESULTS.
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# --- Reorder rows in farms data frame ----
# ---- so they coincide with the order of farm IDs in polygons.

tt1 <- order(farms.DF[, "farmId"])
farms.DF <- farms.DF[tt1, ]

IDs <- as.numeric(sapply(slot(farms.SP, "polygons"), function(x) slot(x, "ID")))

pp1 <- match(IDs, farms.DF$farmId)

farms.DF <- farms.DF[pp1, ]
row.names(farms.DF) <- as.character(farms.DF[, "farmId"])

remove(tt1, pp1, IDs); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Create Spatial Polygons Data Frame objects in various coordinates ----

farms.SPDF <- SpatialPolygonsDataFrame(farms.SP, farms.DF)

# --- Check that polygon IDs and farm IDs in data slot coincide

tile.IDs.this.dept <- sapply(slot(farms.SPDF, "polygons"),
  function(x) {slot(x,"ID")})

uu0 <- as.numeric(tile.IDs.this.dept) == farms.SPDF$farmId
if (any(! uu0)) {
  Log.error('Polygon IDs do not coincide with farm ID in data slot')
}

# --- Convert coordinates from GK (meters) back to lat-lon
# --- and create a second Spatial Polygons Data Frame.

farms.SPDF2 <- spTransform(farms.SPDF, CRSobj	= CRS(ll.crs.string) )

# --- Convert coordinates from Gauss-Kruger (meters) to UTM coordinates
# --- and create a new Spatial Polygons Data Frame.

farms.SPDF3 <- spTransform(farms.SPDF, CRS(utm.crs.string) )

rm(tile.IDs.this.dept, uu0)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Write out created shape files using rgdal library ----

synth.farms.outdir <- "D:/ENSO/Projects/CNH3/data/campos_sinteticos/"
tt0 <- dir.create(synth.farms.outdir, showWarnings = TRUE)

# --- Write out polygons with lat-lon coordinates
synth.farms.outdir.ll <- paste0(synth.farms.outdir, 'latlon_coords/')
tt0 <- dir.create(synth.farms.outdir.ll, showWarnings = TRUE)

oldPath = getwd()
setwd(synth.farms.outdir.ll)

outfile <- paste(dept.names.noblanks[iii],"_synth_farms", sep = "")

rgdal::writeOGR(farms.SPDF2, dsn=".", layer=outfile, driver="ESRI Shapefile",
         overwrite_layer=TRUE, verbose = TRUE)

setwd(oldPath)

# --- Now write out polygons with Gauss-Kruger coordinates

synth.farms.outdir.GK <- paste0(synth.farms.outdir,'GK_coords/')
tt0 <- dir.create(synth.farms.outdir.GK, showWarnings = TRUE)

oldPath = getwd()
setwd(synth.farms.outdir.GK)

outfile <- paste(dept.names.noblanks[iii],"_synth_farms_GK", sep = "")

rgdal::writeOGR(farms.SPDF, dsn=".", layer=outfile, driver="ESRI Shapefile",
    overwrite_layer=TRUE, verbose = TRUE)

setwd(oldPath)

rm(synth.farms.outdir, outfile); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Write out filled-in farm table to output text file ----

synth.farms.tables.outdir <- "D:/ENSO/Projects/CNH3/data/campos_sinteticos/tablas_datos/"
tt0 <- dir.create(synth.farms.tables.outdir, showWarnings = TRUE)

oldPath = getwd()
setwd(synth.farms.tables.outdir)

farm.outfile <- paste(dept.names.noblanks[iii],"_synth_farm_data.txt", sep = "")
write.table(farms.DF, file = farm.outfile,
            append = FALSE, sep = "\t", row.names=FALSE)

setwd(oldPath)
rm(synth.farms.tables.outdir, farm.outfile); gc()
# ------------------------------------------------------------------------------

}   # END of iteration over depts in Salado Basin


# -----------------------------------------------------------------------------#
# --- Write out master tile list to text file brefore ending this run ----

# --- The 'master.tile.list' object is written after a run is finished
# --- and is read before another run begins. In this way, a run that starts
# --- can knw if a tile has been used previously.

master.list.outfile <- paste(master.tile.outdir,"master_tile_data.txt", sep = "")
write.table(master.tile.list, file = master.list.outfile,
  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

rm(master.tile.outdir, master.list.outfile); gc()
# ------------------------------------------------------------------------------






# ------------------------------------------------------------------------------
# --- Clean up all objects ----

remove(list = ls()); gc()
# ------------------------------------------------------------------------------













depts.en.cuenca <- c("9 DE JULIO","ADOLFO ALSINA","BOLIVAR","BRAGADO",
  "CARLOS CASARES","CARLOS TEJEDOR","COLON",
  "DAIREAUX","FLORENTINO AMEGHINO","GENERAL ARENALES",
  "GENERAL PINTO","GENERAL VIAMONTE","GENERAL VILLEGAS","GUAMINI",
  "HIPOLITO YRIGOYEN","JUNIN","LEANDRO N. ALEM", "LINCOLN","PEHUAJO",
  "PELLEGRINI","RIVADAVIA","SALLIQUELO","TRENQUE LAUQUEN","TRES LOMAS",
  "GENERAL ROCA","PRESIDENTE ROQUE SAENZ PENA",
  "ATREUCO","CATRILO","CHAPALEUFU","CONHELO",
  "MARACO","QUEMU QUEMU","TRENEL",
  "GENERAL LOPEZ")




















ccc <- rgeos::gWithin(urban.areas.utm, cuenca.A.UTM, byid = TRUE)
cc2 <- which(ccc)

plot(cuenca.A.UTM, axes = TRUE)
plot(urban.areas.utm[cc2, ], add=TRUE, border = "tomato")

urban.areas.utm@data[cc2, ]

uuu <- (sapply(slot(urban.areas.utm,"polygons"),
               function(x) slot(x,"area"))) / 10000

ccc <- rgeos::gWithin(water.areas, cuenca.A.ll, byid = TRUE)
cc3 <- which(ccc)


plot(water.areas[cc3, ], add=TRUE, border = "steelblue")
plot(road.areas[ , ], add=TRUE, col = "orange")

ppp <- GSIF::getSpatialTiles(cuenca.A.ll, block.x = 0.5, block.y = 0.5, 
                             overlap.percent = 0, limit.bbox = TRUE,
                             return.SpatialPolygons = TRUE)





# --- If the farms do not fit inside the department, because there
# --- not enought available tiles, the number of synthetic farm
# --- needs to be reduced until they fit.
# --- First, we compute 95% of all tiles (not all tiles in a department
# --- are occupied by farms).
# --- Then, original list of synthetic farms is trimmed until the
# --- number of tiles needed reaches the number of available tiles.
# --- The trimming proceeds form the largest farm towards smaller farms.

tt1 <- floor(n.tiles.in.dept * 0.95)  # 95% of tiles in department
tt2 <- cumsum(sim.farms$tiles)        # Cumulative sum of tiles needed
tt3 <- which(tt2 <= tt1)        # List of farms (from largest) that fit inside department

# --- Now sim.farms will be trimmed to eliminate farms that
# --- do not fit inside department (insufficient tiles).
# --- Smallest farms are deleted first until the sum of areas reaches
# --- 95% of tiles inside department.

sim.farms <- sim.farms[tt3, ]   # Delete smallest farms

# --- Adjust these vectors after trimming farms

sim.farm.areas <- sim.farm.areas[tt3] # Adjusted farm size (multiples of point area)
sim.farm.IDs <- sim.farm.IDs[tt3]     # IDs of simulated farms
n.sim.farms <- nrow(sim.farms)    		# number of simulated farms

# --- Plot cumulative number of tiles and tiles in department

if (plot.results) {
  plot(tt2, type="l", col="steelblue",
    xlab="Number of farms",
    ylab="Cumulative number of tiles",
    main = dept.name.lc)
  abline(h = c(n.tiles.in.dept, n.tiles.in.dept * 0.95), col="grey60")
} 
rm(tt1,tt2,tt3)

} #End of check for number of needed tiles > tiles in department

remove(ww1,tiles.required, total.tiles.required, n.tiles.in.dept, uuu); gc()
# ------------------------------------------------------------------------------






# -----------------------------------------------------------------------------#
# --- Function to reduce number of synthetic farms if necessary ----

FilterFarms <- function(synth.farms, max.area.rate = 0.95) {
  
  filtered.synth.farms <- NULL
  
  FilterFarmsPerDept <- function(dept, synth.farms, max.area.rate) {
    dept.farms    <- synth.farms[synth.farms$dept.code == dept$dept.code, ]
    adjusted.area <- sum(dept.farms$adjusted.farm.size)
    
    while (adjusted.area > (max.area.rate * dept$area)) {
      row.number    <- sample(x = seq(from = 1, to = nrow(dept.farms)),
        size=1)
      adjusted.area <- sum(dept.farms$adjusted.farm.size)
      dept.farms    <- dept.farms[-c(row.number),]
    }
    
    if (is.null(filtered.synth.farms)) {
      filtered.synth.farms <<- dept.farms
    } else {
      filtered.synth.farms <<- rbind(filtered.synth.farms, dept.farms)
    }
  }
  
  depts.df  <- unique(synth.farms[c("dept.code", "area")])
  depts.lol <- apply(depts.df, MARGIN=1, FUN=as.list)
  lapply(depts.lol, FUN=FilterFarmsPerDept, synth.farms=synth.farms,
    max.area.rate=max.area.rate)
  return (filtered.synth.farms)
}
# ------------------------------------------------------------------------------


# -----------------------------------------------------------------------------#
# --- Reduce number of farms if their total area is > dept area ----

ddd <- dplyr::inner_join(farm.table, dplyr::tbl_df(basin.depts)) %>%
  dplyr::select(-c(dept.head:cent.lat, area.igm, perim)) %>%
  dplyr::mutate(area = area * 100) %>%
  dplyr::arrange(prov.code, dept.code, farm.id) %>%
  dplyr::arrange(prov.code, dept.code, owner.id, farm.id)

set.seed(128)
farm.table <- FilterFarms(ddd, max.area.rate = 0.95)

dd3 <- dplyr::tbl_df(ddd) %>%
  dplyr::group_by(dept.code, dept.name) %>%
  dplyr::summarise(total.orig = sum(adjusted.farm.size))

dd4 <- dplyr::tbl_df(farm.table) %>%
  dplyr::group_by(dept.code) %>%
  dplyr::summarise(total.adj = sum(adjusted.farm.size))

dplyr::inner_join(dd3, dd4) %>%
  dplyr::filter(total.orig != total.adj)

rm(ddd, dd2, dd3, dd4)
# ------------------------------------------------------------------------------

nacol <- function(spdf) {
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  nunique <- function(x){unique(x[!is.na(x)])}
  
  np = nrow(spdf)
  adjl = spdep::poly2nb(spdf)
  cols = rep(NA, np)
  cols[1]=1
  nextColour = 2
  
  for(k in 2:np){
    adjcolours = nunique(cols[adjl[[k]]])
    if(length(adjcolours)==0){
      cols[k]=resample(cols[!is.na(cols)],1)
    }else{
      avail = setdiff(nunique(cols), nunique(adjcolours))
      if(length(avail)==0){
        cols[k]=nextColour
        nextColour=nextColour+1
      }else{
        cols[k]=resample(avail,size=1)
      }
    }
  }
  return(cols)
}

