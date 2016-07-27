rm(list=ls()); invisible(gc());
# --- Load necessary R packages ----

# --- Define repository from which R packages will be downloaded.

options(repos=c(CRAN="http://mirror.fcaglp.unlp.edu.ar/CRAN/"), error=traceback)

list.of.packages <- c("ncdf4","MASS","fields","geoR","lubridate","sirad","dplyr", "foreach")  

for (pack in list.of.packages) {
    if (!require(pack, character.only = TRUE)) {
        install.packages(pack)
        require(pack, character.only = TRUE)
    }
}

rm(list.of.packages, pack)

# Stations climate.
climate <- read.table('1_setup/data/Salado.dat')

# station metadata
stns = read.table('1_setup/data/Salado_metadata.dat',header=T) %>% arrange(omm_id)
## Lon.Lat && distance matrix
lon.lat = stns[,3:4];

colnames(lon.lat) <- c("lon","lat")
rownames(lon.lat) <- stns$omm_id
# distance matrix of station locations
# converts longitude and latitude to KILOMETERS
dist.mat <- rdist.earth(lon.lat, miles=F)
colnames(dist.mat) <- stns$omm_id
rownames(dist.mat) <- stns$omm_id
# diagonal element is not explicitly equal to zero, so define as such
diag(dist.mat) <- 0

predlocs <- read.table("1_setup/data/grids.txt",header=T)
predloc.GK <- predlocs[,1:2]/1000 # convert from meters to KILOMETERS
predloc <- predlocs[,3:4]
rm(predlocs)


climate$date <- as.Date(climate$date)
climate <- climate %>% arrange(station, date)

unique_stations <- unique(climate$station)




if(require(doMC)) {
    registerDoMC(3)
}

models <- foreach(station=unique_stations) %dopar% {
    # system.time(replicate(500, {
    station_climate <- climate[climate$station == station,] %>%
        mutate(prcp_occ=(prcp > 0.1)+0,
               prcp_intensity=ifelse(prcp < 0.1, NA, prcp),
               prcp_occ_prev=lag(prcp_occ),
               prcp_prev=lag(prcp),
               tx_prev=lag(tx),
               tn_prev=lag(tn),
               year_fraction=2*pi*lubridate::yday(date) / ifelse(leap_year(date), 366, 365),
               ct=cos(year_fraction),
               st=sin(year_fraction),
               row_num=row_number())
    station_climate <- data.frame(station_climate, 
                                  Rt=seq(from=-1, to=1, length.out=nrow(station_climate)),
                                  ST1=season_rainfall_covariats[[1]],
                                  ST2=season_rainfall_covariats[[2]],
                                  ST3=season_rainfall_covariats[[3]],
                                  ST4=season_rainfall_covariats[[4]],
                                  SMX1=season_max_temp_covariats[[1]],
                                  SMX2=season_max_temp_covariats[[2]],
                                  SMX3=season_max_temp_covariats[[3]],
                                  SMX4=season_max_temp_covariats[[4]],
                                  SMN1=season_min_temp_covariats[[1]],
                                  SMN2=season_min_temp_covariats[[2]],
                                  SMN3=season_min_temp_covariats[[3]],
                                  SMN4=season_min_temp_covariats[[4]])
    
    prcp_covariats <- c('ct', 'st', 'ST1', 'ST2', 'ST3', 'ST4')
    # Fit model for precipitation occurrence.
    probit_indexes <- na.omit(station_climate[, c('prcp_occ', 'row_num', 'prcp_occ_prev', prcp_covariats)])$row_num
    
    occ_fit <- glm(formula(paste0('prcp_occ', '~', 'prcp_occ_prev +', paste0(prcp_covariats, collapse='+'))),
                   data=station_climate[probit_indexes, ],
                   family=binomial(probit))
    
    coefocc <- occ_fit$coefficients
    station_climate[probit_indexes, 'probit_residuals'] <- occ_fit$residuals
    
    # Fit model for precipitation amounts.
    gamma_indexes <- na.omit(station_climate[, c('prcp_intensity', 'row_num', prcp_covariats)])$row_num
    # save model in list element "i"
    prcp_fit <- glm(formula(paste0('prcp_intensity', '~', paste0(prcp_covariats, collapse='+'))),
                    data=station_climate[gamma_indexes, ],
                    family=Gamma(link=log))
    
    # coefamt <- prcp_fit$coefficients
    # 
    temps_covariats <- c('tn_prev', 'tx_prev', 'ct', 'st', 'prcp_occ', 'Rt', 'SMN1', 'SMN2', 'SMN3', 'SMN4', 'SMX1', 'SMX2', 'SMX3', 'SMX4')
    tx_indexes <- na.omit(station_climate[, c('tx', 'row_num', temps_covariats)])$row_num
    tn_indexes <- na.omit(station_climate[, c('tn', 'row_num', temps_covariats)])$row_num
    # Fit 
    
    tx_fit <- lm(formula(paste0('tx', '~', paste0(temps_covariats, collapse='+'))),
                 data = station_climate[tx_indexes, ])
    coefmax <- tx_fit$coefficients
    station_climate[tx_indexes, 'tx_residuals'] <- tx_fit$residuals
    
    tn_fit <- lm(formula(paste0('tn', '~', paste0(temps_covariats, collapse='+'))),
                 data = station_climate[tn_indexes, ])
    coefmin <- tn_fit$coefficients
    station_climate[tn_indexes, 'tn_residuals'] <- tn_fit$residuals
    
    return(list(coefficients=list('coefocc'=coefocc,
                                  'coefmin'=coefmin,
                                  'coefmax'=coefmax),
                gamma=list(coef=prcp_fit$coef, alpha=gamma.shape(prcp_fit)$alpha),
                residuals=station_climate[, c('date', 'station', 'probit_residuals', 'tx_residuals', 'tn_residuals')],
                station=station))
    
    # }))
}

first_model <- models[[1]]

# Unwind results.
models_coefficients <- list(
    coefocc=matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefocc),
                   dimnames=list(rownames(dist.mat), names(first_model$coefficients$coefocc))),
    
    coefmin=matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefmin),
                   dimnames=list(rownames(dist.mat), names(first_model$coefficients$coefmin))),
    
    coefmax=matrix(NA, nrow = length(models), ncol = length(first_model$coefficients$coefmax),
                   dimnames=list(rownames(dist.mat), names(first_model$coefficients$coefmax)))
)

models_residuals <- data.frame()

GAMMA <- as.list(rep(NA, times=length(models)))
names(GAMMA) <- rownames(dist.mat)

for(m in models) {
    station <- as.character(m$station)
    models_coefficients$coefocc[station, ] <- m$coefficients$coefocc
    models_coefficients$coefmin[station, ] <- m$coefficients$coefmin
    models_coefficients$coefmax[station, ] <- m$coefficients$coefmax
    GAMMA[[station]] <- m$gamma
    models_residuals <- rbind(models_residuals, m$residuals)
}

rm(first_model, models, m);


# estimate model coefficients on grid using ordinary kriging (OK)
# occurrence
coefocc.sim <- matrix(NA, nrow = nrow(predloc),
                     ncol = ncol(models_coefficients$coefocc),
                     dimnames = list(NULL, colnames(models_coefficients$coefocc)))
for(covariat in 1:ncol(coefocc.sim)){
    coefocc.sim[, covariat] = suppressWarnings(predict(Krig(lon.lat, models_coefficients$coefocc[, covariat]), predloc))
}
# minimum temperature
coefmin.sim <- matrix(NA, nrow = nrow(predloc),
                      ncol = ncol(models_coefficients$coefmin),
                      dimnames = list(NULL, colnames(models_coefficients$coefmin)))
for(covariat in 1:ncol(coefmin.sim)){
    coefmin.sim[, covariat] = suppressWarnings(predict(Krig(lon.lat, models_coefficients$coefmin[, covariat]), predloc))
}
# maximum temperature
coefmax.sim <- matrix(NA, nrow = nrow(predloc),
                      ncol = ncol(models_coefficients$coefmax),
                      dimnames = list(NULL, colnames(models_coefficients$coefmax)))
for(covariat in 1:ncol(coefmax.sim)){
    coefmax.sim[, covariat] = suppressWarnings(predict(Krig(lon.lat, models_coefficients$coefmax[, covariat]), predloc))
}

# rm(models); invisible(gc());

partially_apply_LS <- function(variogram, base_p=c()) {
    ## least squares function for kriging parameters
    LS <- function(p){
        p <- c(base_p, p)
        M <- p[2]*exp((-dist.mat)/p[3])
        diag(M) <- p[1]+p[2]
        return(sum(variogram-M)^2)
    }
    return(LS)
}

month_params <- foreach(m=1:12) %dopar% {
    month_residuals <- models_residuals %>% filter(month(date) == m) %>% arrange(station, date)
    month_climate <- climate %>% filter(month(date) == m)
    
    n_stations <- length(unique_stations)
    
    tx_res_matrix <- matrix(month_residuals$tx_residuals, ncol=n_stations)
    tn_res_matrix <- matrix(month_residuals$tn_residuals, ncol=n_stations)
    probit_res_matrix <- matrix(month_residuals$probit_residuals, ncol=n_stations)
    
    # tx_matrix <- matrix(month_climate$tx, ncol=n_stations)
    # tn_matrix <- matrix()
    
    colnames(tx_res_matrix) <- colnames(tn_res_matrix) <- colnames(probit_res_matrix) <- unique_stations
    
    
    # all(na.omit(tx_matrix[, '87448']) == na.omit(month_climate[month_climate$station == 87448, 'tx']))
    # all(na.omit(tx_res_matrix[,1]) == na.omit(month_residuals[month_residuals$station == 87448, 'tx_residuals']))
    
    prcp_cor <- cor(probit_res_matrix, use="pairwise.complete")
    
    PRCPvario <- var(probit_res_matrix, use="pairwise.complete") * (1 - prcp_cor)
    TMAXvario <- cov(tx_res_matrix, use="pairwise.complete")
    TMINvario <- cov(tn_res_matrix, use="pairwise.complete")
    
    # Assert that the column and row names of the variograms equal the ones of the distance matrix.
    stopifnot(all(colnames(TMAXvario) == colnames(dist.mat)), all(rownames(TMAXvario) == rownames(dist.mat)))
    
    params <- optimize(interval=c(0, max(dist.mat)), f=partially_apply_LS(PRCPvario, base_p=c(0, 1)))$minimum
    params <- c(0, 1, params)
    
    sill_initial_value <- mean(TMAXvario[upper.tri(TMAXvario, diag=T)])
    params.max <- optim(par=c(sill_initial_value, max(dist.mat)), 
                         fn=partially_apply_LS(TMAXvario, base_p = c(0)))$par
    params.max <- c(0, params.max)
    
    sill_initial_value <- mean(TMINvario[upper.tri(TMINvario, diag=T)])
    params.min <- optim(par=c(sill_initial_value, max(dist.mat)),
                         fn=partially_apply_LS(TMINvario, base_p = c(0)))$par
    params.min <- c(0, params.min)
    
    return(list(
        prcp=params,
        tx=params.max,
        tn=params.min
    ))
}

rm(models_residuals)

##########################################
##
## Simulations
##
##########################################
start_date <- '2015-10-01'
end_date <- '2015-12-31'
# number of realizations 
NT=2#00

## if you wish to neglect temperature trends, and be centered on the mean of 1961-2013 
## temperatures, set Rt.sim = 0 for all days

#Rt.sim <- rep(1,length(yr.sim)) # if you want simulated temperatures to be high
#Rt.sim <- rep(-1,length(yr.sim)) # if you want simulated temperatures to be low
# Rt.sim <- rep(0,length(yr.sim))

simulation_dates <- data.frame(date=seq.Date(from=as.Date(start_date), to=as.Date(end_date), by = 'days')) %>%
    mutate(year_fraction=2*pi*lubridate::yday(date) / ifelse(leap_year(date), 366, 365),
           ct=cos(year_fraction),
           st=sin(year_fraction),
           Rt=seq(from=-1, to=1, length.out=n()),
           season=ceiling(month(date)/3),
           ST1=0, ST2=0, ST3=0, ST4=0,
           SMX1=0, SMX2=0, SMX3=0, SMX4=0,
           SMN1=0, SMN2=0, SMN3=0, SMN4=0)


# simulation metadata
x_coords = unique(predloc.GK[,1]) * 1000 # convert back to meters (original form)
y_coords = unique(predloc.GK[,2]) * 1000 # convert back to meters (original form)
n.x_coords <- length(x_coords)
n.y_coords <- length(y_coords)
# elevation at each grid cell 
# elev <- as.matrix(read.table("../1_setup/data/Salado-Abasin-elevation-meters.dat",header=F))



# number of days to simulate
nt.sim <- nrow(simulation_dates)

## years to simulate (just for convenience, they are in fact arbitrary years)
uyr.sim <- unique(yr.sim)
## number of years to simulate (also arbitrary)
nyr.sim <- length(uyr.sim)

# julian days (for .nc file)
# Times <- julian(sim.dates, origin = as.Date("1961-01-01"))
# n.times <- length(Times)    # Number of simulated days

RNum <- 1:NT
n.realizations <- length(RNum)
# id for missing data
mv <- -9999

# bootstrapping OND (i.e., ST4, SMN4, SMX4) values to condition output
# OND 2015 FORECAST FOR PRECIPITATION: 70:20:10
# OND 2015 FORECAST FOR TEMPERATURE:   40:35:25
# see corresponding gif files for IRI forecasts

# sample NT different OND precip totals with IRI probability (B, N, A) = (0.1, 0.2, 0.7)
# st.samp = sample(1:3, NT, prob=c(0.1, 0.2, 0.70), replace=TRUE)

# ordered_rainfall <- d_seatot %>% dplyr::group_by(season) %>% 
#     dplyr::do(data.frame(q = seq(1, length(.$seatot)) %/% ((length(.$seatot)+1)/3),
#                          seatot = sort(.$seatot),
#                          seamax = sort(.$seamax),
#                          seamin = sort(.$seamin)))
# ordered_rainfall <- data.frame(ordered_rainfall)



# get_seasonal_covariat <- function(season_number, season_values) {
#     return(mean(season_values))
# }

get_seasonal_covariat <- function(season_number, season_values, break_quantiles=c(0, 0.33, 0.66, 1), levels_probabilities=c(1/3, 1/3, 1/3)) {
    splitted <- cut(season_values, breaks=quantile(season_values, probs=break_quantiles), include.lowest = T)
    values_level <- sample(levels(splitted), 1, prob=levels_probabilities)
    return(sample(season_values[splitted == values_level], 1))
}

get_rainfall_seasonal_covariat <- get_seasonal_covariat
get_temperatures_seasonal_covariat <- function(season_number, season_tx_values, season_tn_values, break_quantiles=c(0, 0.33, 0.66, 1), levels_probabilities=c(0.25, 0.35, 0.40)) {
    splitted_tx <- cut(season_tx_values, breaks=quantile(season_tx_values, probs=break_quantiles), include.lowest = T)
    splitted_tn <- cut(season_tn_values, breaks=quantile(season_tn_values, probs=break_quantiles), include.lowest = T)
    
    values_level_index <- sample(seq_along(levels(splitted_tx)), 1, prob=levels_probabilities)
    tx_sampled_level <- levels(splitted_tx)[values_level_index]
    tn_sampled_level <- levels(splitted_tn)[values_level_index]
    return(list(tx=sample(season_tx_values[splitted_tx == tx_sampled_level], 1),
             tn=sample(season_tn_values[splitted_tn == tn_sampled_level], 1))
    )
}


# Gamma shape and scale
SH <- c()

for(station in as.character(unique_stations)) SH[station] <- GAMMA[[station]]$alpha
SH.sim <- suppressWarnings(predict(Krig(lon.lat, SH), predloc))


# define arrays for simulated series
# SIMamt.sim <- SIMocc.sim <- SIMmax.sim <- SIMmin.sim <- array(dim=c(nrow(simulation_dates), nrow(predloc), NT))

invisible(gc());
# TODO: remove this.
d <- 1; i <- 1  # For debugging purposes.

# progress bar
pb <- txtProgressBar(min = 0, max = NT, style = 3)
## occurrences
gen_climate <- foreach(i=1:NT, .combine=c) %dopar% {
    simulation_start <- as.Date(start_date) - 1
    
    # initialize with climatology of the previous day
    w2 <- suppressWarnings(grf(nrow(predloc), grid=predloc, cov.model="exponential",
                               cov.pars=month_params[[month(simulation_start)]]$prcp[c(2,3)],
                               nugget=month_params[[month(simulation_start)]]$prcp[1], mean=rep(0, nrow(predloc)), messages=FALSE))
    SIMocc.old <- (w2$data > 0) + 0
    
    start_climatology <- climate %>% filter(month(date) == month(simulation_start), day(date) == day(simulation_start)) %>%
                            group_by(station) %>% summarise(tx = mean(tx, na.rm=T), tn = mean(tn, na.rm=T))
    
    SIMmin.old <- suppressWarnings(as.vector(predict(Krig(lon.lat, start_climatology$tn), predloc)))
    SIMmax.old <- suppressWarnings(as.vector(predict(Krig(lon.lat, start_climatology$tx), predloc)))
    
    SC <- array(data=NA, dim=c(nrow(simulation_dates), nrow(lon.lat)))
    colnames(SC) <- unique_stations
    
    for(season_number in unique(simulation_dates$season)) {
        season_indexes <- simulation_dates$season == season_number
        season_values <- d_seatot[d_seatot$season == season_number, ]
        temp_covariats <- get_temperatures_seasonal_covariat(season_indexes, season_values$seamax, season_values$seamin)
        
        simulation_dates[season_indexes, paste0('ST', season_number)] <- get_seasonal_covariat(season_number, season_values$seatot)
        simulation_dates[season_indexes, paste0('SMX', season_number)] <- temp_covariats$tx
        simulation_dates[season_indexes, paste0('SMN', season_number)] <- temp_covariats$tn
    }
    
    daily_covariats <- as.matrix(t(simulation_dates[, c(3:5, 7:18)]))
    
    daily_covariats <- rbind(1, daily_covariats)
    rownames(daily_covariats)[1] <- '(Intercept)'
    
    for(station in as.character(unique_stations)){
        gamma_coef <- GAMMA[[station]]$coef
        SC[, station] <- exp(apply(daily_covariats[names(gamma_coef), ] * gamma_coef, FUN=sum, MAR=2, na.rm=T)) / SH[station]
    }
    
    simulated_occurrence <- simulated_tx <- simulated_tn <- simulated_prcp <- array(dim=c(nrow(simulation_dates), nrow(predloc)))
    
    temps_retries <- 0
    for(d in 1:nrow(simulation_dates)){
        p <- month_params[[month(simulation_dates[d, 'date'])]]
        
        simulation_matrix <- cbind(SIMocc.old, SIMmin.old, SIMmax.old, matrix(daily_covariats[, d], ncol=nrow(daily_covariats), nrow=length(SIMmax.old), byrow = T))
        colnames(simulation_matrix) <- c('prcp_occ_prev', 'tn_prev', 'tx_prev', rownames(daily_covariats))
        
        
        mu.occ <- apply(simulation_matrix[, colnames(coefocc.sim)]*coefocc.sim, 1, sum, na.rm=T)
        w.occ  <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
                                       cov.pars=c(p$prcp[2],p$prcp[3]),nugget=p$prcp[1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
        
        simulated_occurrence[d,] <- ((mu.occ+w.occ) > 0) + 0
        
        simulation_matrix <- cbind(prcp_occ=simulated_occurrence[d,], simulation_matrix)
        
        mu.min <- apply(simulation_matrix[, colnames(coefmin.sim)]*coefmin.sim, 1, sum, na.rm=T)
        mu.max <- apply(simulation_matrix[, colnames(coefmax.sim)]*coefmax.sim, 1, sum, na.rm=T)
        
        w.min <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
                                      cov.pars=c(p$tn[2],p$tn[3]),nugget=p$tn[1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
        
        w.max <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
                                      cov.pars=c(p$tx[2],p$tx[3]),nugget=p$tx[1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
        
        simulated_tn[d,] <- signif(mu.min+w.min,digits=4)
        simulated_tx[d,] <- signif(mu.max+w.max,digits=4)
        
        # 1/10000
        while(min(simulated_tx[d,] - simulated_tn[d,]) < 0.5) {
            temps_retries <- temps_retries + 1
            
            w.min <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
                                          cov.pars=c(p$tn[2],p$tn[3]),nugget=p$tn[1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
            
            
            w.max <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
                                          cov.pars=c(p$tx[2],p$tx[3]),nugget=p$tx[1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
            
            simulated_tn[d,] <- signif(mu.min+w.min,digits=4)
            simulated_tx[d,] <- signif(mu.max+w.max,digits=4)
        }
        
        SIMocc.old <- simulated_occurrence[d,]
        SIMmax.old <- simulated_tx[d,]
        SIMmin.old <- simulated_tn[d,]
        
        ## amounts
        SC.sim <- suppressWarnings(predict(Krig(lon.lat,SC[d,]), predloc))
        
        w3 <- suppressWarnings(grf(nrow(predloc.GK),grid=predloc.GK,cov.model="exponential",
                                   cov.pars=c(p$prcp[2],p$prcp[3]),nugget=p$prcp[1],mean=rep(0,nrow(predloc.GK)),messages=FALSE)$data)
        
        simulated_prcp[d,] <- signif(qgamma(pnorm(w3), shape=SH.sim, scale=SC.sim), digits=4)
        cat(paste0('\r', NT, ': ', d,'/', ncol(daily_covariats), '. Retries: ', temps_retries))
    }
    
    simulated_prcp[simulated_prcp < 0.1] <- 0
    simulated_occurrence[simulated_prcp == 0] <- 0
    simulated_prcp[simulated_occurrence == 0] <- 0
    
    return(list(
        prcp=simulated_prcp,
        tx=simulated_tx,
        tn=simulated_tn
    ))
}

# estimate potential evapotranspiration
# set up for calculating reference evapotranspiration (ET)
day.of.year <- yday(as.Date(paste(yr.sim,mo.sim,da.sim,sep="-")))
q0 = array(data=NA,dim=dim(SIMamt.sim))
predrad <- numeric(nrow(predloc))

for(kk in 1:nrow(predloc)){ 
    predrad[kk] <- radians(predloc[kk,2])
    q0[,kk,1] <- extrat(i=day.of.year, lat=predrad[kk])$"ExtraTerrestrialSolarRadiationDaily"
}
if(NT > 1) {
    for(i in 2:NT) q0[,,i] <- q0[,,1]
}

# set up for calculating reference evapotranspiration (ET)
coef.a <- 0.001703
coef.b <- 21.967919
coef.c <- 0.083444
coef.d <- 0.541066
# set up for calculating solar radiation (SRAD)
bc.coef=c(A=0.69, B=0.02, C=2.12)

# difference between max and min temperatures
dtr <- SIMmax.sim-SIMmin.sim
# we included a check in the wgen code to ensure max temp is always greater than min temp... so the next line is not necessary
#dtr[dtr<0] <- NA
# average temperature
tavg = (SIMmax.sim+SIMmin.sim)/2
# equation for reference ET
SIMet0.sim <- coef.a * 0.408 * q0 * (tavg + coef.b) * (dtr - (coef.c * SIMamt.sim)) ** coef.d
SIMet0.sim[is.na(SIMet0.sim)]=0 # R cant handle imaginary numbers (i.e. exponentiating negative bases)


# Estimating solar radiation.
# Computate Bristow-Campbell's delta temperature.
# do not consider tomorrow's minimum temperature 
# because this is hourly aggregate Met Service data
# dtr is dtemp
extraT <- array(suppressWarnings(extrat(i=matrix(dayOfYear(matrix(rep(sim.dates,nrow(predloc.GK)),nrow=nt.sim,ncol=nrow(predloc.GK),byrow=F)),nrow=nt.sim,ncol=nrow(predloc.GK),byrow=F), lat=predrad)$ExtraTerrestrialSolarRadiationDaily),dim=dim(dtr))
SIMsrad.sim <- extraT * bc.coef[['A']] * (1 - exp(-bc.coef[['B']] * (dtr^bc.coef[['C']])))


# define dimensions for .nc file 
dimRNum <- ncdim_def(name="rnum", units="number", vals=RNum, unlim=TRUE, create_dimvar=TRUE, longname="Realization number")
dimX <- ncdim_def(name="x_coord", units="meters", vals=x_coords, unlim=FALSE, create_dimvar=TRUE, longname="longitude in Gauss-Krueger coordinates")
dimY <- ncdim_def(name="y_coord", units="meters", vals=y_coords, unlim=FALSE, create_dimvar=TRUE, longname="latitude in Gauss-Krueger coordinates")
dimTime <- ncdim_def(name="time", units="days since 1961-01-01", vals=Times, unlim = FALSE, create_dimvar=TRUE, calendar="standard", longname="Time in days since 1961-01-01")

tmax <- ncvar_def(name="tmax", units="degress_Celsius", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Daily maximum near-surface air temperature", prec="float", verbose=TRUE)
tmin <- ncvar_def(name="tmin", units="degress_Celsius", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Daily minimum near-surface air temperature", prec="float", verbose=TRUE)
prcp <- ncvar_def(name="prcp", units="mm day-1", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Daily total rainfall", prec="float", verbose=TRUE)
srad <- ncvar_def(name="srad", units="Mjoules m-2 day-1", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Daily solar radiation", prec="float", verbose=TRUE)
et0  <- ncvar_def(name="et0",  units="mm day-1", dim=list(dimRNum, dimTime, dimY, dimX), missval=mv, longname="Reference evapotranspiration", prec="float", verbose=TRUE)
dec_lon <- ncvar_def(name="dec_lon", units="degrees_east", missval=mv, dim=list(dimY, dimX), longname="Longitude in decimal degrees", prec="float", verbose=TRUE)
dec_lat <- ncvar_def(name="dec_lat", units="degrees_north", missval=mv, dim=list(dimY, dimX), longname="Latitude in decimal degrees", prec="float", verbose=TRUE)
elevation <- ncvar_def(name="elevation", units="meters from sea level", missval=mv, dim=list(dimY, dimX), longname="Elevation in meters from sea level", prec="float", verbose=TRUE)

STs <- ncvar_def(name="ST", units="mm season-1", dim=list(dimRNum), missval=mv, longname="Seasonal total precipitation covariates", prec="float", verbose=TRUE)
SMXs <- ncvar_def(name="SMX", units="degrees_Celsius", dim=list(dimRNum), missval=mv, longname="Seasonal average maximum temperature covariates", prec="float", verbose=TRUE)
SMNs <- ncvar_def(name="SMN", units="degrees_Celsius", dim=list(dimRNum), missval=mv, longname="Seasonal average minimum temperature covariates", prec="float", verbose=TRUE)

# Create NetCDF file in current working directory
dd <- nc_create(filename="wgen-conditional-seasonal-gridded.nc",
                vars = list(tmax, tmin, prcp, srad, et0, dec_lon, dec_lat, elevation, STs, SMXs, SMNs),
                force_v4 = TRUE, verbose = TRUE)

# define arrays for simulated series
# array(dim=c(NT,nt.sim,n.y_coords,n.x_coords))
# progress bar
pb <- txtProgressBar(min = 0, max = NT, style = 3)
for(i in 1:NT){
    for(d in 1:nt.sim){
        ncvar_put(dd, varid = "tmax", 
                  vals = matrix(SIMmax.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
                  start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords),
                  verbose = TRUE)
        
        ncvar_put(dd, varid = "tmin", 
                  vals = matrix(SIMmin.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
                  start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords),
                  verbose = TRUE)
        
        ncvar_put(dd, varid = "prcp", 
                  vals = matrix(SIMamt.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
                  start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords),
                  verbose = TRUE)
        
        ncvar_put(dd, varid = "et0", 
                  vals = matrix(SIMet0.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
                  start = c(i, d, 1, 1), count = c(1, 1, n.y_coords, n.x_coords),
                  verbose = TRUE)
        
        ncvar_put(dd, varid = "srad",
                  vals = matrix(SIMsrad.sim[d,,i],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
                  start = c(i, d, 1,1 ), count = c(1, 1, n.y_coords, n.x_coords),
                  verbose = TRUE)
    }
    setTxtProgressBar(pb, i)
}

ncvar_put(dd, varid = "dec_lon",
          vals = matrix(predlocs[,"lon"],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
          start = c(1, 1), count = c(n.y_coords, n.x_coords),
          verbose = TRUE)

ncvar_put(dd, varid = "dec_lat",
          vals = matrix(predlocs[,"lat"],nrow=n.y_coords,ncol=n.x_coords,byrow=T),
          start = c(1, 1), count = c(n.y_coords, n.x_coords),
          verbose = TRUE)

ncvar_put(dd, varid = "elevation",
          vals = as.matrix(elev),
          start = c(1, 1), count = c(n.y_coords, n.x_coords),
          verbose = TRUE)

ncvar_put(dd, varid = "ST",
          vals = apply(ST4.sim, 1, unique),
          start = 1, count = NT,
          verbose = TRUE)

ncvar_put(dd, varid = "SMX", 
          vals = apply(SMX4.sim, 1, unique),
          start = 1, count = NT,
          verbose = TRUE)

ncvar_put(dd, varid = "SMN",
          vals = apply(SMN4.sim, 1, unique),
          start = 1, count = NT, 
          verbose = TRUE)

# --- Write global attributes

ncatt_put( dd, varid = 0,
           attname = "title",
           attval = "Synthetic weather data for the Salado River Basin A (Argentina)",
           verbose = TRUE)

ncatt_put( dd, varid = 0,
           attname = "software",
           attval = "Stochastic weather generator version 3.0",
           verbose = TRUE)

ncatt_put( dd, varid = 0,
           attname = "climate driver covariates",
           attval = "Regionally-averaged seasonal total precipitation, average maximum temperature, and average minimum temperature",
           verbose = TRUE)

ncatt_put(dd, varid = 0,
          attname = "Start and end dates",
          attval = paste(sim.start.date, sim.end.date),
          verbose = TRUE)

ncatt_put(dd, varid = "et0",
          attname = "calculation",
          attval = "Hargreaves-Samani modified by Droogers and Allen",
          verbose = TRUE)

ncatt_put(dd, varid = "et0",
          attname = "coefficients",
          attval = "Coefficients for daily data calibrated for Junín",
          verbose = TRUE)  

ncatt_put(dd, varid = "srad",
          attname = "calculation",
          attval = "Bristow-Campbell, not considering tomorrow's min temp because obs data derived from hourly aggregate Met Service data",
          verbose = TRUE)

ncatt_put(dd, varid = "srad",
          attname = "coefficients",
          attval = "Coefficients estimated for Buenos Aires and Pilar",
          verbose = TRUE)

nc_close(dd)
