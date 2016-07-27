  # --- Load necessary R packages ----

# --- Define repository from which R packages will be downloaded.

options(repos=c(CRAN="http://lib.stat.cmu.edu/R/CRAN/"), error=traceback)

# only need ismev package for computing GEV parameters
# and fields package for plotting the spatial maps
# and rgdal to read in the shapefile for A basin
list.of.packages <- c("geoR", "fields", "lubridate", "MASS")  

for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    install.packages(pack)
    require(pack, character.only = TRUE)
  }
}

rm(list.of.packages, pack)

## set the working directory -- this is wherever you put the folder from the zip file
setwd("/Users/andrew/Desktop/workshop-BA/5_metamodel")

load("metamodel-setup.RData")

#n.pred.et <- 11 # intercept, prcp, tmax, tmin, et0, %maize, %girasol, %trigo, %soja1, %trigo/soja2, %pasturas
n.pred.et <- 10 # do not include TRIGOSOJA2 (linearly dependent variable)

z.coefs.et <- array(data=NA,dim=c(n.pred.et,  dim(napa)[2:3], 12))
z.resid.et <- array(data=NA,dim=c(nrow(et.m), dim(napa)[2:3]))

xy.id <- matrix(NA,nrow=700,ncol=2)
cnt = 0

# 
lon.GKs <- lat.GKs <- c()
for(xx in 1:nx){
	for(yy in 1:ny){
		Y.et <- et.m[,xx,yy]

		if((1:nx)[xx] %in% xx.id & (1:ny)[yy] %in% yy.id){
			if(sum(is.na(Y.et)) == (length(Y.et))){
				next 
			}else{
				cnt = cnt + 1
				xy.id[cnt,] <- c(xx,yy)

				for(mo.kk in 1:12){
					Y.et.mo <- Y.et[mo.month==mo.kk]
					#Y.et.mo <- log(Y.et[mo.month==mo.kk])
					X.et <- as.matrix(cbind(prcp.m[mo.month==mo.kk,xx,yy], tmax.m[mo.month==mo.kk,xx,yy], tmin.m[mo.month==mo.kk,xx,yy], et0.m[mo.month==mo.kk,xx,yy], girasol.m[mo.month==mo.kk,xx,yy], maize.m[mo.month==mo.kk,xx,yy], 
							pasturas.m[mo.month==mo.kk,xx,yy], soja1.m[mo.month==mo.kk,xx,yy], trigo.m[mo.month==mo.kk,xx,yy]))#, trigosoja2.m[,xx,yy]))
					#X.et <- as.matrix(cbind(prcp.m[,xx,yy], tmax.m[,xx,yy], tmin.m[,xx,yy], et0.m[,xx,yy], girasol.m[,xx,yy], maize.m[,xx,yy], 
					#		pasturas.m[,xx,yy], soja1.m[,xx,yy], trigo.m[,xx,yy]))#, trigosoja2.m[,xx,yy]))
						
				colnames(X.et) <- c("Prcp", "Tmax", "Tmin", "ET0", "Girasol", "Maize", "Pasturas", "Soja1", "Trigo")#, "TrigoSoja2")

				df.et <- data.frame(cbind(Y.et.mo,X.et))

				zz.et <- lm(Y.et.mo ~ ., data=df.et)
				z.coefs.et[,xx,yy,mo.kk] <- zz.et$coef
				z.resid.et[mo.month==mo.kk,xx,yy] <- zz.et$resid 
				}
				
				lon.GKs <- c(lon.GKs, lon.GK[yy])
				lat.GKs <- c(lat.GKs, lat.GK[xx])
			}
		}else{
			next 
		}
	}
}

z.mean.et <- mean(z.resid.et,na.rm=T)
z.sd.et <- sd(z.resid.et,na.rm=T)

lon.lat.knots <- cbind(as.vector(t(lon[xx.id,yy.id])),as.vector(t(lat[xx.id,yy.id])))
z.coefs.id <- !is.na(as.vector(t(z.coefs.et[1,xx.id,yy.id,1])))
lon.lat.knots <- lon.lat.knots[z.coefs.id,]

lon.lat.knots.GK <- cbind(lon.GKs, lat.GKs)

lon.lat.preds.GK <- expand.grid(lon.GK,lat.GK)#[!is.na(as.vector(t(basin.mask))),]
lon.lat.preds.LL <- cbind(as.vector(lon),as.vector(lat))

n.knots <- nrow(lon.lat.knots.GK)

z.resid.mat.et <- matrix(NA,nrow=(tt.month),ncol=n.knots)
for(k in 1:(tt.month)) z.resid.mat.et[k,] <- as.vector(t(z.resid.et[k,,]))[!is.na(as.vector(t(z.resid.et[k,,])))]

dist.mat <- rdist(lon.lat.knots)
diag(dist.mat) <- 0

##
## least squares function for kriging parameters
## with no nugget 
LS <- function(p){
M <- p[1]*exp((-dist.mat)/p[2])
diag(M) <- p[1]
return(sum(vario-M)^2)
}

params.et <- matrix(0,nrow=12,ncol=3)
for(k in 1:12){
	NAPAcov.et <- cor(z.resid.mat.et[mo.month==k,],use="complete")
	NAPAvario.et <- var(z.resid.mat.et[mo.month==k,])*(1 - NAPAcov.et)
vario <- NAPAvario.et
params.et[k,2:3] <- optim(par=c(max(vario),max(dist.mat)),fn=LS)$par 
}
params.et[params.et<0]=0 

z.coefs.pred.et <- array(data=NA,dim=c(nrow(z.coefs.et),dim(napa)[2:3], 12))

for(mo.kk in 1:12){
for(k in 1:nrow(z.coefs.et)){
	z.coefs.et.tmp <- as.vector(t(z.coefs.et[k,,,mo.kk]))[!is.na(as.vector(t(z.coefs.et[k,,,mo.kk])))]
	z.coefs.pred.et[k,,,mo.kk] <- matrix(predict(Krig(lon.lat.knots, z.coefs.et.tmp,m=1),lon.lat.preds.LL),nrow=70,ncol=74)
}
}

et.m.pred <- array(data=NA,dim=c(tt.month,dim(napa)[2:3]))
for(k in (1:tt.month)){ 
	et.m.pred[k,,] <- 
		z.coefs.pred.et[1,,,mo.month[k]]*matrix(1,nrow=74,ncol=70) + 
		z.coefs.pred.et[2,,,mo.month[k]]*prcp.m[k,,] + 
		z.coefs.pred.et[3,,,mo.month[k]]*tmax.m[k,,] + 
		z.coefs.pred.et[4,,,mo.month[k]]*tmin.m[k,,] + 
		z.coefs.pred.et[5,,,mo.month[k]]*et0.m[k,,] + 
		z.coefs.pred.et[6,,,mo.month[k]]*girasol.m[k,,] + 
		z.coefs.pred.et[7,,,mo.month[k]]*maize.m[k,,] + 
		z.coefs.pred.et[8,,,mo.month[k]]*pasturas.m[k,,] + 
		z.coefs.pred.et[9,,,mo.month[k]]*soja1.m[k,,] + 
		z.coefs.pred.et[10,,,mo.month[k]]*trigo.m[k,,]
}
et.m.pred[et.m.pred < min(et.m,na.rm=T)] <- min(et.m,na.rm=T)
et.m.pred[et.m.pred > max(et.m,na.rm=T)] <- max(et.m,na.rm=T)

# now the napa part of it all 
n.pred.gw <- 6 # intercept, lag1 napa, prcp, tmax, tmin, et

z.coefs.gw <- array(data=NA,dim=c(n.pred.gw, dim(napa)[2:3]))
z.resid.gw <- array(data=NA,dim=c((nrow(napa.m.s)-1), dim(napa)[2:3]))

lon.GKs <- lat.GKs <- c()
for(xx in 1:nx){
	for(yy in 1:ny){
		Y.gw <- napa.m.s[2:tt.month,xx,yy]

		if((1:nx)[xx] %in% xx.id & (1:ny)[yy] %in% yy.id){
			if(sum(is.na(Y.gw)) == (length(Y.gw))){
				next 
			}else{
				#X.gw <- as.matrix(cbind(napa.m.s[1:(tt.month-1),xx,yy], prcp.m[1:(tt.month-1),xx,yy], tmax.m[1:(tt.month-1),xx,yy], tmin.m[1:(tt.month-1),xx,yy], et.m[1:(tt.month-1),xx,yy]))
				X.gw <- as.matrix(cbind(napa.m.s[1:(tt.month-1),xx,yy], prcp.m[2:tt.month,xx,yy], tmax.m[2:tt.month,xx,yy], tmin.m[2:tt.month,xx,yy], et.m.pred[2:tt.month,xx,yy]))
				colnames(X.gw) <- c("Napa","Prcp","Tmax","Tmin","ET")

				df.gw <- data.frame(cbind(Y.gw,X.gw))

				zz.gw <- lm(Y.gw ~ ., data=df.gw)
				z.coefs.gw[,xx,yy] <- zz.gw$coef
				z.resid.gw[,xx,yy] <- zz.gw$resid 

				lon.GKs <- c(lon.GKs, lon.GK[yy])
				lat.GKs <- c(lat.GKs, lat.GK[xx])
			}
		}else{
				next
		}
	}
}


z.mean.gw <- mean(z.resid.gw,na.rm=T)
z.sd.gw <- sd(z.resid.gw,na.rm=T)

lon.lat.knots <- cbind(as.vector(t(lon[xx.id,yy.id])),as.vector(t(lat[xx.id,yy.id])))
z.coefs.id <- !is.na(as.vector(t(z.coefs.gw[1,xx.id,yy.id])))
lon.lat.knots <- lon.lat.knots[z.coefs.id,]

lon.lat.knots.GK <- cbind(lon.GKs, lat.GKs)

lon.lat.preds.GK <- expand.grid(lon.GK,lat.GK)#[!is.na(as.vector(t(basin.mask))),]
lon.lat.preds.LL <- cbind(as.vector(lon),as.vector(lat))

n.knots <- nrow(lon.lat.knots)

z.resid.mat.gw <- matrix(NA,nrow=(tt.month-1),ncol=n.knots)

for(k in 1:(tt.month-1)) z.resid.mat.gw[k,] <- as.vector(t(z.resid.gw[k,,]))[!is.na(as.vector(t(z.resid.gw[k,,])))]

dist.mat <- rdist(lon.lat.knots)
diag(dist.mat) <- 0


mo.month2 <- mo.month[-1]
params.gw <- matrix(0,nrow=12,ncol=3)
for(k in 1:12){
	NAPAcov.gw <- cor(z.resid.mat.gw[mo.month2==k,],use="complete")
	NAPAvario.gw <- var(z.resid.mat.gw[mo.month2==k,])*(1 - NAPAcov.gw)
vario <- NAPAvario.gw
params.gw[k,2:3] <- optim(par=c(max(vario),max(dist.mat)),fn=LS)$par 
}
params.gw[params.gw<0]=0 

z.coefs.pred.gw <- array(data=NA,dim=c(nrow(z.coefs.gw),dim(napa)[2:3]))

for(k in 1:nrow(z.coefs.gw)){
	z.coefs.gw.tmp <- as.vector(t(z.coefs.gw[k,,]))[!is.na(as.vector(t(z.coefs.gw[k,,])))]
	z.coefs.pred.gw[k,,] <- matrix(predict(Krig(lon.lat.knots, z.coefs.gw.tmp,m=1),lon.lat.preds.LL),nrow=70,ncol=74)
}


n.sims = 10#0
tt.month.sim <- tt.month-1 
napa.gw.pred <- array(data=NA,dim=c(tt.month.sim,dim(et.m.pred)[2:3],n.sims))

for(i in 1:n.sims){
print(i)
for(k in 1:(tt.month.sim)){
	#print(k)
	if(k==1){
		#Xo.pred.gw <- rbind(1, as.vector(t(napa.m.s[k,,])), as.vector(t(prcp.m[k+1,,])), as.vector(t(tmax.m[k+1,,])), as.vector(t(tmin.m[k+1,,])), as.vector(t(et.m.pred[k+1,,])))
		mu.gw <- z.coefs.pred.gw[1,,]*matrix(1,nrow=74,ncol=70) + z.coefs.pred.gw[2,,]*napa.m.s[k,,] + z.coefs.pred.gw[3,,]*prcp.m[k,,] + z.coefs.pred.gw[4,,]*tmax.m[k,,] + z.coefs.pred.gw[5,,]*tmin.m[k,,] + z.coefs.pred.gw[6,,]*et.m.pred[k,,]
	}else{
		#Xo.pred.gw <- rbind(1, as.vector(t(napa.gw.pred[k-1,,,i])), as.vector(t(prcp.m[k+1,,])), as.vector(t(tmax.m[k+1,,])), as.vector(t(tmin.m[k+1,,])), as.vector(t(et.m.pred[k+1,,])))
		mu.gw <- z.coefs.pred.gw[1,,]*matrix(1,nrow=74,ncol=70) + z.coefs.pred.gw[2,,]*napa.gw.pred[k-1,,,i] + z.coefs.pred.gw[3,,]*prcp.m[k,,] + z.coefs.pred.gw[4,,]*tmax.m[k,,] + z.coefs.pred.gw[5,,]*tmin.m[k,,] + z.coefs.pred.gw[6,,]*et.m.pred[k,,]
	}

	#mu.gw <- matrix(apply((z.coefs.pred.gw * Xo.pred.gw), 2, sum),nrow=dim(et.m.pred)[2],ncol=dim(et.m.pred)[3],byrow=T)
	
	w.gw  <- matrix(suppressWarnings(grf(nrow(lon.lat.preds.LL),grid=lon.lat.preds.LL,cov.model="exponential",
              cov.pars=c(params.gw[mo.month2[k],2],params.gw[mo.month2[k],3]),nugget=params.gw[mo.month2[k],1],mean=rep(z.mean.gw,nrow(lon.lat.preds.LL)),messages=FALSE)$data),
			  nrow=dim(napa.gw.pred)[2],ncol=dim(napa.gw.pred)[3],byrow=T)

	napa.gw.pred[k,,,i] <- mu.gw + w.gw 
}
}


#plot(apply(napa.gw.pred[,,,1],1,mean,na.rm=T),type='l',ylim=range(c(napa.gw.pred, apply(napa.m.s,1,mean,na.rm=T))))
#for(k in 2:n.sims) lines(apply(napa.gw.pred[,,,k],1,mean,na.rm=T))

boxplot(apply(napa.gw.pred,c(4,1),mean,na.rm=T))
lines(apply(napa.m.s[2:nrow(napa.m.s),,],1,mean,na.rm=T),lwd=2,col='red',lty=2)



library(RColorBrewer)
temp.colors <- colorRampPalette(brewer.pal(9,"RdYlBu"))(100)

wet.id <- which(yr.month2==2000 & mo.month2==12)
dry.id <- which(yr.month2==2012 & mo.month2==2)

par(mfrow=c(2,2))
image.plot(t(napa.m.s[(wet.id+1),,]) - apply(napa.gw.pred[wet.id,,,],c(2,1),mean), col=temp.colors)
image.plot((apply(napa.gw.pred[wet.id,,,],c(2,1),quantile,probs=0.975,na.rm=T) - apply(napa.gw.pred[wet.id,,,],c(2,1),quantile,probs=0.025,na.rm=T)),col=temp.colors[(length(temp.colors)/2):length(temp.colors)])

image.plot(t(napa.m.s[(dry.id+1),,]) - apply(napa.gw.pred[dry.id,,,],c(2,1),mean), col=temp.colors)
image.plot((apply(napa.gw.pred[dry.id,,,],c(2,1),quantile,probs=0.975,na.rm=T) - apply(napa.gw.pred[dry.id,,,],c(2,1),quantile,probs=0.025,na.rm=T)),col=temp.colors[(length(temp.colors)/2):length(temp.colors)]) 


zl.wet <- range(c(napa.m.s[(wet.id+1),,], apply(napa.gw.pred[wet.id,,,],c(2,1),mean,na.rm=T)),na.rm=T)
zl.dry <- range(c(napa.m.s[(dry.id+1),,], apply(napa.gw.pred[dry.id,,,],c(2,1),mean,na.rm=T)),na.rm=T)


library(raster)
library(rgdal)
Salado_A_gk <- readOGR("/Users/andrew/CNH3/Salado/data/",layer="Salado_A_GK")
ext <- extent(min(predloc[,1]),max(predloc[,1]),min(predloc[,2]),max(predloc[,2]))

zl <- range(c(napa.m.s[(wet.id+1),,], apply(napa.gw.pred[wet.id,,,],c(2,1),mean,na.rm=T),napa.m.s[(dry.id+1),,], apply(napa.gw.pred[dry.id,,,],c(2,1),mean,na.rm=T)),na.rm=T)

zl.wet.diff <- range(napa.m.s[(wet.id+1),,] - apply(napa.gw.pred[wet.id,,,],1:2,mean,na.rm=T),na.rm=T)
if(abs(min(zl.wet.diff)) < max(zl.wet.diff)) zl.wet.diff[1] <- -1*max(zl.wet.diff)

zl.dry.diff <- range(napa.m.s[(dry.id+1),,] - apply(napa.gw.pred[dry.id,,,],1:2,mean,na.rm=T),na.rm=T)
if(abs(min(zl.dry.diff)) > max(zl.dry.diff)) zl.dry.diff[2] <- -1*min(zl.dry.diff)


Salado_A_gk <- readOGR("/Users/andrew/CNH3/Salado/data/",layer="Salado_A_GK")
ext <- extent(min(predloc[,1])-2500,max(predloc[,1])+2500,min(predloc[,2])-2500,max(predloc[,2])+2500)

pdf(file="../figs/metamodel-wet-and-dry-months-and-diff-raster.pdf",width=9,height=5,pointsize=10)
par(mfrow=c(2,3),mar=c(2,2,2,2))
napa.wet.r <- raster(napa.m.s[(wet.id+1),,])
extent(napa.wet.r) <- ext 
image.plot(flip(napa.wet.r,2),zlim=zl,xlab="",ylab="",main="(a)",col=temp.colors,axes=F)
plot(Salado_A_gk,add=T,lwd=2)
box()
napa.pred.wet.r <- raster(apply(napa.gw.pred[wet.id,,,],1:2,mean,na.rm=T))
extent(napa.pred.wet.r) <- ext 
image.plot(flip(napa.pred.wet.r,2),zlim=zl,xlab="",ylab="",main="(b)",col=temp.colors,axes=F)
plot(Salado_A_gk,add=T,lwd=2)
box()

napa.wet.diff.r <- raster(napa.m.s[(wet.id+1),,] - apply(napa.gw.pred[wet.id,,,],1:2,mean,na.rm=T))
extent(napa.wet.diff.r) <- ext 
image.plot(flip(napa.wet.diff.r,2),zlim=zl.wet.diff,xlab="",ylab="",main="(c)",col=temp.colors,axes=F)
plot(Salado_A_gk,add=T,lwd=2)
box()

napa.dry.r <- raster(napa.m.s[(dry.id+1),,])
extent(napa.dry.r) <- ext 
image.plot(flip(napa.dry.r,2),zlim=zl,xlab="",ylab="",main="(d)",col=temp.colors,axes=F)
plot(Salado_A_gk,add=T,lwd=2)
box()
napa.pred.dry.r <- raster(apply(napa.gw.pred[dry.id,,,],1:2,mean,na.rm=T))
extent(napa.pred.dry.r) <- ext 
image.plot(flip(napa.pred.dry.r,2),zlim=zl,xlab="",ylab="",main="(e)",col=temp.colors,axes=F)
plot(Salado_A_gk,add=T,lwd=2)
box()

napa.dry.diff.r <- raster(napa.m.s[(dry.id+1),,] - apply(napa.gw.pred[dry.id,,,],1:2,mean,na.rm=T))
extent(napa.dry.diff.r) <- ext 
image.plot(flip(napa.dry.diff.r,2),zlim=zl.dry.diff,xlab="",ylab="",main="(f)",col=temp.colors,axes=F)
plot(Salado_A_gk,add=T,lwd=2)
box()

dev.off()

library(hydroGOF)
rmses.mat <- maes.mat <- pbias.mat <- matrix(NA,nrow=dim(napa.m.s)[2],ncol=dim(napa.m.s)[3])
for(xx in 1:nx){
	for(yy in 1:ny){
		if(!is.na(mean(napa.m.s[,xx,yy]))){
			rmses.mat[xx,yy] <- rmse(apply(napa.gw.pred[,xx,yy,],1,mean), napa.m.s[2:nrow(napa.m.s),xx,yy])
			maes.mat[xx,yy] <- mae(apply(napa.gw.pred[,xx,yy,],1,mean), napa.m.s[2:nrow(napa.m.s),xx,yy])
			pbias.mat[xx,yy] <- pbias(apply(napa.gw.pred[,xx,yy,],1,mean), napa.m.s[2:nrow(napa.m.s),xx,yy])
		}else{
			next 
		}
	}
}

image.plot(t(rmses.mat))
image.plot(t(maes.mat))
image.plot(t(pbias.mat))

zl.rmse <- range(rmses.mat,na.rm=T)
zl.maes <- range(maes.mat,na.rm=T)
zl.bias <- range(pbias.mat,na.rm=T)
if(abs(min(zl.bias)) > max(zl.bias)) zl.bias[2] <- -1*min(zl.bias)

pdf(file="../figs/metamodel-sumstats.pdf",width=9,height=2.5,pointsize=10)
par(mfrow=c(1,3),mar=c(2,2,2,2))
rmses <- raster(rmses.mat)
extent(rmses) <- ext 
image.plot(flip(rmses,2),xlab="",ylab="",main="RMSE",zlim=zl.rmse,col=rev(temp.colors[1:(length(temp.colors)/2)]),axes=F)
plot(Salado_A_gk,add=T,lwd=2)
box()

maes <- raster(maes.mat)
extent(maes) <- ext 
image.plot(flip(maes,2),xlab="",ylab="",main="MAE", zlim=zl.maes, col=rev(temp.colors[1:(length(temp.colors)/2)]),axes=F)
plot(Salado_A_gk, add=T, lwd=2)
box()

pbiass <- raster(pbias.mat)
extent(pbiass) <- ext 
image.plot(flip(pbiass,2),xlab="",ylab="",main="% Bias", zlim=zl.bias, col=temp.colors,axes=F)
plot(Salado_A_gk, add=T, lwd=2)
box() 

dev.off()

