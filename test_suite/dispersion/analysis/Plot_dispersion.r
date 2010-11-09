rm(list=ls())

#Load data
time_step <- 1800
dat.fwd <- read.table("disperse_forwards.txt",col.names=c("time","lon","lat","depth"))
dat.back <- read.table("disperse_backwards.txt",col.names=c("time","lon","lat","depth"))
dat <- rbind(dat.back,dat.fwd)

#Simulation parameters
hdiff <- 10   # m²/s
vdiff <- 0.01 # m²/s
current.dir <- 45
centre <- c(4,55,500)

#Setup pdf
pdf("Dispersion_plots.pdf",width=210/25.4,height=297/25.4,pointsize=14)
par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(5,4,1,1))

#Visualise dispersion from above
plot(dat$lon,dat$lat,xlab="Longitude",ylab="Latitude",asp=1,
        pch=".",col=ifelse(dat$time<0,"blue","red"))
grad <- tan(current.dir*pi/180)*cos(centre[2]/180*pi)
abline(a=centre[2]-grad*centre[1],b=grad,lwd=2,col="black")

#Visualise vertical dispersion
plot(dat$time,dat$depth,xlab="Time(s)",ylab="Depth",ylim=rev(range(pretty(dat$depth))),
        pch=".",col=ifelse(dat$time<0,"blue","red"))
abline(h=centre[3],lwd=2,col="black")

#Check horizontal dispersion coefficients
var.lon <- tapply(dat$lon,dat$time,sd)^2
var.lat <- tapply(dat$lat,dat$time,sd)^2
var.lims <- range(pretty(c(var.lon,var.lat)))
times <- as.numeric(names(var.lon))
plot(NA,xlim=range(pretty(abs(times))),ylim=var.lims,xlab="Time(s)",ylab="Lon/lat variance")
points(abs(times),var.lon,pch=".",col=ifelse(times<0,"blue","red"))
points(abs(times),var.lat,pch="+",col=ifelse(times<0,"blue","red"))
legend("bottomright",legend=c("Lon","Lat"),pch=c(".","+"),col="black")
grid()
abline(a=0,b=2*hdiff*(180/6370000/pi)^2,col="black",lwd=2)
abline(a=0,b=2*hdiff*(180/6370000/cos(centre[2]/180*pi)/pi)^2,col="black",lwd=2)


#Check vertical dispersion coefficient
vars <- tapply(dat$depth,dat$time,sd)^2
var.dat <- data.frame(time=as.numeric(names(vars)),var=vars)
fwd.diff<- abs(mean(diff(subset(var.dat$var,var.dat$time>0)))/time_step/2)
back.diff<- abs(mean(diff(subset(var.dat$var,var.dat$time<0)))/time_step/2)
plot(abs(var.dat$time),var.dat$var,xlab="Time(s)",ylab="Depth Variance (m²)",
        xlim=range(pretty(abs(var.dat$time))),
        pch=19,col=ifelse(var.dat$time<0,"blue","red"))
abline(a=0,b=vdiff*2,col="black",lwd=2)
legend("bottomright",legend=c(sprintf("D fwd  %6.4f",fwd.diff),
      sprintf("D back %6.4f",back.diff)),bty="n")
grid()
dev.off()
