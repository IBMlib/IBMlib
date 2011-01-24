

rm(list=ls())
#library(lattice)

#Load data
dat <- read.table("2D_fields.txt",header=TRUE)
var.name <- colnames(dat)[3]
colnames(dat)[3] <- "var"
mat <- xtabs(var ~ lon+lat,dat)
mat[mat==0] <- NA
dat <- subset(dat,!is.na(dat$var))

if(nrow(dat) == 0) stop("No data in output fields.")


#Setup pdf
pdf(sprintf("%s_2D_fields.pdf",var.name),width=210/25.4,height=210/25.4,pointsize=14)

#Filled contours
filled.contour(as.numeric(rownames(mat)),as.numeric(colnames(mat)),mat,nlevels=25,
  color.palette=terrain.colors,xlab="Longitude",ylab="Latitude",main=var.name,
  zlim=range(pretty(dat$var)))


#Visualise field using lattice
#p <- levelplot(var ~ lon*lat,dat,
#      xlab="Longitude",ylab="Latitude",main=var.name,
#      col.regions = terrain.colors)
#
##End plotting
#print(p)
#
##Visualise using coloured points
#z.rng <- range(pretty(dat$var))
#dat$z.col <- (dat$var-min(z.rng))/(z.rng[2]-z.rng[1])
#plot(dat$lon,dat$lat,pch=19,col=grey(dat$z.col),xlab="Longitude",ylab="Latitude",main=var.name,)
#
dev.off()

cat("----------Plotting of Field complete.--------------\n");flush.console()
