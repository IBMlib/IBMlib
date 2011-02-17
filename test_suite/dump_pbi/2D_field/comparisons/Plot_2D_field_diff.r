

rm(list=ls())
#library(lattice)

#Load data
old.dat <- read.table("one_fields.txt",header=TRUE)
new.dat <- read.table("two_fields.txt",header=TRUE)
var.name <- colnames(old.dat)[3]
colnames(old.dat)[3] <- "var"
var.name <- colnames(new.dat)[3]
colnames(new.dat)[3] <- "var"
new.dat$var[is.na(new.dat$var)] <- -999
old.dat$var[is.na(old.dat$var)] <- -999

new.mat <- xtabs(var ~ lon+lat,new.dat)
old.mat <- xtabs(var ~ lon+lat,old.dat)
new.mat[new.mat==-999] <- NA
old.mat[old.mat==-999] <- NA

mat <- new.mat-old.mat


#Setup pdf
pdf("field_diff.pdf",width=210/25.4,height=210/25.4,pointsize=14)

#Filled contours
filled.contour(as.numeric(rownames(mat)),as.numeric(colnames(mat)),mat,nlevels=50,
  color.palette=rainbow,xlab="Longitude",ylab="Latitude",main=sprintf("%s diff",var.name),
  zlim=range(pretty(mat)))


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
