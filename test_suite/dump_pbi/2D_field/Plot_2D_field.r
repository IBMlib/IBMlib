

rm(list=ls())

#Load data
dat <- read.table("2D_fields.txt",header=TRUE)
var.name <- colnames(dat)[3]
colnames(dat)[3] <- "var"

#Setup pdf
pdf("2D_fields.pdf",width=210/25.4,height=210/25.4,pointsize=14)

#Visualise field
p <- levelplot(var ~ lon*lat,dat,
      xlab="Longitude",ylab="Latitude",main=var.name,
      col.regions = terrain.colors,
      contour=TRUE,labels=FALSE)

#End plotting
print(p)


dev.off()
