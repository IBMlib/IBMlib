dat <- read.table("fort.65");
pdf("bee_test.pdf",width=210/25.4,height=210/25.4,pointsize=14)
plot(dat$V1,dat$V2,asp=1,type="l",xlab="Latitude",ylab="Longitude");
dev.off()
