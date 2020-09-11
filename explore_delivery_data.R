# initialization
rm(list=ls())
dev.off()
cat("\014")  
source("pffunctions_MPJ.R")

data<-read.csv("data/ERYS2009_2019_MPJ NP.csv",header=T)
colnames(data)
colnames(data)<-c("Date","O+","O-","A+","A-","B+","B-","AB+","AB-","total")
colnames(data)

plot(data$total)
max(data$total)
min(data$total)
th<-data[seq(from=1, to=4017, by=7),]
fr<-data[seq(from=2, to=4017, by=7),]
sa<-data[seq(from=3, to=4017, by=7),]
su<-data[seq(from=4, to=4017, by=7),]
mo<-data[seq(from=5, to=4017, by=7),]
tu<-data[seq(from=6, to=4017, by=7),]
we<-data[seq(from=7, to=4017, by=7),]
plot(mo$Date,mo$total)

plot(mo$total, ylim=c(1,2700))
abline(v=seq(from=1, to=600, by=52), col=8)
points(tu$total, col=2)
points(we$total, col=3)
points(th$total, col=4)
points(fr$total, col=5)
points(sa$total, col=6)
points(su$total, col=7)

require(splines)
fit1<-smooth.spline(x=1:nrow(mo),y=mo$total,cv=T) 
fit1
lines(fit1,col=1,lwd=2)


plot(NULL, xlim=c(0,11*52+5), ylim=c(0,2700), ylab="Nr of products delivered", xlab="Weeks since 1/1/2009")
abline(v=seq(from=1, to=600, by= 365.25/7), col=8)
for (i in 1:7) {
  d<-data[seq(from=i, to=nrow(data), by=7),]  
  col=((i+2) %% 7)+1
  print(col)
  points(d$total, col=col)
  fit1<-smooth.spline(x=1:nrow(d),y=d$total,cv=T) 
  lines(fit1,col=col,lwd=4)
  lines(fit1,col=8,lwd=1)
}
legend("topright", c("mo","tu","we","th","fr","sa","su"),col=1:7, lty=rep(1,7))

(nweeks<-4017%/%7)
nweeks*7

# aggregate data per week
weekdata<-colSums(data[1:7,2:ncol(data)])
for (i in 2:nweeks) weekdata<-rbind(weekdata,colSums(data[1:7+7*(i-1),2:ncol(data)]))
weekdata<-as.data.frame(weekdata)
weekdata$week<-1:nrow(weekdata)

plot(NULL, xlim=c(0,nrow(weekdata)), ylim=c(0,max(weekdata)), ylab="Nr of products delivered", xlab="Weeks since 1/1/2009")
abline(v=seq(from=1, to=600, by= 365.25/7), col=8)
for (i in 9:1) {
  d<-weekdata[,i]  
  col=i
  points(d, col=col)
  fit1<-smooth.spline(x=1:length(d),y=d,cv=T) 
  lines(fit1,col=col,lwd=4)
  lines(fit1,col=8,lwd=1)
}
legend("topright", colnames(weekdata)[1:9],col=1:9, lty=rep(1,9))

