library(hash)
library(RColorBrewer)
library(stringr)
library(forecast)
library(ggplot2)
library(gridExtra)
library(knitr)
library(plyr)
library(lubridate)
library(numbers)
library(data.table)
library(R.utils)
library(DT)

# SETTINGS

# working directory
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"
period <- "m"               # m for monthly, w for weakly
rw <- 3
groups <- c("RED", "Ominus", "Oplus", "Aminus", "Aplus", "Bminus", "Bplus", "ABminus", "ABplus", "PLAT")
colors = append(brewer.pal(12, name="Paired"), "#555555")
merge.months = 6

# PRE-PROCESSING

# reshape error data frames and put them together in a hash set with rolling window size as key
df <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_all.txt"), header = FALSE, sep = ";")
df <- df[c(1:140),]
modelnames <- df[c(2:14), "V1"]
# df <- df[c(1:(nrow(df)/10)),]

# abstract year span from original data
d <- read.delim(file = paste0(ROOTDIR, "data/data_processed.csv"), header = TRUE, sep = ",")
d$year <- as.numeric(format(as.Date(d$Date, format="%d/%m/%Y"),"%Y"))

# number of methods tested for each blood group
m <- nrow(df)/length(groups) - 1

# merge months to get a clearer overview
df.sums <- data.frame(df$V1)
df <- df[,-c(1:((ncol(df)-1)%%merge.months + 1))]
n <- ncol(df) / merge.months
for (i in c(0:(n-1))) {
  start <- (i*merge.months)+1
  stop <- start + merge.months - 1
  df.sums[,ncol(df.sums)+1] <- rowSums(df[,c(start:stop)])/merge.months
}

# plot averaged erorrs per method per year in history, 
# for each blood group separately
for (i in c(0:(length(groups)-1))) {
  start <- (i*m) + i + 2
  stop <- start + m - 1
  df.group <- df.sums[c(start:stop),c(2:ncol(df.sums))]
  group <- groups[i+1]
  
  png(file= paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", toString(rw), "_", group, ".png"))
  plot(1, type = "n", main = paste0(group, ", rw", rw, "\naveraged over each ", merge.months, " months"), xlim=c((max(d$year)-(n*(merge.months/12))+1), (max(d$year)+1-(merge.months)/12)), ylim=c(0,max(df.group)), xlab="year", ylab="avg error")
  for (i in c(1:nrow(df.group))) {
    print(df.group[i,],)
    lines(x=seq((max(d$year)-(n*(merge.months/12))+1), (max(d$year)+1-(merge.months)/12), (merge.months/12)), y=df.group[i,], col=colors[i], type="o", lwd=2, pch = 19)
  }
  legend("topleft", legend=modelnames, col=colors, pch = 19)
  dev.off()
}

