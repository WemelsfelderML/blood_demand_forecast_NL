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
library(tidyverse)
library(zoo)

#SETTINGS

#Select are we using NL or FIN data
#NL <- FALSE
NL <- TRUE

# working directory
if (NL) {
  ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"
  groups <- c("RED", "Ominus", "Oplus", "Aminus", "Aplus", "Bminus", "Bplus", "ABminus", "ABplus", "PLAT")
} else {
  ROOTDIR <- "~/Work/proj/OPERATIONAL/blood_demand_forecast_NL/"
  groups <- c("RED")
}

period <- "m"               # m for monthly, w for weakly
rw <- 3
# method.select <- c("snaive", "5-MA", "7-MA", "9-MA", "12-MA", "stl", "ets", "tbats", "stlf", "arimax", "dynreg", "nn", "combined")
method.select <- c("stlf", "7-MA", "ets", "nn", "dynreg", "tbats")        # monthly
# method.select <- c("combined", "snaive", "stl", "5-MA", "ets", "dynreg")    # weekly
merge.months = 6

# reshape error data frames and put them together in a hash set with rolling window size as key
df <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_all.txt"), header = FALSE, sep = ";")
#df <- df[df$V1 %in% method.select,]
#> df$V1 %in% method.select
#[1] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE
#This does not guarantee that they would be in the same order
# arvasmi@tk-arvasmi:~/Work/proj/OPERATIONAL/blood_demand_forecast_NL/rw_testing$ cut -d ";" -f1 m_errors_rwy3_all.txt 
# Red
# snaive
# 5-MA
# 7-MA
# 9-MA
# 12-MA
# stl
# ets
# tbats
# stlf
# arimax
# dynreg
# nn
# combined
# > cbind(df[df$V1 %in% method.select,1:2],method.select)
# V1        V2 method.select
# 4    7-MA 7.2527823          stlf
# 8     ets 6.0121560          7-MA
# 9   tbats 1.2525924           ets
# 10   stlf 0.8383647            nn
# 12 dynreg 3.2047730        dynreg
# 13     nn 5.9707662         tbats
#Put method.select <- c("stlf", "7-MA", "ets", "nn", "dynreg", "tbats")        # monthly
#in the same order that they are in the file. IS THIS THE SAME IN FIN AND NL???
method.select <- c("7-MA", "ets","tbats", "stlf", "dynreg", "nn")        # monthly
#cbind(df[df$V1 %in% method.select,1:2],method.select)
# > cbind(df[df$V1 %in% method.select,1:2],method.select)
# V1        V2 method.select
# 4    7-MA 7.2527823          7-MA
# 8     ets 6.0121560           ets
# 9   tbats 1.2525924         tbats
# 10   stlf 0.8383647          stlf
# 12 dynreg 3.2047730        dynreg
# 13     nn 5.9707662            nn
df <- df[df$V1 %in% method.select,]
rownames(df) <- NULL

# abstract year span from original data
if (NL) {
  d <- read.delim(file = paste0(ROOTDIR, "data/data_processed.csv"), header = TRUE, sep = ",")
  d$year <- as.numeric(format(as.Date(d$Date, format="%d/%m/%Y"),"%Y"))
  #This is only used to get the max value
} else {
  d <- data.frame(year=2020)
  #Finish data extends to 2020.06. Remove last 6 timepoints
  df <- df[,- ((ncol(df) -6) : ncol(df))]
}

# number of methods tested for each blood group
m <- length(method.select)
method.names <- array(df[c(1:m),1])
if (m == 13) {
  colors = append(brewer.pal(12, name="Paired"), "#555555")
} else {
  colors = brewer.pal(m, name="Paired")
}

# merge months to get a clearer overview
df.sums <- data.frame(df$V1)
df <- df[,-c(1:((ncol(df)-1)%%merge.months + 1))]
n <- ncol(df) / merge.months
for (i in c(0:(n-1))) {
  start <- (i*merge.months)+1
  stop <- start + merge.months - 1
  df.sums[,ncol(df.sums)+1] <- rowSums(df[,c(start:stop)])/merge.months
}

# set periodicity for text in plot
if (period == "m") {
  p <- 12
  t <- " months"
  pd <- ", monthly"
} else if (period == "w") {
  p <- 52
  t <- " weeks"
  pd <- ", weekly"
}

# plot averaged errors per method per year in history, 
# for each blood group separately
for (i in c(0:(length(groups)-1))) {
  start <- (i*m) + 1
  stop <- start + m - 1
  df.group <- df.sums[c(start:stop),]
  
  # reorder rows so method order is equal to method.select
  colnames(df.group)[colnames(df.group) == "df.V1"] <- "method"
  df.plot = df.group[,2:ncol(df.group)]
  for (mi in c(1:length(method.select))) {
    method <- method.select[mi]
    df.plot[mi,] = df.group[df.group["method"]==method, 2:ncol(df.group)]
  }
  
  group <- groups[i+1]
  if (NL) {
    png(file= paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", toString(rw), "_", group, ".png"))
  } else {
    png(file= paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", toString(rw), "_", group, "_FIN.png"))
  }
  plot(1, type = "n", main = paste0(group, ", rw", rw, pd, "\naveraged over each ", merge.months, t), xlim=c((max(d$year)-(n*(merge.months/p))+1), (max(d$year)+1-(merge.months)/p)), ylim=c(0,max(df.plot)), xlab="year", ylab="avg error")
  for (i in c(1:nrow(df.group))) {
    lines(x=seq((max(d$year)-(n*(merge.months/p))+1), (max(d$year)+1-(merge.months)/p), (merge.months/p)), y=df.plot[i,], col=colors[i], type="o", lwd=2, pch = 19)
  }
  legend("topleft", legend=method.names, col=colors, pch = 19)
  dev.off()
}

#Alternative plot for RED only as for the moment that is all what we have for FIN data
if (NL) {
  # take only RED section of data
  start <- (0*m) + 1
  stop <- start + m - 1
  df.group <- df.sums[c(start:stop),c(1:ncol(df.sums))]
  colnames(df.group)[colnames(df.group) == "df.V1"] <- "method"
  
  # reorder rows so method order is equal to method.select
  df.plot = df.group[,2:ncol(df.group)]
  for (i in c(1:length(method.select))) {
    method <- method.select[i]
    df.plot[i,] = df.group[df.group["method"]==method,2:ncol(df.group)]
  }
} else {
  df.plot <- as_tibble(df.group)
}
names(colors) <- method.select
colnames(df.plot) <- seq((max(d$year)-(n*(merge.months/p))+1), (max(d$year)+1-(merge.months)/p), (merge.months/p))
df.plot$Method <- method.select
df.plot <- df.plot %>% pivot_longer(cols= !Method) %>% mutate(Year=as.numeric(name))
gr <- ggplot(df.plot)
gr <- gr + geom_point(aes(x=Year,y=value,color=Method,group=Method),size=4)
gr <- gr + geom_line(aes(x=Year,y=value,color=Method,group=Method),size=1.1)
gr <- gr + scale_colour_manual(name = "Method",values = colors) + ylab("Prediction error \n (% of pcs averaged over each 6 months)")
gr <- gr + theme_bw(base_size = 18)
gr <- gr + theme(legend.position = "bottom", legend.direction = "horizontal")
if (NL) {
  filename <- paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", toString(rw), "_", "red", "_gg.pdf")
} else {
  filename <- paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", toString(rw), "_", "red", "_gg_FIN.pdf")
}
ggsave(filename=filename, gr, width = 180,  height = 180,units="mm", dpi=600, scale=1.0)

<<<<<<< HEAD




=======
>>>>>>> d00053341a2b173c09b75264b309059f60350f77
#Another way to count....
# d2 <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_all.txt"), header = FALSE, sep = ";")
# d2 <- d2[d2$V1 %in% method.select,]
# rownames(d2) <- NULL
# as.yearmon(2008 + seq(0, (ncol(d2)-1))/12)

