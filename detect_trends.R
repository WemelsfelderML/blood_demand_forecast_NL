library(forecast)
library(ggplot2)
library(gridExtra)
library(plyr)
library(lubridate)
library(numbers)
library(data.table)
library(R.utils)
library(DT)
library(RColorBrewer)
library(stringr)
library(tsfeatures)
library(nlme)

# SETTINGS

# working directory
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"
source(paste0(ROOTDIR, "pffunctions.R"))

# Read data file
df <- read.delim(file = "data/data.csv", header = FALSE, sep = ";")
df <- df[df$V2 != "",][c(2:10)][-c(1),]
rw_years <- 3

names(df) <- as.matrix(df[1, ]);
df <- df[-1, ];
rownames(df) <- 1:nrow(df);

df1 <- df[0,]
for (i in 2:9){
  df2 <- df[,c(1,i)];
  df2$ABO <- colnames(df)[i];
  names(df2) <- c("date", "quantity", "ABO")
  df1 <- rbind(df1, df2)
}
df <- df1

df$date <- mdy(df$date)
df$quantity <- as.numeric(as.character(df$quantity))
df$quantity[is.na(df$quantity)] <- 0
df$product <- "Erytrocyten"

# Product codes for red cell products
red.distr <- df[,c("date", "product", "quantity", "ABO")] # "volume", "exp"

# Create a full sequence of dates for imputation purposes
all.dates <- (seq.Date(min(red.distr$date),
                       max(red.distr$date),
                       "day"))
all.red <- aggregate(red.distr$quantity, by = list(red.distr$date), sum); colnames(all.red) <- c("date", "pcs")
all.red <- merge(x = data.frame(date = all.dates),
                 y = all.red,
                 all.x = TRUE)

red.monthly <- aggregate(pcs ~ month(date) + year(date), data = all.red, FUN = sum)
months <- seq(from = as.Date("2009-01-01"), to = max(red.distr$date), by = "month")

adj <- as.numeric(bizdays(ts(months, start = decimal_date(as.Date("2009-01-01")), frequency = 12), FinCenter = "Zurich"))
reverse_adj <- as.numeric(bizdays(ts(seq(23), start = decimal_date(months[length(months)]), frequency = 12), FinCenter = "Zurich")) # This is the old implementation that used to be fed into the forecasting function. We'll now repurpose it so it can be used both for tabling and plotting.

# Create a master frame
monthly <- data.frame(date = months,
                      red = red.monthly$pcs/adj)
monthly_real <- data.frame(date = months,
                           red = red.monthly$pcs)

# We will need these stored
modelnames <- c("SNAIVE", "5-MA", "7-MA", "9-MA", "12-MA", "STL", "ETS", "TBATS", "STLF", "ARIMAX", "DYNREG", "NN", "COMBINED")

# Weekly aggregation is a bit more tricky operation as weekly aggregate gives incorrect results, so we'll use
# a function of our own

red.weekly <- aggregate_weekly(all.red)
weekly <- data.frame(week = red.weekly$week,
                     date = red.weekly$date,
                     red = red.weekly$pcs)

daily <- all.red

monthly$year <- as.numeric(format(as.Date(monthly$date, format="%Y/%m/%d"),"%Y"))
weekly$year <- as.numeric(format(as.Date(weekly$date, format="%Y/%m/%d"),"%Y"))
daily$year <- as.numeric(format(as.Date(daily$date, format="%Y/%m/%d"),"%Y"))
daily <- daily[daily$year == tail(daily$year, 1),]
daily <- tail(daily, nrow(daily)-(nrow(daily) %% 7))
rownames(daily) <- NULL
daily$day <- as.numeric(rownames(daily))

ts.m <- ts(data = monthly$red, frequency = 12)
ts.w <- ts(data = weekly$red, frequency = 52)
ts.d <- ts(data = daily$pcs, frequency = 7)

stl.m = stl(ts.m, "periodic")
autoplot(stl.m, main = "STL decomposition of all red blood cell data\nseasonality: month / year")

stl.w = stl(ts.w, "periodic")
autoplot(stl.w, main = "STL decomposition of all red blood cell data\nseasonality: week / year")

# stl.d = stl(ts.d, "periodic")
# autoplot(stl.d, main = paste0("STL decomposition of all red blood cell data\nseasonality: day / week (", tail(daily$year, 1), ")"))

# years <- unique(monthly$year)
# for (y in head(years, length(years) - 2)){
#   m <- monthly[monthly$year %in% c(y:(y+2)),]
#   ts.m <- ts(data = m$red, frequency = 12, start = head(m$year, 1))
#   stl.m = stl(ts.m, "periodic")
#   autoplot(stl.m, main = "STL decomposition of all red blood cell data\nseasonality: month / year")
#   print(paste0(y, "-", y+3))
#   print(paste0("seasonal strength: ", round(stl_features(ts.m)["seasonal_strength"],3)))
#   print(paste0("trend: ", round(stl_features(ts.m)["trend"],5)))
#   print("")
# }


#####################
#####################

# # MANUAL SEASONALITY
# 
# autoplot(as.ts(ts.m), main="red monthly")
# autoplot(as.ts(ts.w), main="red weekly")
# 
# # trend
trend.m = ma(ts.m, order = 12, centre = T)
# trend.w = ma(ts.w, order = 52, centre = T)
# autoplot(trend.m, main="trend monthly")
# autoplot(trend.w, main="trend weekly")
# 
# # detrend
detrend.m = ts.m - trend.m
# detrend.w = ts.w - trend.w
# autoplot(detrend.m, main="detrend monthly")
# autoplot(detrend.w, main="detrend weekly")
# 
# # seasonality
m.m = t(matrix(data = detrend.m, nrow = 12))
seasonal.m = colMeans(m.m, na.rm = T)
autoplot(as.ts(rep(seasonal.m,1)), main="seasonality monthly")
# 
# m.w = t(matrix(data = detrend.w, nrow = 52))
# seasonal.w = colMeans(m.w, na.rm = T)
# autoplot(as.ts(rep(seasonal.w,1)), main="seasonality weekly")
# 
# # random noise
# noise.m = ts.m - trend.m - seasonal.m
# noise.w = ts.w - trend.w - seasonal.w
# autoplot(as.ts(noise.m))
# autoplot(as.ts(noise.w))
# 
# # seasonality of days in a week
# trend.d = ma(ts.d, order = 7, centre = T)
# detrend.d = ts.d - trend.d
# m.d = t(matrix(data = detrend.d, nrow = 7))
# seasonal.d = colMeans(m.d, na.rm = T)
# autoplot(as.ts(rep(seasonal.d,1)), main="seasonality: day / week")
