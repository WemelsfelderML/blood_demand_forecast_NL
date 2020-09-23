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

# SETTINGS

# working directory
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"
source(paste0(ROOTDIR, "pffunctions.R"))


# DATA PRE-PROCESSING

d <- read.delim(file = paste0(ROOTDIR, "data/data.csv"), header = FALSE, sep = ";")
d <- d[d$V2 != "",]

df <- d[c(2:10)][-c(1),]
names(df) <- as.matrix(df[1, ]);
df <- df[-1, ];
rownames(df) <- 1:nrow(df);

df1 <- df[0,]
for (i in 2:9){
  df2 <- df[,c(1,i)];
  names(df2) <- c("date", "quantity")
  df1 <- rbind(df1, df2)
}
df1$date <- mdy(df1$date)
df1$quantity <- as.numeric(df1$quantity)
df <- df1

# Product codes for red cell products
red.distr <- df[,c("date", "quantity")]

# Create a full sequence of dates for imputation purposes
all.dates <- (seq.Date(min(red.distr$date),
                       max(red.distr$date),
                       "day"))

all.red <- aggregate(red.distr$quantity, by = list(red.distr$date), sum); colnames(all.red) <- c("date", "pcs")
all.red <- merge(x = data.frame(date = all.dates),
                 y = all.red,
                 all.x = TRUE)
all.red[is.na(all.red)] <- 0

months <- seq(from = as.Date("2009-01-01"), to = max(df$date), by = "month")
adj <- as.numeric(bizdays(ts(months, start = decimal_date(as.Date("2009-01-01")), frequency = 12), FinCenter = "Zurich"))
reverse_adj <- as.numeric(bizdays(ts(seq(23), start = decimal_date(months[length(months)]), frequency = 12), FinCenter = "Zurich")) 

# Aggregate all by months
red.monthly <- aggregate(pcs ~ month(date) + year(date), data = all.red, FUN = sum)
monthly <- data.frame(date = months,
                      red = red.monthly$pcs/adj)

red.weekly <- aggregate_weekly(all.red)
weekly <- data.frame(week = red.weekly$week,
                     date = red.weekly$date,
                     red = red.weekly$pcs)
weekly <- head(weekly, nrow(weekly)- (nrow(weekly) %% 52))

daily <- all.red


#####################
#####################

# STL SEASONALITY

monthly$year <- as.numeric(format(as.Date(monthly$date, format="%Y/%m/%d"),"%Y"))
weekly$year <- as.numeric(format(as.Date(weekly$date, format="%Y/%m/%d"),"%Y"))
daily$year <- as.numeric(format(as.Date(daily$date, format="%Y/%m/%d"),"%Y"))
daily <- daily[daily$year == tail(daily$year, 1),]
daily <- tail(daily, nrow(daily)-(nrow(daily) %% 7))
rownames(daily) <- NULL
daily$day <- as.numeric(rownames(daily))

ts.m <- ts(data = monthly$red, frequency = 12, start = head(monthly$year, 1))
ts.w <- ts(data = weekly$red, frequency = 52, start = head(weekly$year, 1))
ts.d <- ts(data = daily$pcs, frequency = 7, start = head(daily$day, 1))

stl.m = stl(ts.m, "periodic")
plot(stl.m, main = "STL decomposition of all red blood cell data\nseasonality: month / year")

stl.w = stl(ts.w, "periodic")
plot(stl.w, main = "STL decomposition of all red blood cell data\nseasonality: week / year")

stl.d = stl(ts.d, "periodic")
plot(stl.d, main = paste0("STL decomposition of all red blood cell data\nseasonality: day / week (", tail(daily$year, 1), ")"))


#####################
#####################

# MANUAL SEASONALITY

plot(as.ts(ts.m), main="red monthly")
plot(as.ts(ts.w), main="red weekly")

# trend
trend.m = ma(ts.m, order = 12, centre = T)
trend.w = ma(ts.w, order = 52, centre = T)
plot(as.ts(trend.m), main="trend monthly")
plot(as.ts(trend.w), main="trend weekly")

# detrend
detrend.m = ts.m - trend.m
detrend.w = ts.w - trend.w
plot(as.ts(detrend.m), main="detrend monthly")
plot(as.ts(detrend.w), main="detrend weekly")

# seasonality
m.m = t(matrix(data = detrend.m, nrow = 12))
seasonal.m = colMeans(m.m, na.rm = T)
plot(as.ts(rep(seasonal.m,1)), main="seasonality monthly")

m.w = t(matrix(data = detrend.w, nrow = 52))
seasonal.w = colMeans(m.w, na.rm = T)
plot(as.ts(rep(seasonal.w,1)), main="seasonality weekly")

# random noise
noise.m = ts.m - trend.m - seasonal.m
noise.w = ts.w - trend.w - seasonal.w
plot(as.ts(noise.m))
plot(as.ts(noise.w))

# seasonality of days in a week
trend.d = ma(ts.d, order = 7, centre = T)
detrend.d = ts.d - trend.d
m.d = t(matrix(data = detrend.d, nrow = 7))
seasonal.d = colMeans(m.d, na.rm = T)
plot(as.ts(rep(seasonal.d,1)), main="seasonality: day / week")
