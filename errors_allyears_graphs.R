#cat errors_allyears_graphs.R | R --slave --args 3 6 l
a <- commandArgs()
print(a)
#$ cat errors_allyears_graphs.R | R --slave --args 5 6
#[1] "/usr/lib/R/bin/exec/R" "--slave"               "--args"               
#[4] "5"  "6" "l"

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
NL <- FALSE
#NL <- TRUE

# working directory
if (NL) {
  ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"
  groups <- c("RED", "Ominus", "Oplus", "Aminus", "Aplus", "Bminus", "Bplus", "ABminus", "ABplus", "PLAT")
} else {
  ROOTDIR <- "~/Work/proj/OPERATIONAL/blood_demand_forecast_NL/20201023_error_ts_plots/"
  groups <- c("RED","PLAT")
}

period <- "m"               # m for monthly, w for weakly
rw <- 5
if (!is.na(commandArgs()[4])) {
  rw <- as.numeric(commandArgs()[4])
  cat("Setting rw to",rw,"\n")
}

merge.months = 6
if (!is.na(commandArgs()[5])) {
  merge.months <- as.numeric(commandArgs()[5])
  cat("Setting merge.months to",merge.months,"\n")
}

# method.select <- c("snaive", "5-MA", "7-MA", "9-MA", "12-MA", "stl", "ets", "tbats", "stlf", "arimax", "dynreg", "nn", "combined")
#method.select <- c("stlf", "7-MA", "ets", "nn", "dynreg", "tbats")        # monthly
# method.select <- c("combined", "snaive", "stl", "5-MA", "ets", "dynreg")    # weekly
method.select <- c("7-MA", "ets","tbats", "stlf", "dynreg", "nn")        # monthly

if (!is.na(commandArgs()[6])) {
  merge.param <- commandArgs()[6]
  if (merge.param == "l") {
    method.select <- c("snaive", "5-MA", "7-MA", "9-MA", "12-MA", "stl", "ets", "tbats", "stlf", "arimax", "dynreg", "nn", "combined")
  }
  if (merge.param == "s") {
    method.select <- c("7-MA", "ets","tbats", "stlf", "dynreg", "nn")
  }
}


# DATA PROCESSING

# reshape error data frames and put them together in a hash set with rolling window size as key
file <- paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_all.txt")
if (file.exists(file)) {
df <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_all.txt"), header = FALSE, sep = ";")
cat("Read",nrow(df),"rows from file",file,"\n")
} else {
  stop()
}

#df <- df[df$V1 %in% method.select,]
#> df$V1 %in% method.select
#[1] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE
#This does not guarantee that they would be in the same order
# arvasmi rw_testing$ cut -d ";" -f1 m_errors_rwy3_all.txt 
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


df <- df[as.character(df$V1) %in% method.select,]
rownames(df) <- NULL

# abstract year span from original data
if (NL) {
  d <- read.delim(file = paste0(ROOTDIR, "data/data_processed.csv"), header = TRUE, sep = ",")
  d$year <- as.numeric(format(as.Date(d$Date, format="%d/%m/%Y"),"%Y"))
  #This is only used to get the max value
} else {
  d <- data.frame(year=2020)
  #FIN data extends to 2020.06, predictions extend to 2020.05.
  #NL data extends to 2019.12, predictions extend to 2019.11
  #Hence, remove last 6 timepoints to make error data extend to 2019.11.01 which is the last prediction in NL data
  #df <- df[,- ((ncol(df) -6) : ncol(df))]
  #BUT SHOULD I MAKE THEM MATCH OR NOT!? -> no it is good to show the corona period
}



# number of methods tested for each blood group
m <- length(method.select)
method.names <- array(df[c(1:m),1])
if (m == 13) {
  colors = append(brewer.pal(12, name="Paired"), "#555555")
} else {
  colors = brewer.pal(m, name="Set3")
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

data <- data.frame(matrix(ncol=ncol(df.sums+1),nrow=0))
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
#    png(file= paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", toString(rw), "_", group, ".png"))
  } else {
#    png(file= paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", toString(rw), "_", group, "_FIN.png"))
  }
  plot(1, type = "n", main = paste0(group, ", rw", rw, pd, "\naveraged over each ", merge.months, t), xlim=c((max(d$year)-(n*(merge.months/p))+1), (max(d$year)+1-(merge.months)/p)), ylim=c(0,max(df.plot)), xlab="year", ylab="avg error")
  for (i in c(1:nrow(df.group))) {
    lines(x=seq((max(d$year)-(n*(merge.months/p))+1), (max(d$year)+1-(merge.months)/p), (merge.months/p)), y=df.plot[i,], col=colors[i], type="o", lwd=2, pch = 19)
  }
  legend("topleft", legend=method.names, col=colors, pch = 19)
#  dev.off()
  df.plot$Group <- group
  df.plot$Method <- method.select
  data <- rbind(data,df.plot)
}

####################
# ONCE PLAT PRED DONE UPDATE THIS PART TO WORK WITH RED AND PLAT
####################
data <- as_tibble(data)

names(colors) <- method.select
#colnames(df.plot) <- seq((max(d$year)-(n*(merge.months/p))+1), (max(d$year)+1-(merge.months)/p), (merge.months/p))
colnames(data) <- c(seq((max(d$year)-(n*(merge.months/p))+1), (max(d$year)+1-(merge.months)/p), (merge.months/p)),"Group","Method")

#df.plot$Method <- method.select
data <- data %>% pivot_longer(cols= -c("Method","Group")) %>% mutate(Year=as.numeric(name),
                                                                     Group= recode(Group, RED= "Red cells",PLAT="Platelets"))  
gr <- ggplot(data %>% filter(Group == "Red cells" || Group == "Platelets" ))
gr <- gr + geom_point(aes(x=Year,y=value,color=Method,group=Method),size=4)
gr <- gr + geom_line(aes(x=Year,y=value,color=Method,group=Method),size=1.1)
gr <- gr + scale_colour_manual(name = "Method",values = colors) + ylab("Prediction error \n (% of pcs averaged over each 6 months)")
gr <- gr + theme_bw(base_size = 18)
gr <- gr + theme(legend.position = "bottom", legend.direction = "horizontal")
gr <- gr + facet_wrap(~Group, scales = "free_y") +  theme(axis.text.x = element_text(angle = 90)) 
if (NL) {
  filename <- paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", rw, "_m",merge.months, "_red", "_gg.pdf")
} else {
  filename <- paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", rw, "_m",merge.months, "_red", "_gg_FIN.pdf")
}
ggsave(filename=filename, gr, width = 360,  height = 180,units="mm", dpi=600, scale=1.0)

gr <- gr + facet_grid(Method~Group) +  theme(axis.text.x = element_text(angle = 90)) 

if (NL) {
  filename <- paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", rw, "_m",merge.months,"_l",length(method.select) ,"_red", "_gg_f.pdf")
} else {
  filename <- paste0(ROOTDIR, "rw_testing/img/all_years/", period, "_rwy", rw, "_m",merge.months,"_l",length(method.select) ,"_red", "_gg_f_FIN.pdf")
}
ggsave(filename=filename, gr, width = 360,  height = length(method.select) * 150 ,units="mm", dpi=600, scale=1.0,limitsize = FALSE)



