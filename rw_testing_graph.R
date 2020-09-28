library(hash)
library(RColorBrewer)
library(stringr)

# SETTINGS

# working directory
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"

# "avg" for taking average result of all methods, 
# "best" for taking only result of best performing method
method_avg_best <- "best"
period <- "w"               # m for monthly, w for weakly

# list of all considered rolling window sizes (years)
rolling_windows <- c(3:7)

# names and groups
modelnames <- c("SNAIVE", "5-MA", "7-MA", "9-MA", "12-MA", "STL", "ETS", "TBATS", "STLF", "ARIMAX", "DYNREG", "NN", "COMBINED")
groups <- c("RED", "Ominus", "Oplus", "Aminus", "Aplus", "Bminus", "Bplus", "ABminus", "ABplus", "PLAT")
colors = brewer.pal(length(groups), name="Paired")

# PRE-PROCESSING

# reshape error data frames and put them together in a hash set with rolling window size as key
h <- hash()
for (rw in rolling_windows) {
  df1 <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_mean.txt"), header = FALSE, sep = ";")
  
  # number of methods tested for each blood group
  m <- nrow(df1)/10 - 1
  df <- data.frame(row.names = df1[c(1:m+1),1])
  
  i <- 1
  j <- 1
  colnames <- list()
  while (i < nrow(df1)){
    df[,j] <- df1[c((i+1):(i+m)),2]
    colnames <- append(colnames, toString(df1[i,1]))
    
    i <- i + m + 1
    j <- j + 1
  }
  colnames(df) <- groups
  
  h[toString(rw)] <- df
}

# PLOTTING

# boxplot of average 12-month errors of all methods, for each blood group separately
if (method_avg_best == "avg"){
  for (group in groups){
    df_plot <- data.frame(matrix(nrow=13, ncol=0))
    for (rw in rolling_windows){
      df_plot[c(1:nrow(h[[toString(rw)]])),ncol(df_plot)+1] <- h[[toString(rw)]][,group]
    }
    colnames(df_plot) <- rolling_windows
    boxplot(df_plot, main = paste0("Errors for ", group, "\nfor different rolling window sizes"), xlab = "rolling window (years)", ylab = "errors red blood cells, all methods")
  }
} 

# line plot of 12-month errors resulting from best method for each blood group
if (method_avg_best == "best"){
  df_plot <- data.frame(matrix(nrow=0, ncol=length(rolling_windows)))
  
  for (group in groups){
    for (rw in rolling_windows){
      df_plot[group, rw] <- min(h[[toString(rw)]][,group])
    }
  }
  df_plot <- df_plot[,rolling_windows]
  colnames(df_plot) <- rolling_windows
  
  plot(1, type = "n", main = paste0("Errors for different rolling window sizes"), xlab = "rolling window (years)", ylab = "errors red blood cells", xlim = c(min(rolling_windows), max(rolling_windows)), ylim = c(0,max(df_plot)))
  for (i in c(1:length(groups))) {
    lines(x=rolling_windows, y=df_plot[groups[i],], col=colors[i], type="o", lwd=2, pch = 19)
  }
  legend("bottomright", legend=groups, col=colors, pch = 19, inset = c(0.1, 0.1))
}

