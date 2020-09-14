library(hash)
library(RColorBrewer)
library(stringr)

# settings
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/" # Your directory
HISTORYDIR <- paste0(ROOTDIR, "histories/") # monthly
MONTHLYDIR <- paste0(ROOTDIR, "histories/monthly_") # monthly
WEEKLYDIR <- paste0(ROOTDIR, "histories/weekly_") # weekly

# "avg" for taking average result of all methods, 
# "best" for taking only result of best performing method
method_avg_best <- "best"

# list of all considered rolling window sizes (years)
rolling_windows <- c(3:9)

# names and groups
modelnames <- c("SNAIVE", "5-MA", "7-MA", "9-MA", "12-MA", "STL", "ETS", "TBATS", "STLF", "ARIMAX", "DYNREG", "NN", "COMBINED")
groups <- c("RED", "Ominus", "Oplus", "Aminus", "Aplus", "Bminus", "Bplus", "ABminus", "ABplus", "PLAT")
colors = brewer.pal(length(groups), name="Paired")

# reshape error data frames and put them together in a hash set with rolling window size as key
h <- hash()
for (rw in rolling_windows) {
  df1 <- read.delim(file = paste0(ROOTDIR, "errors/m_errors_rwy", toString(rw), "_df.txt"), header = FALSE, sep = ";")
  df <- data.frame(row.names = modelnames)
  
  i <- 1
  j <- 1
  colnames <- list()
  while (i < nrow(df1)){
    df[,j] <- df1[c((i+1):(i+13)),2]
    colnames <- append(colnames, toString(df1[i,1]))
    
    i <- i + 14
    j <- j + 1
  }
  colnames(df) <- groups
  
  h[toString(rw)] <- df
}


# make a boxplot of average errors of all methods
if (method_avg_best == "avg"){
  for (group in groups){
    df_plot <- data.frame(matrix(nrow=13, ncol=0))
    for (rw in rolling_windows){
      df_plot[,ncol(df_plot)+1] <- h[[toString(rw)]][,group]
    }
    colnames(df_plot) <- rolling_windows
    boxplot(df_plot, main = paste0("Errors for ", group, "\nfor different rolling window sizes"), xlab = "rolling window (years)", ylab = "errors red blood cells, all methods")
  }
} 

# make a line plot of errors for best method for each blood group
if (method_avg_best == "best"){
  df_plot <- data.frame(matrix(nrow=0, ncol=length(rolling_windows)))
  
  for (group in groups){
    for (rw in rolling_windows){
      df_plot[group, rw] <- min(h[[toString(rw)]][,group])
    }
  }
  df_plot[,c(1,2)] <- NULL
  colnames(df_plot) <- rolling_windows
  
  plot(1, type = "n", main = paste0("Errors for different rolling window sizes"), xlab = "rolling window (years)", ylab = "errors red blood cells", xlim = c(min(rolling_windows), max(rolling_windows)), ylim = c(0,max(df_plot)))
  for (i in c(1:length(groups))) {
    lines(x=rolling_windows, y=df_plot[groups[i],], col=colors[i], type="o", lwd=2, pch = 19)
  }
  legend("bottomleft", legend=groups, col=colors, pch = 19, inset = c(0.1, 0.1))
}

