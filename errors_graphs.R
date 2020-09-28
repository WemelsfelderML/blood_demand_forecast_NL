library(RColorBrewer)
library(stringr)
colors = brewer.pal(12, name="Paired")

# Set working directory
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"

types <- c("Red","O-","O+","A-","A+","B-","B+","AB-","AB+","Plat")
period <- "m"   # 4w, 6m
rolling_windows <- c(3:9)

for (rw in rolling_windows) {
  # Read errors file
  df_errors <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_all.txt"), header = FALSE, sep = ";")
  methods <- df_errors[c(2:14),1]
  
  # Read chosen method file
  df_chosen <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_chosen.txt"), header = FALSE, sep = ";")
  
  # Plot all errors
  for (type in types) {
    start <- as.numeric(rownames(df_errors[df_errors$V1==type,])) + 1
    stop <- start + 12
    errors <- df_errors[c(start:stop),]
    rownames(errors) <- NULL
    chosen <- df_chosen[df_chosen$V1 == type,]$V2
    c <- colors
    index <- as.numeric(rownames(errors[str_to_upper(errors$V1) == str_to_upper(chosen),]))
    if (index == 1){
      c <- append("#000000", c)
    } else if (index == length(c) + 1) {
      c <- append(c, "#000000")
    } else {
      c <- append(append(c[1:(index-1)], "#000000"), c[index:length(c)])
    }
    
    # Plot errors
    png(file= paste0(ROOTDIR, "rw_testing/img/elaborate", period, "_rwy", rw, "_", type , ".png"))
    plot(1, type = "n", main = paste0("rw: ", rw, "\ngroup: ", type, "\nchosen: ", chosen), xlab = "month of validation year", ylab = "error", xlim = c(1,11), ylim = c(0,max(errors[,-1])))
    for (i in c(1:13)) {
      lwd = 1
      if (str_to_upper(errors[i,1]) == str_to_upper(chosen)) {
        lwd = 2
      } 
      lines(x=c(1:11), y=errors[i,-1], type="o", col=c[i], lwd = lwd)
    }
    legend("topright", legend=errors[,1], col=c, pch = 19, inset = c(0.1, 0.1))
    dev.off()
  }
}