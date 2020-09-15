library(hash)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggpubr)

# settings
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/" # Your directory

# list of all considered rolling window sizes (years)
rolling_windows <- c(1:9)

# names and groups
modelnames <- c("SNAIVE", "5-MA", "7-MA", "9-MA", "12-MA", "STL", "ETS", "TBATS", "STLF", "ARIMAX", "DYNREG", "NN", "COMBINED")
used <- c("SNAIVE", "5-MA", "7-MA", "9-MA", "ETS", "TBATS", "NN")
unused <- c("12-MA", "STL", "STLF", "ARIMAX", "DYNREG", "COMBINED")

groups <- c("RED", "Ominus", "Oplus", "Aminus", "Aplus", "Bminus", "Bplus", "ABminus", "ABplus", "PLAT")
colors = brewer.pal(length(groups), name="Paired")

# reshape data frame to a list of all monthly validation forecasts for each rolling window,
# and put them together in a hash set with rolling window size as key
all_errors <- data.frame(matrix(nrow=length(used)*11, ncol=0))
for (rw in rolling_windows) {
  df <- read.delim(file = paste0(ROOTDIR, "rw_testing/m_errors_rwy", toString(rw), "_long.txt"), header = FALSE, sep = ";")
  df <- df[c(2:(length(used)+1)),]
  ls <- list()
  for (r in rownames(df)){
    ls <- append(ls, as.vector(t(df[r,c(2:ncol(df))])))
  }
  all_errors[,ncol(all_errors)+1] <- unlist(ls)
}
colnames(all_errors) <- lapply(rolling_windows, function(a) toString(a))

rw_comb <- list()
for (i in rolling_windows[-tail(rolling_windows, 1)]){
  for (j in c((i+1):tail(rolling_windows, 1))){
    rw_comb <- append(rw_comb, list(c(i, j)))
  }
}

for (i in c(1:length(rw_comb))){
  print(rw_comb[i][[1]][1])
  print(rw_comb[i][[1]][2])
  print(t.test(all_errors[,toString(rw_comb[i][[1]][1])], all_errors[,toString(rw_comb[i][[1]][2])], paired = TRUE, alternative = "greater"))
}
