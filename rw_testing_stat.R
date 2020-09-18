library(hash)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggpubr)

# SETTINGS

# working directory
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"

# "avg" for taking average result of all methods, 
# "best" for taking only result of best performing method
method_avg_best <- "best"
period <- "m"               # m for monthly, w for weakly
group <- "RED"              # "RED", "O-", "O+", "A-", "A+", "B-", "B+", "AB-", "AB+", "PLAT"

# list of all considered rolling window sizes (years)
rolling_windows <- c(3:9)

# names and groups
colors = brewer.pal(length(groups), name="Paired")

# PRE-PROCESSING

# read errors of first rolling window to determine the set of used models
df <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rolling_windows[1]), "_all.txt"), header = FALSE, sep = ";")
index_start <- which(grepl(group, toupper(df$V1))) + 1
index_stop <- index_start + nrow(df)/10 - 2
used_models <- df[c(index_start:index_stop),]$V1

# initialize data frame with all prediction errors
if (period == "m"){
  all_errors <- data.frame(matrix(nrow=length(used_models)*11, ncol=0))
} else if (period == "w") {
  all_errors <- data.frame(matrix(nrow=length(used_models)*51, ncol=0))
}

# fill data frame containing prediction errors for all months/weeks for all rolling window sizes
for (rw in rolling_windows) {
  df <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_all.txt"), header = FALSE, sep = ";")
  index_start <- which(grepl(group, toupper(df$V1))) + 1
  index_stop <- index_start + nrow(df)/10 - 2
  df <- df[c(index_start:index_stop),]
  df <- df[df$V1 %in% used_models,]
  
  ls <- list()
  for (r in rownames(df)){
    ls <- append(ls, as.vector(t(df[r,c(2:ncol(df))])))
  }
  all_errors[,ncol(all_errors)+1] <- unlist(ls)
}
colnames(all_errors) <- lapply(rolling_windows, function(a) toString(a))

# list of all combinations of rolling window sizes to be tested
rw_comb <- list()
for (i in rolling_windows[-length(rolling_windows)]){
  for (j in c((i+1):tail(rolling_windows, 1))){
    rw_comb <- append(rw_comb, list(c(i, j)))
  }
}

# T-TESTING

# execute paired t-test on all combinations of rolling window size,
# and write output to rw_testing directory

cat(paste0("DIFFERENT MEAN", "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests_", group, ".txt"), append=TRUE)
for (i in c(1:length(rw_comb))){
  cat(paste0(rw_comb[i][[1]][1],",",rw_comb[i][[1]][2], ": "), file=paste0(ROOTDIR, "rw_testing/t-tests_", group, ".txt"), append=TRUE)
  cat(paste0(t.test(all_errors[,toString(rw_comb[i][[1]][1])], all_errors[,toString(rw_comb[i][[1]][2])], paired = TRUE)[["p.value"]], "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests_", group, ".txt"), append=TRUE)
}
cat("\n", file=paste0(ROOTDIR, "rw_testing/t-tests_", group, ".txt"), append=TRUE)

cat(paste0("FIRST MEAN < SECOND MEAN", "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests_", group, ".txt"), append=TRUE)
for (i in c(1:length(rw_comb))){
  cat(paste0(rw_comb[i][[1]][1],",",rw_comb[i][[1]][2], ": "), file=paste0(ROOTDIR, "rw_testing/t-tests_", group, ".txt"), append=TRUE)
  cat(paste0(t.test(all_errors[,toString(rw_comb[i][[1]][1])], all_errors[,toString(rw_comb[i][[1]][2])], paired = TRUE, alternative = "greater")[["p.value"]], "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests_", group, ".txt"), append=TRUE)
}
cat("\n\n", file=paste0(ROOTDIR, "rw_testing/t-tests_", group, ".txt"), append=TRUE)
