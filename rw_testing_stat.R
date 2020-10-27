library(hash)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# SETTINGS
NL <- FALSE
#NL <- TRUE


# working directory
if (NL) {
  ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"
  groups <- c("RED", "Ominus", "Oplus", "Aminus", "Aplus", "Bminus", "Bplus", "ABminus", "ABplus", "PLAT")
} else {
  ROOTDIR <- "~/Work/proj/OPERATIONAL/blood_demand_forecast_NL/20201026/"
  groups <- c("RED","PLAT")
}


# "avg" for taking average result of all methods, 
# "best" for taking only result of best performing method
period <- "m"               # m for monthly, w for weakly

# list of all considered rolling window sizes (years)
rolling_windows <- c(3:10)

# put everything in a function in order to efficiently execute for each blood group
fun <- function(method_avg_best, period, group, rolling_windows){
  
  # PRE-PROCESSING
  
  # read errors of first rolling window to determine the set of used models
  df <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rolling_windows[1]), "_all.txt"), header = FALSE, sep = ";")
  index_start <- match(groups[1], toupper(df$V1)) + 1
  #index_stop <- index_start + nrow(df)/10 - 2
  index_stop <- index_start + 12
  used_models <- as.character(df[c(index_start:index_stop),]$V1)
  
  
  # initialize data frame with all prediction errors
  if (period == "m"){
    all_errors <- data.frame(matrix(nrow=length(used_models)*11, ncol=0))
  } else if (period == "w") {
    all_errors <- data.frame(matrix(nrow=length(used_models)*51, ncol=0))
  }
  
  # fill data frame containing prediction errors for all months/weeks for all rolling window sizes
  for (rw in rolling_windows) {
    df <- read.delim(file = paste0(ROOTDIR, "rw_testing/", period, "_errors_rwy", toString(rw), "_all.txt"), header = FALSE, sep = ";")
    index_start <- match(group, toupper(df$V1)) + 1
    #index_stop <- index_start + nrow(df)/10 - 2
    index_stop <- index_start + 12
    #df <- df[c(index_start:index_stop),]
    df <- df[c(index_start:index_stop),c(1,(ncol(df) - 10) : ncol(df))]
    df <- df[df$V1 %in% used_models,]
    cat(group,rw,dim(df),"\n")
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
  
  cat(paste0("Ha: MEANS ARE DIFFERENT", "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
  for (i in c(1:length(rw_comb))){
    cat(paste0(rw_comb[i][[1]][1],",",rw_comb[i][[1]][2], ": "), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
    cat(paste0(t.test(all_errors[,toString(rw_comb[i][[1]][1])], all_errors[,toString(rw_comb[i][[1]][2])], paired = TRUE, alternative = "two.sided")[["p.value"]], "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
  }
  cat("\n", file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
  
  cat(paste0("Ha: FIRST MEAN > SECOND MEAN", "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
  for (i in c(1:length(rw_comb))){
    cat(paste0(rw_comb[i][[1]][1],",",rw_comb[i][[1]][2], ": "), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
    cat(paste0(t.test(all_errors[,toString(rw_comb[i][[1]][1])], all_errors[,toString(rw_comb[i][[1]][2])], paired = TRUE, alternative = "greater")[["p.value"]], "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
  }
  cat("\n", file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
  
  cat(paste0("Ha: FIRST MEAN < SECOND MEAN", "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
  for (i in c(1:length(rw_comb))){
    cat(paste0(rw_comb[i][[1]][1],",",rw_comb[i][[1]][2], ": "), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
    cat(paste0(t.test(all_errors[,toString(rw_comb[i][[1]][1])], all_errors[,toString(rw_comb[i][[1]][2])], paired = TRUE, alternative = "less")[["p.value"]], "\n"), file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
  }
  cat("\n\n", file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), append=TRUE)
}

for (group in groups) {
  fun(method_avg_best, period, group, rolling_windows)
}

