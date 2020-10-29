library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(hash)
library(kableExtra)

#First: run forecast_script_FIN_wold.Rmd e.g. echo "rmarkdown::render('forecast_script_FIN_wold.Rmd', clean=TRUE,output_format='html_document',output_file='results/forecast_script_FIN_wold_20201026.html')" | R --slave
#to process data and create predictions namely rw_testing/m_errors* files
#Second: create plots of error with cat errors_allyears_graphs.R | R --slave --args 'rolling_window_for_prediction' 'rolling_window_for_error_averaging' l|s [all methods or just selected]
#Third: use rw_testing_stat.R to create rw comparisons test results
#Fourth: use rw_test_comparing.R to create table of rw test counts

# SETTINGS
NL <- FALSE
#NL <- TRUE
# working directory
if (NL) {
  ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"
  groups <- c("RED", "Ominus", "Oplus", "Aminus", "Aplus", "Bminus", "Bplus", "ABminus", "ABplus", "PLAT")
} else {
  ROOTDIR <- "~/Work/proj/OPERATIONAL/blood_demand_forecast_NL/20201026_ttestcounts//"
  groups <- c("RED","PLAT")
}


period <- "m"               # m for monthly, w for weakly
rw.all <- c(3:10)            # selection of rolling windows to be tested

fun <- function(ROOTDIR, period, group, rw.all){
  # load data
  df <- read.delim(file = paste0(ROOTDIR, "rw_testing/t-tests/", period, "_", group, ".txt"), header = FALSE)
  df <- cbind(df, str_split_fixed(df$V1, "[:,]", 3))[-1]
  colnames(df) <- c("rw1","rw2","pval")
  
  # create data frame with rolling windows and their p-values
  c.names <- df[df$pval == "",]
  c.indices <- append(as.numeric(rownames(c.names)), (nrow(df)+1))
  c.col <- ncol(df)+1
  for (i in c(1:(length(c.indices)-1))){
    if (grepl("DIFFERENT", c.names[i,2])){
      df[c(c.indices[i]:(c.indices[i+1]-1)), c.col] <- "!="
    } else if (grepl(">", c.names[i,2])){
      df[c(c.indices[i]:(c.indices[i+1]-1)), c.col] <- ">"
    } else if (grepl("<", c.names[i,2])){
      df[c(c.indices[i]:(c.indices[i+1]-1)), c.col] <- "<"
    }
    df[c(c.indices[i]:(c.indices[i+1]-1)), (c.col+1)] <- paste0(df[(c.indices[i]+1),1], "-", df[c.indices[i+1]-1,2])
  }
  df <- df[-c.indices,]
  colnames(df)<- c("rw1","rw2","pval", "H.a", "rw.range")
  
  df <- df[df$rw.range == paste0(head(rw.all, 1), "-", tail(rw.all, 1)),]
  
  # filter significant p-values
  df$rw1 <- as.numeric(as.character(df$rw1))
  df$rw2 <- as.numeric(as.character(df$rw2))
  df$pval <- as.numeric(as.character(df$pval))
  df <- df[df$pval < 0.05,]
  
  # make a list of all rolling windows that resulted in significantly smaller errors
  rw.list <- array()
  rw.list <- append(rw.list, array(df[df$H.a == "<", "rw1"]))
  rw.list <- append(rw.list, array(df[df$H.a == ">", "rw2"]))
  rw.list <- rw.list[-1]
  
  return(rw.list)
}

# every time a rolling window is significantly better than
# another, add it to rw.list
rw.hash <- hash()
for (group in groups) {
  rw.hash[group] <- fun(ROOTDIR, period, group, rw.all)
}

# count the number of times for each rw that it was significantly
# better than some other rw, and put the results in a data frame
count <- function(x, n){ length((which(x == n))) }
best.rws <- data.frame(c(1:tail(rw.all, 1)))
for (group in groups){
  rw.list <- rw.hash[[group]]
  group.col <- ncol(best.rws)+1
  for (rw in rw.all) {
    best.rws[rw, group.col] <- count(rw.list, rw)
  }
}
colnames(best.rws) <- append("rw", groups)
best.rws$total = rowSums(best.rws[,c(-1)])
best.rws <- best.rws[complete.cases(best.rws),]

write.csv(best.rws , paste0(ROOTDIR, "rw_testing/t-tests/", period, "_comparison.csv"), row.names = TRUE)
View(best.rws)

#Print latex to file
# library(kableExtra)
# errors <- tibble(
#   Data = "FinDonor",
#   Gender = c("Female", "Male"),
#   Model = "Dynamic linear mixed model",
#   `MAE (g / L)` = c(combined.hfc$mean_errors[1,1], combined.hmc$mean_errors[1,1]),
#   `RMSE (g / L)` = c(combined.hfc$mean_errors[2,1], combined.hmc$mean_errors[2,1]),
#   `MAE (mmol / L)` = c(combined.hfc$mean_errors[1,2], combined.hmc$mean_errors[1,2]),
#   `RMSE (mmol / L)` = c(combined.hfc$mean_errors[2,2], combined.hmc$mean_errors[2,2]),
#   AUC = c(combined.hfc$roc_auc, combined.hmc$roc_auc),
#   AUPR = c(combined.hfc$pr_auc, combined.hmc$pr_auc))
# #FinDonor Prediction errors, average over 4 folds
# as.tibble(best.rws)
# write.csv(errors, file = paste0(table_path, "findonor-errors-both-icp-fix.csv"), row.names = FALSE, quote=FALSE)
if (NL) {
  caption <- "NL: Number of occurrences that a rolling window has lead to significantly lower errors than another rolling window"
  label <- "nl_rw_m_ttesting"
}
if (!NL) {
  caption <- "FIN: Number of occurrences that a rolling window has lead to significantly lower errors than another rolling window"
  label <- "fin_rw_m_ttesting"
}

out <- best.rws[c("rw","RED","PLAT")]
colnames(out) <- c("Rolling window for prediction","Red cells","Platelets")
latex_counts <- kable(out, format="latex", digits=0, caption=caption,label=label,  linesep="",row.names = FALSE,align=c('l','c','c'))
cat(latex_counts,file=paste0(ROOTDIR, "rw_testing/t-tests/", period, "_comparison.latex"),sep="\n")
