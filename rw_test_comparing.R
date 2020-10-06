library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(hash)

# SETTINGS
ROOTDIR <- "/home/merel/Documents/Sanquin/blood_demand_forecast_NL/"
groups <- c("RED", "O-", "O+", "A-", "A+", "B-", "B+", "AB-", "AB+", "PLAT")
period <- "w"               # m for monthly, w for weakly
rw.all <- c(3:5)            # selection of rolling windows to be tested

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
