library(RColorBrewer)
library(stringr)
colors = brewer.pal(12, name="Paired")

# Set working directory
ROOTDIR <- "/home/merel/Documents/Sanquin/Finland" # Your directory
HISTORYDIR <- paste0(ROOTDIR, "/histories/") # monthly
MONTHLYDIR <- paste0(ROOTDIR, "/histories/monthly_") # monthly
WEEKLYDIR <- paste0(ROOTDIR, "/histories/weekly_") # weekly
setwd(ROOTDIR)

type <- "plat"   #"red","Ominus","Oplus","Aminus","Aplus","Bminus","Bplus","ABminus","ABplus","plat" 
period <- "6m"   # 4w, 6m

# Read errors file
df_errors <- read.delim(file = paste0(HISTORYDIR, period, "_errors_", type, ".txt"), header = FALSE, sep = ";")
rownames(df_errors) <- df_errors$V1
df_errors$V1 <- NULL
names(df_errors) <- c(1:12)

# Read chosen method file
df_chosen <- read.delim(file = paste0(HISTORYDIR, period, "_chosen", ".txt"), header = FALSE, sep = ";")
rownames(df_chosen) <- df_chosen$V1
df_chosen$V1 <- NULL

# Plot errors
plot(1, type = "n", main = paste0(type, ": ", df_chosen[c(type),]), xlab = "month", ylab = "error", xlim = c(1,12), ylim = c(0,max(df_errors)))
for (i in c(1:13)) {
  chosen <- df_chosen[c(type),]
  if (str_to_upper(rownames(df_errors[i,])) == str_to_upper(chosen)) {
    color <- "black"
    lwd = 2
  } else {
    color <- colors[i]
    lwd = 1
  }
  lines(x=c(1:12), y=df_errors[i,], type="o", col=color, lwd = lwd)
}

