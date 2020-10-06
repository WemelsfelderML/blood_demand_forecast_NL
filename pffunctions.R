# Helper functions for forecast_script.Rmd

extract_type <- function(red.distr, type){
  # Separate blood types from the master frame.
  # 
  # param:: red.distr: The master frame. Distribution of red cell products.
  # param:: type: 3-length string of the desired blood type, e.g. "AB+", "O -".
  # return:: Dataframe of the desired blood type.
  
  # Create a full sequence of dates for imputation purposes
  all.dates <- (seq.Date(min(red.distr$date),
                         max(red.distr$date),
                         "day"))
  # Find type
  typed <- red.distr[red.distr$ABO == type, ]
  typed.distr <- aggregate(typed$quantity, by = list(typed$date), sum); colnames(typed.distr) <- c("date", "pcs")
  # Merge into a whole set with NAs
  typed.distr <- merge(x = data.frame(date = all.dates),
                       y = typed.distr,
                       all.x = TRUE)
  # Replace with zeroes
  typed.distr[is.na(typed.distr)] <- 0
  # Cut to time after 2014
  # typed.distr <- typed.distr[typed.distr$date >= as.Date("2014-01-01"), ]
  
  return(typed.distr)
}

aggregate_weekly <- function(series){
  # Aggregate distributions into weekly sums.
  # 
  # param:: series: dataframe of a daily series
  # return:: dataframe of a weekly series
  
  pcslist <- list(); weeklist <- list(); datelist <- list()  # Create storage lists
  j <- 0; k <- 0  # j is for resetting week counter, k is for gathering dates
  for(i in seq(to = (length(series$date) - 7), by = 7)){
    if(j == 52){j <- 0}
    j <- j + 1
    k <- k + 1
    pcslist[[k]] <- sum(series$pcs[i : (i + 6)])
    weeklist[[k]] <- j
    datelist[[k]] <- series$date[i]
  }
  weekly <- data.frame(week = unlist(weeklist), date = as.Date.numeric(unlist(datelist), origin = "1970-01-01"), pcs = unlist(pcslist))
  return(weekly)
}


find_errors <- function(method = "none", rw_years, period, train, horizon.m, month.m){
  
  # Fit & Forecast
  if(method == "snaive"){
    fit <- snaive(train)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "5-MA"){
    fit <- ma(train, order = 5)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "7-MA"){
    fit <- ma(train, order = 7)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "9-MA"){
    fit <- ma(train, order = 9)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "12-MA"){
    fit <- ma(train, order = 12)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "stl"){
    fit <- stl(train, s.window = "periodic", t.window = 7)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "ets"){
    fit <- ets(train)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "tbats"){
    fit <- tbats(train)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "stlf"){
    fit <- stlf(train)
    fcast <- forecast(fit, h = 1)$mean
  }
  if(method == "arimax"){
    fit <- auto.arima(train, xreg = month.m)
    fcast <- forecast(fit, xreg = horizon.m, h = 1)$mean
  }
  if(method == "dynreg"){
    fit1 <- tslm(train ~ trend + season)
    fcast1 <- forecast(fit1, h = 1)
    
    fit2 <- auto.arima(fit1$residuals)
    fcast2 <- forecast(fit2, h = 1)
    
    y <- as.numeric(fcast1$mean)
    x <- as.numeric(fcast2$mean)
    fcast <- x + y
  }
  if(method == "nn"){
    fit <- nnetar(train)
    fcast <- forecast(fit, h = 1)$mean
  }
    
  return(fcast)	
}	

select_model <- function(beginning, series.ts, freq, rw_years){	
  
  methods <- c("snaive", "5-MA", "7-MA", "9-MA", "12-MA", "stl", "ets", "tbats", "stlf", "arimax", "dynreg", "nn", "combined")
  apes.hash <- hash()
  
  if (freq == "monthly") {
    period <- "m"
    period.num <- 12
  } else if (freq == "weekly") {
    period <- "w"
    period.num <- 52
  }
  
  for(i in seq(length(series.ts) - (rw_years * period.num + 1))){
    # Define training and testing set as ROLLING WINDOW
    
    train <- ts(series.ts[i:((rw_years * period.num) - 1 + i)], frequency = period.num)
    test <- ts(series.ts[((rw_years * period.num) + i)], frequency = period.num)
    
    # xregs for arimax
    nrow <- length(train)
    onehotyear <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         ncol = 11,
                         byrow = TRUE)
    month.m <- matrix(ncol = 11, nrow = nrow)
    
    if (freq == "monthly") {
      for(j in 1:nrow){
        month.m[j,] <- onehotyear[month(beginning + months(j - 1) + months(i - 1)),]
      }
      horizon.m <- matrix(onehotyear[month(beginning + months(nrow + i - 1)), ], ncol = 11)
    } else if (freq == "weekly") {
      for(j in 1:length(train)){
        month.m[j,] <- onehotyear[month((beginning + weeks(j - 1)) + weeks(i - 1)), ]
      }
      horizon.m <- matrix(onehotyear[month(beginning + weeks(nrow + i - 1)), ], ncol = 11)
    }
    
    fcasts.comb <- c()
    
    for (method in c("snaive", "5-MA", "7-MA", "9-MA", "ets", "tbats", "nn")) {
      fcast <- find_errors(method, rw_years, period, train, horizon.m, month.m)
      fcasts.comb <- append(fcasts.comb, fcast)
      ape <- abs(as.numeric(test) - as.numeric(fcast))/as.numeric(test) * 100
      apes <- c(apes.hash[[method]], ape)
      apes.hash[method] <- apes
    }
    
    for (method in c("12-MA", "arimax", "dynreg")) {
      if (rw_years >= 2) {
        fcast <- find_errors(method, rw_years, period, train, horizon.m, month.m)
        fcasts.comb <- append(fcasts.comb, fcast)
        ape <- abs(as.numeric(test) - as.numeric(fcast))/as.numeric(test) * 100
        apes <- c(apes.hash[[method]], ape)
      } else {
        apes <- c(9000:(9000+length(apes.hash[["snaive"]][1])-1))	
      }
      apes.hash[method] <- apes
    }
    
    for (method in c("stl", "stlf")) {
      if (rw_years >= 3) {
        fcast <- find_errors(method, rw_years, period, train, horizon.m, month.m)
        fcasts.comb <- append(fcasts.comb, fcast)
        ape <- abs(as.numeric(test) - as.numeric(fcast))/as.numeric(test) * 100
        apes <- c(apes.hash[[method]], ape)
      } else {
        apes <- c(9000:(9000+length(apes.hash[["snaive"]])-1))	
      }
      apes.hash[method] <- apes
    }
    
    for (method in c("combined")) {
      if (rw_years >= 3) {
        fcast <- mean(fcasts.comb)
        ape <- abs(as.numeric(test) - as.numeric(fcast))/as.numeric(test) * 100
        apes <- c(apes.hash[[method]], ape)
      } else {
        apes <- c(9000:(9000+length(apes.hash[["snaive"]][1])-1))	
      }
      apes.hash[method] <- apes
    }
  
  }
  
  for (method in methods) {
    apes <- apes.hash[[method]]
    cat(str_replace_all(paste0(method, ";", mean(apes), "\n")," ",""), file=paste0("rw_testing/", period, "_errors_rwy", rw_years, "_mean.txt"), append=TRUE)
    cat(str_replace_all(paste0(method, ";", str_replace_all(toString(apes),", ",";"), "\n")," ",""), file=paste0("rw_testing/", period, "_errors_rwy", rw_years, "_all.txt"), append=TRUE)
  }
  
  all.apes <- c()
  for (method in methods) {
    all.apes <- append(all.apes, apes.hash[[method]])
  }
  
  m <- matrix(all.apes, ncol = 13, byrow = FALSE)
  mdf <- as.data.frame(m)
  chosen.model <- which.min(colMeans(mdf))[[1]]
  
  return(chosen.model)
}

chosen_forecast <- function(chosen.model, series.ts, monthly, freq = "monthly", rw_years){
  
  if(freq == "monthly"){  # Monthly forecast
    train <- ts(tail(series.ts, (rw_years * 12)), start = decimal_date(head(tail(monthly$date, (rw_years * 12)), 1)), frequency = 12)  # 3 year window
    # xregs for arimax
    nrow <- length(train)
    onehotyear <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         ncol = 11,
                         byrow = TRUE)
    
    month.m <- matrix(, nrow = nrow, ncol = 11)
    for(j in 1:nrow){
      month.m[j,] <- onehotyear[month(head(tail(monthly$date, (rw_years * 12)), 1) + months(j - 1)), ]
    }
    horizon.m <- matrix(, nrow = 6, ncol = 11)
    k <- 0
    for(j in nrow:(nrow + 5)){
      k <- k + 1
      horizon.m[k,] <- onehotyear[month(head(tail(monthly$date, (rw_years * 12)), 1) + months(j)), ]
    }
    # Fit model
    if(chosen.model == 1){fit <- snaive(train)}
    if(chosen.model == 2){fit <- ma(train, order = 5)}
    if(chosen.model == 3){fit <- ma(train, order = 7)}
    if(chosen.model == 4){fit <- ma(train, order = 9)}
    if(chosen.model == 5){fit <- ma(train, order = 12)}
    if(chosen.model == 6){fit <- stl(train, s.window = "periodic", t.window = 7)}
    if(chosen.model == 7){fit <- ets(train)}
    if(chosen.model == 8){fit <- tbats(train)}
    if(chosen.model == 9){fit <- stlf(train)}
    if(chosen.model == 10){
      fit <- auto.arima(train, xreg = month.m)
      forecasts <- forecast(fit, xreg = horizon.m, h = 6)$mean
    }
    if(chosen.model == 11){
      fit1 <- tslm(train ~ trend + season)
      fcast1 <- forecast(fit1, h = 6)
      
      fit2 <- auto.arima(fit1$residuals)
      fcast2 <- forecast(fit2, h = 6)
      
      y <- as.numeric(fcast1$mean)
      x <- as.numeric(fcast2$mean)
      forecasts <- (x + y)
      upper80s <- fcast1$upper[[1]]
      upper95s <- fcast1$upper[[2]]
      lower80s <- fcast1$lower[[1]]
      lower95s <- fcast1$lower[[2]]
    }
    if(chosen.model == 12){
      fit <- nnetar(train)
      forecasts <- forecast(fit, h = 6)$mean
    }
    if(chosen.model == 13){
      # Seasonal naive
      snaive.fit <- snaive(train)
      snaive.fcast <- forecast(snaive.fit, h = 6)$mean
      
      # 5-MA
      ma5.fit <- ma(train, order = 5)
      ma5.fcast <- forecast(ma5.fit, h = 6)$mean
      
      # 7-MA
      ma7.fit <- ma(train, order = 7)
      ma7.fcast <- forecast(ma7.fit, h = 6)$mean
      
      # 9-MA
      ma9.fit <- ma(train, order = 9)
      ma9.fcast <- forecast(ma9.fit, h = 6)$mean
      
      # 12-MA
      ma12.fit <- ma(train, order = 12)
      ma12.fcast <- forecast(ma12.fit, h = 6)$mean
      
      # STL
      stl.fit <- stl(train, s.window = "periodic", t.window = 7)
      stl.fcast <- forecast(stl.fit, h = 6)$mean
      
      # ETS
      ets.fit <- ets(train)
      ets.fcast <- forecast(ets.fit, h = 6)$mean
      
      # TBATS
      tbats.fit <- tbats(train)
      tbats.fcast <- forecast(tbats.fit, h = 6)$mean
      
      # STLF
      stlf.fit <- stlf(train)
      stlf.fcast <- forecast(stlf.fit, h = 6)$mean
      
      # Arimax
      arimax.fit <- auto.arima(train, xreg = month.m)
      arimax.fcast <- forecast(arimax.fit, xreg = horizon.m, h = 6)$mean
      
      # Dynamic regression
      fit1 <- tslm(train ~ trend + season)
      fcast1 <- forecast(fit1, h = 6)
      
      fit2 <- auto.arima(fit1$residuals)
      fcast2 <- forecast(fit2, h = 6)
      
      y <- as.numeric(fcast1$mean)
      x <- as.numeric(fcast2$mean)
      dynreg.fcast <- x + y
      
      # NNETAR
      nn.fit <- nnetar(train)
      nn.fcast <- forecast(nn.fit, h = 6)$mean
      
      forecasts <- colMeans(matrix(c(snaive.fcast, ma5.fcast, ma7.fcast, ma9.fcast, ma12.fcast, stl.fcast,
                                     ets.fcast, tbats.fcast, stlf.fcast, arimax.fcast, dynreg.fcast, nn.fcast), 
                                   ncol = 6, 
                                   byrow = TRUE))
    }
    if(!(chosen.model %in% c(10, 11, 12, 13))){
      # If model was not any of the above, forecast here
      fcast <- forecast(fit, h = 6)
      forecasts <- as.numeric(fcast$mean)
      upper80s <- as.numeric(fcast$upper[, 1])
      upper95s <- as.numeric(fcast$upper[, 2])
      lower80s <- as.numeric(fcast$lower[, 1])
      lower95s <- as.numeric(fcast$lower[, 2])
    }
    
    # Save to a returnable dataframe
    if(!(chosen.model %in% c(10, 12, 13))){
      fdf <- data.frame(fcast = forecasts,
                        upper80 = upper80s,
                        upper95 = upper95s,
                        lower80 = lower80s,
                        lower95 = lower95s)
    } else { # Because models 10, 12, and 13 don't have CI
      fdf <- data.frame(fcast = forecasts,
                        upper80 = forecasts,
                        upper95 = forecasts,
                        lower80 = forecasts,
                        lower95 = forecasts)
    }
  }
  if(freq == "weekly"){  # Weekly forecast
    train <- ts(tail(series.ts, (rw_years * 52)), start = decimal_date(head(tail(monthly$date, (rw_years * 52)), 1)), frequency = 52)  # Three year window
    # xregs for arimax
    nrow <- length(train)
    onehotyear <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         ncol = 11,
                         byrow = TRUE)
    
    month.m <- matrix(, nrow = nrow, ncol = 11)
    for(j in 1:nrow){
      month.m[j,] <- onehotyear[month(head(tail(monthly$date, (rw_years * 52)), 1) + weeks(j - 1)), ]
    }
    horizon.m <- matrix(, nrow = 4, ncol = 11)
    k <- 0
    for(j in nrow:(nrow + 3)){
      k <- k + 1
      horizon.m[k,] <- onehotyear[month(head(tail(monthly$date, (rw_years * 52)), 1) + weeks(j)), ]
    }
    # Fit model
    if(chosen.model == 1){fit <- snaive(train)}
    if(chosen.model == 2){fit <- ma(train, order = 5)}
    if(chosen.model == 3){fit <- ma(train, order = 7)}
    if(chosen.model == 4){fit <- ma(train, order = 9)}
    if(chosen.model == 5){fit <- ma(train, order = 12)}
    if(chosen.model == 6){fit <- stl(train, s.window = "periodic", t.window = 7)}
    if(chosen.model == 7){fit <- ets(train)}
    if(chosen.model == 8){fit <- tbats(train)}
    if(chosen.model == 9){fit <- stlf(train)}
    if(chosen.model == 10){
      fit <- auto.arima(train, xreg = month.m)
      forecasts <- forecast(fit, xreg = horizon.m, h = 4)$mean
    }
    if(chosen.model == 11){
      fit1 <- tslm(train ~ trend + season)
      fcast1 <- forecast(fit1, h = 4)
      
      fit2 <- auto.arima(fit1$residuals)
      fcast2 <- forecast(fit2, h = 4)
      
      y <- as.numeric(fcast1$mean)
      x <- as.numeric(fcast2$mean)
      forecasts <- (x + y)
      upper80s <- fcast1$upper[[1]]
      upper95s <- fcast1$upper[[2]]
      lower80s <- fcast1$lower[[1]]
      lower95s <- fcast1$lower[[2]]
    }
    if(chosen.model == 12){
      fit <- nnetar(train)
      forecasts <- forecast(fit, h = 4)$mean
    }
    if(chosen.model == 13){
      # Seasonal naive
      snaive.fit <- snaive(train)
      snaive.fcast <- forecast(snaive.fit, h = 4)$mean
      
      # 5-MA
      ma5.fit <- ma(train, order = 5)
      ma5.fcast <- forecast(ma5.fit, h = 4)$mean
      
      # 7-MA
      ma7.fit <- ma(train, order = 7)
      ma7.fcast <- forecast(ma7.fit, h = 4)$mean
      
      # 9-MA
      ma9.fit <- ma(train, order = 9)
      ma9.fcast <- forecast(ma9.fit, h = 4)$mean
      
      # 12-MA
      ma12.fit <- ma(train, order = 12)
      ma12.fcast <- forecast(ma12.fit, h = 4)$mean
      
      # STL
      stl.fit <- stl(train, s.window = "periodic", t.window = 7)
      stl.fcast <- forecast(stl.fit, h = 4)$mean
      
      # ETS
      ets.fit <- ets(train)
      ets.fcast <- forecast(ets.fit, h = 4)$mean
      
      # TBATS
      tbats.fit <- tbats(train)
      tbats.fcast <- forecast(tbats.fit, h = 4)$mean
      
      # STLF
      stlf.fit <- stlf(train)
      stlf.fcast <- forecast(stlf.fit, h = 4)$mean
      
      # Arimax
      arimax.fit <- auto.arima(train, xreg = month.m)
      arimax.fcast <- forecast(arimax.fit, xreg = horizon.m, h = 4)$mean
      
      # Dynamic regression
      fit1 <- tslm(train ~ trend + season)
      fcast1 <- forecast(fit1, h = 4)
      
      fit2 <- auto.arima(fit1$residuals)
      fcast2 <- forecast(fit2, h = 4)
      
      y <- as.numeric(fcast1$mean)
      x <- as.numeric(fcast2$mean)
      dynreg.fcast <- x + y
      
      # NNETAR
      nn.fit <- nnetar(train)
      nn.fcast <- forecast(nn.fit, h = 4)$mean
      
      print(matrix(c(snaive.fcast, ma5.fcast, ma7.fcast, ma9.fcast, ma12.fcast, stl.fcast,	
                     ets.fcast, tbats.fcast, stlf.fcast, arimax.fcast, dynreg.fcast, nn.fcast), 	
                   ncol = 4, 	
                   byrow = TRUE))
      
      forecasts <- colMeans(matrix(c(snaive.fcast, ma5.fcast, ma7.fcast, ma9.fcast, ma12.fcast, stl.fcast,
                                     ets.fcast, tbats.fcast, stlf.fcast, arimax.fcast, dynreg.fcast, nn.fcast), 
                                   ncol = 4, 
                                   byrow = TRUE))
    }
    if(!(chosen.model %in% c(10, 11, 12, 13))){
      # If model was not any of the above, forecast here
      fcast <- forecast(fit, h = 4)
      forecasts <- as.numeric(fcast$mean)
      upper80s <- as.numeric(fcast$upper[, 1])
      upper95s <- as.numeric(fcast$upper[, 2])
      lower80s <- as.numeric(fcast$lower[, 1])
      lower95s <- as.numeric(fcast$lower[, 2])
    }
    # Save to a returnable dataframe
    if(!(chosen.model %in% c(10, 12, 13))){
      fdf <- data.frame(fcast = forecasts,
                        upper80 = upper80s,
                        upper95 = upper95s,
                        lower80 = lower80s,
                        lower95 = lower95s)
    } else {
      fdf <- data.frame(fcast = forecasts,
                        upper80 = forecasts,
                        upper95 = forecasts,
                        lower80 = forecasts,
                        lower95 = forecasts)
    }
  }
  return(fdf)
}

# An extended forecasting function is needed to safely create longer forecast sequences. It's not efficient per se,
# but keeps everything else intact.
# Forecasted periods are changed from 6 to 23 (to accommodate two whole years). The actual required length is determined in
# forecast_script.Rmd
chosen_forecast_extended <- function(chosen.model, series.ts, monthly, freq = "monthly", rw_years){	
  if(freq == "monthly"){  # Monthly forecast	
    train <- ts(tail(series.ts, (rw_years * 12)), start = decimal_date(head(tail(monthly$date, (rw_years * 12)), 1)), frequency = 12)  # 3 year window
    # xregs for arimax
    nrow <- length(train)
    onehotyear <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         ncol = 11,
                         byrow = TRUE)
    
    month.m <- matrix(, nrow = nrow, ncol = 11)
    for(j in 1:length(train)){
      month.m[j,] <- onehotyear[month(head(tail(monthly$date, (rw_years * 12)), 1) + months(j - 1)), ]
    }
    horizon.m <- matrix(, nrow = 23, ncol = 11)
    k <- 0
    for(j in length(train):(length(train) + 22)){
      k <- k + 1
      horizon.m[k,] <- onehotyear[month(head(tail(monthly$date, (rw_years * 12)), 1) + months(j)), ]
    }
    # Fit model
    if(chosen.model == 1){fit <- snaive(train)}
    if(chosen.model == 2){fit <- ma(train, order = 5)}
    if(chosen.model == 3){fit <- ma(train, order = 7)}
    if(chosen.model == 4){fit <- ma(train, order = 9)}
    if(chosen.model == 5){fit <- ma(train, order = 12)}
    if(chosen.model == 6){fit <- stl(train, s.window = "periodic", t.window = 7)}
    if(chosen.model == 7){fit <- ets(train)}
    if(chosen.model == 8){fit <- tbats(train)}
    if(chosen.model == 9){fit <- stlf(train)}
    if(chosen.model == 10){
      fit <- auto.arima(train, xreg = month.m)
      forecasts <- forecast(fit, xreg = horizon.m, h = 23)$mean
    }
    if(chosen.model == 11){
      fit1 <- tslm(train ~ trend + season)
      fcast1 <- forecast(fit1, h = 23)
      
      fit2 <- auto.arima(fit1$residuals)
      fcast2 <- forecast(fit2, h = 23)
      
      y <- as.numeric(fcast1$mean)
      x <- as.numeric(fcast2$mean)
      forecasts <- (x + y)
      upper80s <- fcast1$upper[[1]]
      upper95s <- fcast1$upper[[2]]
      lower80s <- fcast1$lower[[1]]
      lower95s <- fcast1$lower[[2]]
    }
    if(chosen.model == 12){
      fit <- nnetar(train)
      forecasts <- forecast(fit, h = 23)$mean
    }
    if(chosen.model == 13){
      # Seasonal naive
      snaive.fit <- snaive(train)
      snaive.fcast <- forecast(snaive.fit, h = 23)$mean
      
      # 5-MA
      ma5.fit <- ma(train, order = 5)
      ma5.fcast <- forecast(ma5.fit, h = 23)$mean
      
      # 7-MA
      ma7.fit <- ma(train, order = 7)
      ma7.fcast <- forecast(ma7.fit, h = 23)$mean
      
      # 9-MA
      ma9.fit <- ma(train, order = 9)
      ma9.fcast <- forecast(ma9.fit, h = 23)$mean
      
      # 12-MA
      ma12.fit <- ma(train, order = 12)
      ma12.fcast <- forecast(ma12.fit, h = 23)$mean
      
      # STL
      stl.fit <- stl(train, s.window = "periodic", t.window = 7)
      stl.fcast <- forecast(stl.fit, h = 23)$mean
      
      # ETS
      ets.fit <- ets(train)
      ets.fcast <- forecast(ets.fit, h = 23)$mean
      
      # TBATS
      tbats.fit <- tbats(train)
      tbats.fcast <- forecast(tbats.fit, h = 23)$mean
      
      # STLF
      stlf.fit <- stlf(train)
      stlf.fcast <- forecast(stlf.fit, h = 23)$mean
      
      # Arimax
      arimax.fit <- auto.arima(train, xreg = month.m)
      arimax.fcast <- forecast(arimax.fit, xreg = horizon.m, h = 23)$mean
      
      # Dynamic regression
      fit1 <- tslm(train ~ trend + season)
      fcast1 <- forecast(fit1, h = 23)
      
      fit2 <- auto.arima(fit1$residuals)
      fcast2 <- forecast(fit2, h = 23)
      
      y <- as.numeric(fcast1$mean)
      x <- as.numeric(fcast2$mean)
      dynreg.fcast <- x + y
      
      # NNETAR
      nn.fit <- nnetar(train)
      nn.fcast <- forecast(nn.fit, h = 23)$mean
      
      forecasts <- colMeans(matrix(c(snaive.fcast, ma5.fcast, ma7.fcast, ma9.fcast, ma12.fcast, stl.fcast,
                                     ets.fcast, tbats.fcast, stlf.fcast, arimax.fcast, dynreg.fcast, nn.fcast), 
                                   ncol = 23, 
                                   byrow = TRUE))
    }
    if(!(chosen.model %in% c(10, 11, 12, 13))){
      # If model was not any of the above, forecast here
      fcast <- forecast(fit, h = 23)
      
      forecasts <- as.numeric(fcast$mean)
      upper80s <- as.numeric(fcast$upper[, 1])
      upper95s <- as.numeric(fcast$upper[, 2])
      lower80s <- as.numeric(fcast$lower[, 1])
      lower95s <- as.numeric(fcast$lower[, 2])
    }
    
    # Save to a returnable dataframe
    if(!(chosen.model %in% c(10, 12, 13))){
      fdf <- data.frame(fcast = forecasts,
                        upper80 = upper80s,
                        upper95 = upper95s,
                        lower80 = lower80s,
                        lower95 = lower95s)
    } else { # Because some models don't have CI
      fdf <- data.frame(fcast = forecasts,
                        upper80 = forecasts,
                        upper95 = forecasts,
                        lower80 = forecasts,
                        lower95 = forecasts)
    }
  }
  
  return(fdf)
}

draw_forecast <- function(forecast_dataframe, freq, history, type, modelname, palette){
  if(freq == "monthly"){
    ggplot(data = forecast_dataframe, aes(x = seq.Date(from = tail(history$date, 1), to = tail(history$date, 1) + months(5), by = "month"), y = fcast)) + 
      geom_segment(x = tail(history$date, 2)[1], y = tail(history[,type], 2)[1], xend = tail(history$date, 1), yend = head(forecast_dataframe$fcast, 1), color = "#D55E00", linetype = "dashed", alpha = palette[["alphaSeg"]]) +
      geom_ribbon(data = forecast_dataframe, aes(ymin = lower95, ymax = upper95), alpha = palette[["alpha95"]], fill = palette[["fill95"]]) + 
      geom_ribbon(data = forecast_dataframe, aes(ymin = lower80, ymax = upper80), alpha = palette[["alpha80"]], fill = palette[["fill80"]]) +
      geom_line(aes(colour = "Forecast"), size = 0.8) +
      geom_line(data = tail(head(history, -1), 24), aes(x = date, y = get(type), colour = palette[["colData"]]), size = 0.8, alpha = palette[["alphaSeg"]] ) +
      
      labs(title = paste0("6 month forecast: ", type),
           subtitle = paste0("Best performing model of the previous year: ", modelname),
           x = "",
           y = "Units per month",
           colour = "Series: ") +
      theme_minimal() + 
      theme(legend.position = "bottom", legend.margin = margin(t = -20, b = 0)) + 
      scale_colour_manual(values = c(palette[["colData"]], palette[["colPred"]]), labels = c("Data", "Forecast"))
  } else {
    ggplot(data = forecast_dataframe, aes(x = seq.Date(from = tail(history$date, 1), to = tail(history$date, 1) + weeks(3), by = "week"), y = fcast)) + 
      geom_segment(x = tail(history$date, 2)[1], y = tail(history[,type], 2)[1], xend = tail(history$date, 1), yend = head(forecast_dataframe$fcast, 1), color = "#D55E00", linetype = "dashed", alpha = palette[["alphaSeg"]]) +
      geom_ribbon(data = forecast_dataframe, aes(ymin = lower95, ymax = upper95), alpha = palette[["alpha95"]], fill = palette[["fill95"]]) + 
      geom_ribbon(data = forecast_dataframe, aes(ymin = lower80, ymax = upper80), alpha = palette[["alpha80"]], fill = palette[["fill80"]]) +
      geom_line(aes(colour = "Forecast"), size = 0.8) +
      geom_line(data = tail(head(history, -1), 16), aes(x = date, y = get(type), colour = palette[["colData"]]), size = 0.8, alpha = palette[["alphaSeg"]] ) +
      
      labs(title = paste0("4 week forecast: ", type),
           subtitle = paste0("Best performing model of the previous year: ", modelname),
           x = "",
           y = "Units per week",
           colour = "Series: ") +
      theme_minimal() + 
      theme(legend.position = "bottom", legend.margin = margin(t = -20, b = 0)) + 
      scale_colour_manual(values = c(palette[["colData"]], palette[["colPred"]]), labels = c("Data", "Forecast"))
  }
  
}

save_forecast <- function(fdf, months = TRUE, modelname, history, file, reverse_adj, rw_years){
  # We will need these stored
  modelnames <- c("SNAIVE", "5-MA", "7-MA", "9-MA", "12-MA", "STL", "ETS", "TBATS", "STLF", "ARIMAX", "DYNREG", "NN", "COMBINED")
  
  # Extract type from filename
  if(months){
    type <- regmatches(file, regexec("monthly_(.*)_rwy[0-9].csv", file))[[1]][2]
  } else{
    type <- regmatches(file, regexec("weekly_(.*)_rwy[0-9].csv", file))[[1]][2]
  }
  
  
  if(months){
    # Prepare to save the new forecast into file
    save.this <- data.frame(time = tail(history$date, 1),
                            model = modelname, 
                            forecast = fdf$fcast[1] * reverse_adj[1],
                            upper80 = fdf$upper80[1] * reverse_adj[1],
                            upper95 = fdf$upper95[1] * reverse_adj[1],
                            lower80 = fdf$lower80[1] * reverse_adj[1],
                            lower95 = fdf$lower95[1] * reverse_adj[1], stringsAsFactors = F) 
  } else{
    save.this <- data.frame(time = tail(history$date, 1), 
                            model = modelname, 
                            forecast = fdf$fcast[1],
                            upper80 = fdf$upper80[1],
                            upper95 = fdf$upper95[1],
                            lower80 = fdf$lower80[1],
                            lower95 = fdf$lower95[1], stringsAsFactors = F) 
  }
  
  # Does forecast history file exist?
  if(file.exists(file)){
    fcast.history <- read.csv(file, header = TRUE, stringsAsFactors = F)
    fcast.history$time <- as.Date(fcast.history$time)
    
    if(months){
      # Check if missing dates in forecast history
      all.dates <- seq.Date(min(fcast.history$time, na.rm=T), max(fcast.history$time, na.rm=T), "month") # Build a sequence of all dates
      missing <- all.dates[!all.dates %in% fcast.history$time] # Check if some are missing
    } else{
      # Check if missing dates in forecast history
      all.dates <- seq.Date(min(fcast.history$time, na.rm=T), max(fcast.history$time, na.rm=T), "week") # Build a sequence of all dates
      missing <- all.dates[!all.dates %in% fcast.history$time] # Check if some are missing
    }
    
    if(length(missing) == 0){
      # Save newest forecast as usual
      if(save.this$time %in% fcast.history$time){
        warning("Data point already in file!") # Do not save if a forecast exists already!
      } else{
        write.table(save.this, file = file, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    } else{
      # Fix missing dates
      dat2 <- data.frame(time = all.dates, stringsAsFactors = F)
      fixed.fhistory <- merge(fcast.history, dat2, all = TRUE)
      
      # Iterate through missing days and compute forecasts for them, then insert to fixed.fhistory
      if(months){
        for(mdate in missing){
          cutreal<- history[history$date <= mdate, ] # Cut real history at the missing date
          beginning <- head(tail(cutreal$date, (((rw_years+1)*12) + 1)), 1) # Define beginning of 4 year window 
          segment <- head(tail(cutreal[, type], (((rw_years+1)*12) + 1)), ((rw_years+1)*12)) # Define 4 year window
          series.ts <- ts(segment, start = decimal_date(beginning), frequency = 12) # Transform into a ts object
          
          scaleback <- as.numeric(bizdays(ts(seq(1), start = decimal_date(cutreal$date[length(cutreal$date)]), frequency = 12), FinCenter = "Zurich")) # Scaler for saving
          
          # Forecast
          chosen.model <- select_model(beginning, series.ts, freq = "monthly", rw_years) # Choose model
          mcast <- chosen_forecast(chosen.model, series.ts, history, freq = "monthly") # Output a forecast
          missing.fcast <- data.frame(time = tail(cutreal$date, 1), 
                                      model = modelnames[chosen.model], 
                                      forecast = mcast$fcast[1] * scaleback,
                                      upper80 = mcast$upper80[1] * scaleback,
                                      upper95 = mcast$upper95[1] * scaleback,
                                      lower80 = mcast$lower80[1] * scaleback,
                                      lower95 = mcast$lower95[1] * scaleback, stringsAsFactors = F) 
          
          fixed.fhistory[fixed.fhistory$time == mdate, ] <- missing.fcast
          
        }
      } else{
        for(mdate in missing){
          cutreal<- history[history$date <= mdate, ] # Cut real history at the missing date
          wbeginning <- head(tail(cutreal$date, (((rw_years+1)*52) + 1)), 1) # Define beginning of 4 year window 
          segment <- head(tail(cutreal[, type], (((rw_years+1)*52) + 1)), ((rw_years+1)*52)) # Define 4 year window
          series.ts <- ts(segment, start = decimal_date(wbeginning), frequency = 52) # Transform into a ts object
          
          # Forecast
          chosen.model <- select_model(beginning, series.ts, freq = "weekly", rw_years) # Choose model
          mcast <- chosen_forecast(chosen.model, series.ts, history, freq = "weekly", rw_years) # Output a forecast
          missing.fcast <- data.frame(time = tail(cutreal$date, 1), 
                                      model = modelnames[chosen.model], 
                                      forecast = mcast$fcast[1],
                                      upper80 = mcast$upper80[1],
                                      upper95 = mcast$upper95[1],
                                      lower80 = mcast$lower80[1],
                                      lower95 = mcast$lower95[1], stringsAsFactors = F) 
          
          fixed.fhistory[fixed.fhistory$time == mdate, ] <- missing.fcast
        }
      }
      
      # First save the fixed forecast history into a file
      write.csv(fixed.fhistory, file = file, row.names = FALSE)
      
      # Then save newest forecast as usual
      if(save.this$time %in% fcast.history$time){
        warning("Data point already in file!") # Do not save if a forecast exists already!
      } else{
        write.table(save.this, file = file, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
    
  } else{
    # If history is missing entirely, we choose to generate 2 years for monthly fcasts and 6 months for weekly fcasts (required for plotting)
    if(months){
      N <- 24 # Preallocate 24 rows
      forecast.history <- data.frame(time = rep(as.Date(head(tail(history$date, (rw_years+1)*12), 1)), N),
                                     model = rep("", N),
                                     forecast = rep(NA, N),
                                     upper80 = rep(NA, N),
                                     upper95 = rep(NA, N),
                                     lower80 = rep(NA, N),
                                     lower95 = rep(NA, N),
                                     stringsAsFactors = F)
      
      for(i in seq(0, 23)){
        cutreal <- head(history, -(24-i)) # Erase (24-x) months from real history
        beginning <- head(tail(cutreal$date, (((rw_years+1)*12) + 1)), 1) # Define beginning of 4 year window 
        segment <- head(tail(cutreal[, type], (((rw_years+1)*12) + 1)), ((rw_years+1)*12)) # Define 4 year window
        series.ts <- ts(segment, start = decimal_date(beginning), frequency = 12) # Transform into a ts object
        
        scaleback <- as.numeric(bizdays(ts(seq(1), start = decimal_date(cutreal$date[length(cutreal$date)]), frequency = 12), FinCenter = "Zurich")) # Scaler for saving
        
        # Forecast
        chosen.model <- select_model(beginning, series.ts, freq = "monthly", rw_years) # Choose model
        mcast <- chosen_forecast(chosen.model, series.ts, history, freq = "monthly", rw_years) # Output a forecast
        missing.fcast <- data.frame(time = tail(cutreal$date, 1), 
                                    model = modelnames[chosen.model], 
                                    forecast = mcast$fcast[1] * scaleback,
                                    upper80 = mcast$upper80[1] * scaleback,
                                    upper95 = mcast$upper95[1] * scaleback,
                                    lower80 = mcast$lower80[1] * scaleback,
                                    lower95 = mcast$lower95[1] * scaleback, stringsAsFactors = FALSE)
        forecast.history[(i+1), ] <- missing.fcast
      }
      
      # Save into a csv
      write.csv(forecast.history, file = file, row.names = FALSE)
      
    } else{
      N <- 24 # Preallocate 24 rows
      forecast.history <- data.frame(time = rep(as.Date(head(tail(history$date, (rw_years+1)*12), 1)), N),
                                     model = rep("", N),
                                     forecast = rep(NA, N),
                                     upper80 = rep(NA, N),
                                     upper95 = rep(NA, N),
                                     lower80 = rep(NA, N),
                                     lower95 = rep(NA, N),
                                     stringsAsFactors = F)
      
      for(i in seq(0, 23)){
        cutreal <- head(history, -(24-i)) # Erase (24-x) weeks from real history
        wbeginning <- head(tail(cutreal$date, (((rw_years+1)*52) + 1)), 1) # Define beginning of 4 year window 
        segment <- head(tail(cutreal[, type], (((rw_years+1)*52) + 1)), ((rw_years+1)*52))
        series.ts <- ts(segment, start = decimal_date(wbeginning), frequency = 52)
        
        # Forecast
        chosen.model <- select_model(wbeginning, series.ts, "weekly", rw_years) # Choose model
        mcast <- chosen_forecast(chosen.model, series.ts, history, freq = "weekly", rw_years) # Output a forecast
        missing.fcast <- data.frame(time = tail(cutreal$date, 1), 
                                    model = modelnames[chosen.model], 
                                    forecast = mcast$fcast[1],
                                    upper80 = mcast$upper80[1],
                                    upper95 = mcast$upper95[1],
                                    lower80 = mcast$lower80[1],
                                    lower95 = mcast$lower95[1], stringsAsFactors = FALSE)
        forecast.history[(i+1), ] <- missing.fcast
      }
      
      # Save into a csv
      write.csv(forecast.history, file = file, row.names = FALSE)
    }
    
    # Then save newest forecast as usual
    write.table(save.this, file = file, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
    
  }
  
}

draw_history <- function(forecasthistory, history, freq = "monthly", type, palette){
  # Monthly forecast history plot: 2 years
  # Weekly forecast history plot: 6 months
  # Tables: full history for both
  
  fcast.history.full <- head(read.csv(forecasthistory, header = TRUE, stringsAsFactors = FALSE), -1)
  fcast.history.full$err <- abs((tail(head(history[[type]], -1), length(fcast.history.full$forecast)) - fcast.history.full$forecast) / tail(head(history[[type]], -1), length(fcast.history.full$forecast))) * 100
  
  # Plot
  if(freq == "monthly"){
    # Cut for plots: 2 years
    fcast.history <- tail(fcast.history.full, 24)
    
    p <- ggplot(data = fcast.history, aes(x = as.Date(time), y = forecast, colour = "Forecast"), label = model) +
      geom_line(aes(colour = "Forecast"), size = 0.9, alpha = palette[["alphaSeg"]]) +
      geom_line(data = head(history, -1), aes(x = date, y = get(type), colour = "Data"), size = 0.8, alpha = palette[["alphaSeg"]]) +
      labs(title = "Monthly history (2 years)",
           subtitle = paste0("Mean error: ", round(mean(fcast.history$err), 2), "%"),
           x = "",
           y = "Units per month",
           colour = "Series: ") +
      theme_minimal() + 
      theme(legend.position = "bottom", legend.margin = margin(t = -20, b = 0)) + 
      scale_colour_manual(values = c(palette[["colData"]], palette[["colPred"]]), labels = c("Data", "Forecast"))
  } else{
    # Cut for plots: 6 months
    fcast.history <- tail(fcast.history.full, 24)
    
    p <- ggplot(data = fcast.history, aes(x = as.Date(time), y = forecast, colour = "Forecast"), label = model) +
      geom_line(aes(colour = "Forecast"), size = 0.9, alpha = palette[["alphaSeg"]]) +
      geom_line(data = head(history, -1), aes(x = date, y = get(type), colour = "Data"), size = 0.8, alpha = palette[["alphaSeg"]]) +
      labs(title = "Weekly history (6 months)",
           subtitle = paste0("Mean error: ", round(mean(fcast.history$err), 2), "%"),
           x = "",
           y = "Units per week",
           colour = "Series: ") +
      theme_minimal() + 
      theme(legend.position = "bottom", legend.margin = margin(t = -20, b = 0)) + 
      scale_colour_manual(values = c(palette[["colData"]], palette[["colPred"]]), labels = c("Data", "Forecast"))
  }
  print(p)
  model.freq <- as.data.frame(table(fcast.history.full$model), stringsAsFactors = FALSE); colnames(model.freq) <- c("model", "frequency")
  model.freq$model <- as.character(model.freq$model)
  model.freq$error <- round(aggregate(fcast.history.full$err, by = list(fcast.history.full$model), mean)$x, 2)
  for(name in modelnames){
    if(!(name %in% model.freq$model)){
      model.freq <- rbind(model.freq, list(name, 0, 0))
    }
  }
  model.freq$frequency <- round(model.freq$frequency / sum(model.freq$frequency) * 100, 2)
  model.freq <- data.frame(model.freq, row.names = 1)
  datatable(model.freq, rownames = TRUE, filter="top", options = list(pageLength = 9))
  
}