benchmarkCVForecastingMethods <- function(elementName, index, frequency, 
                                        dateFrom, dateTo, maxHorizonCV,
                                        windowSizeCVPercent, error_index) {
  # TIME SERIES PARAMETERS AND SEARCH PARAMETERS #
  dateFrom <<- dateFrom   # TIME WINDOW TO SUBSET FROM THE INDEX
  dateTo <<- dateTo
  cvHorizon <- 5
  # Time Series, model fitting and forecasts --------------------------------
  error_happened <- tryCatch(
    {
      metric <- searchElement(name = elementName, index = index, tsfrequency = frequency, returnDF = FALSE)
      windowSize <- length(metric) * windowSizeCVPercent
      if (windowSize + 2 * maxHorizonCV > length(metric)) {
        TRUE #return
      } else {
        FALSE #return
      }
    },
    error = function(e) {
      print(e)
      return(TRUE)
    }
  )
  
  if (error_happened == TRUE) {
    return(stringMethods[4]) # ETS (DEFAULT)
  }
  
  windowSizeCV = floor(length(metric)*windowSizeCVPercent)
  # FIT MODELS, COMPUTE FORECASTS #
  a_a <- tryCatch(
    expr = {
      arima_cv <- cvts(metric, FUN = function(x) auto.arima(x, stepwise = FALSE, approximation = FALSE),
                       rolling = TRUE, maxHorizon = maxHorizonCV,
                       windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                       verbose = TRUE, horizonAverage = TRUE)
      cat("ARIMA OK\n")
      a_a<-as.numeric(accuracy(arima_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_a<- NA
    })
  
  a_a_fs <- tryCatch(
    expr = {
      arima_fs_cv <- cvts(metric, FUN = function(x) auto.arima(x, D=1, stepwise = FALSE, approximation = FALSE),
                          rolling = TRUE, maxHorizon = maxHorizonCV,
                          windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                          verbose = TRUE, horizonAverage = TRUE)
      cat("ARIMA_FS OK\n")
      a_a_fs<-as.numeric(accuracy(arima_fs_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_a_fs <- NA
    })
  
  a_t <- tryCatch(
    expr = {
      theta_cv <- cvts(metric, thetaf, rolling = TRUE, maxHorizon = maxHorizonCV,
                       windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                       verbose = TRUE, horizonAverage = TRUE)
      cat("THETA OK\n")
      a_t<-as.numeric(accuracy(theta_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_t <- NA
    })
  
  a_ets <- tryCatch(
    expr = {
      ets_cv <- cvts(metric, ets, rolling = TRUE, maxHorizon = maxHorizonCV,
                     windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                     verbose = TRUE, horizonAverage = TRUE)
      cat("ETS OK\n")
      a_ets<-as.numeric(accuracy(ets_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_ets <- NA
    })
  
  a_ets_fc <- tryCatch(
    expr = {
      ets_fs_cv <- cvts(metric, FUN = function(x) ets(x, damped=TRUE), rolling = TRUE, maxHorizon = maxHorizonCV,
                        windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                        verbose = TRUE, horizonAverage = TRUE)
      cat("ETS_FD OK\n")
      a_ets_fc<-as.numeric(accuracy(ets_fs_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_ets_fc <- NA
    })
  
  a_bets <- tryCatch(
    expr = {
      bagged_ets_cv <- cvts(metric, baggedETS, rolling = TRUE, maxHorizon = maxHorizonCV,
                            windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                            verbose = TRUE, horizonAverage = TRUE)
      cat("BAGGED_ETS OK\n")
      a_bets<-as.numeric(accuracy(bagged_ets_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_bets <- NA
    })
  
  a_stl <- tryCatch(
    expr = {
      stl_cv <- cvts(metric, mstl, rolling = TRUE, maxHorizon = maxHorizonCV,
                     windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                     verbose = TRUE, horizonAverage = TRUE)
      cat("STL OK\n")
      a_stl<-as.numeric(accuracy(stl_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_stl <- NA
    })
  
  a_nn <- tryCatch(
    expr = {
      nn_cv <- cvts(metric, nnetar, rolling = TRUE, maxHorizon = maxHorizonCV,
                    windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                    verbose = TRUE, horizonAverage = TRUE)
      cat("NN OK\n")
      a_nn<-as.numeric(accuracy(nn_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_nn <- NA
    })
  
  a_p <- tryCatch(
    expr = {
      NA # not implemented
    },
    error = function(e) {
      print(e)
      a_p <- NA
    })
  
  a_tb <- tryCatch(
    expr = {
      tbats_cv <- cvts(metric, tbats, rolling = TRUE, maxHorizon = maxHorizonCV,
                       windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                       verbose = TRUE, horizonAverage = TRUE)
      cat("TBATS OK\n")
      a_tb<-as.numeric(accuracy(tbats_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_tb <- NA
    })
  
  a_h <- tryCatch(
    expr = {
      hybrid_cv <- cvts(metric, FUN = function(x) hybridModel(x, models="aenstz"),
                        rolling = TRUE, maxHorizon = maxHorizonCV,
                        windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                        verbose = TRUE, horizonAverage = TRUE)
      cat("HYBRID OK\n")
      a_h<-as.numeric(accuracy(hybrid_cv)[error_index])
    },
    error = function(e) {
      print(e)
      a_h <- NA
    })

  a_naive <- tryCatch(
    expr = {
      snaive_cv <- cvts(metric, naive, rolling = TRUE, maxHorizon = maxHorizonCV,
                        windowSize = windowSizeCV, saveModels = FALSE, saveForecasts = FALSE,
                        verbose = TRUE, horizonAverage = TRUE)
      cat("NAIVE OK\n")
      a_naive<-as.numeric(accuracy(snaive_cv)[2])
    },
    error = function(e) {
      print(e)
      a_naive <- NA
    })
  
  error_matrix <- c(a_a,a_a_fs,a_t,a_ets,a_ets_fc,
                    a_bets,a_stl,a_nn,a_h,a_p,a_tb,a_naive)
  mdat <- matrix(error_matrix, nrow=length(error_matrix), ncol = 1,
                 byrow = TRUE, dimnames = 
                   list(stringMethods))
  
  best_method <- rownames(mdat)[which.min(abs(mdat))]
  cat(mdat)
  
  return(best_method)
}
