benchmarkForecastingMethods <- function(elementName, index, frequency, 
                                        dateTrainFrom, dateTrainTo, 
                                        horizon, error_index) {
  # TIME SERIES PARAMETERS AND SEARCH PARAMETERS #
  dateFrom <<- dateTrainFrom   # TIME WINDOW TO SUBSET FROM THE INDEX (TRAIN+TEST)
  cat(dateFrom)
  dateTo <<- as.Date(dateTrainTo) + horizon
  cat(dateTo)
  cvHorizon <- 5
  # Time Series, model fitting and forecasts --------------------------------
  error_happened <- tryCatch(
    {
      metric <- searchElement(name = elementName, index = index, tsfrequency = frequency, returnDF = FALSE)
      metric_prophet <- searchElement(name = elementName, index = index, tsfrequency = frequency, returnDF = TRUE)
      submetric <- subset(metric,end=length(metric)-horizon)
      submetric_prophet <- metric_prophet[1:(length(metric_prophet$y)-horizon),]
      test_data <- subset(metric, start = length(metric)-horizon+1)
      FALSE #return
    },
    error = function(e) {
      print(e)
      return(TRUE)
    }
  )
  
  if (error_happened == TRUE) {
    return(stringMethods[4]) # ETS (DEFAULT)
  }

  # FIT MODELS, COMPUTE FORECASTS #
  a_a <- tryCatch(
    expr = {
      model_arima <- auto.arima(submetric, stepwise = FALSE, approximation = FALSE)
      f_arima <- forecast(model_arima, h = horizon)
      cat("ARIMA OK\n")
      a_a<-accuracy(f_arima,test_data)
    },
    error = function(e) {
      print(e)
      a_a<- NA
    })
  
  a_a_fs <- tryCatch(
    expr = {
      model_arima_fs <- auto.arima(submetric, D=1, stepwise = FALSE, approximation = FALSE)
      f_arima_fs <- forecast(model_arima_fs, h = horizon)
      cat("ARIMA_FS OK\n")
      a_a_fs<-accuracy(f_arima_fs,test_data)
    },
    error = function(e) {
      print(e)
      a_a_fs <- NA
    })
  
  a_t <- tryCatch(
    expr = {
      model_theta <- thetam(submetric)
      f_theta <- forecast(model_theta, h = horizon)
      cat("THETA OK\n")
      a_t<-accuracy(f_theta,test_data)
    },
    error = function(e) {
      print(e)
      a_t <- NA
    })
  
  a_ets <- tryCatch(
    expr = {
      model_ets <- ets(submetric)
      f_ets <- forecast(model_ets, h = horizon)
      cat("ETS OK\n")
      a_ets<-accuracy(f_ets,test_data)
    },
    error = function(e) {
      print(e)
      a_ets <- NA
    })
  
  a_ets_fc <- tryCatch(
    expr = {
      model_ets_fc <- ets(submetric, damped=TRUE)
      f_ets_fc <- forecast(model_ets_fc, h = horizon)
      cat("ETS_FD OK\n")
      a_ets_fc<-accuracy(f_ets_fc,test_data)
    },
    error = function(e) {
      print(e)
      a_ets_fc <- NA
    })
  
  a_bets <- tryCatch(
    expr = {
      model_baggedETS <- baggedETS(submetric)
      f_baggedETS <- forecast(model_baggedETS, h = horizon)
      cat("BAGGED_ETS OK\n")
      a_bets<-accuracy(f_baggedETS,test_data)
    },
    error = function(e) {
      print(e)
      a_bets <- NA
    })
  
  a_stl <- tryCatch(
    expr = {
      model_STL <- mstl(submetric)
      f_STL <- forecast(model_STL, h = horizon)
      cat("STL OK\n")
      a_stl<-accuracy(f_STL,test_data)
    },
    error = function(e) {
      print(e)
      a_stl <- NA
    })
  
  a_nn <- tryCatch(
    expr = {
      model_nn <- nnetar(submetric)
      f_nn <- forecast(model_nn, h = horizon, PI = TRUE)
      cat("NN OK\n")
      a_nn<-accuracy(f_nn,test_data)
    },
    error = function(e) {
      print(e)
      a_nn <- NA
    })
  
  a_p <- tryCatch(
    expr = {
      model_ets <- ets(submetric)
      f_ets <- forecast(model_ets, h = horizon)
      model_prophet <- prophet(submetric_prophet, daily.seasonality = 'auto', weekly.seasonality = 'auto')
      f_prophet <- predict(model_prophet, make_future_dataframe(model_prophet, periods = horizon, freq = 'day', include_history = FALSE), future)
      f_prophet_hack <- duplicate(f_ets, shallow = FALSE)
      f_prophet_hack$mean <- ts(data = f_prophet$yhat, start = start(f_ets$mean), end = end(f_ets$mean), frequency = tsFrequency)
      cat("PROPHET OK\n")
      a_p<-accuracy(f_prophet_hack,test_data)
    },
    error = function(e) {
      print(e)
      a_p <- NA
    })
  
  a_tb <- tryCatch(
    expr = {
      model_tbats <- tbats(submetric)
      f_tbats <- forecast(model_tbats, h = horizon)
      cat("TBATS OK\n")
      a_tb<-accuracy(f_tbats,test_data)
    },
    error = function(e) {
      print(e)
      a_tb <- NA
    })
  
  a_h <- tryCatch(
    expr = {
      model_hybridCVModel <- hybridModel(submetric,
                                         lambda = "auto", 
                                         windowSize = (length(submetric)-cvHorizon*2),
                                         weights = "cv.errors", cvHorizon = cvHorizon,
                                         horizonAverage = TRUE, 
                                         a.args = list(stepwise = FALSE, trace = FALSE),
                                         e.args = list(allow.multiplicative.trend = TRUE),
                                         parallel = TRUE,
                                         num.cores = 2)
      f_hybrid <- forecast(model_hybridCVModel, h = horizon)
      cat("HYBRID OK\n")
      a_h<-accuracy(f_hybrid,test_data)
    },
    error = function(e) {
      print(e)
      a_h <- NA
    })
  
  # a_stlm <- tryCatch(
  #   expr = {
  #     model_stlm <- stlm(submetric)
  #     f_stlm <- forecast(model_stlm, h = horizon)
  #     cat("STLM OK\n")
  #     a_stlm<-accuracy(f_stlm, test_data)
  #   },
  #   error = function(e) {
  #     print(e)
  #     a_stlm <- NA
  #   })
  
  a_naive <- tryCatch(
    expr = {
      f_naive <- naive(submetric, h = horizon)
      cat("NAIVE OK\n")
      a_naive<-accuracy(f_naive, test_data)
    },
    error = function(e) {
      print(e)
      a_naive <- NA
    })
  
  error_matrix <- c(a_a[error_index],a_a_fs[error_index],a_t[error_index],a_ets[error_index],
           a_ets_fc[error_index],a_bets[error_index],a_stl[error_index],
           a_nn[error_index],a_h[error_index],a_p[error_index],a_tb[error_index],
           a_naive[error_index])
  mdat <- matrix(error_matrix, nrow=length(error_matrix), ncol = 1, byrow = TRUE, dimnames = 
                   list(stringMethods))
  
  best_method <- rownames(mdat)[which.min(abs(mdat))]
  cat(mdat)

  return(best_method)
}
