## functions for KCDE estimation

#' Estimate KCDE model
#'
#' @param y a univariate time series or numeric vector.
#' @param ts_frequency frequency of time series.  Must be provided if y is not
#'   of class "ts".  See the help for stats::ts for more.
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", "forecast-box-cox", or "none".  See details for more.
#' @param bc_gamma numeric offset used in Box-Cox transformation; gamma is
#'   added to all observations before transforming.  Default value of 0.5
#'   allows us to use the Box-Cox transform (which requires positive inputs)
#'   in case of observations of 0, and also ensures that the de-transformed
#'   values will always be at least -0.5, so that they round up to non-negative
#'   values.
#' @param seasonal_difference boolean; take a seasonal difference before passing
#'   to auto.arima?

#' @return a KCDE model fit
#'
#' @details Fit KCDE using np
#'
#' @export
fit_kcde <- function(
  y,
  t,
  ts_frequency,
  transformation = "none",
  bc_gamma = 0.5,
  seasonal_difference = FALSE,
  h=4,
  d = NA,
  D = NA) {

  library(quantmod)
  library(np)
  library(tsfknn)
  ## Validate arguments
  if(!(is.numeric(y) || is.ts(y))) {
    stop("The argument y must be a numeric vector or object of class 'ts'.")
  }

  if(!is.ts(y) && missing(ts_frequency)) {
    stop("If y is not an object of class 'ts', the ts_frequency argument must be supplied.")
  }

  if(is.ts(y)) {
    ts_frequency <- frequency(y)
  }

  ## Initial transformation, if necessary
  if(identical(transformation, "box-cox")) {
    est_bc_params <- car::powerTransform(y + bc_gamma, family = "bcPower")
    est_bc_params <- list(
      lambda = est_bc_params$lambda,
      gamma = bc_gamma)
  }
  transformed_y <- do_initial_transform(
    y = y,
    transformation = transformation,
    bc_params = est_bc_params)

  ## Initial seasonal differencing, if necessary
  if(seasonal_difference) {
    differenced_y <- do_seasonal_difference(
      y = transformed_y,
      ts_frequency = ts_frequency)
  } else {
    differenced_y <- ts(transformed_y, frequency = ts_frequency)
  }
  kcde_fit <- list()
  ## Get KCDE fit
  pred <- knn_forecasting(differenced_y, h = 4, lags = 1:12, k = 10)
  kcde_fit[[1]] <- nearest_neighbors(pred)
  kcde_fit$kcde_call <- match.call()
  for(param_name in c("y", "ts_frequency", "transformation", "seasonal_difference", "d", "D")) {
    kcde_fit[[paste0("KCDE_used_", param_name)]] <- get(param_name)
  }
  if(identical(transformation, "box-cox")) {
    kcde_fit$KCDE_est_bc_params <- est_bc_params
  }

  class(kcde_fit) <- c("KCDE", class(kcde_fit))

  return(kcde_fit)
}



## TEST

#library(np)

#fit <-fit_kcde(1:10,1,transformation = "none",seasonal_difference = FALSE)
