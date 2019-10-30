
#' Simulate predictive trajectories from a KCDE model with class KCDE
#'
#' This function does a few things worth noting.  It linearly interpolates
#' missing values on the interior of the time series, which is necessary for
#' some model fits with moving average components.  It also handles
#' transformations and possible seasonal differencing that may have been done
#' in the fit_sarima function.
#'
#' @param object a KCDE fit of class "sarima_td", as returned by fit_sarima
#' @param nsim number of sample trajectories to simulate
#' @param seed either `NULL` or an integer that will be used in a call to
#'   `set.seed` before simulating the response vectors.  If set, the value is
#'   saved as the "seed" attribute of the returned value.  The default, `Null`,
#'   will not change the random generator state, and return `.Random.seed`
#'   as the "seed" attribute
#' @param newdata numeric vector of new data to simulate forward from
#' @param h number of time steps forwards to simulate
#'
#' @return an nsim by h matrix with simulated values
#'
#' @export
simulate.KCDE <- function(
  object,
  nsim = 1,
  seed = NULL,
  newX,
  ts_frequency = 52,
  h = 1,
  seasonal_difference=FALSE
) {
  library(truncnorm)
  newdata <- newX
  if(is.null(seed)) {
    set.seed(1)
    seed <- .Random.seed
  } else {
    set.seed(seed)
  }

  if(is.null(object$sarimaTD_call)) {
    transformation <- "none"
    seasonal_difference <- seasonal_difference
    ts_frequency <- ts_frequency
  } else {
    transformation <- object$sarimaTD_used_transformation
    seasonal_difference <-
      object$sarimaTD_used_seasonal_difference
    ts_frequency <- object$arma[5]
  }

  ## Update SARIMA fit object with transformed and seasonally differenced data
  ## Initial transformation, if necessary
  if(identical(transformation, "box-cox")) {
    est_bc_params <- object$KCDE_est_bc_params
  } else {
    est_bc_params <- NULL
  }
  transformed_y <- do_initial_transform(
    y = newdata,
    transformation = transformation,
    bc_params = est_bc_params)

  ## Initial seasonal differencing, if necessary
  if(seasonal_difference) {
    differenced_y <- do_seasonal_difference(
      y = transformed_y,
      ts_frequency = ts_frequency)
  } else {
    differenced_y <- ts(transformed_y, frequency = 52)
  }

  ## Drop leading missing values, fill in internal missing values via linear
  ## interpolation.  This is necessary to ensure non-missing predictions if the
  ## sarima model has a moving average component.
  ## deal with internal missing or infinite values in y that can
  ## result in simulated trajectories of all NAs if the model has a moving
  ## average component.  Here we do this by linear interpolation.
  ##
  ## Another (better?) solution would be to write a version of stats::filter
  ## that does the "right thing" with NAs if filter coefficients are 0 and then
  ## use that function in forecast:::myarima.sim
  interpolated_y <- interpolate_and_clean_missing(differenced_y)

  kcde_fit <- object[[1]]$nneighbors
  raw_trajectory_samples <- matrix(NA,nrow=nsim,ncol=h)

  for (h_itr in 1:h){
    raw_trajectory_samples[,h_itr] <- rtruncnorm(nsim,a=0,b=Inf,mean=mean(kcde_fit[,ncol(kcde_fit) - h+ h_itr]),sd=3*h_itr*var(kcde_fit[,ncol(kcde_fit) - h+ h_itr]))
  }


  ## Sampled trajectories are of seasonally differenced transformed time series
  ## Get to trajectories for originally observed time series ("orig") by
  ## adding seasonal lag of incidence and inverting the transformation
  orig_trajectory_samples <- raw_trajectory_samples
  if(seasonal_difference) {
    for(i in seq_len(nsim)) {
      orig_trajectory_samples[i, ] <-
        invert_seasonal_difference(
          dy = raw_trajectory_samples[i, ],
          y = transformed_y,
          ts_frequency = ts_frequency)
      orig_trajectory_samples[i, ] <-
        invert_initial_transform(
          y = orig_trajectory_samples[i, ],
          transformation = transformation,
          bc_params = est_bc_params)
    }
  } else {
    for(i in seq_len(nsim)) {
      orig_trajectory_samples[i, ] <-
        invert_initial_transform(
          y = raw_trajectory_samples[i, ],
          transformation = transformation,
          bc_params = est_bc_params)
    }
  }

  #attr(orig_trajectory_samples, "seed") <- seed

  return(orig_trajectory_samples)
}



## TEST

#sims <-simulate(fit,newdata = 1:10,nsim=100)
