rbf <- function(x, y, sigma = 1)
{
  exp(- sigma * (x - y) ^ 2)
}
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
  epiweeks,
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
    #raw_trajectory_samples[,h_itr] <- rtruncnorm(nsim,a=0,b=max(kcde_fit[,ncol(kcde_fit) - h+ h_itr]),mean=mean(kcde_fit[,ncol(kcde_fit) - h+ h_itr]),sd=var(kcde_fit[,ncol(kcde_fit) - h+ h_itr]))
    kcde_samples <- kcde_fit[,ncol(kcde_fit) - h+ h_itr]
    tmp_traj_samples <- rep(NA,nsim)

    X_to_sim<- kcde_fit[,1:(ncol(kcde_fit) - h)]
    similarities <- rowSums(rbf(X_to_sim,c(interpolated_y),sigma = 1))
    similarities <- similarities/(sum(similarities))
    for (samp_idx in 1:nsim){
      obs_to_sample <- sample(kcde_samples,1,prob =similarities )
     if (tail(epiweeks,1) >= 40 & tail(epiweeks,1) <= 50 ){
        tmp_traj_samples[samp_idx] <- rnorm(1,mean = obs_to_sample,sd=.005)
     } else if (tail(epiweeks,1) >= 50 & tail(epiweeks,1) <= 52 ){
       tmp_traj_samples[samp_idx] <- rnorm(1,mean = obs_to_sample,sd=.01)
     } else if (tail(epiweeks,1) <= 10 ){
       tmp_traj_samples[samp_idx] <- rnorm(1,mean = obs_to_sample,sd=.01)
     } else if (tail(epiweeks,1) > 10 & tail(epiweeks,1) <= 20 ){
       tmp_traj_samples[samp_idx] <- rnorm(1,mean = obs_to_sample,sd=0.00000000001)
     }
    }
    raw_trajectory_samples[,h_itr] <- .99*tmp_traj_samples + .01*runif(nsim,0,100)
  }

  ### smooth the trajectories
  for ( i in 1:nrow(raw_trajectory_samples)){
    y <- raw_trajectory_samples[i,]
    loes_fit <- loess(y~., data=data.frame(y=y),n=.05)
    raw_trajectory_samples[i,] <- predict(loes_fit)
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
# library(sarimaTD)
#
# state_data <-read.csv("data/state_data.csv")
# unique_states <- unique(state_data$region)
# nsim <- 10000
#
# fs_mat_kcde <- matrix(NA,nrow=100*length(unique_states),ncol=4)
# fs_mat_sarima <- matrix(NA,nrow=100*length(unique_states),ncol=4)
#
# mse_mat_kcde <- matrix(NA,nrow=100*length(unique_states),ncol=4)
# mse_mat_sarima <- matrix(NA,nrow=100*length(unique_states),ncol=4)
#
#
# ls_index <- 1
# for (state in sample(unique_states,10)){
#   ma_data <- state_data[state_data$region == state,]
#   for (first_test_index in 200:220){
#     kcde_fit <- fit_kcde(ma_data$unweighted_ili[1:(first_test_index-1)],
#
#                          ts_frequency = 52,
#                          h=4,transformation = "none")
#
#     preds <- simulate.KCDE(kcde_fit,nsim,newX = ma_data$unweighted_ili[1:(first_test_index-1)],ts_frequency = 52,
#                            h=4,seed = 1,seasonal_difference = F)
#
#
#     sarima_fit <- fit_sarima(tail(ma_data$unweighted_ili[1:(first_test_index-1)],100),52,transformation = "box-cox",seasonal_difference = TRUE)
#
#     sarima_pred <-
#       simulate(
#         object = sarima_fit,
#         nsim = 1000,
#         seed = 1,
#         newdata = ma_data$unweighted_ili[1:(first_test_index-1)],
#         h = 4
#       )
#     truth <- ma_data$unweighted_ili[first_test_index:(first_test_index+3)]
#     for (h in 1:4){
#       fs_mat_kcde[ls_index,h] <- sum(round(preds[,ncol(preds) -4 + h],2) == round(truth[h],2))/1000
#       fs_mat_sarima[ls_index,h] <- sum(round(sarima_pred[,ncol(sarima_pred) -4 + h],2) == round(truth[h],2))/1000
#       mse_mat_kcde[ls_index,h] <- (mean(preds[,ncol(preds)-4+h],na.rm=T) - truth[h])^2
#       mse_mat_sarima[ls_index,h] <- (mean(sarima_pred[,ncol(preds)-4+h],na.rm=T) - truth[h])^2
#
#     }
#     ls_index <- ls_index+1
#   }
# }
# #ggplot(data=data.frame(y=c(t(preds)),x=rep(1:ncol(preds),1000),group=rep(1:1000,each=ncol(preds))),aes(x=x,y=y,group=group)) + geom_line() + geom_line(data=ma_data[first_test_index:(first_test_index+3),],aes(x=ncol(preds)-4+ 1:4,y=unweighted_ili,group=1,col='truth'))
#
# exp(mean(pmax(-10,log(fs_mat_sarima)),na.rm=T))
# exp(mean(pmax(-10,log(fs_mat_kcde)),na.rm=T))
# mean(mse_mat_kcde,na.rm=T)
# mean(mse_mat_sarima,na.rm=T)
