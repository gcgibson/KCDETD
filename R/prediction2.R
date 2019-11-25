
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
simulate.KCDE2 <- function(
  object,
  nsim = 1,
  seed = NULL,
  newX,
  ts_frequency = 52,
  h = 1
) {
  library(truncnorm)
  rbf <- function(x, y, sigma = bandwidth)
  {
    exp(- sigma * (x - y) ^ 2)
  }

  newdata <- newX

  kcde_fit <-knn_forecasting(newX,h=h,lags=1:6,k=20)
  kcde_fit <- nearest_neighbors(kcde_fit)$nneighbors

  raw_trajectory_samples <- matrix(NA,nrow=nsim,ncol=h)

  for (h_itr in 1:h){
    #raw_trajectory_samples[,h_itr] <- rtruncnorm(nsim,a=0,b=max(kcde_fit[,ncol(kcde_fit) - h+ h_itr]),mean=mean(kcde_fit[,ncol(kcde_fit) - h+ h_itr]),sd=var(kcde_fit[,ncol(kcde_fit) - h+ h_itr]))
    kcde_samples <- kcde_fit[,ncol(kcde_fit) - h+ h_itr]
    errors <- matrix(NA,nrow=1000,ncol=3)
    err_idx <- 1
    for (bandwith1 in c(.1,1,10)){
      for (bandwith2 in c(.1,1,10)){
        for (sample_idx in 1:length(kcde_samples)){
          kcde_samples_without <- kcde_samples[-sample_idx]
          loo_truth <- kcde_samples[sample_idx]
          X_to_sim<- kcde_fit[,1:(ncol(kcde_fit) - h)]
          similarities <- rowSums(rbf(X_to_sim,c(newX),sigma = bandwith1))
          similarities <- similarities/(sum(similarities))
          tmp_traj_samples <- rep(NA,nsim)
          for (samp_idx in 1:nsim){
            obs_to_sample <- sample(kcde_samples,1,prob =similarities )
            tmp_traj_samples[samp_idx] <- rnorm(1,mean = obs_to_sample,sd=bandwith2)
          }
          ls <- sum(loo_truth <= tmp_traj_samples +.5 & loo_truth>= tmp_traj_samples-.5)/1000
          errors[err_idx,] <-c(log(ls),bandwith1,bandwith2)
          err_idx <- err_idx + 1
        }
      }
    }
    errors <- as.data.frame(errors)
    colnames(errors) <- c("ls","bw1","bw2")
    avg_ls <-errors %>% dplyr::group_by(bw1,bw2) %>% summarize(mls=mean(ls,na.rm=T))
    bw1_optim <- avg_ls[avg_ls$mls ==max(avg_ls$mls,na.rm=T),]$bw1[1]
    bw2_optim <- avg_ls[avg_ls$mls ==max(avg_ls$mls,na.rm=T),]$bw2[1]


    for (sample_idx in 1:length(kcde_samples)){
      kcde_samples_without <- kcde_samples[-sample_idx]
      loo_truth <- kcde_samples[sample_idx]
      X_to_sim<- kcde_fit[,1:(ncol(kcde_fit) - h)]
      similarities <- rowSums(rbf(X_to_sim,c(newX),sigma = bw1_optim))
      similarities <- similarities/(sum(similarities))
      tmp_traj_samples <- rep(NA,nsim)
      for (samp_idx in 1:nsim){
        obs_to_sample <- sample(kcde_samples,1,prob =similarities )
        tmp_traj_samples[samp_idx] <- rnorm(1,mean = obs_to_sample,sd=bw2_optim)
      }
      ls <- sum(loo_truth <= tmp_traj_samples +.5 & loo_truth>= tmp_traj_samples-.5)/1000
      errors[err_idx,] <-c(log(ls),bandwith1,bandwith2)
      err_idx <- err_idx + 1
    }

    raw_trajectory_samples[,h_itr] <- .99*tmp_traj_samples + .01*runif(nsim,0,100)
  }

  ### smooth the trajectories
  for ( i in 1:nrow(raw_trajectory_samples)){
    y <- raw_trajectory_samples[i,]
    loes_fit <- loess(y~x, data=data.frame(y=y,x=1:length(y)),span=.4)
    raw_trajectory_samples[i,] <- pmax(predict(loes_fit),0)
  }
  return(raw_trajectory_samples)
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
