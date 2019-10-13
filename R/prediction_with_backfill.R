#'
#'
#' @param object a KCDE fit of class "sarima_td", as returned by fit_sarima
#' @param nsim number of sample trajectories to simulate
#' @param seed either `NULL` or an integer that will be used in a call to
#'   `set.seed` before simulating the response vectors.  If set, the value is
#'   saved as the "seed" attribute of the returned value.  The default, `Null`,
#'   will not change the random generator state, and return `.Random.seed`
#'   as the "seed" attribute
#' @param newX numeric vector of new data to simulate forward from
#' @param newt numeric vector of new time points to use (if applicable)
#' @param h number of time steps forwards to simulate
#' @param epiweek the current epiweek
#' @param season the current season
#' @return an nsim by h + current time in season matrix with simulated values
#'
simulate_with_backfill <- function(
  object,
  nsim = 1,
  seed = NULL,
  newX,
  newt,
  ts_frequency = 52,
  h = 1,
  epiweek=40,
  season=2016,
  region="nat"
)
{
  kcde_preds <- simulate(kcde_fit,newX = newX,newt=newt,nsim=nsim,h=h)

  if (epiweek > 20){
    time_in <- epiweek - 40 +1
  } else{
    time_in <- 12 + epiweek
  }
  print (time_in)
  backfill_sim <-cdcfluutils::rRevisedILI(n = nsim,observed_inc = tail(newX,time_in),region = region,epiweek_idx = epiweek,season = season,add_nowcast = TRUE)
  print (dim(backfill_sim))
  return (cbind(backfill_sim[,1:(ncol(backfill_sim-2))],
                .5*backfill_sim[,ncol(backfill_sim)-1] + .5*kcde_preds[,ncol(kcde_preds)-1],
                .5*backfill_sim[,ncol(backfill_sim)] + .5*kcde_preds[,ncol(kcde_preds)],
                kcde_preds[,3:ncol(kcde_preds)]))


}
