
#'Compute historical variance
#'

compute_historical_variance <- function(){
    library(dplyr)
    state_data_with_backfill <- readRDS("../flu_data_with_backfill.rds")
    state_data_with_backfill_max <- state_data_with_backfill %>%
            group_by(region,epiweek) %>% filter(lag==max(lag))

    regions <- unique(state_data_with_backfill$region)
    vars <- matrix(NA,nrow=length(regions),ncol=11)
    rgn_idx_num <- 1
    for (rgn_idx in regions){
      for (l in 0:10){
        state_data_with_backfill_current_lag <- state_data_with_backfill %>%
          group_by(region,epiweek) %>% filter(lag==l)

        common_dates <- intersect(state_data_with_backfill_current_lag[state_data_with_backfill_current_lag$region == rgn_idx,]$epiweek,
                                  state_data_with_backfill_max[state_data_with_backfill_max$region == rgn_idx,]$epiweek)
        vars[rgn_idx_num,l] <- var(state_data_with_backfill_current_lag[state_data_with_backfill_current_lag$region == rgn_idx & state_data_with_backfill_current_lag$epiweek %in% common_dates,]$ili-state_data_with_backfill_max[state_data_with_backfill_max$region == rgn_idx & state_data_with_backfill_max$epiweek %in% common_dates,]$ili)
      }
      rgn_idx_num <- rgn_idx_num+1

    }

    return (vars)
}



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
  ts_frequency = 52,
  h = 1,
  epiweek=40,
  season=2016,
  region="nat",
  season_start_epiweek = 40
)
{
  library(jsonlite)
  epiweek_idx <- epiweek
  kcde_preds <- simulate(kcde_fit,newX = newX,newt=newt,nsim=nsim,h=h)
  if (nchar(epiweek) == 1){
    current_season_epiweek <- paste0(substr(season,6,10),"0",epiweek)

  } else{
    current_season_epiweek <- paste0(substr(season,1,4),epiweek)

  }

 # req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,1)))
 # nowcast_json <- jsonlite::prettify(rawToChar(req$content))
 # nowcast_obj_1_wk_ahead <- fromJSON(nowcast_json)
 # oneweek_ahead_nowcast <- nowcast_obj_1_wk_ahead$epidata$value

 # req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,2)))
 # nowcast_json <- jsonlite::prettify(rawToChar(req$content))
 ## nowcast_obj_2_wk_ahead <- fromJSON(nowcast_json)
 # twoweek_ahead_nowcast <- nowcast_obj_2_wk_ahead$epidata$value

  #kcde_preds[,1] <- .75*kcde_preds[,1] + .25*oneweek_ahead_nowcast
 # kcde_preds[,2] <- .75*kcde_preds[,2] + .25*twoweek_ahead_nowcast

  if (epiweek_idx <= 20){
    time_in <- cdcfluutils::get_num_MMWR_weeks_in_first_season_year(season) - season_start_epiweek + epiweek_idx
  } else{
    time_in <- epiweek_idx - season_start_epiweek + 1
  }


  sampled_historical <- matrix(NA,nrow=nsim,ncol=time_in)
  for (i in 1:nsim){
    sampled_historical[i,] <- rnorm(time_in,tail(newX,time_in),.001*1:time_in)
  }


  return (cbind(sampled_historical,kcde_preds))


}
