#' sample revised ILI trajectories using normal approximation
#'
#' @param n number of revised trajectories to sample
#' @param observed_inc observed incidence so far this season, starting at EW40 and going to the most recent report
#' @param epiweek_idx most recent epidemic week (equals 40 + length(observed_inc) - # weeks in season)
#' @param region region trajectory was made from
#' @param season current season in 20xx/20xx+1 format
#' @param add_nowcast logical; add nowcast based on delphi epicast?
#' @param min_value minimum value for revised wILI; any lower values are truncated here.
#' @param return_sampled_id logical; return a data frame with identifiers of sampled region, season, and epiweek?
#' @return n by length(observed_inc) matrix of samples for possible revised ili.
#'
#' @export
rRevisedILI_fast <- function(
  n,
  observed_inc,
  epiweek_idx,
  region,
  season,
  season_start_epiweek = 40,
  add_nowcast = FALSE,
  min_value = 0.05,
  return_sampled_id = FALSE) {
  library(jsonlite)
  library(dplyr)

  if(region %in% c('nat', paste0('hhs', 1:10))) {
    flu_data_with_backfill <- cdcfluutils::nat_reg_flu_data_with_backfill
    historical_vars <- readRDS("./data/historical_vars_regional.RDS")
    regions <- c(paste0("hhs",1:10),"nat")
    region_idx <- which(regions==region)
  } else if(region %in% c(
    'al', 'ak', 'az', 'ar', 'ca', 'co', 'ct', 'de', 'fl', 'ga', 'hi', 'id', 'il',
    'in', 'ia', 'ks', 'ky', 'la', 'me', 'md', 'ma', 'mi', 'mn', 'ms', 'mo', 'mt',
    'ne', 'nv', 'nh', 'nj', 'nm', 'ny_minus_jfk', 'nc', 'nd', 'oh', 'ok', 'or',
    'pa', 'ri', 'sc', 'sd', 'tn', 'tx', 'ut', 'vt', 'va', 'wa', 'wv', 'wi', 'wy',
    'as', 'mp', 'dc', 'gu', 'pr', 'vi', 'ord', 'lax', 'jfk')) {
    flu_data_with_backfill <- cdcfluutils::state_local_flu_data_with_backfill
    historical_vars <- readRDS("/Users/gcgibson/historical_vars_state.RDS")
    regions <- c(paste0("hhs",1:10),"nat")
    region_idx <- which(regions==region)
  } else {
    stop("Invalid region provided to rRevisedILI")
  }



  if (epiweek_idx <= 20){
    time_in <- cdcfluutils::get_num_MMWR_weeks_in_first_season_year(season) - season_start_epiweek + epiweek_idx + 1
  } else{
    time_in <- epiweek_idx - season_start_epiweek + 1
  }

  # fully observed data



  total_traj <-  matrix(nrow = n, ncol= time_in)


  total_traj <- matrix(rnorm(n*time_in,tail(observed_inc,time_in),rev(historical_vars[,region_idx][1:time_in])),nrow=n,byrow=TRUE)

  total_traj[total_traj < min_value] <- min_value

  ## add nowcast
  if(add_nowcast) {
    current_season_epiweek <- ifelse(epiweek_idx <=20,paste0(substr(season,6,12),epiweek_idx),paste0(substr(season,1,4),epiweek_idx))

    req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,1)))
    nowcast_json <- jsonlite::prettify(rawToChar(req$content))
    nowcast_obj_1_wk_ahead <- fromJSON(nowcast_json)

    req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,2)))
    nowcast_json <- jsonlite::prettify(rawToChar(req$content))
    nowcast_obj_2_wk_ahead <- fromJSON(nowcast_json)

    total_traj <- cbind(total_traj, rep(nowcast_obj_1_wk_ahead$epidata$value,n),rep(nowcast_obj_2_wk_ahead$epidata$value,n))
  }

  if(return_sampled_id) {
    return(list(total_traj = total_traj, sampled_id = sampled_id))
  } else {
    return (total_traj)
  }
}

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
  library(cdcfluutils)
  epiweek_idx <- epiweek
  kcde_preds <- simulate(kcde_fit,newX = newX,nsim=nsim,h=h)
  if (nchar(epiweek) == 1){
    current_season_epiweek <- paste0(substr(season,6,10),"0",epiweek)

  } else{
    current_season_epiweek <- paste0(substr(season,1,4),epiweek)

  }
  req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,1)))
  nowcast_json <- jsonlite::prettify(rawToChar(req$content))
  nowcast_obj_1_wk_ahead <- fromJSON(nowcast_json)
  oneweek_ahead_nowcast <- nowcast_obj_1_wk_ahead$epidata$value

  req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,2)))
  nowcast_json <- jsonlite::prettify(rawToChar(req$content))
  nowcast_obj_2_wk_ahead <- fromJSON(nowcast_json)
  twoweek_ahead_nowcast <- nowcast_obj_2_wk_ahead$epidata$value

  kcde_preds[,1] <- .5*kcde_preds[,1] + .5*oneweek_ahead_nowcast
  kcde_preds[,2] <- .5*kcde_preds[,2] + .5*twoweek_ahead_nowcast

  if (epiweek_idx <= 20){
    time_in <- cdcfluutils::get_num_MMWR_weeks_in_first_season_year(season) - season_start_epiweek + epiweek_idx
  } else{
    time_in <- epiweek_idx - season_start_epiweek + 1
  }


  sampled_historical <-rRevisedILI_fast(n = 1000,tail(newX,time_in),region = region,add_nowcast = FALSE,epiweek_idx = epiweek)



  return (cbind(sampled_historical,kcde_preds))


}
