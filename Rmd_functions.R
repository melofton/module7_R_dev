# Functions needed to complete Module 7: Using Data to Improve Ecological Forecasts
# Author: Mary Lofton
# Date: 06APR23

# Load packages
library(tidyverse)
library(lubridate)
library(zoo)

## Define functions----

### AR model related functions

# can be used to calculate the significance threshold of an ACF or PACF object
calculate_pacf_threshold <- function(x, ci=0.95, ci.type="white"){
  #' Gets confidence limit data from acf object `x`
  if (!ci.type %in% c("white", "ma")) stop('`ci.type` must be "white" or "ma"')
  if (class(x) != "acf") stop('pass in object of class "acf"')
  clim0 <- qnorm((1 + ci)/2) / sqrt(x$n.used)
  if (ci.type == "ma") {
    clim <- clim0 * sqrt(cumsum(c(1, 2 * x$acf[-1]^2))) 
    return(clim[-length(clim)])
  } else {
    return(clim0)
  }
}


### ENKF functions from Jake Zwart GLEON workshop---

#### Function to create vector to hold states and parameters for updating---
#' vector for holding states and parameters for updating
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_params_est number of parameters we're calibrating
#' @param n_step number of model timesteps
#' @param n_en number of ensembles
get_Y_vector = function(n_states, n_step, n_en){
  
  Y_ic = array(dim = c(n_states, n_step, n_en))
  Y_pred = array(dim = c(n_states, n_step, n_en))
  
  return(list(Y_ic = Y_ic, Y_pred = Y_pred))
}

#### Function to create observation error matrix----
#' observation error matrix, should be a square matrix where
#'   col & row = the number of states and params for which you have observations
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_param_obs number of parameters for which we have observations
#' @param n_step number of model timesteps
#' @param state_sd vector of state observation standard deviation; assuming sd is constant through time
#' @param param_sd vector of parameter observation standard deviation; assuming sd is constant through time
get_obs_error_matrix = function(n_states, n_step, state_sd){
  
  R = array(0, dim = c(n_states, n_states, n_step))
  
  state_var = state_sd^2 #variance of chl-a observations
  
  for(i in 1:n_step){
    # variance is the same for each depth and time step; could make dynamic or varying by time step if we have good reason to do so
    R[,,i] = diag(state_var, n_states, n_states)
  }
  
  return(R)
}

#### Function to create matrix that identifies when observations are available----
#' Measurement operator matrix saying 1 if there is observation data available, 0 otherwise
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_param_obs number of parameters for which we have observations
#' @param n_params_est number of parameters we're calibrating
#' @param n_step number of model timesteps
#' @param obs observation matrix created with get_obs_matrix function
get_obs_id_matrix = function(n_states, n_step, obs){
  
  H = array(0, dim=c(n_states, n_states, n_step))
  
  # order goes 1) states, 2)params for which we have obs, 3) params for which we're estimating but don't have obs
  
  for(t in 1:n_step){
    H[1:(n_states), 1:(n_states), t] = diag(ifelse(is.na(obs[,,t]),0, 1), n_states, n_states)
  }
  
  return(H)
}

#### Function to turn observation dataframe into matrix----
#' turn observation dataframe into matrix
#'
#' @param obs_df observation data frame
#' @param model_dates dates over which you're modeling
#' @param n_step number of model time steps
#' @param n_states number of states we're updating in data assimilation routine
#' @param states character string vector of state names in obs_file
get_obs_matrix = function(obs_df, model_dates, n_step, n_states, states){
  
  # need to know location and time of observation
  
  obs_df_filtered = obs_df %>%
    dplyr::filter(as.Date(datetime) %in% model_dates) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, chla) %>%
    mutate(date_step = which(model_dates %in% date))
  obs_matrix = array(NA, dim = c(n_states, 1, n_step))
  for(i in 1:n_states){
    for(j in obs_df_filtered$date_step){
      obs_matrix[i, 1, j] = dplyr::filter(obs_df_filtered,
                                          date_step == j) %>%
        pull(states[i])
    }}
  
  return(obs_matrix)
}

#### Kalman filter function----
##' @param Y vector for holding states and parameters you're estimating
##' @param R observation error matrix
##' @param obs observations at current timestep
##' @param H observation identity matrix
##' @param n_en number of ensembles
##' @param cur_step current model timestep
kalman_filter = function(Y, R, obs, H, n_en, cur_step){
  
  cur_obs = obs[ , , cur_step]
  
  cur_obs = ifelse(is.na(cur_obs), 0, cur_obs) # setting NA's to zero so there is no 'error' when compared to estimated states
  
  ###### estimate the spread of your ensembles #####
  Y_mean = matrix(mean(Y[ , cur_step, ], na.rm = TRUE), nrow = length(Y[ , 1, 1])) # calculating the mean of each temp and parameter estimate
  delta_Y = Y[ , cur_step, ] - matrix(rep(Y_mean, n_en), nrow = length(Y[ , 1, 1])) # difference in ensemble state/parameter and mean of all ensemble states/parameters
  
  ###### estimate Kalman gain #########
  K = ((1 / (n_en - 1)) * delta_Y %*% t(delta_Y) %*% t(H[, , cur_step])) %*%
    qr.solve(((1 / (n_en - 1)) * H[, , cur_step] %*% delta_Y %*% t(delta_Y) %*% t(H[, , cur_step]) + R[, , cur_step]))
  
  ###### update Y vector ######
  for(q in 1:n_en){
    Y[, cur_step, q] = Y[, cur_step, q] + K %*% (cur_obs - H[, , cur_step] %*% Y[, cur_step, q]) # adjusting each ensemble using kalman gain and observations
  }
  return(Y)
}

#### Function to initialize state and parameter vector----
#' initialize Y vector with draws from distribution of obs
#'
#' @param Y Y vector
#' @param obs observation matrix
initialize_Y = function(Y, obs, n_states_est, n_step, n_en, state_sd, yini){
  
  # initializing states with earliest observations and parameters
  first_obs = yini 
  
  Y$Y_ic[ , 1, ] = array(abs(rnorm(n = n_en * (n_states_est),
                                   mean = c(first_obs),
                                   sd = c(state_sd))),
                         dim = c(c(n_states_est), n_en))
  
  Y$Y_pred[ , 1, ] = array(abs(rnorm(n = n_en * (n_states_est),
                                     mean = c(first_obs),
                                     sd = c(state_sd))),
                           dim = c(c(n_states_est), n_en))
  
  return(Y)
}

#### EnKF wrapper----
#' wrapper for running EnKF 
#' 
#' @param n_en number of model ensembles 
#' @param start start date of model run 
#' @param stop date of model run
#' @param forecast_data observation file 
#' @param ic_sd coefficient of variation of observations 
#' @param ic vector of initial conditions for states (chla)
#' @param model R object of fitted autoregressive model
#' @param residuals vector of residuals from model fit
run_forecasts = function(n_en = 30, 
                start = '2020-09-25', # start date 
                stop = '2020-10-29', 
                forecast_data = lake_data,
                ic_sd = c(0.1),
                ic = yini,
                model = ar_model,
                residuals = residuals){
  
  
  n_en = n_en
  start = as.Date(start)
  stop = as.Date(stop)
  dates = seq.Date(from = as.Date(start), to = as.Date(stop), by = "days")
  n_step = length(dates)
  ar_model = model
  residuals = residuals
  
  # get observation matrix
  obs_df = forecast_data %>% 
    select(datetime, chla) 
  
  yini <- c( #initial estimate of chl-a
    chla = ic[1]) #ug/L
  
  #define sds
  state_sd = ic_sd
  init_cond_sd = ic_sd

  # setting up matrices
  # observations as matrix
  obs = get_obs_matrix(obs_df = obs_df,
                       model_dates = dates,
                       n_step = n_step,
                       n_states = 1,
                       states = "chla")
  
  # Y vector for storing state estimates and updates
  Y = get_Y_vector(n_states = 1,
                   n_step = n_step,
                   n_en = n_en)
  Y_ic = Y$Y_ic
  Y_pred = Y$Y_pred
  
  # observation error matrix
  R = get_obs_error_matrix(n_states = 1,
                           n_step = n_step,
                           state_sd = state_sd)
  
  # observation identity matrix
  H = get_obs_id_matrix(n_states = 1,
                        n_step = n_step,
                        obs = obs)
  
  # initialize Y vector
  Y = initialize_Y(Y = Y, obs = obs, n_states_est = 1,
                   n_step = n_step, n_en = n_en, state_sd = init_cond_sd, yini = yini)
  Y_ic = Y$Y_ic
  Y_pred = Y$Y_pred
  
  if(any(!is.na(obs[ , , 1]))){
    Y_ic = kalman_filter(Y = Y_ic,
                         R = R,
                         obs = obs,
                         H = H,
                         n_en = n_en,
                         cur_step = 1) # updating params / states if obs available
  }
  
  
  # start modeling
  for(t in 2:n_step){
    
    for(n in 1:n_en){

      # pull parameter values from a distribution
      ar1 = rnorm(n = 1, mean = ar_model$ar, sd = ar_model$asy.se.coef$ar)
      chla_mean = rnorm(n = 1, mean = ar_model$x.mean, sd = ar_model$asy.se.coef$x.mean)
      intercept = c(ar_model$x.intercept)
      
      # define sigma
      sigma = sd(residuals, na.rm = TRUE)
      
      # draw process error value from a distribution
      W = rnorm(n = 1, mean = 0, sd = sigma)

      # run model; 
      chla_pred <- intercept + ar1 * (Y_ic[1, t-1, n] - chla_mean) + chla_mean + W
      
      # put prediction in Y matrices
      Y_ic[1 , t, n] = chla_pred
      Y_pred[1 , t, n] = chla_pred
      
      }
      
      
    # check if there are any observations to assimilate 
    if(any(!is.na(obs[ , , t]))){
      Y_ic = kalman_filter(Y = Y_ic,
                           R = R,
                           obs = obs,
                           H = H,
                           n_en = n_en,
                           cur_step = t) # updating params / states if obs available
    }
  }
  out = list(Y_ic = Y_ic, Y_pred = Y_pred, dates = dates, R = R, obs_file = forecast_data, state_sd = state_sd)
  
  return(out)
}


### Plotting functions ----
#### Function to plot NEON chl-a data----

# Plot chl-a data 
#' @param lake_data NEON lake dataset 
plot_chla_obs <- function(lake_data){
  p <- ggplot(data = lake_data, aes(x = datetime, y = chla))+
    geom_line(aes(color = "Chl-a"))+
    xlab("")+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    scale_color_manual(values = c("Chl-a" = "chartreuse4"), name = "")+
    theme_bw()
  return(p)
}

# Plot chl-a data + 1 day lag: timeseries
#' @param plot_data chlorophyll-a dataframe with datetime, chla and chla_lag

plot_chla_lag <- function(plot_data){
  p <- ggplot(data = plot_data)+
    geom_line(aes(x = datetime, y = chla_lag, color = "1 day lag of chlorophyll"))+
    geom_line(aes(x = datetime, y = chla, color = "chlorophyll"))+
    xlab("2018")+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    scale_color_manual(values = c("chlorophyll" = "darkgreen","1 day lag of chlorophyll" = "lightgreen"))+
    labs(color = NULL)+
    theme_bw()+
    theme(legend.position = "bottom")
  return(p)
}

#### Function to plot AR model fit ----
#' @param model_fit_plot_data data frame of lake observations and model predictions
#'  
plot_mod_predictions_chla <- function(model_fit_plot_data){
  cols <- RColorBrewer::brewer.pal(8, "Dark2") # Set custom color palette for our plot - ooh we are fancy!! :-)
  
  ggplot(data = model_fit_plot_data) +
    geom_point(aes(date, chla, color = "Observed")) +
    geom_line(aes(date, model, color = "Modeled")) +
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")"))) +
    xlab("Time") +
    scale_color_manual(values = c( "Observed" = "black", "Modeled" = cols[4]),
                       name = "",
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid","blank"),
                         shape = c(NA,16)))) +
    theme_bw(base_size = 12) 
}

#### Function to plot distribution of initial conditions ----
#'@param curr_chla mean initial condition of chla
#'@param ic_uc vector draws from a distribution of initial conditions
#'
plot_ic_dist <- function(curr_chla, ic_uc){
  
  #Set colors
  l.cols <- RColorBrewer::brewer.pal(8, "Set2")[-c(1, 2)] # Defining another custom color palette :-)
  
  #Build plot
  ggplot() +
    # geom_vline(data = df, aes(xintercept = x, color = label)) +
    geom_vline(xintercept = curr_chla) +
    geom_density(aes(ic_uc), fill = l.cols[2], alpha = 0.3) +
    xlab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")"))) +
    ylab("Density") +
    theme_bw(base_size = 18)+
    ggtitle("Initial condition distribution")
}



#### Functions to plot chl-a forecast----
#' plot chlorophyll forecast and observations 
#' 
#' @param est_out forecast output from EnKF wrapper 
#' @param lake_data NEON data for selected site formatted using format_enkf_inputs function
#' @param obs_file data frame of observations
#' @param start start date
#' @param stop stop date
#' @param n_en number of ensemble members
#' 
plot_chla = function(est_out, lake_data, obs_file, start, stop, n_en){
  
  plot_dates <- seq.Date(from = as.Date(start), to = as.Date(stop), by = "days")
  
  assim_data <- obs_file %>%
    rename(obs_assimilated = chla) %>%
    # mutate(assim_key1 = ifelse(!is.na(obs_assimilated) | !is.na(lag(obs_assimilated)), 1, 0),
    #        assim_key2 = ifelse(!is.na(obs_assimilated) | !is.na(lead(obs_assimilated)), 1, 0),
    #        assim_key3 = ifelse(!is.na(obs_assimilated), 1, 0),
    #        obs_assimilated = exp(obs_assimilated))
    mutate(assim_key1 = ifelse(!is.na(obs_assimilated), 1, 0),
           assim_key2 = ifelse(!is.na(lag(obs_assimilated)), 1, 0),
           obs_assimilated = obs_assimilated)
  
  chla_pred = NULL
  chla_ic = NULL
  for(i in 1:length(plot_dates)){
    chla_pred <- c(chla_pred, est_out$Y_pred[1,i,])
    chla_ic <- c(chla_ic, est_out$Y_ic[1,i,])
  }
  
  plot_data <- tibble(datetime = rep(plot_dates, each = n_en),
                      ensemble_member = rep(1:n_en, times = length(plot_dates)),
                      chla_pred = chla_pred,
                      chla_ic = chla_ic) %>%
    left_join(., lake_data, by = "datetime") %>%
    left_join(., assim_data, by = "datetime")
  
  plot_data2 <- plot_data %>%
    add_column(datefactor = as.factor(format(as.Date(plot_data$datetime), "%m-%d"))) %>%
    mutate(chla_fc = ifelse(datefactor == "09-25", NA, chla_pred))
  
  ens <- plot_data2 %>%
    select(datetime, ensemble_member)
  
  for(d in 1:(length(plot_dates)-1)){
    
    obs <- obs_file %>% filter(datetime == plot_dates[d])
    
    if(!is.na(obs$chla)){
      
      ic <- plot_data2 %>%
        filter(datetime == ymd(plot_dates[d])) %>%
        pull(chla_ic)
      
      pred <- plot_data2 %>%
        filter(datetime == ymd(plot_dates[d])+1) %>%
        pull(chla_pred)
      
      temp <- tibble(datetime = rep(seq(from = ymd(plot_dates[d]), to = ymd(plot_dates[d])+1, by = "day"), each = n_en),
                     ensemble_member = rep(1:n_en, times = 2),
                     chla_ens_today = c(ic, pred))
      
      newName <- setNames("chla_ens_today", paste0("chla_ens_",plot_dates[d]))
      
      ens <- ens %>%
        left_join(., temp, by = c("datetime","ensemble_member")) %>%
        rename(all_of(newName))
    } else {
      ic <- plot_data2 %>%
        filter(datetime == ymd(plot_dates[d])) %>%
        pull(chla_pred)
      
      pred <- plot_data2 %>%
        filter(datetime == ymd(plot_dates[d])+1) %>%
        pull(chla_pred)
      
      temp <- tibble(datetime = rep(seq(from = ymd(plot_dates[d]), to = ymd(plot_dates[d])+1, by = "day"), each = n_en),
                     ensemble_member = rep(1:n_en, times = 2),
                     chla_ens_today = c(ic, pred))
      
      newName <- setNames("chla_ens_today", paste0("chla_ens_",plot_dates[d]))
      
      ens <- ens %>%
        left_join(., temp, by = c("datetime","ensemble_member")) %>%
        rename(all_of(newName))
    }
    
  }
  
  ens <- ens %>%
    add_column(datefactor = as.factor(format(as.Date(plot_data$datetime), "%m-%d"))) %>%
    pivot_longer(cols = -c(datetime, ensemble_member, datefactor), names_to = "ensemble_date",
                 values_to = "ens")  %>%
    filter(complete.cases(.)) 
  
  p <- ggplot()+
    geom_line(data = ens, aes(x = datefactor, y = ens, group = interaction(ensemble_date, ensemble_member), color = "Ensemble members"))+
    geom_violin(data = plot_data2, aes(x = datefactor, y = chla_fc, fill = "Forecast for today"), color = "black",
                scale = "width", width = 0.7)+
    geom_violin(data = plot_data2, aes(x = datefactor, y = chla_ic, fill = "Initial condition for tomorrow's forecast"), color = "cornflowerblue", alpha = 0.4, scale = "width", width = 0.7)+
    geom_point(data = plot_data2, aes(x = datefactor, y = chla, color = "Non-assimilated observations"))+
    geom_point(data = plot_data2, aes(x = datefactor, y = obs_assimilated, color = "Assimilated observations"))+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values = c("Ensemble members" = "lightgray",
                                  "Non-assimilated observations" = "black",
                                  "Assimilated observations" = "orange"), 
                       name = "",
                       guide = guide_legend(override.aes = list(
                         linetype = c("blank","solid","blank"),
                         shape = c(16,NA, 16))))+
    scale_fill_manual(values = c("Forecast for today" = "white", "Initial condition for tomorrow's forecast" = "cornflowerblue"),
                      name = "",
                      guide = guide_legend(override.aes = list(
                        color = c("black","cornflowerblue"))))+
    ggtitle("1-day-ahead forecasts generated every day")
  
  return(p)
}

#### Function to plot predictions vs. observations for chl-a----
#' plot chl-a mean forecast and observations 
#' 
#' @param est_out forecast output from EnKF wrapper 
#' @param lake_data NEON data for selected site 
#' 
pred_v_obs_chla = function(forecasts, lake_data){
  mean_chla_est = apply(forecasts$Y_pred[1,,] , 1, FUN = mean)
  
  # this could be used to show a 95% confidence error bar on predicted
  # but I think the error bars make the plot a bit hard to read
  # top_din_est = apply(forecasts$Y[2,,] , 1, FUN = quantile, probs=c(0.975))
  # bottom_din_est = apply(forecasts$Y[2,,] , 1, FUN = quantile, probs=c(0.025))
  lake_data <- lake_data %>%
    mutate(datetime = as.Date(datetime)) %>%
    filter(datetime %in% forecasts$dates) 

  plot(mean_chla_est ~ lake_data$chla, type ='p', 
       ylim = c(min(c(range(mean_chla_est,na.rm = TRUE),range(lake_data$chla,na.rm = TRUE))),c(max(c(range(mean_chla_est,na.rm = TRUE),range(lake_data$chla,na.rm = TRUE))))),
       xlim = c(min(c(range(mean_chla_est,na.rm = TRUE),range(lake_data$chla,na.rm = TRUE))),c(max(c(range(mean_chla_est,na.rm = TRUE),range(lake_data$chla,na.rm = TRUE))))),
       col = 'black', ylab = 'predicted Chl-a (ug L-1)', xlab = 'observed Chl-a (ug L-1)')
  abline(a=0, b=1)
  
  #this code could be used to show error bars representing uncertainty but 
  #I think it makes the plot kinda hard to read
  # arrows(lake_data$din[1:35]  - forecasts$state_sd[2], mean_din_est,
  #        lake_data$din[1:35]  + forecasts$state_sd[2], mean_din_est, 
  #        code = 3, length = 0.1, angle = 90, col = 'black')
  # arrows(lake_data$din, mean_din_est  - forecasts$state_sd[2], 
  #        lake_data$din, mean_din_est  + forecasts$state_sd[2], 
  #        code = 3, length = 0.1, angle = 90, col = 'black')
}


####Function to plot DA frequency experiment results----
#'@param da_frequency_experiment_output list of forecast output with different DA frequencies
#'@param chla_assimilation_frequencies vector of chla assimilation frequencies in days
#'
plot_da_frequency_experiment_results <- function(da_frequency_experiment_output,
                                                 chla_assimilation_frequencies){
  
  rmse <- rep(NA,length(chla_assimilation_frequencies))
  
  for(i in 1:length(da_frequency_experiment_output)){
    
  forecast = apply(da_frequency_experiment_output[[i]]$Y_pred[1,,] , 1, FUN = mean)
  
  #limit obs to forecast dates
  forecast_obs <- lake_data %>%
    mutate(datetime = as.Date(datetime)) %>%
    filter(datetime %in% da_frequency_experiment_output[[i]]$dates) 
  
  #calculate RMSE
  rmse[i] <- sqrt(mean((forecast_obs$chla - forecast)^2, na.rm = TRUE))
  }
  
  #create plot data
  plot_data <- tibble(da_freq = chla_assimilation_frequencies,
                      rmse = rmse)
  
  #plot
  p <- ggplot(data = plot_data, aes(x = da_freq, y = rmse)) +
    geom_point(shape = 19, size = 2)+
    xlab("Data assimilation frequency (days)")+
    ylab(expression(paste("RMSE (",mu,g,~L^-1,")")))+
    scale_x_continuous(breaks = chla_assimilation_frequencies)+
    theme_bw()
  
  return(p)
}

####Function to plot obs uncertainty experiment results----
#'@param obs_uncertainty_experiment_output list of forecast output with different levels of observation uncertainty
#'@param obs_uncertainty vector of observation uncertaintys in log ugl

plot_obs_uncertainty_experiment_results <- function(obs_uncertainty_experiment_output,
                                                    obs_uncertainty){
  
  rmse <- rep(NA,length(obs_uncertainty))
  
  for(i in 1:length(obs_uncertainty_experiment_output)){
    
    forecast = apply(obs_uncertainty_experiment_output[[i]]$Y_pred[1,,] , 1, FUN = mean)
    
    #limit obs to forecast dates
    forecast_obs <- lake_data %>%
      mutate(datetime = as.Date(datetime)) %>%
      filter(datetime %in% obs_uncertainty_experiment_output[[i]]$dates) 
    
    #calculate RMSE
    rmse[i] <- sqrt(mean((forecast_obs$chla - forecast)^2, na.rm = TRUE))
  }
  
  #create plot data
  plot_data <- tibble(obs_uc = obs_uncertainty,
                      rmse = rmse)
  
  #plot
  p <- ggplot(data = plot_data, aes(x = obs_uc, y = rmse)) +
    geom_point(shape = 19, size = 2)+
    xlab(expression(paste("Observation uncertainty (log ",mu,g,~L^-1,")")))+
    ylab(expression(paste("RMSE (",mu,g,~L^-1,")")))+
    scale_x_continuous(breaks = obs_uncertainty)+
    theme_bw()
  
  return(p)
}

