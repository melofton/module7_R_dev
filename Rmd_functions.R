# Functions needed to complete Module 7: Using Data to Improve Ecological Forecasts
# Author: Mary Lofton
# Date: 06APR23

# Load packages
library(tidyverse)
library(lubridate)
library(zoo)
library(mvtnorm)
library(see) #see annotation below on how to avoid using this package if it becomes an issue

##Notes on see package:
#' This is literally only used for one plotting command: geom_violinhalf().
#' To completely avoid use of this package, ctrl+F for geom_violinhalf() and
#' replace with geom_violin() - VOILA!

## Define functions----

EnKF <- function(forecast, new_observation, ic_sd){
  
  #Allocate matrices
  x_corr <- matrix(forecast)
  y <- matrix(new_observation)
  h_matrix <- matrix(0, nrow = 1, ncol = 1)
  R_matrix <- matrix(0, nrow = 1, ncol = 1)
  dit <- matrix(NA, nrow = length(x_corr[,1]), ncol = 1) 
  y_corr <- matrix(NA, nrow =  length(x_corr[,1]), ncol = length(y))
  x_update <- matrix(NA, nrow = length(x_corr[,1]), ncol = 1)
  
  #Only do EnKF if observations are present that day
  #there has to be at least 1 non-NA observation.
  if(length(which(!is.na(y))) > 0){
    
    #Assign observations to depths
    h_matrix[1, 1] <- 1

    #Create observational uncertainty matrix
    R_matrix[1,1] <- ic_sd^2

    #Calculate mean prediction for each depth
    ens_mean <- colMeans(x_corr)
    
    #Loop through ensemble members
    for(m in 1:length(x_corr[,1])){  
      #Ensemble specific deviation
      dit[m, ] <- x_corr[m, ] - ens_mean
      
      #if the first ensemble then create the matrix that is then averaged
      if(m == 1){
        p_it <- dit[m, ] %*% t(dit[m, ]) 
      }else{
        #if not the first ensemble then add the matrix to the previous matrix
        p_it <- dit[m, ] %*% t(dit[m, ]) +  p_it 
      }
    }
    
    #Calculate Cxx matrix
    Cxx_matrix <- p_it / (length(x_corr[,1]) - 1)
    
    #Add noise to observations
    for(m in 1:length(x_corr[,1])){
      y_corr[m, ] <- y + t(rmvnorm(n = 1, mean = c(0), sigma = R_matrix))
    }
    
    #Calculate Kalman Gain
    K <- Cxx_matrix %*% t(h_matrix) %*% solve(h_matrix %*% Cxx_matrix %*% t(h_matrix) + R_matrix)
    
    #Update model states based on Kalman Gain and devivations
    for(m in 1:length(x_corr[,1])){
      x_update[m, ] <- x_corr[m,] + K %*% (y_corr[m,] - h_matrix %*% x_corr[m,])
    }
  }else{
    #Only add noise if observations are missing
    x_update <- x_corr
  }
  
  ic_update <- c(x_update[,1])
  return(ic_update)
}

### Plotting functions ----
#### Function to plot NEON chl-a data----

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
plot_mod_predictions <- function(model_fit_plot_data, variable_name){
  cols <- RColorBrewer::brewer.pal(8, "Dark2") # Set custom color palette for our plot - ooh we are fancy!! :-)
  
  ggplot(data = model_fit_plot_data) +
    geom_point(aes(date, chla, color = "Observed")) +
    geom_line(aes(date, model, color = "Modeled")) +
    ylab(variable_name) +
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

#### Function to plot distribution of process uncertainty ----
#'@param proc_uc vector draws from a distribution of proc uc
#'
plot_process_dist <- function(proc_uc){
  
  #Set colors
  l.cols <- RColorBrewer::brewer.pal(8, "Set2")[-c(1, 2)] # Defining another custom color palette :-)
  
  #Build plot
  ggplot() +
    geom_vline(xintercept = 0) +
    geom_density(aes(proc_uc), fill = l.cols[1], alpha = 0.3) +
    xlab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")"))) +
    ylab("Density") +
    theme_bw(base_size = 18)+
    ggtitle("Process uncertainty distribution")
}

#### Function to plot distribution of a forecast ----
#'@param forecast_dist vector of forecast distribution
#'
plot_fc_dist <- function(forecast_dist){
  
  #Set colors
  l.cols <- RColorBrewer::brewer.pal(8, "Set2")[-c(1, 2)] # Defining another custom color palette :-)
  #Build plot
  ggplot() +
    geom_density(aes(forecast_dist), fill = l.cols[3], alpha = 0.3) +
    xlab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")"))) +
    ylab("Density") +
    theme_bw(base_size = 18)+
    ggtitle("Chl-a forecast distribution")
}



#### Functions to plot chl-a forecast----

#' One-day forecast plot
#' 

plot_fc_1day <- function(curr_chla, start_date, forecast_date, ic_distribution, forecast_chla, n_members){
  
  ens <- tibble(date = c(rep(start_date, times = length(ic_distribution)),
                         rep(forecast_date, times = length(forecast_chla))),
                ens = c(ic_distribution, forecast_chla),
                ensemble_member = rep(1:n_members, times = 2))
  ic <- ens %>%
    filter(date == start_date)
  fc <- ens %>%
    filter(date == forecast_date)
  obs <- tibble(date = start_date,
                obs = curr_chla)
  
  p <- ggplot()+
    geom_line(data = ens, aes(x = date, y = ens, group = ensemble_member, color = "Ensemble members"))+
    geom_violinhalf(data = fc, aes(x = date, y = ens, fill = "Forecast"), color = "black",
                scale = "width", width = 0.7)+
    geom_violinhalf(data = ic, aes(x = date, y = ens, fill = "Initial condition"), color = "cornflowerblue", alpha = 0.4, scale = "width", width = 0.7)+
    geom_point(data = obs, aes(x = date, y = obs, color = "Observation"), size = 3)+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values = c("Ensemble members" = "lightgray",
                                  "Observation" = "orange"), 
                       name = "",
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid","blank"),
                         shape = c(NA, 16))))+
    scale_fill_manual(values = c("Forecast" = "white", "Initial condition" = "cornflowerblue"),
                      name = "",
                      guide = guide_legend(override.aes = list(
                        color = c("black","cornflowerblue"))))+
    ggtitle("1-day-ahead forecast")
  
  return(p)
}

#' One-day forecast plot
#' 

plot_fc_update <- function(chla_obs, start_date, forecast_date, ic_distribution, ic_update, forecast_chla, n_members){
  
  ens <- tibble(date = c(rep(start_date, times = length(ic_distribution)),
                         rep(forecast_date, times = length(forecast_chla)*2)),
                ens = c(ic_distribution, forecast_chla, ic_update),
                ensemble_member = rep(1:n_members, times = 3),
                data_type = c(rep("ic", times = length(ic_distribution)),
                              rep("fc", times = length(forecast_chla)),
                              rep("ic", times = length(ic_update))))
  ic <- ens %>%
    filter(data_type == "ic")
  fc <- ens %>%
    filter(data_type == "fc")
  obs <- tibble(date = c(start_date, forecast_date),
                obs = chla_obs)
  
  p <- ggplot()+
    geom_line(data = ens, aes(x = date, y = ens, group = ensemble_member, color = "Ensemble members"))+
    geom_violinhalf(data = fc, aes(x = date, y = ens, fill = "Forecast"), color = "black",
                scale = "width", width = 0.7)+
    geom_violinhalf(data = ic, aes(x = date, y = ens, fill = "Initial condition"), color = "cornflowerblue", alpha = 0.4, scale = "width", width = 0.7)+
    geom_point(data = obs, aes(x = date, y = obs, color = "Observation"), size = 3)+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values = c("Ensemble members" = "lightgray",
                                  "Observation" = "orange"), 
                       name = "",
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid","blank"),
                         shape = c(NA, 16))))+
    scale_fill_manual(values = c("Forecast" = "white", "Initial condition" = "cornflowerblue"),
                      name = "",
                      guide = guide_legend(override.aes = list(
                        color = c("black","cornflowerblue"))))+
    ggtitle("1-day-ahead forecast")
  
  return(p)
}

plot_second_forecast <- function(chla_obs, start_date, forecast_dates, ic_distribution, ic_update, forecast_chla, second_forecast, n_members){
  
  ens <- tibble(date = c(rep(start_date, times = length(ic_distribution)),
                         rep(forecast_dates[1], times = length(forecast_chla)*2),
                         rep(forecast_dates[2], times = length(second_forecast))),
                ens = c(ic_distribution, forecast_chla, ic_update, second_forecast),
                ensemble_member = rep(1:n_members, times = 4),
                data_type = c(rep("ic", times = length(ic_distribution)),
                              rep("fc", times = length(forecast_chla)),
                              rep("ic", times = length(ic_update)),
                              rep("fc", times = length(second_forecast))))
  ic <- ens %>%
    filter(data_type == "ic")
  fc <- ens %>%
    filter(data_type == "fc")
  obs <- tibble(date = c(start_date, forecast_dates[1]),
                obs = chla_obs)
  
  p <- ggplot()+
    geom_line(data = ens, aes(x = date, y = ens, group = ensemble_member, color = "Ensemble members"))+
    geom_violinhalf(data = fc, aes(x = date, y = ens, fill = "Forecast"), color = "black",
                    scale = "width", width = 0.7)+
    geom_violinhalf(data = ic, aes(x = date, y = ens, fill = "Initial condition"), color = "cornflowerblue", alpha = 0.4, scale = "width", width = 0.7)+
    geom_point(data = obs, aes(x = date, y = obs, color = "Observation"), size = 3)+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values = c("Ensemble members" = "lightgray",
                                  "Observation" = "orange"), 
                       name = "",
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid","blank"),
                         shape = c(NA, 16))))+
    scale_fill_manual(values = c("Forecast" = "white", "Initial condition" = "cornflowerblue"),
                      name = "",
                      guide = guide_legend(override.aes = list(
                        color = c("black","cornflowerblue"))))+
    ggtitle("Two 1-day-ahead forecasts")
  
  return(p)
}

plot_many_forecasts <- function(forecast_data, forecast_series){

  forecast_dates <- unique(forecast_series$date)
  forecast_series <- forecast_series %>%
    mutate(datefactor = as.factor(format(date, "%m-%d")),
           chla = ifelse((date == last(forecast_dates) & data_type == "ic"),NA,chla))
  ic <- forecast_series %>%
    filter(data_type == "ic")
  fc <- forecast_series %>%
    filter(data_type == "fc")
  obs <- tibble(date = forecast_data$datetime,
                obs = forecast_data$chla) %>%
    mutate(datefactor = as.factor(format(date, "%m-%d")),
           obs = ifelse(date == last(forecast_dates),NA,obs))
  
  p <- ggplot()+
    geom_line(data = forecast_series, aes(x = datefactor, y = chla, group = ensemble_member, color = "Ensemble members"))+
    geom_violinhalf(data = fc, aes(x = datefactor, y = chla, fill = "Forecast"), color = "black",
                    scale = "width", width = 0.7)+
    geom_violinhalf(data = ic, aes(x = datefactor, y = chla, fill = "Initial condition"), color = "cornflowerblue", alpha = 0.4, scale = "width", width = 0.7)+
    geom_point(data = obs, aes(x = datefactor, y = obs, color = "Observations"), size = 3)+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values = c("Ensemble members" = "lightgray",
                                  "Observations" = "orange"), 
                       name = "",
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid","blank"),
                         shape = c(NA, 16))))+
    scale_fill_manual(values = c("Forecast" = "white", "Initial condition" = "cornflowerblue"),
                      name = "",
                      guide = guide_legend(override.aes = list(
                        color = c("black","cornflowerblue"))))+
    ggtitle("A series of 1-day-ahead forecasts")
  
  return(suppressWarnings(print(p)))
}

plot_many_forecasts_with_obs <- function(forecast_data, forecast_series, chla_observations){
  
  forecast_dates <- unique(forecast_series$date)
  forecast_series <- forecast_series %>%
    mutate(datefactor = as.factor(format(date, "%m-%d")),
           chla = ifelse((date == last(forecast_dates) & data_type == "ic"),NA,chla))
  ic <- forecast_series %>%
    filter(data_type == "ic")
  fc <- forecast_series %>%
    filter(data_type == "fc")
  obs_assim <- tibble(date = forecast_data$datetime,
                obs = forecast_data$chla) %>%
    mutate(datefactor = as.factor(format(date, "%m-%d")),
           obs = ifelse(date == last(forecast_dates),NA,obs))
  obs_not_assim <- chla_observations %>%
    mutate(datefactor = as.factor(format(datetime, "%m-%d")))
  
  p <- ggplot()+
    geom_line(data = forecast_series, aes(x = datefactor, y = chla, group = ensemble_member, color = "Ensemble members"))+
    geom_violinhalf(data = fc, aes(x = datefactor, y = chla, fill = "Forecast"), color = "black",
                    scale = "width", width = 0.7)+
    geom_violinhalf(data = ic, aes(x = datefactor, y = chla, fill = "Initial condition"), color = "cornflowerblue", alpha = 0.4, scale = "width", width = 0.7)+
    geom_point(data = obs_not_assim, aes(x = datefactor, y = chla, color = "Obs. - not assimilated"), size = 3, shape = 21)+
    geom_point(data = obs_assim, aes(x = datefactor, y = obs, color = "Obs. - assimilated"), size = 3)+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values = c("Ensemble members" = "lightgray",
                                  "Obs. - not assimilated" = "black",
                                  "Obs. - assimilated" = "orange"), 
                       name = "",
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid","blank", "blank"),
                         shape = c(NA, 16, 21))))+
    scale_fill_manual(values = c("Forecast" = "white", "Initial condition" = "cornflowerblue"),
                      name = "",
                      guide = guide_legend(override.aes = list(
                        color = c("black","cornflowerblue"))))+
    ggtitle("A series of 1-day-ahead forecasts")
  
  return(suppressWarnings(print(p)))
}

plot_scenario_forecasts <- function(forecast_data, ens, show_final_obs=NULL){
  
  forecast_dates <- unique(ens$date)
  ens <- ens %>%
    mutate(datefactor = as.factor(format(date, "%m-%d")),
           ens = ifelse((date == last(forecast_dates) & data_type == "ic"),NA,ens))
  ic <- ens %>%
    filter(data_type == "ic")
  fc <- ens %>%
    filter(data_type == "fc")
  obs <- tibble(date = forecast_data$datetime,
                obs = forecast_data$chla) %>%
    mutate(datefactor = as.factor(format(date, "%m-%d")))
  
  if(is.null(show_final_obs)){
    obs <- obs %>%
      mutate(obs = ifelse(date == last(forecast_dates),NA,obs))
  }
  
  p <- ggplot()+
    geom_line(data = ens, aes(x = datefactor, y = ens, group = ensemble_member, color = "Ensemble members"))+
    geom_violinhalf(data = fc, aes(x = datefactor, y = ens, fill = "Forecast"), color = "black",
                    scale = "width", width = 0.7)+
    geom_violinhalf(data = ic, aes(x = datefactor, y = ens, fill = "Initial condition"), color = "cornflowerblue", alpha = 0.4, scale = "width", width = 0.7)+
    geom_point(data = obs, aes(x = datefactor, y = obs, color = "Observations"), size = 3)+
    geom_hline(aes(yintercept = 10, color = "Water quality threshold"))+
    ylab(expression(paste("Chlorophyll-a (",mu,g,~L^-1,")")))+
    xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values = c("Ensemble members" = "lightgray",
                                  "Observations" = "orange",
                                  "Water quality threshold" ="red"), 
                       name = "",
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid","blank","solid"),
                         shape = c(NA, 16, NA))))+
    scale_fill_manual(values = c("Forecast" = "white", "Initial condition" = "cornflowerblue"),
                      name = "",
                      guide = guide_legend(override.aes = list(
                        color = c("black","cornflowerblue"))))+
    ggtitle("A series of 1-day-ahead forecasts")
  
  return(p)
}

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

