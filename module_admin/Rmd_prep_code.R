# Prep code for Rmd version of Module 7
# Author: Mary Lofton
# Date: 05Jan24

# Purpose: prepare data files and other ancillary code for RMarkdown version of Mod 7

# Load packages
pacman::p_load(tidyverse, lubridate, data.table, zoo, sparklyr, neonUtilities, hms)

#define NEON token
source("./neon_token_source.R")

# Get high-frequency data to inform IC distribution
#read in site info
sites <- readr::read_csv("./data/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(field_site_id == "BARC")
lake_sites <- sites$field_site_id[(which(sites$field_site_subtype == "Lake"))]

#read in wq data

neon <- arrow::s3_bucket("neon4cast-targets/neon",
                         endpoint_override = "data.ecoforecast.org",
                         anonymous = TRUE)

####Chl-a ===============================================#######

wq_portal <- purrr::map_dfr(sites$field_site_id, function(site){
  message(site)
  arrow::open_dataset(neon$path("waq_instantaneous-basic-DP1.20288.001"), partitioning = "siteID") %>%   # waq_instantaneous
    #wq_portal <- neonstore::neon_table("waq_instantaneous", site = sites$field_site_id, lazy = TRUE) %>%   # waq_instantaneous
    dplyr::filter(siteID %in% site) %>%
    dplyr::select(siteID, startDateTime, sensorDepth,
                  chlorophyll, chlorophyllExpUncert,chlorophyllFinalQF,
                  chlaRelativeFluorescence, chlaRelFluoroFinalQF) %>%
    dplyr::mutate(sensorDepth = as.numeric(sensorDepth),
                  chla = as.numeric(chlorophyll),
                  RFU = as.numeric(chlaRelativeFluorescence),
                  chlorophyllExpUncert = as.numeric(chlorophyllExpUncert)) %>%
    dplyr::rename(site_id = siteID) %>% 
    dplyr::filter(((sensorDepth > 0 & sensorDepth < 1)| is.na(sensorDepth))) |> 
    collect() %>% # sensor depth of NA == surface?
    dplyr::mutate(startDateTime = as_datetime(startDateTime)) %>%
    # QF (quality flag) == 0, is a pass (1 == fail), 
    # make NA so these values are not used in the mean summary
    dplyr::mutate(chla = ifelse(chlorophyllFinalQF == 0, chla, NA),
                  RFU = ifelse(chlaRelFluoroFinalQF == 0, RFU, NA)) %>% 
    dplyr::select(startDateTime, site_id, chla, RFU) %>%
    dplyr::filter(site_id == "BARC" & month(startDateTime) %in% c(9:10) & year(startDateTime) == 2019)
}
)

#QC water quality data

wq_full <- wq_portal

# chlorophyll ranges
chla_max <- 200
chla_min <- 0

# GR flag will be true if either the DO concentration or the chlorophyll are 
# outside the ranges specified about

wq_cleaned <- wq_full %>%
  dplyr::mutate(chla = ifelse(is.na(chla),chla, ifelse(chla >= chla_min & chla <= chla_max, chla, NA)),
                RFU = ifelse(is.na(RFU),RFU, ifelse(RFU >= chla_min & RFU <= chla_max, RFU, NA))) 
# manual cleaning based on visual inspection

##############################################################################
ggplot(data = wq_cleaned, aes(x = startDateTime, y = chla))+
  geom_point()+
  facet_wrap(vars(site_id), scales = "free", nrow = 2)+
  theme_bw()
##############################################################################

chla_final <- wq_cleaned %>%
  rename(datetime = startDateTime) %>%
  select(datetime, chla) %>%
  filter(complete.cases(.))

write.csv(chla_final, "./data/BARC_chla_microgramsPerLiter_highFrequency.csv",row.names = FALSE)

##########WATER TEMP===============================================#####

#define function for temperature QC
#### Temperature data QC ####
# Function to QC temperature profiles 
# option to group by site and depth or just by site
# need to specify the 1) absolute range (vector - length 2) and the 2) the acceptable rate of change (absolute value)
QC.temp <- function(df, range, spike, by.depth = T) {
  if (by.depth == F) {
    df_QC <- df %>%
      arrange(site_id, time) %>%
      group_by(site_id) %>%
      mutate(temp_change = observation - dplyr::lag(observation), 
             
             flagged_abs_val = ifelse(between(observation, min(range), max(range)), F, T),
             flagged_spike = ifelse(abs(temp_change) > spike & time != first(time),
                                    T, F),
             flagged_flat = ifelse(temp_change == 0  & time != first(time),
                                   T, F),
             final_flag = ifelse(flagged_spike != F | 
                                   flagged_abs_val != F |
                                   flagged_flat != F,
                                 T, F),
             observation = ifelse(final_flag == T, 
                                  NA, observation)) %>%
      select(-c(contains('flag'), temp_change))
  } else {
    df_QC <- df %>%
      arrange(site_id, time) %>%
      group_by(site_id, depth) %>%
      mutate(temp_change = observation - dplyr::lag(observation), 
             
             flagged_abs_val = ifelse(between(observation, min(range), max(range)), F, T),
             flagged_spike = ifelse(abs(temp_change) > spike & time != first(time),
                                    T, F),
             flagged_flat = ifelse(temp_change == 0 & time != first(time),
                                   T, F),
             final_flag = ifelse((flagged_spike == T | 
                                    flagged_abs_val == T |
                                    flagged_flat == T),
                                 T, F),
             observation = ifelse((final_flag == F | is.na(final_flag)), 
                                  observation, NA)) %>%
      select(-c(contains('flag'), temp_change))
  }
  return(df_QC)
}

message("#### Generate hourly temperature profiles for lake #############")
message("##### NEON portal data #####")
hourly_temp_profile_portal <- arrow::open_dataset(neon$path("TSD_30_min-basic-DP1.20264.001"), partitioning = "siteID") %>%
  #hourly_temp_profile_portal <- neonstore::neon_table("TSD_30_min", site = sites$field_site_id) %>%
  dplyr::filter(siteID %in% lake_sites) %>%
  rename(site_id = siteID,
         depth = thermistorDepth) %>%
  dplyr::select(startDateTime, site_id, tsdWaterTempMean, depth, tsdWaterTempExpUncert, tsdWaterTempFinalQF, verticalPosition) %>%
  dplyr::collect() %>%
  # errors in the sensor depths reported - see "https://www.neonscience.org/impact/observatory-blog/incorrect-depths-associated-lake-and-river-temperature-profiles"
  # sensor depths are manually assigned based on "vertical position" variable as per table on webpage
  dplyr::mutate(depth = ifelse(site_id == "CRAM" & 
                                 lubridate::as_date(startDateTime) < lubridate::as_date("2020-11-01") & 
                                 verticalPosition == 502, 1.75, depth),
                depth = ifelse(site_id == "CRAM" & 
                                 lubridate::as_date(startDateTime) < lubridate::as_date("2020-11-01") & 
                                 verticalPosition == 503, 3.45, depth),
                depth = ifelse(site_id == "CRAM" & 
                                 lubridate::as_date(startDateTime) < lubridate::as_date("2020-11-01") & 
                                 verticalPosition == 504, 5.15, depth),
                depth = ifelse(site_id == "CRAM" & 
                                 lubridate::as_date(startDateTime) < lubridate::as_date("2020-11-01") & 
                                 verticalPosition == 505, 6.85, depth),
                depth = ifelse(site_id == "CRAM" & 
                                 lubridate::as_date(startDateTime) < lubridate::as_date("2020-11-01") & 
                                 verticalPosition == 506, 8.55, depth),
                depth = ifelse(site_id == "CRAM" & 
                                 lubridate::as_date(startDateTime) < lubridate::as_date("2020-11-01") & 
                                 verticalPosition == 505, 10.25, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) < lubridate::as_date("2021-06-08") | 
                                    lubridate::as_date(startDateTime) > lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 502, 0.55, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) < lubridate::as_date("2021-06-08") | 
                                    lubridate::as_date(startDateTime) > lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 503, 1.05, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) < lubridate::as_date("2021-06-08") | 
                                    lubridate::as_date(startDateTime) > lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 504, 1.55, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) < lubridate::as_date("2021-06-08") | 
                                    lubridate::as_date(startDateTime) > lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 505, 2.05, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) < lubridate::as_date("2021-06-08") | 
                                    lubridate::as_date(startDateTime) > lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 506, 2.55, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) < lubridate::as_date("2021-06-08") | 
                                    lubridate::as_date(startDateTime) > lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 507, 3.05, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") & 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 502, 0.3, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") & 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 503, 0.55, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") & 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 504, 0.8, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") & 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 505, 1.05, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") & 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 506, 1.3, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") & 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 507, 1.55, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") & 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 508, 2.05, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") & 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 509, 2.55, depth),
                depth = ifelse(site_id == "BARC" & 
                                 (lubridate::as_date(startDateTime) >= lubridate::as_date("2021-06-08") | 
                                    lubridate::as_date(startDateTime) <= lubridate::as_date("2022-01-07") )& 
                                 verticalPosition == 510, 3.05, depth)) %>%
  dplyr::filter((tsdWaterTempFinalQF == 0 | (tsdWaterTempFinalQF == 1 & as_date(startDateTime) > as_date("2020-07-01")))) %>% 
  dplyr::mutate(time = ymd_h(format(startDateTime, "%y-%m-%d %H")),
                depth = round(depth, 1)) %>% # round to the nearest 0.1 m
  group_by(site_id, depth, time) %>%
  dplyr::summarize(temperature = median(tsdWaterTempMean, na.rm = TRUE),.groups = "drop") %>%
  dplyr::select(time, site_id, temperature, depth) |> 
  rename(observation = temperature) |> 
  mutate(variable = "temperature") |> 
  QC.temp(range = c(-5, 40), spike = 5, by.depth = T) %>%
  mutate(data_source = 'NEON_portal')

#### Temp QC protocol=================
temp_full <- hourly_temp_profile_portal %>%
  dplyr::filter(((depth > 0 & depth < 1))) |> 
  dplyr::filter(site_id == "BARC" & month(time) %in% c(9:10) & year(time) == 2019)
# additional QC steps implemented (FO, 2022-07-13)
##### check 1 Gross range tests on temperature
# temperature ranges 
T_max <- 32 # gross max
T_min <- -2 # gross min

# GR flag will be true if the temperature is outside the range specified 
temp_cleaned <-  temp_full %>%
  dplyr::mutate(observation =ifelse(observation >= T_min & observation <= T_max , 
                                    observation, NA))  %>%
  rename(wtemp = observation, datetime = time) %>%
  ungroup() %>%
  select(-variable, -depth, -data_source)

#plot cleaned data
ggplot(data = temp_cleaned, aes(x = datetime, y = wtemp))+
  geom_point()+
  facet_wrap(vars(site_id), scales = "free", nrow = 2)+
  theme_bw()

#final dataset
temp_final <- temp_cleaned %>%
  select(datetime, wtemp) %>%
  filter(complete.cases(.))

#write to file
write.csv(temp_final, "./data/BARC_wtemp_celsius_highFrequency.csv",row.names = FALSE)

######################DISSOLVED OXYGEN and CHL-A##############################

wq <- loadByProduct(dpID="DP1.20288.001", site=c("BARC"),
                     startdate="2018-01", enddate="2022-12", 
                     package="expanded", 
                     token = NEON_TOKEN,
                     check.size = F,
                     release = "current",
                     include.provisional = T)

# unlist the variables and add to the global environment
list2env(wq, .GlobalEnv)
colnames(waq_instantaneous)

wq_full <- waq_instantaneous %>%
  dplyr::select(siteID, startDateTime, sensorDepth,
         dissolvedOxygen) %>%
  dplyr::mutate(sensorDepth = as.numeric(sensorDepth),
                DO = as.numeric(dissolvedOxygen)) %>%
  dplyr::rename(site_id = siteID) %>% 
  dplyr::filter(((sensorDepth > 0 & sensorDepth < 1)| is.na(sensorDepth))) |> 
  dplyr::mutate(startDateTime = as_datetime(startDateTime)) %>%
  dplyr::select(startDateTime, site_id, DO) %>%
  dplyr::filter(site_id == "BARC" & month(startDateTime) %in% c(9:10) & year(startDateTime) == 2019)

# dissolved oxygen ranges
do_max <- 20
do_min <- 0

wq_cleaned <- wq_full %>%
  dplyr::mutate(DO = ifelse(is.na(DO),DO, ifelse(DO >= do_min & DO <= do_max, DO, NA))) 
# manual cleaning based on visual inspection

##############################################################################
gplot(data = wq_cleaned, aes(x = startDateTime, y = DO))+
  geom_point()+
  facet_wrap(vars(site_id), scales = "free", nrow = 2)+
  theme_bw()
##############################################################################

do_final <- wq_cleaned %>%
  rename(datetime = startDateTime) %>%
  select(datetime, DO) %>%
  filter(complete.cases(.))

write.csv(do_final, "./data/BARC_DO_milligramsPerLiter_highFrequency.csv",row.names = FALSE)

######################NITRATE##############################

nit <- loadByProduct(dpID="DP1.20033.001", site=c("BARC"),
                    startdate="2018-01", enddate="2022-12", 
                    package="expanded", 
                    token = NEON_TOKEN,
                    check.size = F,
                    release = "current",
                    include.provisional = T)

# unlist the variables and add to the global environment
list2env(nit, .GlobalEnv)
colnames(NSW_15_minute)

nit_full <- NSW_15_minute %>%
  dplyr::select(siteID, startDateTime, 
                surfWaterNitrateMean) %>%
  dplyr::mutate(nitrate = as.numeric(surfWaterNitrateMean)) %>%
  dplyr::rename(site_id = siteID) %>% 
  dplyr::mutate(startDateTime = as_datetime(startDateTime)) %>%
  dplyr::select(startDateTime, site_id, nitrate) %>%
  dplyr::filter(site_id == "BARC" & month(startDateTime) %in% c(9:10) & year(startDateTime) == 2019)

# nitrate ranges
nitrate_max <- 200
nitrate_min <- 0

nit_cleaned <- nit_full %>%
  dplyr::mutate(nitrate = ifelse(is.na(nitrate),nitrate, ifelse(nitrate >= nitrate_min & nitrate <= nitrate_max, nitrate, NA))) 
# manual cleaning based on visual inspection

##############################################################################
ggplot(data = nit_cleaned, aes(x = startDateTime, y = nitrate))+
  geom_point()+
  facet_wrap(vars(site_id), scales = "free", nrow = 2)+
  theme_bw()
##############################################################################

nit_final <- nit_cleaned %>%
  rename(datetime = startDateTime) %>%
  select(datetime, nitrate) %>%
  filter(complete.cases(.))

write.csv(nit_final, "./data/BARC_surfN_micromolesPerLiter_highFrequency.csv",row.names = FALSE)

##### Re-wrangle high-frequency data to be "lab experiment" ##################
chla <- read_csv("./module_admin/data/BARC_chla_microgramsPerLiter_highFrequency.csv") %>%
  mutate(time = as_hms(ymd_hms(datetime)),
         date = date(datetime)) %>%
  filter(date >= "2019-10-09" & date <= "2019-10-12") %>%
  group_by(date) %>%
  mutate(day = cur_group_id()) %>%
  ungroup() %>%
  select(-date, -datetime)
write.csv(chla,"./assignment/data/chla_microgramsPerLiter_highFrequency.csv",row.names = FALSE)

do <- read_csv("./module_admin/data/BARC_DO_milligramsPerLiter_highFrequency.csv") %>%
  mutate(time = as_hms(ymd_hms(datetime)),
         date = date(datetime)) %>%
  filter(date >= "2019-10-09" & date <= "2019-10-12") %>%
  group_by(date) %>%
  mutate(day = cur_group_id()) %>%
  ungroup() %>%
  select(-date, -datetime)
write.csv(do,"./assignment/data/DO_milligramsPerLiter_highFrequency.csv",row.names = FALSE)

surfn <- read_csv("./module_admin/data/BARC_surfN_micromolesPerLiter_highFrequency.csv") %>%
  mutate(time = as_hms(ymd_hms(datetime)),
         date = date(datetime)) %>%
  filter(date >= "2019-10-09" & date <= "2019-10-12") %>%
  group_by(date) %>%
  mutate(day = cur_group_id()) %>%
  ungroup() %>%
  select(-date, -datetime)
write.csv(surfn,"./assignment/data/surfN_micromolesPerLiter_highFrequency.csv",row.names = FALSE)

wtemp <- read_csv("./module_admin/data/BARC_wtemp_celsius_highFrequency.csv") %>%
  mutate(time = as_hms(ymd_hms(datetime)),
         date = date(datetime)) %>%
  filter(date >= "2019-10-09" & date <= "2019-10-12") %>%
  group_by(date) %>%
  mutate(day = cur_group_id()) %>%
  ungroup() %>%
  select(-date, -datetime)
write.csv(wtemp,"./assignment/data/wtemp_celsius_highFrequency.csv",row.names = FALSE)
