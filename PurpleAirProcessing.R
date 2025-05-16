#########################################################################
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 05/16/2025
###
#########################################################################

#' Process PurpleAir Sensor Data
#' 
#' @description 
#' This function processes raw PurpleAir sensor data by:
#' 1. Converting timestamps and adjusting to local time zone
#' 2. Calculating meteorological parameters (mixing ratio, etc.)
#' 3. Averaging particulate matter readings between dual sensors
#' 4. Applying EPA correction formulas to PM2.5 measurements
#' 5. Performing quality control checks on sensor readings
#' 
#' @param df A data frame containing raw PurpleAir sensor data
#' @param LocalTimeZone Character string specifying the local time zone (e.g., "America/Los_Angeles")
#' @param platform_name Character string identifying the sensor platform name
#' 
#' @return A processed data frame with additional calculated columns
#' 
#' @details
#' The function expects specific column names in the input data frame:
#' - UTCDateTime: Timestamp in UTC
#' - current_temp_f: Temperature in Fahrenheit
#' - current_dewpoint_f: Dew point in Fahrenheit
#' - pressure: Atmospheric pressure
#' - current_humidity: Relative humidity
#' - pm1_0_atm, pm2_5_atm, pm10_0_atm: Particulate matter readings from primary sensor
#' - pm1_0_atm_b, pm2_5_atm_b, pm10_0_atm_b: Particulate matter readings from secondary sensor
#' 
#' The PM2.5 correction formula is based on EPA research:
#' https://cfpub.epa.gov/si/si_public_record_report.cfm?dirEntryId=353088&Lab=CEMM
#' 
#' @examples
#' # Example usage:
#' # processed_data <- process_PA(raw_data, "America/New_York", "PurpleAir-A1234")
#' 
#' @import dplyr
#' @import lubridate
process_PA <- function(df, LocalTimeZone, platform_name) {
  df %>%
    # Convert timestamps and adjust to local timezone
    mutate(UTCDateTime = ymd_hms(UTCDateTime, tz = 'UTC'),
           local_time = with_tz(UTCDateTime, tzone = LocalTimeZone),
           
           # Calculate meteorological parameters
           # Convert temperature and dew point from Fahrenheit to Celsius
           temp_c = (5/9 * (current_temp_f - 32)),
           dp_c = (5/9 * (current_dewpoint_f - 32)),
           
           # Calculate vapor pressure and saturated vapor pressure using the Magnus-Tetens formula
           Vapor_pressure = (6.11*10^((7.5*dp_c)/(237.7+dp_c))),
           sat_vapor = (6.11*10^((7.5*temp_c)/(237.7+temp_c))),
           
           # Calculate mixing ratio and saturated mixing ratio
           # The constant 621.97 converts to g/kg
           mixing_ratio = (621.97*(Vapor_pressure/(pressure - Vapor_pressure))),
           saturated_mr = (621.97*(sat_vapor/(pressure - sat_vapor))),
           
           # Add platform identifier
           platform = platform_name) %>%
    
    # Calculate averages of particulate matter between the two PMS5003 sensors
    # Apply EPA correction formulas to PM2.5 measurements
    mutate(
      # Calculate averages for PM1.0, PM2.5, and PM10.0
      across(
        c(pm1_0_atm, pm2_5_atm, pm10_0_atm),
        ~ rowMeans(data.frame(., get(paste0(cur_column(), "_b"))), na.rm = TRUE),
        .names = "{.col}_avg"),
      
      # Updated PM2.5 correction formula based on EPA research:
      # https://cfpub.epa.gov/si/si_public_record_report.cfm?dirEntryId=353088&Lab=CEMM
      # The correction depends on the PM2.5 concentration range
      pm2_5_corr = as.numeric(case_when(
        # Range 1: 0 ≤ x < 30 μg/m³
        pm2_5_atm_avg >= 0 & pm2_5_atm_avg < 30 ~ 
          0.524 * pm2_5_atm_avg - 0.0862 * current_humidity + 5.75,
        
        # Range 2: 30 ≤ x < 50 μg/m³ (transitional blending of formulas)
        pm2_5_atm_avg >= 30 & pm2_5_atm_avg < 50 ~ 
          (0.786 * (pm2_5_atm_avg/20 - 3/2) + 0.524 * (1 - (pm2_5_atm_avg/20 - 3/2))) * pm2_5_atm_avg - 
          0.0862 * current_humidity + 5.75,
        
        # Range 3: 50 ≤ x < 210 μg/m³
        pm2_5_atm_avg >= 50 & pm2_5_atm_avg < 210 ~ 
          0.786 * pm2_5_atm_avg - 0.0862 * current_humidity + 5.75,
        
        # Range 4: 210 ≤ x < 260 μg/m³ (transitional blending of formulas)
        pm2_5_atm_avg >= 210 & pm2_5_atm_avg < 260 ~ {
          factor = (pm2_5_atm_avg/50 - 21/5)
          (0.69 * factor + 0.786 * (1 - factor)) * pm2_5_atm_avg - 
            0.0862 * current_humidity * (1 - factor) + 
            2.966 * factor + 5.75 * (1 - factor) + 
            8.84 * (10^(-4)) * pm2_5_atm_avg^2 * factor},
        
        # Range 5: 260 ≤ x μg/m³
        pm2_5_atm_avg >= 260 ~ 
          2.966 + 0.69 * pm2_5_atm_avg + 8.84 * 10^(-4) * pm2_5_atm_avg^2
      ))
    ) %>%
    
    # Perform quality control by comparing readings between the two identical sensors
    mutate(
      # Calculate absolute difference between sensors
      difference = abs(pm2_5_atm - pm2_5_atm_b), 
      
      # Calculate Relative Percent Difference (RPD)
      RPD = abs((pm2_5_atm - pm2_5_atm_b)/((pm2_5_atm + pm2_5_atm_b)/2)*100),
      
      # Flag data points where sensors disagree significantly
      # - RPD > 70% indicates poor agreement
      # - Absolute difference > 5 μg/m³ indicates significant deviation
      RPD_diff = if_else(RPD > 70, 1, 0) + if_else(difference > 5, 1, 0),
      
      # Mark data quality as "bad" if both conditions are met
      QC = case_when(RPD_diff == 2 ~ "bad", TRUE ~ "good"),
      
      # Set PM values to NA if quality control check failed
      pm2_5_corr = ifelse(QC == "bad", NA_character_, pm2_5_corr),
      pm2_5_atm_avg = ifelse(QC == "bad", NA_character_, pm2_5_atm_avg),
      pm10_0_atm_avg = ifelse(QC == "bad", NA_character_, pm10_0_atm_avg),
      pm1_0_atm_avg = ifelse(QC == "bad", NA_character_, pm1_0_atm_avg))
}