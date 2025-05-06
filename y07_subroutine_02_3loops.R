


library(dplyr)
library(tidyr)

# Define model configs
woody_configs <- list(
      list(size = 1,   pool_cols = c("BR2_C1",   "BR2_C2",   "BR2_C3",   "BR2_C4",   "BR2_C5")),
      list(size = 4.5, pool_cols = c("BR2_7_C1", "BR2_7_C2", "BR2_7_C3", "BR2_7_C4", "BR2_7_C5")),
      list(size = 10,  pool_cols = c("BR7_C1",   "BR7_C2",   "BR7_C3",   "BR7_C4",   "BR7_C5"))
)

woody_configs <- list(
      list(size = 1,   pool_cols = c(0.4, 0.25, 0.15, 0.15, 0.5)),
      list(size = 4.5, pool_cols = c(0.45, 0.25, 0.15, 0.10, 0.5)),
      list(size = 10,  pool_cols = c(0.55, 0.25, 0.1, 0.1, 0.5))
)

# Create an empty list to collect results
results_list <- list()

# Loop through each config
for (cfg in woody_configs) {
      for (j in 1:nrow(input_data)) {
            for (i in 1:length(sim_time)) {
                  year <- sim_time[i]
                  
                  # Get carbon pool input
                  c_pool <- input_data[j, ][[cfg$pool_col]] * rep(0.2, 5)
                  
                  # Run the model
                  out <- yasso07.light(
                        MeanTemperature = input_data[j, ]$MeanTemp_avg,
                        TemperatureAmplitude = input_data[j, ]$TempAmp_avg,
                        Precipitation = input_data[j, ]$Precip_annual_avg,
                        InitialCPool = c_pool,
                        LitterInput = c(0, 0, 0, 0, 0),
                        WoodySize = cfg$size,
                        Yasso07Parameters = Yasso07Parameters_load$value,
                        SimulationTime = year
                  )
                  
                  # Store results in long format
                  results_list[[length(results_list) + 1]] <- data.frame(
                        Id_Inventari = input_data[j, ]$Id_Inventari,
                        Species = input_data[j, ]$Species,
                        WoodySize = cfg$size,
                        Year = year,
                        C_A = out[1],
                        C_W = out[2],
                        C_E = out[3],
                        C_N = out[4],
                        C_H = out[5],
                        TotalC = sum(out)
                  )
            }
      }
}

# Combine all into a single data frame
model_outputs <- bind_rows(results_list)

model_outputs[model_outputs$Id_Inventari == 164 & model_outputs$Species == "Quercus ilex", ]


#_____________________________________________________________________________

library(dplyr)
library(tidyr)

# Define model configs
woody_configs <- list(
      list(size = 1,   pool_col = "BR2_tha_diff"),
      list(size = 4.5, pool_col = "BR2_7_tha_diff"),
      list(size = 10,  pool_col = "BR7_tha_diff")
)

# Create an empty list to collect results
results_list <- list()

# Loop through each config
for (cfg in woody_configs) {
      for (j in 1:nrow(input_data)) {
            for (i in 1:length(sim_time)) {
                  year <- sim_time[i]
                  
                  # Get carbon pool input
                  c_pool <- input_data[j, ][[cfg$pool_col]] * rep(0.2, 5)
                  
                  # Run the model
                  out <- yasso07.light(
                        MeanTemperature = input_data[j, ]$MeanTemp_avg,
                        TemperatureAmplitude = input_data[j, ]$TempAmp_avg,
                        Precipitation = input_data[j, ]$Precip_annual_avg,
                        InitialCPool = c_pool,
                        LitterInput = c(0, 0, 0, 0, 0),
                        WoodySize = cfg$size,
                        Yasso07Parameters = Yasso07Parameters_load$value,
                        SimulationTime = year
                  )
                  
                  # Store results in long format
                  results_list[[length(results_list) + 1]] <- data.frame(
                        Id_Inventari = input_data[j, ]$Id_Inventari,
                        Species = input_data[j, ]$Species,
                        WoodySize = cfg$size,
                        Year = year,
                        C_A = out[1],
                        C_W = out[2],
                        C_E = out[3],
                        C_N = out[4],
                        C_H = out[5],
                        TotalC = sum(out)
                  )
            }
      }
}

# Combine all into a single data frame
model_outputs <- bind_rows(results_list)

