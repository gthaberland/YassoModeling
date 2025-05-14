# 1) Input forest debris and residues
# 2) Input climatic information (Mean, max and min temperature, Temp amplitude, and precipitation)

# Input Debris DB
input_data <- readRDS("../../../2024/03_Forest_debris_decomposition/01_Std_Restes_CPF/DB_processades_rds/DB_debris_species.rds")
input_data$Id_Inventari <- as.character(input_data$Id_Inventari)
input_data$Species <- as.factor(input_data$Species)

# Input climatic data
library(readxl)  # to read Excel
library(lubridate)  # to extract year and month from date
library(dplyr)

input_clima <- read_excel("../../../2024/Meteo/Historical_meteo/Historical_meteo.xlsx")

# Add year and month to the climate data
input_clima$Year <- year(input_clima$dates)
input_clima$Month <- month(input_clima$dates)

# Calculate monthly mean temperature, temperature amplitude, and other necessary values
input_clima <- input_clima %>%
      group_by(Id_Inventari, Year) %>%
      summarise(
            MeanTemp = mean(MeanTemperature, na.rm = TRUE),
            MinTemp = min(MinTemperature, na.rm = TRUE),
            MaxTemp = max(MaxTemperature, na.rm = TRUE),
            Precip = sum(Precipitation, na.rm = TRUE),
            .groups = "drop"
      ) %>%
      mutate(
            # Temperature Amplitude: Difference between max and min mean temperature for the year
            TemperatureAmplitude = MaxTemp - MinTemp
      )

# Now let's join with inventory data
input_clima <- input_clima %>%
      mutate(Id_Inventari = as.character(Id_Inventari))

input_data <- input_data %>%
      mutate(Id_Inventari = as.character(Id_Inventari))

input_clima <- input_clima %>%
      left_join(input_data[, c("Id_Inventari", "Inv_Pre_Year")], by = "Id_Inventari") %>%
      distinct()  # Remove duplicate rows

# Summarize the climate data by Id_Inventari
climate_summary <- input_clima %>%
      filter(Year >= Inv_Pre_Year) %>%
      group_by(Id_Inventari) %>%
      summarise(
            MeanTemp_avg = round(mean(MeanTemp, na.rm = TRUE), 2),
            TempAmp_avg = round(mean(TemperatureAmplitude, na.rm = TRUE), 2),
            Precip_total = sum(Precip, na.rm = TRUE),
            Years_count = n_distinct(Year),
            Precip_annual_avg = round(Precip_total / Years_count, 2),
            .groups = "drop"
      ) %>%
      mutate(
            MeanTemp_avg = ifelse(!is.finite(MeanTemp_avg), 10, MeanTemp_avg),
            TempAmp_avg = ifelse(!is.finite(TempAmp_avg), 5, TempAmp_avg),
            Precip_annual_avg = ifelse(!is.finite(Precip_annual_avg) | Precip_annual_avg == 0, 650, Precip_annual_avg)
      ) %>%
      select(-Precip_total)


# The final dataset `climate_summary` contains:
# - MeanTemp_avg
# - TempAmp_avg
# - Precip_annual_avg


# Join input_data
input_data <- climate_summary %>%
      left_join(input_data, by = "Id_Inventari") %>%
      filter(
            !is.na(BR7_tha_diff) & !is.na(BR2_7_tha_diff) & !is.na(BR2_tha_diff),  # Remove rows with any NA
            (BR7_tha_diff + BR2_7_tha_diff + BR2_tha_diff) != 0                    # Remove rows where the sum is 0
      ) %>%
      select(-(8:16))

names(input_data)

# Biochemical per species and wood size
unique_species <- read_excel("../../../2025/Exchange period - Finland/YassoModeling/unique_species.xlsx")
