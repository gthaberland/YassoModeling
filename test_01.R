
source("y07_subroutine_01.r")
Yasso07Parameters_load = read.csv("y07par_gui.csv")

# we just follow one single pulse
InitialCPool_load = c(0,0,0,0,0)

# simulation runs for long time to understand how long it takes to go to zero
SimulationTime_load = seq(1:200)

# Input Debris DB
input_data = readRDS("../../../2024/03_Forest_debris_decomposition/01_Std_Restes_CPF/DB_processades_rds/DB_debris_species.rds")
input_data$Id_Inventari = as.factor(input_data$Id_Inventari)
input_data$Species = as.factor(input_data$Species)

plots_ID = levels(input_data$Id_Inventari)

# Input climatic data
library(readxl) #importar de excel
library(lubridate) #extract the year and month of a date
library(dplyr)

input_clima = read_excel("../../../2024/Meteo/Historical_meteo/Historical_meteo.xlsx")
str(input_clima)
input_clima$Year <- year(input_clima$dates)
input_clima$Month <- month(input_clima$dates)
input_clima$TemperatureAmplitude <- abs(input_clima$MinTemperature - input_clima$MaxTemperature)/2
names(input_clima)

input_clima <- input_clima %>%
      mutate(Id_Inventari = as.character(Id_Inventari))

input_data <- input_data %>%
      mutate(Id_Inventari = as.character(Id_Inventari))

# Now do the join
input_clima <- input_clima %>%
      left_join(input_data[, c("Id_Inventari", "Inv_Pre_Year")], by = "Id_Inventari")

library(dplyr)

input_clima <- input_clima %>% 
      distinct()

climate_summary <- input_clima %>%
      filter(Year >= Inv_Pre_Year) %>%
      group_by(Id_Inventari) %>%
      summarise(
            MeanTemp_avg = round(mean(MeanTemperature, na.rm = TRUE), 2),
            TempAmp_avg = round(mean(TemperatureAmplitude, na.rm = TRUE), 2),
            Precip_total = sum(Precipitation, na.rm = TRUE),
            Years_count = n_distinct(Year),
            Precip_annual_avg = round(Precip_total / Years_count, 2)
      ) %>%
      mutate(
            MeanTemp_avg = ifelse(is.na(MeanTemp_avg), 10, MeanTemp_avg),
            TempAmp_avg = ifelse(is.na(TempAmp_avg), 5, TempAmp_avg),
            Precip_annual_avg = ifelse(is.na(Precip_annual_avg) | Precip_annual_avg == 0, 650, Precip_annual_avg)
      ) %>%
      select(-Precip_total)


# Join input_data
input_data <- climate_summary %>%
      left_join(input_data, by = "Id_Inventari") %>%
      filter(
            !is.na(BR7_tha_diff) & !is.na(BR2_7_tha_diff) & !is.na(BR2_tha_diff),  # Remove rows with any NA
            (BR7_tha_diff + BR2_7_tha_diff + BR2_tha_diff) != 0                    # Remove rows where the sum is 0
      )

# save
saveRDS(input_data, "/../../../DB_processades_rds/DB_allometrics.rds")

rm(input_clima, input_data)


      # draft of a nested loop
# #outer loop by plot id
# for(i in 1:length(plots_ID)){
#       selected_plot = plots_ID[i]
#       
#       subset_data = input_data[input_data$Id_Inventari == selected_plot, ]
#       subset_data$Species = droplevels(subset_data$Species)
#      
#       subset_species = levels(subset_data$Species)
#      
#       #inner loop by species inside eah plot
#       for (j in 1:length(subset_species)){
#             subset_data[subset_data$Species == subset_species[j],]
#             
#             
#             
#             
#             
#             
#       } # end inner loop (for j, species)
#      
#      } #end outer loop (for i, plot)
# 
# 

for(i in 1:dim(input_data)[1]){
      input_data[i,]
   
      #yasso instance for >7 cm
      yasso07.light(MeanTemperature = 20, # time series, should probably be a vector of lenth N
                    TemperatureAmplitude = 5, # time series, should probably be a vector of lenth N
                    Precipitation = 0.600 ,# time series, should probably be a vector of lenth N
                    InitialCPool = InitialCPool_load,
                    LitterInput = input_data[i,]$BR7_tha_diff,
                    WoodySize = 20,
                    Yasso07Parameters = Yasso07Parameters_load,
                    SimulationTime = SimulationTime_load)   
      
      #yasso instance for 2 to 7cm
      yasso07.light(MeanTemperature = 20, # time series, should probably be a vector of lenth N
                    TemperatureAmplitude = 5, # time series, should probably be a vector of lenth N
                    Precipitation = 0.6,# time series, should probably be a vector of lenth N
                    InitialCPool = InitialCPool_load,
                    LitterInput = input_data[i,]$BR2_7_tha_diff,
                    WoodySize = 4.5,
                    Yasso07Parameters = Yasso07Parameters_load,
                    SimulationTime = SimulationTime_load)   
      
      #yasso instance for >2 cm
      yasso07.light(MeanTemperature, # time series, should probably be a vector of lenth N
                    TemperatureAmplitude, # time series, should probably be a vector of lenth N
                    Precipitation,# time series, should probably be a vector of lenth N
                    InitialCPool = InitialCPool_load,
                    LitterInput,
                    WoodySize=1,
                    Yasso07Parameters = Yasso07Parameters_load,
                    SimulationTime=SimulationTime_load)   
}


WoodySize = 

#N is the number of years of the simulation

yasso07.light(MeanTemperature, # time series, should probably be a vector of lenth N
              TemperatureAmplitude, # time series, should probably be a vector of lenth N
              Precipitation,# time series, should probably be a vector of lenth N
              InitialCPool = InitialCPool_load,
              LitterInput,
              WoodySize,
              Yasso07Parameters = Yasso07Parameters_load,
              SimulationTime)




      #    MeanTemperature,...                % Mean annual temperature, C
      #    TemperatureAmplitude,...           % (T_max-T_min)/2, C
      #    Precipitation,...                  % Annual rainfall, mm
      #    InitialCPool,...                   % AWENH, kg
      #    LitterInput,...                    % AWENH, kg/a
      
      MT=MeanTe