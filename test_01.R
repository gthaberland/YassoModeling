
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
names(input_clima)
input_clima <- input_clima %>%
      select(2, 8, 12, 23)

na_ids_year_93 <- input_clima %>%
      filter(is.na(MeanTemperature), Year > 1993) %>%
      select(Id_Inventari) %>%
      distinct()

Clima <- input_clima %>%
      select(Id_Inventari, MeanTemperature, Year) %>%
      group_by(Id_Inventari, Year) %>%
      summarise(
            MeanTemp = round(mean(MeanTemperature, na.rm = TRUE), 2),
            N = sum(!is.na(MeanTemperature)),  # or N = n() if no NAs
            .groups = "drop"
      )



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