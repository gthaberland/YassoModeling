
# This file contains an R-function of Yasso07. The version is # based on matrix-version created by Jaakko Heikkinen with Matlab and # Yasso07 description by Tuomi & Liski 17.3.2008  (Yasso07.pdf) # Created by Taru Palosuo, Jaakko Heikkinen & Anu Akuj?rvi in December 2011

#  Instructions  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

# 1) first run the source code for the function with "source(R.Yasso_[DATE].r)"
# 2) then you can use the yasso07-function by just calling it yasso07(..) # 3) Input needed for the function:
#        1. MeanTemperature - vector of mean annual temperatures [C]
#        2. TemperatureAmplitude - vector of temperature amplitudes [C] 
#        3. Precipitation - vector of annual precipiations [mm]
#        4. InitialCPool - vector of initial C pools of model compartments, length 5, [whatever]
#        5. LitterInput - litter input matrix, 5 columns x simulation years as rows, [whatever]
#        6. WoodySize - size of woody litter (for non-woody litter this is 0)
#        7. Yasso07Parameters - these in the format applied in the fortran version, length 44

# NOTE that this function eats only one type of material at the time. So, non-woody and different woody litter # materials needs to be calculated separately.

# The output of the function is the matrix AWENH compartments as columns x rows as simulation years


# Basics  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

# additional R libraries (as needed)

# install.packages("Matrix", dependencies = TRUE)

library(Matrix)  # tai Matrix tms  

rm(list = ls())

# Function definition   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

# The input for this function is provided as vectors, row = year

yasso07 =  function(MeanTemperature, TemperatureAmplitude,Precipitation,InitialCPool,LitterInput,WoodySize,Yasso07Parameters) {
      
      #    MeanTemperature,...                % Mean annual temperature, C
      #    TemperatureAmplitude,...           % (T_max-T_min)/2, C
      #    Precipitation,...                  % Annual rainfall, mm
      #    InitialCPool,...                   % AWENH, kg
      #    LitterInput,...                    % AWENH, kg/a
      
      MT=MeanTemperature
      TA=TemperatureAmplitude
      PR=Precipitation
      PR=PR/1000;               # conversion from mm to meters
      
      LC=InitialCPool
      LI=LitterInput
      
      PA=Yasso07Parameters
      WS=WoodySize;              
      
      YR=length(MT);            # =length of simulation in years
      
      alfa=c(-PA[1], -PA[2], -PA[3], -PA[4], -PA[35])   # Vector of decomposition rates
      
      # Creating the matrix A_p (here called p)
      
      row1 = c(-1, PA[5], PA[6], PA[7], 0)
      row2 = c(PA[8], -1, PA[9], PA[10], 0)
      row3 = c(PA[11], PA[12], -1, PA[13], 0)
      row4 = c(PA[14], PA[15], PA[16], -1, 0)
      row5 = c(PA[36], PA[36], PA[36], PA[36], -1)
      
      p = matrix(c(row1, row2, row3, row4, row5), 5, 5, byrow=T)  
      
      # temperature dependence parameters
      beta1 = PA[17]
      beta2 = PA[18]
      gamma = PA[26]
      
      # Woody litter size dependence parameters
      delta1 = PA[39]
      delta2 = PA[40]
      r = PA[41]
      
      LC = matrix(InitialCPool, nrow=YR+1, ncol=5, byrow=TRUE)  # byrow added 30.10.2012! /TP
      
      for (h in 1:YR) {
            T1=MT[h]+4*TA[h]/pi*(1/sqrt(2)-1)          # Eq. 2.4 in model description
            T2=MT[h]-4*TA[h]/(sqrt(2)*pi)              # Eq. 2.5 in model description
            T3=MT[h]+4*TA[h]/pi*(1-1/sqrt(2))          # Eq. 2.6 in model description
            T4=MT[h]+4*TA[h]/(sqrt(2)*pi)              # Eq. 2.7 in model description 
            
            # k following Eq. 3 in Tuomi et al. 2009. Eco.Mod. 220: 3362-3371
            k=alfa*mean(exp(beta1*c(T1,T2,T3,T4)+beta2*(c(T1,T2,T3,T4)^2))*(1-exp(gamma*PR[h])))     
            
            # the effect of wl size as in Eq. 3.1 in model description
            k = c(k[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),k[5]) # If 0 it will just multiply '1 * k'
            
            A=p%*%diag(k)                             # Matrix multiplication in R: %*%
            
            #analytical solution as in Eq. 1.3 in model description
            LC[h+1,] = as.array(solve(A)%*% (expm(A)%*%(A%*%LC[h,]+LI[h,])-LI[h,]))
            
      }  # end of for h
      
      LC  # prints this as a result of a function
      
}  # end of yasso07 function

##############################
#### Here couple of examples
##############################


# Global parameters by Tuomi et al. (2011) 
# Boris code - 09.05.2025

Yasso07Parameters <- c(-0.7300,-5.8000,-0.2900,-0.031,0.4800,0.0100,0.8300,0.9900,0.0000,0.0100,0.0000,0.0000,0.0300,0.0000,0.0100,0.9200,0.0960,-0.0014,0.000000,0.000000,0.0000,0.0000,0.0000,0.0000,0.000000,-1.2100,0.000000,0.000000,0.000000,0.00000,0.00000,0.00000,0.000000,0.000000,-0.001700,0.004500,0.000000,0.000000,-1.7000,0.8600,-0.3060,0.0000,0.0000,0.0000)

# Years of simulation
j=100 # how many years

MeanTemperature <- rep(10,j)  # mean T 10 degrees C
TemperatureAmplitude <- rep(15,j) # temperature amplitude 
Precipitation <- rep(800,j) # precipitation 800mm
#InitialCPool<- c(0,0,0,0,0)
InitialCPool<- c(5.9099265,0.8196918,0.8101619,25.5382060,26.5636798)*0
LitterInput <- 4 * c(0.52,0.18,0.08,0.2,0) # 4 ton c /ha litter times AWEN fractions 
WoodySize <- c(10) # leaf litter so s = 0, for woody it would be else

for (i in seq(1,j, by=1)) { 
      LitterInput <- rbind(LitterInput, c(0,0,0,0,0)) } # looping input litter for i (time)

res1 <- yasso07(MeanTemperature, TemperatureAmplitude,Precipitation,InitialCPool,LitterInput,WoodySize,Yasso07Parameters)
View(LitterInput)
# degradation, litter goes down to 60%
LitterInput <- 0.6 * LitterInput
res2 <- yasso07(MeanTemperature, TemperatureAmplitude,Precipitation,InitialCPool,LitterInput,WoodySize,Yasso07Parameters)
# degrdatation + climate change 4 C
MeanTemperature <- MeanTemperature+4 # we add 2 degrees (Climate change)
res3 <- yasso07(MeanTemperature, TemperatureAmplitude,Precipitation,InitialCPool,LitterInput,WoodySize,Yasso07Parameters)
#
plot(rowSums(res1), type="l", col="darkgreen", 
     ylab="Carbon Mg/ha", xlab="Years", xlim=c(0,100), #ylim=c(0,60), 
     bty="n")
lines(rowSums(res2), col="green")
lines(rowSums(res3), col="blue") 
legend("topright", cex=0.6,c("business as usual","degradation","degradation + climate change"), pch=1, col=c("darkgreen","green","blue"))
abline(h=60, col="gray")
abline(h=35, col="gray")
abline(h=30, col="gray")

degchange <- diff(rowSums(res2))
degchangecc <- diff(rowSums(res3))

plot(degchange, type="l", col="darkgreen", 
     ylab="Change in Carbon Mg/ha", xlab="Years", ylim=c(-2,0), xlim=c(0,20), bty="n")
lines(degchangecc, col="green")





#________________________________________________________________________________
# My data and code ####

library(dplyr)
library(tidyr)


source("y07_subroutine_01.r")
Yasso07Parameters_load = read.csv("y07par_gui.csv")

# simulation runs for long time to understand how long it takes to go to zero
SimulationTime_load = seq(1:200)

# Years of simulation
j=100 # how many years

# Define model configs
woody_configs <- list(
      list(size = 0,   pool_col = "BR2_tha_diff"),
      list(size = 4, pool_col = "BR2_7_tha_diff"),
      list(size = 10,  pool_col = "BR7_tha_diff")
)

LitterInput <- input_data[j, ]

# Create an empty list to collect results
results_list <- list()

# Loop through each config
for (cfg in woody_configs) {
      for (i in seq(1,j, by=1)) { 
            LitterInput <- rbind(LitterInput, c(0,0,0,0,0)) }
                  
                  # Get carbon pool input
                  c_pool <- input_data[j, ][[cfg$pool_col]] * rep(0.2, 5)
                  
                  # Run the model
                  out <- yasso07(
                        MeanTemperature = input_data[j, ]$MeanTemp_avg,
                        TemperatureAmplitude = input_data[j, ]$TempAmp_avg,
                        Precipitation = input_data[j, ]$Precip_annual_avg,
                        InitialCPool = c(0, 0, 0, 0, 0),
                        LitterInput = input_data[j,]$BR2_tha_diff*c(0.2,0.2,0.2,0.2,0.2),
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
      

# Combine all into a single data frame
model_outputs <- bind_rows(results_list)

print(paste("Woody Size (WS): ", cfg$size))


#_____________#####


# Number of years to simulate
simulation_years <- 20

# Create an empty list to store results
results_list <- list()

# Loop through input_data
for (j in 1:nrow(input_data)) {
      
      row_data <- input_data[j, ]
      
      # Determine the wood size and corresponding biomass
      if (!is.na(row_data$BR2_tha_diff)) {
            ws_size <- 0
            biomass <- row_data$BR2_tha_diff
      } else if (!is.na(row_data$BR2_7_tha_diff)) {
            ws_size <- 4
            biomass <- row_data$BR2_7_tha_diff
      } else if (!is.na(row_data$BR7_tha_diff)) {
            ws_size <- 10
            biomass <- row_data$BR7_tha_diff
      } else {
            next  # Skip if no biomass data
      }
      
      # Get species-specific carbon shares for this wood size
      species_row <- unique_species[
            unique_species$Species == row_data$Species & unique_species$WS == ws_size, 
      ]
      
      if (nrow(species_row) == 0) next
      
      # Extract carbon share values and ensure they are numeric
      carbon_shares <- as.numeric(species_row[1, c("C1", "C2", "C3", "C4", "C5")])
      
      # Build LitterInput: only year 0 has input
      litter_input <- matrix(0, nrow = simulation_years +1, ncol = 5)
      litter_input[1, ] <- biomass * carbon_shares  # kg/ha per year for year 0
      
      # Replicate climate variables for all years (no litter input in year 1+)
      MeanTemperature <- rep(row_data$MeanTemp_avg, simulation_years)
      TemperatureAmplitude <- rep(row_data$TempAmp_avg, simulation_years)
      Precipitation <- rep(row_data$Precip_annual_avg, simulation_years)
      
      # Run the Yasso model
      out <- yasso07(
            MeanTemperature = MeanTemperature,
            TemperatureAmplitude = TemperatureAmplitude,
            Precipitation = Precipitation,
            InitialCPool = rep(0, 5),
            LitterInput = litter_input,
            WoodySize = ws_size,
            Yasso07Parameters = Yasso07Parameters_load$value
      )
      
      # Store results with corrected Year range
      results_list[[length(results_list) +1]] <- data.frame(
            Id_Inventari = row_data$Id_Inventari,
            Species = row_data$Species,
            WoodySize = ws_size,
            Year = 0:simulation_years,  # Adjusted to start from 0
            C_A = out[, 1],
            C_W = out[, 2],
            C_E = out[, 3],
            C_N = out[, 4],
            C_H = out[, 5],
            TotalC = rowSums(out)
      )
}

# Combine all into a single data frame
model_outputs <- bind_rows(results_list)
