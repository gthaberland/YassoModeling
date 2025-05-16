
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


run_yasso_and_format


#______________________________________________________________________________
############################################################

# My data and code ####
library(Matrix)  # tai Matrix tms  
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

Yasso07Parameters_load = read.csv("y07par_gui.csv")

# Define number of simulation years
simulation_years <- 29  # or whatever you set

# Then Year vector must be length 22 (0 to 21)
Year_vector <- 0:(simulation_years +1)
length(Year_vector)
      
final_results <- list()

for (i in 1:nrow(input_data)) {
      
      row_data <- input_data[i, ]
      
      # Filter species-specific Cpools
      species_row <- unique_species %>% filter(Species == row_data$Species)
      
      for (ws_size in c(0, 4, 10)) {
            
            col_name <- case_when(
                  ws_size == 0 ~ "BR2_tha_diff",
                  ws_size == 4 ~ "BR2_7_tha_diff",
                  ws_size == 10 ~ "BR7_tha_diff"
            )
            
            if (!col_name %in% colnames(row_data) || is.na(row_data[[col_name]])) next
            
            input_tha <- row_data[[col_name]]
            
            cpool <- species_row %>% filter(WS == ws_size)
            if (nrow(cpool) == 0) next
            
            litter_input <- matrix(0, nrow = simulation_years + 1, ncol = 5)
            litter_input[1, ] <- c(
                  input_tha * cpool$C1,
                  input_tha * cpool$C2,
                  input_tha * cpool$C3,
                  input_tha * cpool$C4,
                  input_tha * cpool$C5
            )
            
            # Climate variables for both scenarios
            MeanTemperature_base <- rep(row_data$MeanTemp_avg, simulation_years + 1)
            MeanTemperature_cc <- MeanTemperature_base + 2.5
            TemperatureAmplitude <- rep(row_data$TempAmp_avg / 2, simulation_years + 1)
            Precipitation_base <- rep(row_data$Precip_annual_avg, simulation_years + 1)
            Precipitation_cc <- Precipitation_base * 0.8
            
            # Function to run Yasso and return dataframe
            run_yasso_and_format <- function(temp, precip, scenario_label) {
                  out <- yasso07(
                        MeanTemperature = temp,
                        TemperatureAmplitude = TemperatureAmplitude,
                        Precipitation = precip,
                        InitialCPool = rep(0, 5),
                        LitterInput = litter_input,
                        WoodySize = ws_size,
                        Yasso07Parameters = Yasso07Parameters_load$value
                  )
                  
                  out_df <- as.data.frame(out)
                  colnames(out_df) <- c("C_A", "C_W", "C_E", "C_N", "C_H")
                  
                  # Add Year 0 manually
                  year0 <- data.frame(
                        C_A = litter_input[1, 1],
                        C_W = litter_input[1, 2],
                        C_E = litter_input[1, 3],
                        C_N = litter_input[1, 4],
                        C_H = litter_input[1, 5]
                  )
                  combined <- rbind(year0, out_df[-1, ])
                  combined$TotalC <- rowSums(combined)
                  
                  data.frame(
                        Id_Inventari = rep(row_data$Id_Inventari, length(Year_vector)),
                        Species = rep(row_data$Species, length(Year_vector)),
                        WoodySize = rep(ws_size, length(Year_vector)),
                        Year = Year_vector,
                        C_A = combined$C_A,
                        C_W = combined$C_W,
                        C_E = combined$C_E,
                        C_N = combined$C_N,
                        C_H = combined$C_H,
                        TotalC = combined$TotalC,
                        Scenario = scenario_label,
                        stringsAsFactors = FALSE
                  )
            }
            
            # Run both scenarios
            res_base  <- run_yasso_and_format(MeanTemperature_base, Precipitation_base, "Baseline")
            res_cc_t  <- run_yasso_and_format(MeanTemperature_cc, Precipitation_base, "CC_T2.5")
            res_cc_p  <- run_yasso_and_format(MeanTemperature_base, Precipitation_cc, "CC_P0.8")
            res_cc_tp <- run_yasso_and_format(MeanTemperature_cc, Precipitation_cc, "CC_T2.5_P0.8")
            
            # Append both to results
            final_results[[length(final_results) + 1]] <- res_base
            final_results[[length(final_results) + 1]] <- res_cc_t
            final_results[[length(final_results) + 1]] <- res_cc_p
            final_results[[length(final_results) + 1]] <- res_cc_tp
      }
}

# Combine all results into one dataframe
final_output <- bind_rows(final_results)

final_output <- final_output %>%
      mutate(across(c(C_A, C_W, C_E, C_N, C_H, TotalC),
                    ~ ifelse(. == 0, 1e-6, .)))





final_output <- final_output %>%
      group_by(Id_Inventari, WoodySize, Scenario, Species) %>%
      mutate(
            initial_sum = first(C_A + C_W + C_E + C_N + C_H),
            
            rel_C_A = if_else(initial_sum > 0.001, C_A / first(C_A), NA_real_),
            rel_C_W = if_else(initial_sum > 0.001, C_W / first(C_W), NA_real_),
            rel_C_E = if_else(initial_sum > 0.001, C_E / first(C_E), NA_real_),
            rel_C_N = if_else(initial_sum > 0.001, C_N / first(C_N), NA_real_),
            rel_C_H = if_else(initial_sum > 0.001, C_H / first(C_H), NA_real_),
            rel_TotalC = if_else(initial_sum > 0.001, TotalC / first(TotalC), NA_real_)
      ) %>%
      select(-initial_sum) %>%
      ungroup()












final_output <- final_output %>%
      group_by(Id_Inventari, WoodySize, Scenario, Species) %>%
      mutate(
            rel_C_A = C_A / first(C_A),
            rel_C_W = C_W / first(C_W),
            rel_C_E = C_E / first(C_E),
            rel_C_N = C_N / first(C_N),
            rel_C_H = C_H / first(C_H),
            rel_TotalC = TotalC / first(TotalC)
      ) %>%
      ungroup()

#_______________________________________________________________________________
#  PLOTS #####

##### Scenarios
final_output %>%
      group_by(Year, Scenario) %>%
      summarise(TotalC = sum(TotalC, na.rm = TRUE), .groups = "drop") %>%
      ggplot(aes(x = Year, y = TotalC, color = Scenario)) +
      geom_line(size = 1) +
      labs(
            title = "Simulated Carbon Stocks under Climate Scenarios",
            x = "Year",
            y = "Total Carbon (Mg/ha)"
      ) +
      theme_minimal()

##### By C pool_______________________________________

##### By SCENARIO
plot_pool_scenario <- function(scenario_label) {
  final_output %>%
    filter(Scenario == scenario_label) %>%
    group_by(Year) %>%
    summarise(
      C_A = sum(C_A, na.rm = TRUE),
      C_W = sum(C_W, na.rm = TRUE),
      C_E = sum(C_E, na.rm = TRUE),
      C_N = sum(C_N, na.rm = TRUE),
      C_H = sum(C_H, na.rm = TRUE)
    ) %>%
    pivot_longer(cols = starts_with("C_"), names_to = "Pool", values_to = "Carbon") %>%
    mutate(Pool = factor(Pool, levels = c("C_A", "C_W", "C_E", "C_N", "C_H"))) %>% 
    ggplot(aes(x = Year, y = Carbon, color = Pool)) +
    geom_line() +
    labs(
      title = paste("Carbon Pools -", scenario_label),
      y = "Carbon Mg/ha",
      x = "Years",
      color = "Pool"
    ) +
    theme_minimal()
}

p1 <- plot_pool_scenario("Baseline")
p2 <- plot_pool_scenario("CC_T2.5")
p3 <- plot_pool_scenario("CC_P0.8")
p4 <- plot_pool_scenario("CC_T2.5_P0.8")

(p1 | p2) / (p3 | p4) + 
      plot_layout(guides = "collect") & 
      theme(legend.position = "bottom")

##### By SPECIES
plot_pool_species <- function(species_name) {
      final_output %>%
            filter(Species == species_name) %>%
            group_by(Year) %>%
            summarise(
                  C_A = sum(C_A, na.rm = TRUE),
                  C_W = sum(C_W, na.rm = TRUE),
                  C_E = sum(C_E, na.rm = TRUE),
                  C_N = sum(C_N, na.rm = TRUE),
                  C_H = sum(C_H, na.rm = TRUE)
            ) %>%
            pivot_longer(cols = starts_with("C_"), names_to = "Pool", values_to = "Carbon") %>%
            mutate(Pool = factor(Pool, levels = c("C_A", "C_W", "C_E", "C_N", "C_H"))) %>% 
            ggplot(aes(x = Year, y = Carbon, color = Pool)) +
            geom_line() +
            labs(
                  title = paste("Carbon Pools -", species_name),
                  y = "Carbon Mg/ha",
                  x = "Years",
                  color = "Pool"
            ) +
            theme_minimal()
}

s1 <- plot_pool_species("Quercus ilex")
s2 <- plot_pool_species("Pinus halepensis")
s3 <- plot_pool_species("Quercus humilis")
s4 <- plot_pool_species("Pinus nigra")

(s1 | s2) / (s3 | s4) + 
      plot_layout(guides = "collect") & 
      theme(legend.position = "bottom")

# Species AND WoodSize
plot_pool_species <- function(species_name) {
      final_output %>%
            filter(Species == species_name) %>%
            group_by(WoodySize, Year) %>%
            summarise(mean_rel_TotalC = mean(rel_TotalC, na.rm = TRUE), .groups = "drop") %>%
            ggplot(aes(x = Year, y = mean_rel_TotalC, color = factor(WoodySize))) +
            geom_line() +
            labs(
                  title = paste("Debris Type -", species_name),
                  y = "Remaining weight",
                  x = "Years",
                  color = "Debris Class"
            ) +
            theme_minimal()
}

s1 <- plot_pool_species("Quercus ilex")
s2 <- plot_pool_species("Pinus halepensis")
s3 <- plot_pool_species("Quercus humilis")
s4 <- plot_pool_species("Pinus nigra")

(s1 | s2) / (s3 | s4) + 
      plot_layout(guides = "collect") & 
      theme(legend.position = "bottom")

# Boxplot
final_output %>%
      filter(Year == 30, 
             rel_TotalC < 0.99,
             Species %in% c("Quercus ilex", "Pinus halepensis", "Pinus nigra", "Quercus humilis", 
                            "Quercus suber", "Pinus sylvestris")) %>%
      ggplot(aes(x = Species, y = rel_TotalC * 100, fill = Species)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 1.2) +
      labs(
            title = "% Mass Remaining - 30 Years",
            y = "% Mass Remaining"
      ) +
      facet_wrap(~WoodySize) +
      theme_minimal() +
      theme(
            legend.position = "bottom",
            axis.text.x = element_blank()  # Hide species names on x-axis
      )


#______________________________________________________________________________
# Adding Rodal information to the final output

# Convert Rodal$Id_Inventari to character
Rodal <- Rodal %>%
      mutate(Id_Inventari = as.character(Id_Inventari))

# Now do the join
final_output <- final_output %>%
      left_join(
            Rodal %>% select(Id_Inventari, Alt, Class_Actuacio),
            by = "Id_Inventari"
      )
#______________________________________________________________________________

##### By ACTUACIO
final_output %>%
      group_by(Class_Actuacio, Year) %>%
      summarise(mean_rel_TotalC = mean(rel_TotalC, na.rm = TRUE), .groups = "drop") %>%
      ggplot(aes(x = Year, y = mean_rel_TotalC, color = Class_Actuacio)) +
      geom_line() +
      labs(
            title = "Mean relative Total Carbon by Management",
            y = "Mean Relative Total Carbon",
            x = "Year"
      ) +
      theme_minimal()

##### By CLASS OF DEBRIS
final_output %>%
      group_by(WoodySize, Year) %>%
      summarise(mean_rel_TotalC = mean(rel_TotalC, na.rm = TRUE), .groups = "drop") %>%
      ggplot(aes(x = Year, y = mean_rel_TotalC, color = factor(WoodySize))) +
      geom_line() +
      labs(
            title = "Mean Relative Total Carbon by Debris Size",
            y = "Mean Relative Total Carbon",
            x = "Year",
            color = "Woody Size"
      ) +
      theme_minimal()

# Boxplot
final_output %>%
      filter(Year == 30, rel_TotalC < 0.99, !is.na(Alt)) %>%
      mutate(Alt_class = case_when(
            Alt < 300 ~ "<300 m",
            Alt >= 300 & Alt <= 600 ~ "300–600 m",
            Alt > 600 ~ ">600 m"
      ),
      Alt_class = factor(Alt_class, levels = c("<300 m", "300–600 m", ">600 m"))
      ) %>%
      ggplot(aes(x = Alt_class, y = rel_TotalC * 100, fill = Alt_class)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 1.2) +
      labs(
            title = "% Mass Remaining After 30 Years by Debris Class",
            x = "Altitude Class",
            y = "% Mass Remaining"
      ) +
      facet_wrap(~WoodySize) +
      theme_minimal() +
      theme(legend.position = "none")
