
Yasso07Parameters_load = read.csv("y07par_gui.csv")

# the temperature functions are four, one for each of the AWEN pools:

# T1 = MT + 4 * TA / pi * (1 / sqrt(2) - 1)          # Eq. 2.4 in model description
# T2 = MT - 4 * TA / (sqrt(2) * pi)                  # Eq. 2.5 in model description
# T3 = MT + 4 * TA / pi * (1 - 1 / sqrt(2))          # Eq. 2.6 in model description
# T4 = MT + 4 * TA / (sqrt(2) * pi)   
# mean(exp(beta1 * c(T1, T2, T3, T4) + beta2 * (c(T1, T2, T3, T4)^2)))
                
#Build exploratory grid of the temp functions

# Define the T1 function
T1 <- function(MT, TA) {
      MT + 4 * TA / pi * (1 / sqrt(2) - 1)
}
T2 <- function(MT, TA) {
      MT - 4 * TA / (sqrt(2) * pi) 
}
T3 <- function(MT, TA) {
      MT + 4 * TA / pi * (1 - 1 / sqrt(2))       
}
T4 <- function(MT, TA) {
      MT + 4 * TA / (sqrt(2) * pi)   
}

beta1 = Yasso07Parameters_load$value[17]
beta2 = Yasso07Parameters_load$value[18]

Temperature_function = function(MT, TA, beta1, beta2) {
      
      response = mat.or.vec(length(MT), length(TA))
      colnames(response) = TA
      rownames(response) = MT
      
      for(i in 1:length(MT)) {
            print(paste("Processing MT:", MT[i], "of", length(MT)))
            for(j in 1:length(TA)) {
                  print(paste("Processing TA:", TA[j], "of", length(TA)))
                  # Calculate T1, T2, T3, and T4
                  T1_val <- T1(MT[i], TA[j])
                  T2_val <- T2(MT[i], TA[j])
                  T3_val <- T3(MT[i], TA[j])
                  T4_val <- T4(MT[i], TA[j])
                  
                  # Calculate the mean of the exponentiated values
                  # using the provided beta1 and beta2
                  response[i, j] <- mean(exp(beta1 * c(T1_val, T2_val, T3_val, T4_val) + beta2 * (c(T1_val, T2_val, T3_val, T4_val)^2)))
                  print(paste("Response value:", response[i, j]))
                  }
      }
      return(response)
}        
               
               
# Create the grid of MT and TA values
MT_values <- seq(0, 20, length.out = 100)
TA_values <- seq(5, 15, length.out = 100)

T_response_matrix  = Temperature_function(MT_values, TA_values, beta1, beta2)



# Set up a nice color palette
library(viridis)


png("./study_functions/T_response_contour_plot.png", width = 800, height = 600, res = 100)
par(mar = c(5, 5, 4, 6))  # Adjust margins for better layout

# Create filled contour plot with custom levels
filled.contour(
      x = MT_values,
      y = TA_values,
      z = T_response_matrix,
      color.palette = viridis,
      levels = pretty(range(T_response_matrix), 20),
      plot.title = title(
            main = "T response surface",
            xlab = "Mean temp",
            ylab = "Temp amplitude",
            cex.main = 1.5,
            cex.lab = 1.2
      ),
      key.title = title(main = "response temperature", cex.main = 1),
      plot.axes = {
            axis(1, cex.axis = 1.1)
            axis(2, cex.axis = 1.1)
            # Add gridlines
            grid(lty = "dotted", col = "gray70")
      }
)
dev.off()


#### study moisture function

gamma = Yasso07Parameters_load$value[26]

moisture <- function(PR) {
      1 - exp(gamma * (PR/ 1000))
}


png("./study_functions/Moisture_function_plot.png", width = 800, height = 600, res = 100)
PR_seq = seq(0, 800, length.out = 100)
plot(PR_seq, moisture(PR = PR_seq), type = "l", col = "blue", lwd = 2,
     xlab = "Precipitation (mm)", ylab = "response moisture",
     main = "Moisture response function"
)
dev.off()




# plot the size scling, applied only to AWEN and not H pool

#(1 + delta1 * WS + delta2 * (WS^2))^(r), k[5]

delta1 = Yasso07Parameters_load$value[39]
delta2 = Yasso07Parameters_load$value[40]
r = Yasso07Parameters_load$value[41]

size <- function(WS) {
   (1 + delta1 * WS + delta2 * (WS^2))^(r)
   }

png("./study_functions/Size_function_plot.png", width = 800, height = 600, res = 100)
WS_seq = seq(0, 50, length.out = 100)
plot(WS_seq, size(WS = WS_seq), type = "l", col = "red", lwd = 2,
     xlab = "Diameter (cm)", ylab = "Response size",
     main = "Size response function"
)
dev.off()
