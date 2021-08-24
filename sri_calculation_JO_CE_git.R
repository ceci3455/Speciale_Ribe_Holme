## Solar radiation Index calculation

## script by J. von Oppen and N. Chardon
## Start: April 2020 / January 2021

## Used by Cecilie Kristensen for thesis 2020 - 2021

# Installing packages
install.packages("GeoLight")
install.packages("insol")
install.packages("writexl")

rm(list=ls()) #clear environment
showTmpFiles() #show temporary files
removeTmpFiles(h = 0) #clear all temporary files


## DATA INPUT ##

#setwd("/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data") 
#PLOT_DATA <- read.csv("miljo_data.csv", header = TRUE, sep = ';') #plot metadata = PLOT_DATA - Note variable names: lat, long, SLOPE, ASPECT

# Softcoding working directory
Dir.Base <- getwd() # to find the project folder
Dir.Data <- file.path(Dir.Base, 'miljo_data.csv') # index the data folder


view(PLOT_DATA)
str(PLOT_DATA)

####################################################################################################

# # SOLAR RADIATION INDEX # # 

####################################################################################################

# Code based on Keating et al. 2007. Journal of Wildlife Management, 71:1344-1348

library(GeoLight) # zenith angles
library(insol) # declination
library(writexl) # exporting output back into excel format

# Calculate solar time and declination
doy <- 263 # Julian day in middle of fieldwork interval; 20st September chosen here
time <- as.POSIXct("2020-09-20 12:00:00 UTC") #POSIXct time; REPLACE WITH TIME OF INTEREST

sun <- GeoLight::solar(time)

# Calculate solar zenith angle
zen <- GeoLight::zenith(sun, PLOT_DATA$long, PLOT_DATA$lat)

# Calculate declination
declin <- insol::declination(doy)

# Radial angle of earth to sun at day of year (for calculation of ecventricity factor)
day_angle <- 2*pi*(doy-1)/365

# Calculate excentricity correction factor following Eq. 1.2.1 in Iqbal 1983
# https://doi.org/10.1016/B978-0-12-373750-2.50006-9
# https://www.sciencedirect.com/topics/physics-and-astronomy/elliptical-orbits
E0 <- 1.000110 + 0.034221 * cos(day_angle) + 0.00128 * sin(day_angle) + 0.000719 * 
  cos(2 * day_angle) + 0.000077 * sin(2 * day_angle) 

# calculate SRI (here creating sri column in PLOTS_DATA dataframe, adjust as needed)
sinDeg <- function(angle) sin(angle*pi/180)
cosDeg <- function(angle) cos(angle*pi/180)

for (i in 1:nrow(PLOT_DATA)) {
  PLOT_DATA$srii[i] <- E0 * 
    ((sinDeg(PLOT_DATA$lat[i]) * cosDeg(PLOT_DATA$haeldning[i]) - cosDeg(PLOT_DATA$lat[i]) *
        sinDeg(PLOT_DATA$haeldning[i]) * cosDeg(180 - PLOT_DATA$retning[i] - declin)) *
       sinDeg(zen[i]) + 
       (cosDeg(PLOT_DATA$lat[i]) * cosDeg(PLOT_DATA$haeldning[i]) + sinDeg(PLOT_DATA$lat[i]) *
          sinDeg(PLOT_DATA$haeldning[i]) * cosDeg(180 - PLOT_DATA$retning[i] - declin)) *
       (cosDeg(zen[i]) * cosDeg(0)) +
       cosDeg(zen[i]) * sinDeg(PLOT_DATA$haeldning[i]) *
       sinDeg(180 - PLOT_DATA$retning[i] - declin) * sinDeg(0))
}

# Check data distribution
summary(PLOT_DATA$srii)


# Writing output: write_xlsx(the dataframe name,"path to store the Excel file\\file name.xlsx")
setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/')

write.xlsx(PLOT_DATA, "sri_data.xlxs")


write.csv(PLOT_DATA, "sri.csv", row.names = FALSE)

# end of script ----

