

library(devtools)
library(roxygen2)


#create("C:\\RPAKET_multPower\\multPower")
create("C:\\multPower_17_07/multPower")

#dann R Datei in Unterordner R kopieren
document("C:\\multPower_17_07/multPower")

#Dann DESCRIPTION und NAMESPACE bearbeiten
check("C:\\multPower_17_07/multPower")

build("C:\\multPower_17_07/multPower")

