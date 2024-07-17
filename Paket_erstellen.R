

library(devtools)
library(roxygen2)


#create("C:\\RPAKET_multPower\\multPower")
create("C:\\multPower")

#dann R Datei in Unterordner R kopieren
document("C:\\multPower")

#Dann DESCRIPTION und NAMESPACE bearbeiten
check("C:\\multPower")

build("C:\\multPower")

