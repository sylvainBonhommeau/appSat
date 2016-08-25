### LOAD everything
## don't forget to do that if you don't have the last version
## install.packages('raster', repos = 'http://r-forge.r-project.org/', type = 'source')
library(raster)
library(R.utils)
library(ncdf4)
library(fields)
library(mapdata)
library(shiny)
library(highcharter)

load("data/globaDataFrame.Rdata")
timeLag <- as.character(seq(0,105, 15))
for (i in 1:length(timeLag)){
  fileName1 <- paste("data/corNrtSLA",    timeLag[i], ".Rdata", sep="")
  fileName2 <- paste("data/pValueNrtSLA", timeLag[i], ".Rdata", sep="")
  fileName3 <- paste("data/corRepSLA",    timeLag[i], ".Rdata", sep="")
  fileName4 <- paste("data/pValueRepSLA", timeLag[i], ".Rdata", sep="")
  load("fileName1")
  load("fileName2")
  load("fileName3")
  load("fileName4")
}
load('data/latLon.Rdata')
## For the Dropbox
token     <- readRDS("data/droptoken.rds")

####################
#### FUNCTIONS   ###
####################
### Function to plot the time series of NRT and REP SLA and SOnel data
quant10SLA   <- quantile(globaDataFrame$SLANrt, prob=0.1)
quant33SLA   <- quantile(globaDataFrame$SLANrt, prob=0.3333)
currentSLA   <- globaDataFrame$SLANrt[length(globaDataFrame$SLANrt)]
minSLA       <- min(globaDataFrame$SLANrt, globaDataFrame$SLARep)
maxSLA       <- max(globaDataFrame$SLANrt, globaDataFrame$SLARep)

plotSonelSLA <- function(data){
  hc_opts <- list()
  hc_opts$title <- list(text="")
  hc <- highchart(hc_opts) %>%
    hc_chart(zoomType="x")%>%
    hc_yAxis(
      list(
        title = list(text = "SLA"),
        align = "left",
        showFirstLabel = FALSE,
        showLastLabel = FALSE,
        labels = list(format = "{value} m"),
        plotBands = list(
          list(from = minSLA, to = as.numeric(quant10SLA), color = "rgba(100, 0, 0, 0.1)",
               label = list(text = "Risk of low level", verticalAlign=c("bottom"))))
      ),
      list(
        title = list(text = "Sea level (Sonel data)"),
        align = "right",
        showFirstLabel = FALSE,
        showLastLabel = FALSE,
        labels = list(format = "{value} mm", useHTML = TRUE),
        opposite = TRUE
      )
    )  %>%
    hc_xAxis(categories = data$date) %>% 
    hc_add_series(data  = data$SLANrt, name="NRT SLA") %>% 
    hc_add_series(data  = data$SLARep, name="NRT REP") %>% 
    hc_add_series(data  = data$Sonel, name="Sonel data", yAxis = 1)
  return(hc)
}

### Function to plot the correlation map
lat       <- latLon$y[which(latLon$y>= -60 & latLon$y<= 10)]
lon       <- latLon$x[which(latLon$x>=  20 & latLon$x<= 90)]
colorTable<- designer.colors(100, c( "blue","white", "red") )
plotCor <- function(data){
  image.plot(x=lon, y=lat, data, graphics.reset = FALSE, col=colorTable, ylab="", xlab="")
  map("worldHires", add=T, fill=T)
}

