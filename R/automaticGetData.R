library(raster)
library(R.utils)
library(ncdf4)
library(fields)
library(mapdata)
library(reshape)
library(plyr)

setwd('/media/sbonhomm/Caparmor/syl/boulot/Article/BonhMouq2017/data/NRT/')
print("#######################################################################################")
print(paste("#####################", Sys.Date(), "############################################"))
print("#######################################################################################")
### Function to smooth the data for zoom area (Indian and Reunion)
smoothData <- function(data){
  tmpData           <- data$z
  colnames(tmpData) <- data$x
  rownames(tmpData) <- data$y
  tmpData           <- melt(tmpData)
  tmpData           <- smooth.2d(tmpData$value,
                                 x=data.frame(lon=tmpData$X2, lat=tmpData$X1),
                                 theta=0.25, nrow=1000, ncol=1000)
  tmpData$z[tmpData$z > 0.3] <- 0.3
  tmpData$z[tmpData$z< -0.3] <- -0.3
  tmpData$z                  <- t(tmpData$z)
  return(tmpData)
}

### Function to get the last SLA files. It checks whether the CRON was not run during some days before...
updateSLADatabase <- function(day){
  ### Download the last archive
  dateStart <- format(Sys.Date()-6-day, "%Y%m%d")
  dateEnd   <- format(Sys.Date()-day,   "%Y%m%d")
  datePlot  <- format(Sys.Date()-6-day, "%Y-%m-%d") 
  fileName  <- paste("nrt_global_allsat_msla_h_", dateStart, "_", dateEnd, ".nc.gz", sep="")
  address   <- paste("ftp://ftp.sltac.cls.fr:21/Core/SEALEVEL_GLO_SLA_MAP_L4_NRT_OBSERVATIONS_008_026/dataset-duacs-nrt-global-merged-allsat-msla-l4/", fileName, sep="")
  download.file(dest=fileName, 
                address,
                method="wget",
                extra=paste("--user=", loginpasswdSyl[1], " --password=", loginpasswdSyl[2], sep=""))
  nc1       <- gunzip(fileName, overwrite=T)
  nc1       <- as.character(nc1)
  nc        <- nc_open(nc1, write = T)
  lon       <- ncvar_get(nc, "lon")
  lat       <- ncvar_get(nc, "lat")
  sla       <- ncvar_get(nc, "sla")
  nc_close(nc)
  gzip(paste(nc1))
  lon2      <- c((lon[which(lon>180)]-360),lon[lon<=180])
  sla2      <- rbind(sla[which(lon>180),], sla[which(lon<=180),])
  sshdata   <- list(x=lon2, y=lat, z=sla2)
  
  iwhlon            <- which(sshdata$x>=55 & sshdata$x<=55.5)
  iwhlat            <- which(sshdata$y>= -21.25 & sshdata$y<= -20.75)
  timeSLAReunionTmp <- data.frame(date=dateStart, SLA=mean(sshdata$z[iwhlon,iwhlat], na.rm=T))
  
  sshdata$z[sshdata$z>= 0.3] <-  0.3
  sshdata$z[sshdata$z<=-0.3] <- -0.3
  
  iwhlon <- which(sshdata$x>=  20 & sshdata$x<= 90)
  iwhlat <- which(sshdata$y>= -60 & sshdata$y<= 10)
  sshdataIndian <- list(x=sshdata$x[iwhlon], 
                        y=sshdata$y[iwhlat], 
                        z=sshdata$z[iwhlon,iwhlat], 
                        date=datePlot)
  sshdataIndian2 <- smoothData(sshdataIndian)
  
  iwhlon <- which(sshdata$x>=  45 & sshdata$x<=  65)
  iwhlat <- which(sshdata$y>= -30 & sshdata$y<= -10)
  sshdatareunion <- list(x=sshdata$x[iwhlon], 
                         y=sshdata$y[iwhlat], 
                         z=sshdata$z[iwhlon,iwhlat])
  sshdatareunion <- smoothData(sshdatareunion)
  
  iwhlon <- which(sshdata$x>=  55   & sshdata$x<=  56)
  iwhlat <- which(sshdata$y>= -21.5 & sshdata$y<= -20.75)
  sshReunion <- list(x=sshdata$x[iwhlon], 
                     y=sshdata$y[iwhlat], 
                     z=sshdata$z[iwhlon,iwhlat])
  
  colorTable<- designer.colors(100, c( "blue","white", "red") )
  ### PLOT of the world SLA
  fileName <- paste("/home/sbonhomm/dropboxScientific/Dropbox/SSH/NRT/world_", dateStart, ".png", sep="")
  png(fileName, width=1400, height=1000)
    image.plot(sshdata, graphics.reset=FALSE, main=datePlot, zlim=c(-0.3,0.3), col=colorTable)
    map("world", add=T, fill=T)
  dev.off()
  ### PLOT of the Indian Ocean SLA
  fileName <- paste("/home/sbonhomm/dropboxScientific/Dropbox/SSH/NRT/indian_", dateStart, ".png", sep="")
  png(fileName, width=1400, height=1000)
    image.plot(sshdataIndian2, graphics.reset=FALSE, main=datePlot, zlim=c(-0.3,0.3), col=colorTable)
    map("worldHires", add=T, fill=T)
  dev.off()
  fileName <- paste("/home/sbonhomm/dropboxScientific/Dropbox/SSH/NRT/reunion_", dateStart, ".png", sep="")
  png(fileName, width=1400, height=1000)
    image.plot(sshdatareunion, graphics.reset=FALSE, main=datePlot, zlim=c(-0.3,0.3), col=colorTable)
    map("worldHires", add=T, fill=T)
  dev.off()
  fileName <- paste("/home/sbonhomm/dropboxScientific/Dropbox/SSH/NRT/reunionIsland_", dateStart, ".png", sep="")
  png(fileName, width=1400, height=1000)
    image.plot(sshReunion, graphics.reset=FALSE, zlim=c(-0.3,0.3), col=colorTable)
    map("worldHires", add=T, fill=T)
  dev.off()
  
  timeSeriesSLANrt <- rbind(timeSeriesSLANrt, timeSLAReunionTmp)
  save(timeSeriesSLANrt, file='/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/timeSeriesSLANrt.Rdata' )
  output <- sshdataIndian$z
  return(output)
}

### Open the SLA time series for Reunion Island pixels
load('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/timeSeriesSLANrt.Rdata')
load('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/sshNrtData.Rdata')
load("~/.loginpasswdSyl.Rdata")
lastDownload   <- as.POSIXct(timeSeriesSLANrt$date[length(timeSeriesSLANrt$date)], format="%Y%m%d")
missingDate    <- as.numeric(as.POSIXct(format(Sys.Date()-6, "%Y-%m-%d")) - lastDownload)


if (missingDate>0){
  sshNrtDataTmp <- array(NA, dim=c(280, 280, (dim(sshNrtData)[3]+missingDate)))
  for (i in 1:missingDate){
    sshNrtDataTmp[,,(dim(sshNrtData)[3]+i)] <- updateSLADatabase(i)
  }
  sshNrtData <- sshNrtDataTmp
  save(sshNrtData, file='/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/sshNrtData.Rdata')
}

### Now you can automatically launch this script using the command: 
### crontab -u sbonhomm -e
 

print("#################################################################")