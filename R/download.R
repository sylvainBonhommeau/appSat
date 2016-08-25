################################################################
#### AIM: SCRIPT to download the SLA data from Copernicus    ###
#### create the images, and make the time series for Reunion ###
#### INPUT: Copernicus ftp server is given in the script     ###
#### OUTPUT: - NetCDF files                                  ###
####         - images for world, Indian Ocean, Reunion area, ###
####           and Reunion Island                            ###
####         - time series of SLA for                        ###
#### FUNCTIONS: smoothData, imageMaker                       ###
################################################################

##################
### LIBRARIES  ###
##################
library(pacman)
p_load(R.utils, ncdf4, fields, mapdata, raster, reshape, plyr, progress)

###########################
### SET the directory   ###
### to save the data    ###
###########################
setwd('/media/sbonhomm/Caparmor/syl/boulot/Article/BonhMouq2017/data')

#################
### FUNCTIONS ###
#################
### Function to smooth the data for zoom area (Indian and Reunion)
smoothData <- function(data){
  tmpData <- data$z
  colnames(tmpData) <- data$x
  rownames(tmpData) <- data$y
  tmpData           <- melt(tmpData)
  tmpData           <- smooth.2d(tmpData$value,
                                 x=data.frame(lon=tmpData$X2, lat=tmpData$X1),
                                 theta=0.25, nrow=1000, ncol=1000)
  tmpData$z[tmpData$z>0.3]   <- 0.3
  tmpData$z[tmpData$z< -0.3] <- -0.3
  tmpData$z <- t(tmpData$z)
  return(tmpData)
}

### Function to make an image of the different areas and extract the Indian Ocean, Reunion Island vicinity,
### and extract the mean value of the pixels around Reunion Island

imageMaker <- function(fileName, slaType="NRT", destFile="/home/sbonhomm/dropboxScientific/Dropbox/SSH/"){
  ### Unzip and read the netCDF file
  fileUnzip <- gunzip(fileName, overwrite=T)
  fileUnzip <- as.character(fileUnzip)
  nc        <- nc_open(fileUnzip, write = T)
  lon       <- ncvar_get(nc, "lon")
  lat       <- ncvar_get(nc, "lat")
  sla       <- ncvar_get(nc, "sla")
  nc_close(nc)
  gzip(paste(fileUnzip))
  datePlot <- 
  ### transform the longitudes because they are on a 0°-360° range
  lon2      <- c((lon[which(lon>180)]-360),lon[lon<=180])
  sla2      <- rbind(sla[which(lon>180),], sla[which(lon<=180),])
  sshdata   <- list(x=lon2, y=lat, z=sla2)
  ### Extract the mean value for the pixels West of Reunion Island
  iwhlon         <- which(sshdata$x>=  55    & sshdata$x<=  55.5)
  iwhlat         <- which(sshdata$y>= -21.25 & sshdata$y<= -20.75)
  dateStart      <- unlist(strsplit(fileUnzip, split="_"))[6]
  timeSLAReunion <- data.frame(date=dateStart, SLA=mean(sshdata$z[iwhlon,iwhlat], na.rm=T))
  
  ### Limit the SLA values for images
  sshdata$z[sshdata$z>=  0.3] <-  0.3
  sshdata$z[sshdata$z<= -0.3] <- -0.3
  ### Extract and smooth the Indian Ocean area
  iwhlon         <- which(sshdata$x>=  20 & sshdata$x<= 90)
  iwhlat         <- which(sshdata$y>= -60 & sshdata$y<= 10)
  sshdataIndian  <- list(x=sshdata$x[iwhlon], y=sshdata$y[iwhlat], z=sshdata$z[iwhlon,iwhlat])
  sshdataIndian2 <- smoothData(sshdataIndian)
  ### Extract and smooth the Reunion Island area
  iwhlon         <- which(sshdata$x>=  45 & sshdata$x<=  65)
  iwhlat         <- which(sshdata$y>= -30 & sshdata$y<= -10)
  sshdatareunion <- list(x=sshdata$x[iwhlon], y=sshdata$y[iwhlat], z=sshdata$z[iwhlon,iwhlat])
  sshdatareunion <- smoothData(sshdatareunion)
  ### Extract the Reunion Island pixels
  iwhlon         <- which(sshdata$x>=  55   & sshdata$x<=  56)
  iwhlat         <- which(sshdata$y>= -21.5 & sshdata$y<= -20.75)
  sshReunion     <- list(x=sshdata$x[iwhlon], y=sshdata$y[iwhlat], z=sshdata$z[iwhlon,iwhlat])
  
  colorTable  <- designer.colors(100, c( "blue","white", "red") )
  fileName    <- paste(destFile, slaType, "/world_", dateStart, ".png", sep="")
  png(fileName, width=1400, height=1000)
  image.plot(sshdata, graphics.reset=FALSE, zlim=c(-0.3,0.3), col=colorTable)
  map("world", add=T, fill=T)
  dev.off()
  fileName <- paste(destFile, slaType, "/indian_", dateStart, ".png", sep="")
  png(fileName, width=1400, height=1000)
  image.plot(sshdataIndian2, graphics.reset=FALSE, zlim=c(-0.3,0.3), col=colorTable)
  map("worldHires", add=T, fill=T)
  dev.off()
  fileName <- paste(destFile, slaType, "/reunion_", dateStart, ".png", sep="")
  png(fileName, width=1400, height=1000)
  image.plot(sshdatareunion, graphics.reset=FALSE, zlim=c(-0.3,0.3), col=colorTable)
  map("worldHires", add=T, fill=T)
  dev.off()
  fileName <- paste(destFile, slaType, "/reunionIsland_", dateStart, ".png", sep="")
  png(fileName, width=1400, height=1000)
  image.plot(sshReunion, graphics.reset=FALSE, zlim=c(-0.3,0.3), col=colorTable)
  map("worldHires", add=T, fill=T)
  dev.off()
  output <- list(date=dateStart, timeSLAReunion=timeSLAReunion$SLA, sshdataIndian=sshdataIndian$z)
  return(output)
}

############################################
### FTP server for NRT and REP SLA data  ###
############################################
ftpRep <- "ftp://ftp.sltac.cls.fr/Core/SEALEVEL_GLO_SLA_MAP_L4_REP_OBSERVATIONS_008_027/dataset-duacs-rep-global-merged-allsat-msla-l4/"
ftpNrt <- "ftp://ftp.sltac.cls.fr:21/Core/SEALEVEL_GLO_SLA_MAP_L4_NRT_OBSERVATIONS_008_026/dataset-duacs-nrt-global-merged-allsat-msla-l4/"

#####################
### GET the data  ###
#####################
load("~/.loginpasswdSyl.Rdata")
### REP SLA data
dir.create("REP")
setwd('REP')
system(paste("wget -r --ftp-user=", loginpasswdSyl[1], " --ftp-password=", loginpasswdSyl[2], " ", ftpRep,"*", sep=""))

### NRT SLA data
setwd('..')
dir.create("NRT")
setwd('NRT')
system(paste("wget -r --ftp-user=", loginpasswdSyl[1], " --ftp-password=", loginpasswdSyl[2], " ", ftpNRT,"*", sep=""))

#########################
### CLEAN UP folders  ###
#########################
### Clean up the NRT folder
fileList    <- list.files(".", recursive = T)
ncFilesOnly <- grep(pattern = "nc.gz", x = fileList)
for (i in 1:length(ncFilesOnly)){
  file.copy(fileList[ncFilesOnly[i]], ".", copy.date = TRUE)
  file.remove(fileList[ncFilesOnly[i]])
}
### Clean up the REP folder
setwd('../REP')
fileList    <- list.files(".", recursive = T)
ncFilesOnly <- grep(pattern = "nc.gz", x = fileList)
for (i in 1:length(ncFilesOnly)){
  file.copy(fileList[ncFilesOnly[i]], ".", copy.date = TRUE)
  file.remove(fileList[ncFilesOnly[i]])
}

#########################################################################################################
### Make the images, save the time series for Reunion Island pixels and save the Indian Ocean data    ###
#########################################################################################################
### For the REP data
destFile1 <- "/home/sbonhomm/dropboxScientific/Dropbox/SSH/"
slaType1  <- "REP"
fileNames <- list.files(path = paste("/media/sbonhomm/Caparmor/syl/boulot/Article/BonhMouq2017/data/", slaType1, sep=""), pattern="nc.gz")
# number of nc files
nbNc             <- length(fileNames)
sshRepData       <- array(NA, dim=c(280, 280, nbNc))
timeSeriesSLARep <- data.frame(date=NA, SLA=NA)
pb <- progress_bar$new(total = nbNc)
for (i in 1:nbNc){
  pb$tick()
  tmp                 <- imageMaker(fileName = fileNames[i], slaType=slaType1, destFile=destFile1)
  sshRepData[,,i]     <- tmp$sshdataIndian
  timeSeriesSLARep    <- rbind(timeSeriesSLARep, data.frame(date=tmp$date, SLA=tmp$timeSLAReunion))
  
}
timeSeriesSLARep <- timeSeriesSLARep[-1,]
save(sshRepData,       file="/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/sshRepData.Rdata")
save(timeSeriesSLARep, file="/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/timeSeriesSLARep.Rdata")

### For the NRT data
setwd('../NRT')
destFile1 <- "/home/sbonhomm/dropboxScientific/Dropbox/SSH/"
slaType1  <- "NRT"
fileNames <- list.files(path = paste("/media/sbonhomm/Caparmor/syl/boulot/Article/BonhMouq2017/data/", slaType1, sep=""), pattern="nc.gz")
# number of nc files
nbNc             <- length(fileNames)
sshNrtData       <- array(NA, dim=c(280, 280, nbNc))
timeSeriesSLANrt <- data.frame(date=NA, SLA=NA)
pb <- progress_bar$new(total = nbNc)
for (i in 1:nbNc){
  print(i)
  pb$tick()
  tmp                 <- imageMaker(fileName = fileNames[i], slaType=slaType1, destFile=destFile1)
  sshNrtData[,,i]     <- tmp$sshdataIndian
  timeSeriesSLANrt    <- rbind(timeSeriesSLANrt, data.frame(date=tmp$date, SLA=tmp$timeSLAReunion))
}
timeSeriesSLANrt <- timeSeriesSLANrt[-1,]
save(sshNrtData,       file="/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/sshNrtData.Rdata")
save(timeSeriesSLANrt, file="/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/timeSeriesSLANrt.Rdata")
