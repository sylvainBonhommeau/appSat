### LOAD everything
## don't forget to do that if you don't have the last version
## install.packages('raster', repos = 'http://r-forge.r-project.org/', type = 'source')
library(raster)
library(R.utils)
library(ncdf4)
library(fields)
library(mapdata)
library(shiny)

date1 <- paste(format(Sys.Date()-1, format="%Y%m%d"), format(Sys.Date()-1, format="%Y%m%d"), sep="_")
extractSSH <- function(date=date1){
  # download.file(dest="titi.txt", 
  #               "ftp://ftp.sltac.cls.fr/Core/SEALEVEL_GLO_SLA_MAP_L4_NRT_OBSERVATIONS_008_026/dataset-duacs-nrt-global-merged-allsat-msla-l4/",
  #               method="wget",
  #               extra="--user=sbonhommeau1 --password=Mundaka35!")
  # tata <- readLines("titi.txt")
  # name <- grep(patt=date, x=tata)
  # tata               <- tata[name]
  # startAddress       <- regexpr('<a href=\"', tata)[1]+attr(regexpr('<a href=\"', tata), "match.length")
  # endAddress         <- regexpr('\">', tata)[1]-1
  # address            <- substr(tata,startAddress, endAddress)
  # fileName           <- unlist(strsplit(address, split="/"))[7]
  date <- paste(format(as.POSIXct(date, format="%Y%m%d")-6*3600*24, format="%Y%m%d"),"_",format(as.POSIXct(date, format="%Y%m%d"), format="%Y%m%d"), sep="")

  fileName <- paste("nrt_global_allsat_msla_h_", date, ".nc.gz", sep="")
  address  <- paste("ftp://ftp.sltac.cls.fr:21/Core/SEALEVEL_GLO_SLA_MAP_L4_NRT_OBSERVATIONS_008_026/dataset-duacs-nrt-global-merged-allsat-msla-l4/", fileName, sep="")
  download.file(dest=fileName, 
                address,
                method="wget",
                extra="--user=sbonhommeau1 --password=Mundaka35!")
  fileUnzip <- paste(unlist(strsplit(fileName, split="[.]"))[1], ".nc", sep="")
  nc1      <- gunzip(fileName, overwrite=T)
  nc       <- nc_open(fileUnzip, write = T)
  lat      <- ncvar_get(nc, "lat")
  lon      <- ncvar_get(nc, "lon")
  lonBnds  <- ncvar_get(nc, "lon_bnds")
  lon2 <- c((lon[which(lon>180)]-360),lon[lon<=180])
  lonBnds2 <- lonBnds
  lonBnds2[1,] <- c((lonBnds[1,][which(lon>180)]-360),lonBnds[1,][lon<=180])
  lonBnds2[2,] <- c((lonBnds[2,][which(lon>180)]-360),lonBnds[2,][lon<=180])
  sla  <- ncvar_get(nc, "sla")
  sla2 <- rbind(sla[which(lon>180),], sla[which(lon<=180),])
  nc_close(nc)
  file.remove(c(dir(pattern="gz"), dir(pattern="nc")))
  sshdata <- list(x=lon2, y=lat, z=sla2)
  sshdata$z[sshdata$z>=0.2] <- 0.2
  sshdata$z[sshdata$z<=-0.2] <- -0.2
  return(sshdata)
  
  
  

  
  # ### create the nc
  # nx <- length(lon)
  # ny <- length(lat)
  # dimX   <- dim.def.ncdf( "Longitude", "degrees", lon2 )
  # dimY   <- dim.def.ncdf( "Latitude", "degrees", lat )
  # mv=-2147483647
  # SLA <- var.def.ncdf( "SLA", "meters", list(dimX,dimY), mv )
  # nc <- create.ncdf( "essai.nc", list(SLA) )
  # put.var.ncdf( nc, SLA, sla2 )
  # close.ncdf(nc)
}  
  
# r <- raster("essai.nc", band=1)
# crs(r) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# pal    <- colorNumeric(c("#377EB8", "#4DAF4A", "#E41A1C"), seq(-0.7, 0.5, by=0.4),
#                     na.color = "transparent")


plotSSHworld <- function(sshdata=sshdata){
  image.plot(sshdata, graphics.reset=FALSE, main="Sea Level Anomaly (in m)", zlim=c(-0.2,0.2))
  map("world", add=T, fill=T)
}

plotSSHreunion <- function(sshdata=sshdata){
  iwhlon <- which(sshdata$x>=45 & sshdata$x<=70)
  iwhlat <- which(sshdata$y>= -30 & sshdata$y<= 0)
  sshdatareunion <- list(x=sshdata$x[iwhlon], y=sshdata$y[iwhlat], z=sshdata$z[iwhlon,iwhlat])
  image.plot(sshdatareunion, graphics.reset=FALSE, main="Sea Level Anomaly (in m)", xlim=c(45,70), ylim=c(-30,00),zlim=c(-0.2,0.2))
  map("worldHires", add=T, fill=T)
}

# plotSSH(date1)
# 
# 
# 
# image.plot(x=lon, y=lat, z=sla, zlim=c(-0.7,1), xlim=c(20,80), ylim=c(-45,5), graphics.reset=FALSE)
# map("worldHires", add=T, fill=T)
# 
# 
# world <- readShapePoly(fn="/home/sbonhomm/ownCloudIfremer/ThonRouge/ICCAT/Statistics/Donnees2015/data/TM_WORLD_BORDERS-0.3")
# lat <- get.var.ncdf(nc, "lat")
# lon <- get.var.ncdf(nc, "lon")-180
# ssha <- get.var.ncdf(nc, "ssha")
# sshaTrans<- data.frame(lon=lon, lat=lat, ssh=ssha)
# sshaTrans <- sshaTrans[-which(is.na(sshaTrans$ssh)==TRUE),]
# coordinates(sshaTrans) <- ~lon+lat
# r <- raster(ncol=length(lon), nrow=length(lat))
# 
# ras <- rasterize(sshaTrans, r, "ssh", fun=sum)
# 
# esshaTrans <- cast(sshTrans, lat~lon, value="value")
# 
# sshaTrans <- matrix(ssha, nrow=length(lat), ncol=length(lon))
# image.plot(sshaTrans, zlim=range(sshaTrans, na.rm=T))
# 
# wget ftp://sbonhommeau1:Mundaka35!@ftp.sltac.cls.fr/Core/SEALEVEL_GLO_SLA_MAP_L4_NRT_OBSERVATIONS_008_026/dataset-duacs-nrt-global-merged-allsat-msla-l4 