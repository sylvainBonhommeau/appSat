
### LOAD everything
## don't forget to do that if you don't have the last version
## install.packages('raster', repos = 'http://r-forge.r-project.org/', type = 'source')
if (!require("pacman")) install.packages("pacman")
pacman::p_load(xlsx, tm, stringr, plyr, reshape, progress,
               ggplot2, xtable, knitr, pander, grid, pander,
               mapdata, sp, rgdal, maptools, rgeos, mailR, XML, ncdf,
               raster, fields, RColorBrewer, ncdf4, leaflet)

data.frame(shortName="SSH")

download.file(dest="fileList.txt", 
              "ftp://ftp.sltac.cls.fr/Core/SEALEVEL_GLO_SLA_MAP_L4_NRT_OBSERVATIONS_008_026/dataset-duacs-nrt-global-merged-allsat-msla-l4/",
              method="wget",
              extra="--user=sbonhommeau1 --password=Mundaka35!")
tata <- readLines("fileList.txt")
name <- grep(patt="20160817_20160817", x=tata)

tata               <- tata[name]
startAddress       <- regexpr('<a href=\"', tata)[1]+attr(regexpr('<a href=\"', tata), "match.length")
endAddress         <- regexpr('\">', tata)[1]-1

address            <- substr(tata,startAddress, endAddress)
fileName           <- unlist(strsplit(address, split="/"))[7]
download.file(dest=fileName, 
              address,
              method="wget",
              extra="--user=sbonhommeau1 --password=Mundaka35!")


library(raster)
nc <- open.ncdf('~/Downloads/nrt_global_allsat_msla_h_20160603_20160603.nc', write = T)
lat <- get.var.ncdf(nc, "lat")
lon <- get.var.ncdf(nc, "lon")
lonBnds <- get.var.ncdf(nc, "lon_bnds")
lon2 <- c((lon[which(lon>180)]-360),lon[lon<=180])
lonBnds2 <- lonBnds
lonBnds2[1,] <- c((lonBnds[1,][which(lon>180)]-360),lonBnds[1,][lon<=180])
lonBnds2[2,] <- c((lonBnds[2,][which(lon>180)]-360),lonBnds[2,][lon<=180])
sla  <- get.var.ncdf(nc, "sla")
sla2 <- rbind(sla[which(lon>180),], sla[which(lon<=180),])
put.var.ncdf(nc, varid = "lon", vals = lon2 )
put.var.ncdf(nc, varid = "sla", vals = sla2 )
put.var.ncdf(nc, varid = "lon_bnds", vals = lonBnds2 )
close.ncdf(nc)


### create the nc
nx <- length(lon)
ny <- length(lat)
dimX   <- dim.def.ncdf( "Longitude", "degrees", lon2 )
dimY   <- dim.def.ncdf( "Latitude", "degrees", lat )

mv=-2147483647
SLA <- var.def.ncdf( "SLA", "meters", list(dimX,dimY), mv )
nc <- create.ncdf( "essai.nc", list(SLA) )
put.var.ncdf( nc, SLA, sla2 )
close.ncdf(nc)
r <- raster("essai.nc", band=1)
crs(r) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
pal <- colorNumeric(c("#377EB8", "#4DAF4A", "#E41A1C"), seq(-0.7, 0.5, by=0.4),
                    na.color = "transparent")

leaflet() %>% addTiles() %>%
  addRasterImage(r, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = seq(-0.7, 0.5, by=0.4),
            title = "Sea Level Anomaly")


image.plot(x=lon2, y=lat, z=sla2, graphics.reset=FALSE)
map("worldHires", add=T, fill=T)


image.plot(x=lon, y=lat, z=sla, zlim=c(-0.7,1), xlim=c(20,80), ylim=c(-45,5), graphics.reset=FALSE)
map("worldHires", add=T, fill=T)


world <- readShapePoly(fn="/home/sbonhomm/ownCloudIfremer/ThonRouge/ICCAT/Statistics/Donnees2015/data/TM_WORLD_BORDERS-0.3")
lat <- get.var.ncdf(nc, "lat")
lon <- get.var.ncdf(nc, "lon")-180
ssha <- get.var.ncdf(nc, "ssha")
sshaTrans<- data.frame(lon=lon, lat=lat, ssh=ssha)
sshaTrans <- sshaTrans[-which(is.na(sshaTrans$ssh)==TRUE),]
coordinates(sshaTrans) <- ~lon+lat
r <- raster(ncol=length(lon), nrow=length(lat))

ras <- rasterize(sshaTrans, r, "ssh", fun=sum)

esshaTrans <- cast(sshTrans, lat~lon, value="value")

sshaTrans <- matrix(ssha, nrow=length(lat), ncol=length(lon))
image.plot(sshaTrans, zlim=range(sshaTrans, na.rm=T))

wget ftp://sbonhommeau1:Mundaka35!@ftp.sltac.cls.fr/Core/SEALEVEL_GLO_SLA_MAP_L4_NRT_OBSERVATIONS_008_026/dataset-duacs-nrt-global-merged-allsat-msla-l4 
