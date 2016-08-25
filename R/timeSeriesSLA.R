
setwd('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data')

sonelData           <- read.csv("sonelData.csv", sep=";", header=F)
sonelData$V1        <- gsub("-", "", sonelData$V1)
colnames(sonelData) <- c("date", "Sonel")
load('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/timeSeriesSLARep.Rdata')
load('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/timeSeriesSLANrt.Rdata')

dateTot        <- data.frame(date=unique(c(timeSeriesSLARep$date, timeSeriesSLANrt$date)))
globaDataFrame <- merge(dateTot, sonelData, by.x="date", by.y="date", all.x=T, sort=T)
globaDataFrame <- merge(globaDataFrame, timeSeriesSLARep, by.x="date", by.y="date", all.x=T, sort=T)
globaDataFrame <- merge(globaDataFrame, timeSeriesSLANrt, by.x="date", by.y="date", all.x=T, sort=T)
colnames(globaDataFrame) <- c("date", "Sonel", "SLARep", "SLANrt")

save(globaDataFrame, file="/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/globaDataFrame.Rdata")