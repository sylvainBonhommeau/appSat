####################################################################
#### AIM: SCRIPT to calculate the correlation maps between       ###
#### Reunion Island SLA and all other pixels in the Indian Ocean ###
#### INPUT: - time series of REP and NRT SLA,                    ###
####        - array of REP and NRT SLA in the Indian Ocean       ###
#### OUTPUT: - correlation maps for different time lags          ###
#### FUNCTIONS: pyperPete,          ###
####           and Reunion Island                            ###
####         - time series of SLA for                        ###
################################################################
##################
### LIBRARIES  ###
##################
library(ncdf4)
library(progress)

###################
#### FUNCTIONS  ###
###################
## Function to calculate the p-value accounting for the autocorrelation
pyperPeter <- function(X,Y){
  # N = length of time-series
  N <- length(X)
  # if you want to use Pyper and Peterman method and not choose the lag by yourself, uncomment the line below
  l <- N%/%5
  # Computing the autocorrelation of the X time series
  # (Box et Jenkins (1976)
  # Function acf
  # ------------------------------------
  #Here, it opens a window where you can see the autocorrelation function
  # and hence the lag you should use to compute autocorrelation (where it is under the significant value: blue dashed line)
  #windows()
  rxx <- acf(X,na.action=na.pass, plot = FALSE)$acf[seq(2,l+1)]
  # Computing the autocorrelation of the Y time series
  # (Box et Jenkins (1976)
  # Function acf
  # ------------------------------------
  #same as above
  #windows()
  ryy <- acf(Y,na.action=na.pass, plot = FALSE)$acf[seq(2,l+1)]
  # Computing the Nstars
  # 1 for Chelton ; 2 for Chelton modified by Chatfield
  # ------------------------------------------------------
  # Computing the autocorrelation of the X time series following Chatfield (1989)
  # Weightened by N/(N-lag_j)
  # --------------------------------------
  weigth <- rep(N,l)
  for (k in 1:l)
  {
    weigth[k]<-N/(N-k)
  }
  rxx1 <- rxx * weigth
  ryy1 <- ryy * weigth
  
  # Computing the Nstars
  # --------------------------------------
  #somme = sum in French :-)
  somme1 <- rxx*ryy
  somme2 <- rxx1*ryy1
  
  # We only use the positive terms
  somme1 <- (somme1>0)*somme1
  somme2 <- (somme2>0)*somme2
  inv.Nstar1 <- (1/N) + (2/N)*sum(somme1, na.rm=T)
  Nstar1 <- min(1/inv.Nstar1, N)
  inv.Nstar2 <- (1/N) + (2/N)*sum(somme2, na.rm=T)
  Nstar2 <- min(1/inv.Nstar2, N)

  # Significance Tests
  # --------------------------------------------------------
  RXY <- cor(X,Y, use="pairwise.complete.obs")
  T0 <- RXY*(N-2)^0.5/(1-RXY^2)^0.5
  T1 <- RXY*(Nstar1-2)^0.5/(1-RXY^2)^0.5
  T2 <- RXY*(Nstar2-2)^0.5/(1-RXY^2)^0.5
  
  # p-value
  # --------------------------------------------------------
  p.crit0 <- 1-pt(q=abs(T0), df=(N-2), ncp=0, lower.tail = TRUE, log.p = FALSE)
  p.crit1 <- 1-pt(q=abs(T1), df=(Nstar1-2), ncp=0, lower.tail = TRUE, log.p = FALSE)
  p.crit2 <- 1-pt(q=abs(T2), df=(Nstar2-2), ncp=0, lower.tail = TRUE, log.p = FALSE)
  
  return(list(
    coef.cor = RXY,
    p.crit0=p.crit0,
    p.crit1=p.crit1,
    p.crit2=p.crit2))
}

######################
### Load the data  ###
######################
### For REP SLA data
load('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/sshRepData.Rdata')
load('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/timeSeriesSLARep.Rdata')

corRepSLA0      <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
pValueRepSLA0   <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
corRepSLA15     <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
pValueRepSLA15  <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
corRepSLA30     <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
pValueRepSLA30  <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
corRepSLA45     <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
pValueRepSLA45  <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
corRepSLA60     <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
pValueRepSLA60  <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
corRepSLA75     <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
pValueRepSLA75  <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
corRepSLA90     <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
pValueRepSLA90  <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
corRepSLA105    <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))
pValueRepSLA105 <- matrix(NA, nrow=nrow(sshRepData), ncol=ncol(sshRepData))

for (i in 1:nrow(sshRepData)){
  print(i)
  for (j in 1:ncol(sshRepData)){
    corRepSLA0[i,j]      <- pyperPeter(timeSeriesSLARep$SLA, sshRepData[i,j,])$coef.cor
    pValueRepSLA0[i,j]   <- pyperPeter(timeSeriesSLARep$SLA, sshRepData[i,j,])$p.crit2
    corRepSLA15[i,j]     <- pyperPeter(timeSeriesSLARep$SLA[16:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-15)])$coef.cor
    pValueRepSLA15[i,j]  <- pyperPeter(timeSeriesSLARep$SLA[16:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-15)])$p.crit2
    corRepSLA30[i,j]     <- pyperPeter(timeSeriesSLARep$SLA[31:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-30)])$coef.cor
    pValueRepSLA30[i,j]  <- pyperPeter(timeSeriesSLARep$SLA[31:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-30)])$p.crit2
    corRepSLA45[i,j]     <- pyperPeter(timeSeriesSLARep$SLA[46:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-45)])$coef.cor
    pValueRepSLA45[i,j]  <- pyperPeter(timeSeriesSLARep$SLA[46:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-45)])$p.crit2
    corRepSLA60[i,j]     <- pyperPeter(timeSeriesSLARep$SLA[61:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-60)])$coef.cor
    pValueRepSLA60[i,j]  <- pyperPeter(timeSeriesSLARep$SLA[61:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-60)])$p.crit2
    corRepSLA75[i,j]     <- pyperPeter(timeSeriesSLARep$SLA[76:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-75)])$coef.cor
    pValueRepSLA75[i,j]  <- pyperPeter(timeSeriesSLARep$SLA[76:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-75)])$p.crit2
    corRepSLA90[i,j]     <- pyperPeter(timeSeriesSLARep$SLA[91:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-90)])$coef.cor
    pValueRepSLA90[i,j]  <- pyperPeter(timeSeriesSLARep$SLA[91:length(timeSeriesSLARep)],  sshRepData[i,j,1:(dim(sshRepData)[3]-90)])$p.crit2
    corRepSLA105[i,j]    <- pyperPeter(timeSeriesSLARep$SLA[106:length(timeSeriesSLARep)], sshRepData[i,j,1:(dim(sshRepData)[3]-105)])$coef.cor
    pValueRepSLA105[i,j] <- pyperPeter(timeSeriesSLARep$SLA[196:length(timeSeriesSLARep)], sshRepData[i,j,1:(dim(sshRepData)[3]-105)])$p.crit2
      }
}

for (i in 1:length(timeLag)){
  timeLag <- as.character(seq(0,105, 15))
  save(eval(parse(text(paste("corRepSLA", timeLag[i], sep="")))), file=paste("/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/corRepSLA", timeLag[i], ".Rdata", sep=""))
  save(eval(parse(text(paste("pValueRepSLA", timeLag[i], sep="")))), file=paste("/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/pValueRepSLA", timeLag[i], ".Rdata", sep=""))
}

### For REP SLA data
load('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/sshNrtData.Rdata')
load('/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/timeSeriesSLANrt.Rdata')

corNrtSLA0      <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
pValueNrtSLA0   <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
corNrtSLA15     <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
pValueNrtSLA15  <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
corNrtSLA30     <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
pValueNrtSLA30  <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
corNrtSLA45     <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
pValueNrtSLA45  <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
corNrtSLA60     <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
pValueNrtSLA60  <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
corNrtSLA75     <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
pValueNrtSLA75  <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
corNrtSLA90     <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
pValueNrtSLA90  <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))
corNrtSLA105    <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshRepData))
pValueNrtSLA105 <- matrix(NA, nrow=nrow(sshNrtData), ncol=ncol(sshNrtData))

for (i in 1:nrow(sshNrtData)){
  print(i)
  for (j in 1:ncol(sshNrtData)){
    corNrtSLA0[i,j]      <- pyperPeter(timeSeriesSLANrt$SLA, sshNrtData[i,j,])$coef.cor
    pValueNrtSLA0[i,j]   <- pyperPeter(timeSeriesSLANrt$SLA, sshNrtData[i,j,])$p.crit2
    corNrtSLA15[i,j]     <- pyperPeter(timeSeriesSLANrt$SLA[16:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-15)])$coef.cor
    pValueNrtSLA15[i,j]  <- pyperPeter(timeSeriesSLANrt$SLA[16:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-15)])$p.crit2
    corNrtSLA30[i,j]     <- pyperPeter(timeSeriesSLANrt$SLA[31:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-30)])$coef.cor
    pValueNrtSLA30[i,j]  <- pyperPeter(timeSeriesSLANrt$SLA[31:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-30)])$p.crit2
    corRNrtSLA45[i,j]    <- pyperPeter(timeSeriesSLANrt$SLA[46:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-45)])$coef.cor
    pValueNrtSLA45[i,j]  <- pyperPeter(timeSeriesSLANrt$SLA[46:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-45)])$p.crit2
    corNrtSLA60[i,j]     <- pyperPeter(timeSeriesSLANrt$SLA[61:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-60)])$coef.cor
    pValueNrtSLA60[i,j]  <- pyperPeter(timeSeriesSLANrt$SLA[61:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-60)])$p.crit2
    corNrtSLA75[i,j]     <- pyperPeter(timeSeriesSLANrt$SLA[76:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-75)])$coef.cor
    pValueNrtSLA75[i,j]  <- pyperPeter(timeSeriesSLANrt$SLA[76:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-75)])$p.crit2
    corNrtSLA90[i,j]     <- pyperPeter(timeSeriesSLANrt$SLA[91:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-90)])$coef.cor
    pValueNrtSLA90[i,j]  <- pyperPeter(timeSeriesSLANrt$SLA[91:length(timeSeriesSLANrt)],  sshNrtData[i,j,1:(dim(sshNrtData)[3]-90)])$p.crit2
    corNrtSLA105[i,j]    <- pyperPeter(timeSeriesSLANrt$SLA[106:length(timeSeriesSLANrt)], sshNrtData[i,j,1:(dim(sshNrtData)[3]-105)])$coef.cor
    pValueNrtSLA105[i,j] <- pyperPeter(timeSeriesSLANrt$SLA[196:length(timeSeriesSLANrt)], sshNrtData[i,j,1:(dim(sshNrtData)[3]-105)])$p.crit2
  }
}
timeLag <- as.character(seq(0,105, 15))
for (i in 1:length(timeLag)){
  save(eval(parse(text(paste("corNrtSLA", timeLag[i], sep="")))), file=paste("/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/corNrtSLA", timeLag[i], ".Rdata", sep=""))
  save(eval(parse(text(paste("pValueNrtSLA", timeLag[i], sep="")))), file=paste("/home/sbonhomm/ownCloudIfremer/Recherche/dataDiffusion/satData/shiny/data/pValueNrtSLA", timeLag[i], ".Rdata", sep=""))
}

