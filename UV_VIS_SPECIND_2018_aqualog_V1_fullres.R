library(data.table)
library(ggplot2)
library(reshape2)
library(kriging)
library(gstat)
library(geoR)
library(sp)



#You need to read in the files for emcolhead and ex.axis in order to properly handle EEM data inside of loop

#######################################################################################################################################################################
#read in .csv containing Excitation wavlength values (as column)
ex.axis <-  fread("E:/Aqualog/R_scripts/ex_vector_120919.csv",header=TRUE)


#clears any pre-existing spec.final data tables
spec.final = NULL # Setup abs spec indices final table

emcolhead <- fread("E:/Aqualog/R_scripts/em_header.csv")
emcolhead <- as.character(data.frame(emcolhead))


#clears any pre-existing spec.final data tables
spec.final = NULL # Setup abs spec indices final table

########################################################################################################################################################################
# UPDATE location of folder (path) where your Aqualog and MATLAB-exported ASCII are located in next line
wd = 'E:/Aqualog/Data/ASCII_exported_data_processed_outputs_in_matlab_norm_data_subfolders/TREEPEAT' 

setwd(wd)

dirs = dir(getwd())

# Loop through sampling folders from given year

for (i in 1:length(dirs)) {
    
  setwd(paste(wd,dirs[i],"matlab_norm_outputs",sep = "/")) # move to sub folder with processed EEMs and abs
  
  absfiles = list.files(pattern = "*abs_coefficients.csv") # list abs files
  
  EEMfiles = list.files(pattern = "*Units.csv") # list EEM files
  
  spec.merge = NULL

absfiles = list.files(pattern = "*abs_coefficients.csv") # list abs files

EEMfiles = list.files(pattern = "*Units.csv") # list EEM files

spec.merge = NULL

par(mfrow = c(3,3))

# Loop through abs files

  for (j in 1:length(absfiles)) {

  abs <- fread(absfiles[j]) 

  u_IDs = as.data.frame(strsplit(absfiles[j], '[.]'))
  
  unique_ID = as.character(u_IDs[1,])

#log (natural) transform absorbance data, and fit a linear model to calculate spectral slope (240-600, 275,295, 350-400 nm).  Then use predicted lm response surface to calculate absorbance indices (Sr, E2:E3, a254, a300)
setnames(abs, c("Wavelength_nm",  "Absoprtion_coeff_m_1"), c("wl_nm","abs_coef"))
abs<-abs[,abs.log:=log(abs_coef)]

#Spectral Slope for all available data (240:600 nm)
ss.lm<-lm(abs.log~wl_nm,data=abs)
ss.full<-summary(ss.lm)
ss240.600<-ss.full$coefficients[2,1]

#Spectral Slope for (275:295 nm)
ss.275<-summary(lm(abs.log~wl_nm,data=abs[wl_nm>275 & wl_nm<295]))
ss275.295<-ss.275$coefficients[2,1]

#Spectral Slope for (350:400 nm)
ss.350<-summary(lm(abs.log~wl_nm,data=abs[wl_nm>350 & wl_nm<400]))
ss350.400<-ss.350$coefficients[2,1]

#Spectral Ratio (Sr: Helms etal 2008)
Sr<-ss275.295/ss350.400

#---------------------------------------------------------------------------
###predict absorbance at 1 nm intervals using log tranformed linear model
nm1<- seq(240,800,1)
pd.abs = predict(ss.lm,data.frame(wl_nm=nm1),interval="none")
nm1<-data.table(nm1)
pd.abs<-data.table(pd.abs)
nm1<-nm1[,abs_coef_p:=exp(pd.abs)]

#----------------------------------------------------------------------------
#use lm model to calculate absorbance indices
#E2:E3
e2<-nm1[nm1==250]
e3<-nm1[nm1==365]
e2e3<-e2[,abs_coef_p]/e3[,abs_coef_p]

#e4:e6
e4<-nm1[nm1==465]
e6<-nm1[nm1==665]
e4e6<-e4[,abs_coef_p]/e6[,abs_coef_p]

#Absorbance @ 254 
a254<-nm1[nm1==254]
a254<-a254[,abs_coef_p]/2.303

#Absorbance @ 300
a300<-nm1[nm1==300]
a300<-a300[,abs_coef_p]/2.303

#Absorbance @ 350
a350<-nm1[nm1==350]
a350<-a350[,abs_coef_p]/2.303

#Absorbance @ 370
a370<-nm1[nm1==370]
a370<-a370[,abs_coef_p]/2.303


#Absorbance @ 440
a440<-nm1[nm1==440]
a440<-a440[,abs_coef_p]/2.303

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#EEM SPEC INDEX CALCULATIONS

#Read in Corrected EEM file
eem <- fread(EEMfiles[j]) 

#rename EEM columns with emission wavelengths (emcolhead) and paste excitation wavelenghths (ex.axis) as new column
setnames(eem, c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20","V21","V22","V23","V24","V25","V26","V27","V28","V29","V30","V31","V32","V33","V34","V35","V36","V37","V38","V39","V40","V41","V42","V43","V44","V45","V46","V47","V48","V49","V50","V51","V52","V53","V54","V55","V56","V57","V58","V59","V60","V61","V62","V63","V64","V65","V66","V67","V68","V69","V70","V71","V72","V73","V74","V75","V76","V77","V78","V79","V80","V81","V82","V83","V84","V85","V86","V87","V88","V89","V90","V91","V92","V93","V94","V95","V96","V97","V98","V99","V100","V101","V102","V103","V104","V105","V106","V107","V108","V109","V110","V111","V112","V113","V114","V115","V116","V117","V118","V119","V120","V121","V122","V123","V124","V125"), emcolhead)
eem <- eem[,ex.wl:=ex.axis[,wl_nm]]

#melt data into long format (x, y, z), and use melted data table to calculate EEM indices
eem3d<-melt(eem,id=c("ex.wl"))
setnames(eem3d, c("ex.wl","variable",	"value"), c("ex.wl","em.wl","Icor"))
eem3d[,ex.wl:=as.numeric(ex.wl)]
eem3d[,em.wl:=as.character(em.wl)]
eem3d[,em.wl:=as.numeric(em.wl)]

#calculate FI using nearest measured exported pixel
e450<-eem3d[ex.wl==369]
e450<-e450[em.wl==449.123]
e450<-e450[,Icor]

e500<-eem3d[ex.wl==369]
e500<-e500[em.wl==498.898]
e500<-e500[,Icor]

FI<- e450/e500

#calculate Humification Index (HIX) and normalized HIX (HIX.norm) according to T. Ohno 2002 
eem3d435.480 <- eem3d[em.wl>435 & em.wl<480]
emsum435.480<- sum(eem3d435.480[,Icor],na.rm=T)

eem3d300.345 <- eem3d[em.wl>300 & em.wl<345]
emsum300.345<- sum(eem3d300.345[,Icor],na.rm=T)

HIX<-emsum435.480/emsum300.345
HIX.norm<-emsum435.480/(emsum300.345+emsum435.480)

# Calculate Peaks
# A

peakA = eem3d[ex.wl == 261]
peakA = peakA[em.wl == 449.123]
peakA = peakA[,Icor]

# C

peakC = eem3d[ex.wl == 339]
peakC = peakC[em.wl == 439.205]
peakC = peakC[,Icor]

# M

peakM = eem3d[ex.wl == 300]
peakM = peakM[em.wl == 389.828]
peakM = peakM[,Icor]

# T
peakT = eem3d[ex.wl == 261]
peakT = peakT[em.wl == 305.222]
peakT = peakT[,Icor]

#Ct(Kothawala 2012)
CT = eem3d[ex.wl == 282]
CT = CT[em.wl == 363.654]
CT = CT[,Icor]


A.T = peakA/peakT

C.A = peakC/peakA

C.M = peakC/peakM

C.T = peakC/peakT

# Freshness index

FI.380 = eem3d[ex.wl == 309]
FI.380 = FI.380[em.wl == 379.999]
FI.380 = FI.380[,Icor]

FI.420.435 = eem3d[ex.wl == 309]
FI.420.435 = FI.420.435[em.wl>420 & em.wl<435]
FI.420.435 = FI.420.435[,Icor]

Fresh = FI.380/max(FI.420.435, na.rm = T)

# BIX

BIX.380 = eem3d[ex.wl == 309]
BIX.380 = BIX.380[em.wl == 379.999]
BIX.380 = BIX.380[,Icor]


BIX.430 = eem3d[ex.wl == 309]
BIX.430 = BIX.430[em.wl == 429.3]
BIX.430 = BIX.430[,Icor]

BIX = BIX.380/BIX.430

# Relative flourescence efficiency

RFE.FL = eem3d[ex.wl == 369]
RFE.FL = RFE.FL[em.wl == 459.054]
RFE.FL = RFE.FL[,Icor]

RFE.ABS = a370

RFE = RFE.FL/RFE.ABS

# Tryptophan Index

TI = eem3d[ex.wl>275 & ex.wl<285]
TI = TI[em.wl>346 & em.wl<354]
TI = TI[,Icor]

TI = mean(TI)

#generate blank table for Inverse Distance Weighted (IDW) interpolation at 1 nm intervals (em and ex).  Perform IDW on eem3d, then re-calculate FI with IDW table
eem.grid <- data.table(ex.wl = rep(240:800, each=length(212:620)), em.wl = rep(212:620, times=length(240:800)))
eem.grid[,Icor:=NA_real_]
idw.out <- idw(Icor~1, locations = ~ex.wl+em.wl, data = eem3d[!is.na(Icor)], newdata = eem.grid, idp = 2.5, maxdist = 10)
#ggplot(idw.out, aes(x=ex.wl, y=em.wl, z = var1.pred)) + stat_contour() + stat_contour(aes( z = Icor), data= eem3d,col="red")

idw.out<-data.table(idw.out)
e450.idw<-idw.out[ex.wl==370]
e450.idw<-e450.idw[em.wl==450]
e450.idw<-e450.idw[,var1.pred]

e500.idw<-idw.out[ex.wl==370]
e500.idw<-e500.idw[em.wl==500]
e500.idw<-e500.idw[,var1.pred]

FI.idw<- e450.idw/e500.idw

# Rose Cory Redox Index

SQ11 = eem3d[ex.wl>290 & ex.wl<300]
SQ11 = SQ11[em.wl>255 & em.wl<265]
SQ11 = SQ11[,Icor]
SQ11 = mean(SQ11)

SQ12 = eem3d[ex.wl>457 & ex.wl<467]
SQ12 = SQ12[em.wl>255 & em.wl<265]
SQ12 = SQ12[,Icor]
SQ12 = mean(SQ12)

SQ1 = SQ11 + SQ12

SQ21 = eem3d[ex.wl>375 & ex.wl<385]
SQ21 = SQ21[em.wl>265 & em.wl<275]
SQ21 = SQ21[,Icor]
SQ21 = mean(SQ21)

SQ22 = eem3d[ex.wl>515 & ex.wl<525]
SQ22 = SQ22[em.wl>265 & em.wl<275]
SQ22 = SQ22[,Icor]
SQ22 = mean(SQ22)

SQ2 = SQ21+SQ22

SQ31 = eem3d[ex.wl>340 & ex.wl<350]
SQ31 = SQ31[em.wl>260 & em.wl<270]
SQ31 = SQ31[,Icor]
SQ31 = mean(SQ31)

SQ32 = eem3d[ex.wl>407 & ex.wl<417]
SQ32 = SQ32[em.wl>260 & em.wl<270]
SQ32 = SQ32[,Icor]
SQ32 = mean(SQ32)

SQ3 = SQ31+SQ32

HQ = eem3d[ex.wl>245 & ex.wl<255]
HQ = HQ[em.wl>545 & em.wl<555]
HQ = HQ[,Icor]
HQ = mean(HQ)

Q1 = eem3d[ex.wl>250 & ex.wl<260]
Q1 = Q1[em.wl>435 & em.wl<445]
Q1 = Q1[,Icor]
Q1 = mean(Q1)

Q2 = eem3d[ex.wl>250 & ex.wl<260]
Q2 = Q2[em.wl>454 & em.wl<460]
Q2 = Q2[,Icor]
Q2 = mean(Q2)

Q3 = eem3d[ex.wl>250 & ex.wl<260]
Q3 = Q3[em.wl>383 & em.wl<393]
Q3 = Q3[,Icor]
Q3 = mean(Q3)

RI_cory = sum(SQ1 + SQ2 + SQ3 + HQ)/(sum(Q1 + Q2 + Q3)+sum(SQ1 + SQ2 + SQ3 + HQ))


#END OF INDICE CALCULATIONS


#populate summary table with spec indices for each unique_ID
spec.tmp = cbind(unique_ID,ss240.600,ss275.295,ss350.400,Sr,e2e3,e4e6,a254,a300,a350,a440,e450,e500,FI,e450.idw,e500.idw,FI.idw,HIX,HIX.norm, Fresh, A.T, C.A, C.M, C.T,BIX, RFE, TI, RI_cory)

spec.merge = rbind(spec.merge,spec.tmp)

}

spec.final = rbind(spec.final, spec.merge)

rm(spec.merge)

}
##################################################################################################################################
#Set the location (path) you want your spec_indice file saved
setwd(wd)

write.csv(spec.final, file = "spec_indices0810.csv", row.names =F)


#END


