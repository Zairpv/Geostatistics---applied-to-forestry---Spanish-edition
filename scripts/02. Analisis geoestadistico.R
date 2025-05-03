##### 1. PREPARAR CONSOLA


options(max.print=1000000)

#Librerias
library(geoR)
library(scatterplot3d)
library(sm)
library(gstat)
library(lattice)
library(maptools)
library(rgdal)
library(raster)
library(sp)
library(automap)

#Bases de datos
Base_Ho13<-read.csv("Ho_2013.csv")
names(Base_Ho13)
Base_F13<-read.csv("F_2013.csv")
names(Base_F13)
Base_Ho18<-read.csv("Ho_2018.csv")
names(Base_Ho18)
Base_F18<-read.csv("F_2018.csv")
names(Base_F18)


##### 1. DEFINIR GEOBASE
TA<-CRS("+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

coordinates(Base_Ho13)<-~UTMX+UTMY
proj4string(Base_Ho13)<-CRS("+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Base_Ho13_geo<-spTransform(Base_Ho13,TA)

coordinates(Base_F13)<-~UTMX+UTMY
proj4string(Base_F13)<-CRS("+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Base_F13_geo<-spTransform(Base_F13,TA)

coordinates(Base_Ho18)<-~UTMX+UTMY
proj4string(Base_Ho18)<-CRS("+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Base_Ho18_geo<-spTransform(Base_Ho18,TA)

coordinates(Base_F18)<-~UTMX+UTMY
proj4string(Base_F18)<-CRS("+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Base_F18_geo<-spTransform(Base_F18,TA)


##### 2. HACER RASTER CONTENEDOR
ref_raster<-CRS("+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
ext<-extent(540795, 543796, 2277951, 2280951)
raster_white<-raster(ext,res=1,nrow=3000,ncol=3001,ref_raster)
rastergrid<-as(raster_white,"SpatialPixelsDataFrame")
plot(rastergrid) #malla contenedora
points(Base_Ho13_geo,pch=20,col="red")

###Nota. En este enfoque, se realiza el analisis en cada una de las capas por separado debido a que es univariado


##### 3. Ho: HOJARASCA 2013

### Variograma empirico y ajuste a modelo teorico
gs_ho13<-gstat(formula = CC~1,locations = Base_Ho13)
Var1_ho13<-variogram(gs_ho13,width = 50)
plot(Var1_ho13) #variograma empirico

?autofitVariogram
fitvar1_ho13<-autofitVariogram(CC~1,Base_Ho13,model = c("Sph", "Exp", "Gau", "Ste","Exc","Mat","Cir","Lin","Bes","Pen","Per","Wav","Hol","Log","Leg"))
summary(fitvar1_ho13) #Ajuste automático
plot(fitvar1_ho13)
#Nota. El ajuste automático indica que el mejor modelo es Matern
#Parametros: psill=0.1722; range=137.1554; kappa=1.9

#Realizando un ajuste manual
fitvar2_ho13<-variogram(CC~1,Base_Ho13)
plot(fitvar2_ho13,pch=20,cex=1,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "CC Ho 2013")
fitvar3_ho13<-vgm(psill = 0.1722,model = "Mat",range = 137.1554,kappa = 1.9,nugget = 0.39)
varmat_Ho13<-fit.variogram(object = fitvar2_ho13,model =fitvar3_ho13)

#Variograma final
plot(fitvar2_ho13,pch=20,cex=1,col="black",ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "a",model=varmat_Ho13)

### Interpolación kriging
krig_Ho13<-krige(CC~1,Base_Ho13,rastergrid,model=varmat_Ho13)
names(krig_Ho13)  
spplot(krig_Ho13) #Esta función tarda bastante (~1h) debido a la cantidad de pixeles
raster_predHo13<-rasterFromXYZ(krig_Ho13)
#writeRaster(raster_predHo13, "ho13_krig.img") #Tarda alrededor de 3 a 5 horas.... El archivo img ya esta generado en la carpeta
raster_sdHo13<-raster(krig_Ho13,layer="var1.var")
#writeRaster(raster_sdHo13, "ho13_krig_sd.img") #Tarda alrededor de 3 a 5 horas.... El archivo img ya esta generado en la carpeta


##### 4. FERMENTACION 2013
#Crear variograma empirico
gs_F13<-gstat(formula = CC~1,locations = Base_F13)
Var1_F13<-variogram(gs_F13,width = 50)
plot(Var1_F13)

#Ajuste a modelo teórico
fitvar1_F13<-autofitVariogram(CC~1,Base_F13,model = c("Sph", "Exp", "Gau", "Ste","Exc","Mat","Cir","Lin","Bes","Pen","Per","Wav","Hol","Log","Leg"))
summary(fitvar1_F13)
plot(fitvar1_F13)
#Nota. El ajuste automático indica que el mejor modelo es Ste
#Parametros: psill=29.2324 range=49.7256; kappa=1.1, nugget=0

#ajuste manual
fitvar2_F13<-variogram(CC~1,Base_F13)
plot(fitvar2_F13,pch=20,cex=1,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "CC F 2013")
fitvar3_F13<-vgm(psill = 29.2324,model = "Ste",range = 49.7256,nugget = 0,kappa = 1.1)
varste_F13<-fit.variogram(object = fitvar2_F13,model =fitvar3_F13)

#Variograma final
plot(fitvar2_F13,pch=20,cex=1,col="black",ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "b",model=varste_F13)

#Interpolacion kriging
krig_F13<-krige(CC~1,Base_F13,rastergrid,model=varste_F13)
names(krig_F13)  
spplot(krig_F13) #Esta funcion tarda aprox 1 hora
raster_predF13<-rasterFromXYZ(krig_F13)
#writeRaster(raster_predF13, "F13_krig.img") #Tarda alrededor de 3 a 5 horas.... El archivo img ya esta generado en la carpeta
raster_sdF13<-raster(krig_F13,layer="var1.var")
#writeRaster(raster_sdF13, "F13_krig_sd.img") #Tarda alrededor de 3 a 5 horas.... El archivo img ya esta generado en la carpeta



##### 5. Ho: HOJARASCA 2018

### Variograma empirico y ajuste a modelo teorico
gs_ho18<-gstat(formula = CC~1,locations = Base_Ho18)
Var1_ho18<-variogram(gs_ho18,width = 50)
plot(Var1_ho18) #variograma empirico

?autofitVariogram
fitvar1_ho18<-autofitVariogram(CC~1,Base_Ho18,model = c("Sph", "Exp", "Gau", "Ste","Exc","Mat","Cir","Lin","Bes","Pen","Per","Wav","Hol","Log","Leg"))
summary(fitvar1_ho18) #Ajuste automático
plot(fitvar1_ho18)
#Nota. El ajuste automático indica que el mejor modelo es Per
#Parametros: psill=0.03709; range=58.1450; nugget=0.3732

#Realizando un ajuste manual
fitvar2_ho18<-variogram(CC~1,Base_Ho18)
plot(fitvar2_ho18,pch=20,cex=1,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "CC Ho 2018")
fitvar3_ho18<-vgm(psill = 0.03709,model = "Per",range = 58.1450,nugget = 0.3732)
varper_Ho18<-fit.variogram(object = fitvar2_ho18,model =fitvar3_ho18)

#Variograma final
plot(fitvar2_ho18,pch=20,cex=1,col="black",ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "a",model=varper_Ho18)

### Interpolación kriging
krig_Ho18<-krige(CC~1,Base_Ho18,rastergrid,model=varper_Ho18)
names(krig_Ho18)  
spplot(krig_Ho18) #Esta función tarda bastante (~1h) debido a la cantidad de pixeles
raster_predHo18<-rasterFromXYZ(krig_Ho18)
#writeRaster(raster_predHo18, "ho18_krig.img") #Tarda alrededor de 3 a 5 horas.... El archivo img ya esta generado en la carpeta
raster_sdHo18<-raster(krig_Ho18,layer="var1.var")
#writeRaster(raster_sdHo18, "ho18_krig_sd.img") #Tarda alrededor de 3 a 5 horas.... El archivo img ya esta generado en la carpeta


##### 6. FERMENTACION 2018
#Crear variograma empirico
gs_F18<-gstat(formula = CC~1,locations = Base_F18)
Var1_F18<-variogram(gs_F18,width = 50)
plot(Var1_F18)

#Ajuste a modelo teórico
fitvar1_F18<-autofitVariogram(CC~1,Base_F18,model = c("Sph", "Exp", "Gau", "Ste","Exc","Mat","Cir","Lin","Bes","Pen","Per","Wav","Hol","Log","Leg"))
summary(fitvar1_F18)
plot(fitvar1_F18)
#Nota. El ajuste automático indica que el mejor modelo es Ste
#Parametros: psill=2.3423 range=116.0586; kappa=0.2, nugget=1.5059

#ajuste manual
fitvar2_F18<-variogram(CC~1,Base_F18)
plot(fitvar2_F18,pch=20,cex=1,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "CC F 2018")
fitvar3_F18<-vgm(psill = 2.3423,model = "Ste",range = 116.0586,nugget = 1.5059,kappa = 0.2)
varste_F18<-fit.variogram(object = fitvar2_F18,model =fitvar3_F18)

#Variograma final
plot(fitvar2_F18,pch=20,cex=1,col="black",ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "b",model=varste_F18)

#Interpolacion kriging
krig_F18<-krige(CC~1,Base_F18,rastergrid,model=varste_F18)
names(krig_F18)  
spplot(krig_F18) #Esta funcion tarda aprox 1 hora
raster_predF18<-rasterFromXYZ(krig_F18)
#writeRaster(raster_predF18, "F18_krig.img")  #Tarda alrededor de 3 a 5 horas.... El archivo img ya esta generado en la carpeta
raster_sdF18<-raster(krig_F18,layer="var1.var")
#writeRaster(raster_sdF18, "F18_krig_sd.img") #Tarda alrededor de 3 a 5 horas.... El archivo img ya esta generado en la carpeta



##### 7. GRAFICAR VARIOGRAMAS
opts <- trellis.par.get()
opts$add.text$col <- "transparent"

#Variogramas 2013
VARIOG_HO13<-plot(fitvar1_ho13)
V1<-update(VARIOG_HO13, par.settings = opts)
VARIOG_F13<-plot(fitvar1_F13)
V2<-update(VARIOG_F13,par.settings=opts)

#Variogramas 2018
VARIOG_HO18<-plot(fitvar1_ho18)
V3<-update(VARIOG_HO18, par.settings = opts)
VARIOG_F18<-plot(fitvar1_F18)
V4<-update(VARIOG_F18,par.settings=opts)

library(gridExtra)
grid.arrange(V1,V2,V3,V4,nrow=2)



#Segundo tipo de variogramas (Reportados en tesis)
VG1<-plot(fitvar2_ho13,pch=20,cex=0.5,col="black",ylab=expression("Semivariance ("*gamma*")"),
          xlab="Distance (m)", main = "Ho",model=varmat_Ho13,sub="2013")
VG2<-plot(fitvar2_F13,pch=20,cex=0.5,col="black",ylab=expression("Semivariance ("*gamma*")"),
          xlab="Distance (m)", main = "F",model=varste_F13,sub="2013")
VG3<-plot(fitvar2_ho18,pch=20,cex=0.5,col="black",ylab=expression("Semivariance ("*gamma*")"),
          xlab="Distance (m)", main = "",model=varper_Ho18,sub="2018")
VG4<-plot(fitvar2_F18,pch=20,cex=0.5,col="black",ylab=expression("Semivariance ("*gamma*")"),
          xlab="Distance (m)", main = "",model=varste_F18,sub="2018")
dev.off()
tiff("Semivariogramas.tiff",width = 190,height = 145,units = "mm",res = 1000)
opar<-par(mar=c(0,0,0,0),
          oma=c(0,0,0,0),
          cex.lab=0.1,font.lab=2,font=1,
          cex.main=0.8,col="black",adj = 0)
grid.arrange(VG1,VG2,VG3,VG4,nrow=2)
dev.off()



##### 8. VALIDACION CRUZADA
#Ho13
cvOK_Ho13<-krige.cv(CC~1,Base_Ho13,nfold=10,model=varmat_Ho13)
names(cvOK_Ho13)
write.csv(cvOK_Ho13,"CVOK_Ho13.csv")
plot(cvOK_Ho13$var1.pred~cvOK_Ho13$observed,cvOK_Ho13, main="CV Ho 13",xlab="Observed",ylab="Predicted")
cor(cvOK_Ho13$var1.pred,cvOK_Ho13$observed) #R2
mean(cvOK_Ho13$residual) #ME
sqrt(sum(cvOK_Ho13$residual^2)/length(cvOK_Ho13$residual)) #RMSE
sum(abs(cvOK_Ho13$residual))/length(cvOK_Ho13$residual) #MAE

#F13
cvOK_F13<-krige.cv(CC~1,Base_F13,nfold=10,model=varste_F13)
names(cvOK_F13)
write.csv(cvOK_F13,"CVOK_F13.csv")
plot(cvOK_F13$var1.pred~cvOK_F13$observed,cvOK_F13, main="CV F 13",xlab="Observed",ylab="Predicted")
cor(cvOK_F13$var1.pred,cvOK_F13$observed) #R2
mean(cvOK_F13$residual) #ME
sqrt(sum(cvOK_F13$residual^2)/length(cvOK_F13$residual)) #RMSE
sum(abs(cvOK_F13$residual))/length(cvOK_F13$residual) #MAE

#Ho18
cvOK_Ho18<-krige.cv(CC~1,Base_Ho18,nfold=10,model=varper_Ho18)
names(cvOK_Ho18)
write.csv(cvOK_Ho18,"CVOK_Ho18.csv")
plot(cvOK_Ho18$var1.pred~cvOK_Ho18$observed,cvOK_Ho18, main="CV Ho 18",xlab="Observed",ylab="Predicted")
cor(cvOK_Ho18$var1.pred,cvOK_Ho18$observed) #R2
mean(cvOK_Ho18$residual) #ME
sqrt(sum(cvOK_Ho18$residual^2)/length(cvOK_Ho18$residual)) #RMSE

#F18
cvOK_F18<-krige.cv(CC~1,Base_F18,nfold=10,model=varste_F18)
names(cvOK_F18)
write.csv(cvOK_F18,"CVOK_F18.csv")
plot(cvOK_F18$var1.pred~cvOK_F18$observed,cvOK_F18, main="CV F 18",xlab="Observed",ylab="Predicted")
cor(cvOK_F18$var1.pred,cvOK_F18$observed) #R2
mean(cvOK_F18$residual) #ME
sqrt(sum(cvOK_F18$residual^2)/length(cvOK_F18$residual)) #RMSE



##### 9. GRAFICAS DE VALIDACI?N

CVOK_1<-plot(cvOK_Ho13$var1.pred~cvOK_Ho13$observed,cvOK_Ho13, main="Ho",xlab="Observed",ylab="Predicted",sub="2013")
CVOK_2<-plot(cvOK_F13$var1.pred~cvOK_F13$observed,cvOK_F13, main="F",xlab="Observed",ylab="Predicted",sub="2013")
CVOK_3<-plot(cvOK_Ho18$var1.pred~cvOK_Ho18$observed,cvOK_Ho18, main="",xlab="Observed",ylab="Predicted",sub="2018")
CVOK_4<-plot(cvOK_F18$var1.pred~cvOK_F18$observed,cvOK_F18, main="",xlab="Observed",ylab="Predicted",sub="2018")

library(ggplot2)
#Default settings
theme_set(theme_bw()) # A ggplot2 theme with white background

CVHO13_OK<-as.data.frame(cvOK_Ho13)
VALID_CVOK1<-ggplot(CVHO13_OK,aes(x=observed,y=var1.pred))+
  geom_point(size=.2)+
  scale_y_continuous(limits=c(0, 5))+
  scale_x_continuous(limits=c(0, 5))+
  xlab("")+
  ylab(expression("Predicted CC  "("Mg ha"^-1)))+
  labs(title="Ho",subtitle="2013")+
  stat_smooth(method = "lm",se=FALSE,color="black",size=.4)+
  geom_abline(intercept = 0,slope = 1, linetype="dashed",size=.3)+ 
  theme(text = element_text(size=10),plot.subtitle = element_text(hjust = 1),axis.text = element_text(size = 10),plot.title = element_text(hjust = 0.5))
VALID_CVOK1

CVF13_OK<-as.data.frame(cvOK_F13)
VALID_CVOK2<-ggplot(CVF13_OK,aes(x=observed,y=var1.pred))+
  geom_point(size=.2)+
  scale_y_continuous(limits=c(0, 30))+
  scale_x_continuous(limits=c(0, 30))+
  xlab("")+
  ylab("")+
  labs(title="F",subtitle="2013")+
  stat_smooth(method = "lm",se=FALSE,color="black",size=.4)+
  geom_abline(intercept = 0,slope = 1, linetype="dashed",size=.3)+ 
  theme(text = element_text(size=10),plot.subtitle = element_text(hjust = 1),axis.text = element_text(size = 10),plot.title = element_text(hjust = 0.5))
VALID_CVOK2

CVHO18_OK<-as.data.frame(cvOK_Ho18)
VALID_CVOK3<-ggplot(CVHO18_OK,aes(x=observed,y=var1.pred))+
  geom_point(size=.2)+
  scale_y_continuous(limits=c(0, 6))+
  scale_x_continuous(limits=c(0, 6))+
  xlab(expression("Observed CC "("Mg ha"^-1)))+
  ylab(expression("Predicted CC "("Mg ha"^-1)))+
  labs(title="",subtitle="2018") +
  stat_smooth(method = "lm",se=FALSE,color="black",size=.4)+
  geom_abline(intercept = 0,slope = 1, linetype="dashed",size=.3)+ 
  theme(text = element_text(size=10),plot.subtitle = element_text(hjust = 1),axis.text = element_text(size = 10),plot.title = element_text(hjust = 0.5))
VALID_CVOK3

CVF18_OK<-as.data.frame(cvOK_F18)
VALID_CVOK4<-ggplot(CVF18_OK,aes(x=observed,y=var1.pred))+
  geom_point(size=.2)+
  scale_y_continuous(limits=c(0, 10))+
  scale_x_continuous(limits=c(0, 10))+
  xlab(expression("Observed CC "("Mg ha"^-1)))+
  ylab("")+
  labs(title="",subtitle="2018") +
  stat_smooth(method = "lm",se=FALSE,color="black",size=.4)+
  geom_abline(intercept = 0,slope = 1, linetype="dashed",size=.3)+ 
  theme(text = element_text(size=10),plot.subtitle = element_text(hjust = 1),axis.text = element_text(size = 10),plot.title = element_text(hjust = 0.5))
VALID_CVOK4

library(ggpubr)
ggarrange(VALID_CVOK1,VALID_CVOK2,VALID_CVOK3,VALID_CVOK4)
ggsave("Validacion cruzada.tiff",width = 190,height = 140,units = "mm",dpi = 500,compression="lzw")

