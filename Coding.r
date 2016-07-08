library(Matrix)
library(Kendall)
library(lattice)
library(DAAG)
library(car)
library(locfit)
library(boot)
################################ Import data and Export data ######################
ing95 = read.table(file.choose(),head=T)
head(ing95)
ing95[ing95 == -1.23456e+25]  <- NA
head(ing95)
ing96 = read.table(file.choose(),head=T)
ing96[ing96 == -1.23456e+25]  <- NA
ing96[270:320,]
ing97 = read.table(file.choose(),head=T)
ing97[ing97 == -1.23456e+25]  <- NA
head(ing95)
head(ing96)
head(ing97)
ing = rbind(ing95,ing96,ing97)
dim(ing)

write.csv(ing,file="ing.csv")

NESRS95 = read.table(file.choose(),head=T)
head(NESRS95)
NESRS95[NESRS95 == -1.23456e+25]  <- NA
head(NESRS95)
NESRS96 = read.table(file.choose(),head=T)
NESRS96[NESRS96 == -1.23456e+25]  <- NA
NESRS96[270:320,]
NESRS97 = read.table(file.choose(),head=T)
NESRS97[NESRS97 == -1.23456e+25]  <- NA
head(NESRS95)
head(NESRS96)
head(NESRS97)
NESRS = rbind(NESRS95,NESRS96,NESRS97)
dim(NESRS)

write.csv(NESRS,file="NESRS.csv")

p95 = read.table(file.choose(),head=T)
head(p95)
p95[p95 == -1.23456e+25]  <- NA
head(p95)
p96 = read.table(file.choose(),head=T)
p96[p96 == -1.23456e+25]  <- NA
p96[270:320,]
p97 = read.table(file.choose(),head=T)
p97[p97 == -1.23456e+25]  <- NA
head(p95)
head(p96)
head(p97)
p = rbind(p95,p96,p97)
dim(p)

write.csv(p,file="p.csv")

################################### Data from spot Ing ############################
ing = read.csv(file.choose(),head=T)

ing = ing[1:63046,]
plot(ing$VAP.PRES.GRAD~ing$DATE,type="n",xlab="Time Spot",ylab="vap.press",main="Plot of Ing spot data")
lines(ing$VAP.PRES.GRAD~ing$DATE,lt=2)
hist(ing$VAP.PRES.GRAD,main="Histogram of Vap.Press")
acf(ing$VAP.PRES.GRAD,na.action=na.pass,main="ACF of Vap.Press")
mean(ing$VAP.PRES.GRAD,na.rm = T)
var(ing$VAP.PRES.GRAD,na.rm = T)


date1 <-  as.Date(ing$DATE,'%Y-%m-%d')
year = as.numeric(format(date1,'%Y'))
month = as.numeric(format(date1,'%m'))
day = as.numeric(format(date1,'%d'))
time = as.POSIXct(ing$TIME,format='%H:%M')
hour = as.numeric(format(time,'%H'))
minute = as.numeric(format(time,"%M"))
head(minute)

fulling = cbind(ing,year,month,day,hour,minute)
boxplot(VAP.PRES.GRAD~hour,data=fulling,xlab="Hour", ylab="vap.press",main="Boxplot of vap per hour")
boxplot(VAP.PRES.GRAD~month,data=fulling,xlab="Month", ylab="vap.press",main="Boxplot of vap per month")
boxplot(VAP.PRES.GRAD~day,data=fulling,xlab="Day spot", ylab="vap.press",main="Boxplot of vap per day")


trend= MannKendall(fulling$VAP.PRES.GRAD)
trend

dayaverage = vector()
dayaverage[1]=mean(fulling$VAP.PRES.GRAD[1:59])
for (i in 2:length(unique(fulling$DATE))){
  dayaverage[i] = mean(fulling$VAP.PRES.GRAD[fulling$DATE == fulling$DATE[(i-2)*96+60]])
}
plot(dayaverage,main = "Day average",xlab="Days")
lines(lowess(dayaverage),lw=2,col=4)
dayvar = vector()
dayvar[1]=sqrt(var(fulling$VAP.PRES.GRAD[1:59]))
for (i in 2:length(unique(fulling$DATE))){
  dayvar[i] = sqrt(var(fulling$VAP.PRES.GRAD[fulling$DATE == fulling$DATE[(i-2)*96+60]]))
}
plot(dayvar,main = "Day average variance",xlab="Days")
lines(lowess(dayvar,f=0.3),lw=2,col=4)

daydata = as.data.frame(cbind(as.vector(unique(fulling$DATE)),as.numeric(round(dayaverage,4)),as.numeric(round(dayvar,4))))
head(daydata)
dayvar.trend = MannKendall(dayvar)
dayvar.trend

houraverage = vector()
for (i in 1:24){
  houraverage[i] = median(fulling$VAP.PRES.GRAD[fulling$hour==(i-1)],na.rm=TRUE)
}
plot(houraverage,main = "Hour average ",xlab="Hour")


############  linear regression ########## 
newing = na.omit(fulling)
dim(newing)
model = lm(VAP.PRES.GRAD~.,data=newing[,4:15])
vif(model)

## Refined model 
model1 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP, data=newing[,4:15])
vif(model1)
model2 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newing[,4:15])
vif(model2)

error = model2$residuals
plot(error,main="Resuduals from linear regression",xlab="Time Spots",ylab = "Residuals")
boxplot(error~newing$month+newing$year,main="Boxplot of Resuduals from linear regression",xlab="Month",ylab="Residual")


## Check for the error variance in each month
var = vector()
for(i in (1995:1997)){
  varmonth = vector()
  for(j in (1:12)){
    varmonth[j] = var(error[newing$month==j & newing$year==i])
  }
  var = cbind(var,varmonth)
}
var = as.data.frame(var)
names(var) = c("1995","1996","1997")

## Check for the error variance for each hour 
houraverage.errorvar = vector()
for (i in 1:24){
  houraverage.errorvar[i] = var(error[newing$hour==(i-1)],na.rm=TRUE)
}
plot(houraverage.errorvar,main="Error variance of Hour average",xlab="Hour",ylab="Error Variance")
lines(houraverage.errorvar,lt=2,col=2)


## Calculate variance for each point as weight 
error.var = vector()
for (i in 1:length(error)){
  error.var[i] = var(error[newing$year==newing$year[i]
                           & newing$month == newing$month[i]])
}
unique(error.var)



## iterative weight least square ## 
IRLS1 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newing[,4:15],weights = 1/error.var )
IRLS1$coefficient
error.var2 = vector()
for (i in 1:length(IRLS1$residuals)){
  error.var2[i] = var(IRLS1$residuals[newing$year==newing$year[i]
                                      & newing$month == newing$month[i]])
}



i=1
coefficient = matrix(NA,10,10)
coefficient[,i] = as.vector(model2$coefficient)
coefficient[,i+1] = as.vector(IRLS1$coefficient)

while (any(coefficient[,i] !=coefficient[,i+1])) {
  IRLS = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newing[,4:15],weights = 1/error.var2 )
  coefficient[,i+2] = as.vector(IRLS$coefficient)
  print(coefficient)
  if (i>=5) stop
  error.var2 = vector()
  for (j in 1:length(IRLS$residuals)){
    error.var2[j] = var(IRLS$residuals[newing$year==newing$year[j]
                                       & newing$month == newing$month[j]])
  }
  i = i+1
}
final.model.ing = IRLS
final.model.ing$coefficient




################################ Data from spot P #####################################
P = read.csv(file.choose(),head=T)
names(P)
P = P[1:63776,]
dim(P)
##### Analyzing 
plot(P$VAP.PRES.GRAD,type="n",xlab="time",ylab="vap.press")
lines(P$VAP.PRES.GRAD,lt=2)
hist(P$VAP.PRES.GRAD)
acf(P$VAP.PRES.GRAD,na.action=na.pass)
mean(P$VAP.PRES.GRAD,na.rm = T)
var(P$VAP.PRES.GRAD,na.rm = T)


#### converting date values
date2 <-  as.Date(P$DATE,'%Y-%m-%d')
year = as.numeric(format(date2,'%Y'))
month = as.numeric(format(date2,'%m'))
day = as.numeric(format(date2,'%d'))
time2 = as.POSIXct(P$TIME,format='%H:%M')
hour = as.numeric(format(time2,'%H'))
minute = as.numeric(format(time2,"%M"))
head(minute)

fullP = cbind(P,year,month,day,hour,minute)
boxplot(VAP.PRES.GRAD~hour,data=fullP,xlab="hour", ylab="vap.press")
boxplot(VAP.PRES.GRAD~month,data=fullP,xlab="month", ylab="vap.press")
boxplot(VAP.PRES.GRAD~day,data=fullP,xlab="day", ylab="vap.press")

### test of trend ### 
library(Kendall)
trend= MannKendall(fullP$VAP.PRES.GRAD)
trend

length(unique(fullP$DATE))
Pdayaverage = vector()
Pdayaverage[1]=mean(fullP$VAP.PRES.GRAD[1:17])
for (i in 2:length(unique(fullP$DATE))){
  Pdayaverage[i] = mean(fullP$VAP.PRES.GRAD[fullP$DATE == fullP$DATE[(i-2)*96+18]])
}
plot(Pdayaverage,main = "Day average of spot P", ylab = "Day average",xlab="Day")
Pdayvar = vector()
Pdayvar[1]=sqrt(var(fullP$VAP.PRES.GRAD[1:17]))
for (i in 2:length(unique(fullP$DATE))){
  Pdayvar[i] = sqrt(var(fullP$VAP.PRES.GRAD[fullP$DATE == fullP$DATE[(i-2)*96+18]]))
}
plot(Pdayvar,main = "Day variance of spot P", ylab = "Day variance",xlab="Day")
lines(lowess(Pdayvar,f=0.3),lw=2,col=4)
dayvar.trend = MannKendall(Pdayvar)
dayvar.trend

Phouraverage = vector()
for (i in 1:24){
  Phouraverage[i] = mean(fullP$VAP.PRES.GRAD[fullP$hour==(i-1)],na.rm=TRUE)
}
plot(Phouraverage,main = "Houe mean of spot P", ylab = "Hour mean",xlab="Hour")

acf(P$VAP.PRES.GRAD,na.action = na.pass,main="ACF for vap at ing spot ")

############### linear regression ######################
newP = na.omit(fullP)
dim(newP)


modelP = lm(VAP.PRES.GRAD~.,data=newP[,4:15])
vif(modelP)

## Refined model 
modelP1 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP, data=newP[,4:15])
vif(modelP1)
modelP2 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newP[,4:15])
vif(modelP2)

error.P = modelP2$residuals
plot(error.P,main="Resuduals from linear regression for spot P",xlab="Time Spots",ylab = "Residuals")
boxplot(error.P~newP$month+newP$year,main="Boxplot of Resuduals from linear regression for P",xlab="Month",ylab="Residual")


### Check for the error variance in each month
var.p = vector()
for(i in (1995:1997)){
  varmonth.p = vector()
  for(j in (1:12)){
    varmonth.p[j] = var(error.P[newP$month==j & newP$year==i])
  }
  var.p = cbind(var.p,varmonth.p)
}
var.p = as.data.frame(var.p)
names(var.p) = c("1995","1996","1997")

### Check for the error variance for each hour 
houraverage.errorvar.p = vector()
for (i in 1:24){
  houraverage.errorvar.p[i] = var(error.P[newP$hour==(i-1)],na.rm=TRUE)
}
plot(houraverage.errorvar.p,main="Error variance of Hour average for spot P",xlab="Hour",ylab="Error Variance")
lines(houraverage.errorvar.p,lt=2,col=2)


## Calculate variance for each point as weight 
error.var.p = vector()
for (i in 1:length(error.P)){
  error.var.p[i] = var(error.P[newP$year==newP$year[i]
                               & newP$month == newP$month[i]])
}
unique(error.var.p)

#### WLS 
IRLS1.p = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newP[,4:15],weights = 1/error.var.p )
IRLS1.p$coefficient
error.var2.p = vector()
for (i in 1:length(IRLS1.p$residuals)){
  error.var2.p[i] = var(IRLS1.p$residuals[newP$year==newP$year[i]
                                          & newP$month == newP$month[i]])
}

### Iteration WLS  
i=1
coefficient.p = matrix(NA,10,10)
coefficient.p[,i] = as.vector(modelP2$coefficient)
coefficient.p[,i+1] = as.vector(IRLS1.p$coefficient)

while (any(coefficient.p[,i] !=coefficient.p[,i+1])) {
  IRLS.p = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newP[,4:15],weights = 1/error.var2.p )
  coefficient.p[,i+2] = as.vector(IRLS.p$coefficient)
  print(coefficient.p)
  if(i >= 5) stop("Exceeded maximum number of iterations")
  error.var2.p = vector()
  for (j in 1:length(IRLS.p$residuals)){
    error.var2.p[j] = var(IRLS.p$residuals[newP$year==newP$year[j]
                                           & newP$month == newP$month[j]])
  }
  i = i+1
}

## Final model 
final.model.p = IRLS.p
final.model.p$coefficient




###################################### Data from spot nesre ############################################
################# Import data 
nesre = read.csv(file.choose(),head=T)
names(nesre)
dim(nesre)
##### Analyzing 
plot(nesre$VAP.PRES.GRAD,type="n",xlab="time",ylab="vap.press",main = "Plot for nesre spot")
lines(nesre$VAP.PRES.GRAD,lt=2)
hist(nesre$VAP.PRES.GRAD)
acf(nesre$VAP.PRES.GRAD,na.action=na.pass,main="ACF for nesre spot")
mean(nesre$VAP.PRES.GRAD,na.rm = T)
var(nesre$VAP.PRES.GRAD,na.rm = T)



#### converting date values
date3 <-  as.Date(nesre$DATE,'%Y-%m-%d')
year = as.numeric(format(date3,'%Y'))
month = as.numeric(format(date3,'%m'))
day = as.numeric(format(date3,'%d'))
time3 = as.POSIXct(nesre$TIME,format='%H:%M')
hour = as.numeric(format(time3,'%H'))
minute = as.numeric(format(time3,"%M"))
head(minute)

fullnesre = cbind(nesre,year,month,day,hour,minute)
boxplot(VAP.PRES.GRAD~hour,data=fullnesre,xlab="hour", ylab="vap.press")
boxplot(VAP.PRES.GRAD~month,data=fullnesre,xlab="month", ylab="vap.press")
boxplot(VAP.PRES.GRAD~day,data=fullnesre,xlab="day", ylab="vap.press")

### test of trend ### 
library(Kendall)
trend= MannKendall(fullnesre$VAP.PRES.GRAD)
trend

length(unique(fullnesre$DATE))
nesre.dayaverage = vector()
nesre.dayaverage[1]=mean(fullnesre$VAP.PRES.GRAD[1:52])
for (i in 2:length(unique(fullnesre$DATE))){
  nesre.dayaverage[i] = mean(fullnesre$VAP.PRES.GRAD[fullnesre$DATE == fullnesre$DATE[(i-2)*96+53]])
}
plot(nesre.dayaverage,main="Dayaverage for spot Nesre",xlab = "Day",ylab= "Day averge")
nesre.dayvar = vector()
nesre.dayvar[1]=sqrt(var(fullnesre$VAP.PRES.GRAD[1:52]))
for (i in 2:length(unique(fullnesre$DATE))){
  nesre.dayvar[i] = sqrt(var(fullnesre$VAP.PRES.GRAD[fullnesre$DATE == fullnesre$DATE[(i-2)*96+53]]))
}
plot(nesre.dayvar,main="Day variance for spot Nesre",xlab = "Day",ylab= "Day variance")
lines(lowess(nesre.dayvar,f=0.3),lw=2,col=4)
dayvar.trend = MannKendall(nesre.dayvar)
dayvar.trend

nesre.houraverage = vector()
for (i in 1:24){
  nesre.houraverage[i] = mean(fullnesre$VAP.PRES.GRAD[fullnesre$hour==(i-1)],na.rm=TRUE)
}
plot(nesre.houraverage,main="Hour variance for spot Nesre",xlab = "Hour",ylab= "Hour variance")

acf(nesre$VAP.PRES.GRAD,na.action = na.pass,main="ACF for vap at ing spot ")

################### linear regression ###################
new.nesre = na.omit(fullnesre)
dim(new.nesre)
model.nesre = lm(VAP.PRES.GRAD~.,data=new.nesre[,4:15])
vif(model.nesre)

## Refined model 
model.nesre1 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP, data=new.nesre[,4:15])
vif(model.nesre1)
model.nesre2 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=new.nesre[,4:15])
vif(model.nesre2)

error.N = model.nesre2$residuals
plot(error.N,main="Resuduals from linear regression for spot Nesre",xlab="Time Spots",ylab = "Residuals")
boxplot(error.N~new.nesre$month+new.nesre$year,main="Boxplot of Resuduals from linear regression for Nesre",xlab="Month",ylab="Residual")


### Check for the error variance in each month
var.n = vector()
for(i in (1995:1997)){
  varmonth.n = vector()
  for(j in (1:12)){
    varmonth.n[j] = var(error.N[new.nesre$month==j & new.nesre$year==i])
  }
  var.n = cbind(var.n,varmonth.n)
}
var.n = as.data.frame(var.n)
names(var.n) = c("1995","1996","1997")
var.n   #### variance for month 1,2,3 are larger, variance for the rest of month are smaller

### Check for the error variance for each hour 
houraverage.errorvar.n = vector()
for (i in 1:24){
  houraverage.errorvar.n[i] = var(error.N[new.nesre$hour==(i-1)],na.rm=TRUE)
}
plot(houraverage.errorvar.n,main="Error variance of Hour average for spot N",xlab="Hour",ylab="Error Variance")
lines(houraverage.errorvar.n,lt=2,col=2)


## Calculate variance for each point as weight 
error.var.n = vector()
for (i in 1:length(error.N)){
  error.var.n[i] = var(error.N[new.nesre$year==new.nesre$year[i]
                               & new.nesre$month == new.nesre$month[i]])
}
unique(error.var.n)

#### WLS 
IRLS1.n = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=new.nesre[,4:15],weights = 1/error.var.n )
IRLS1.n$coefficient
error.var2.n = vector()
for (i in 1:length(IRLS1.n$residuals)){
  error.var2.n[i] = var(IRLS1.n$residuals[new.nesre$year==new.nesre$year[i]
                                          & new.nesre$month == new.nesre$month[i]])
}

### Iteration WLS  
i=1
coefficient.n = matrix(NA,10,10)
coefficient.n[,i] = as.vector(model.nesre2$coefficient)
coefficient.n[,i+1] = as.vector(IRLS1.n$coefficient)

while (any(coefficient.n[,i] !=coefficient.n[,i+1])) {
  IRLS.n = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=new.nesre[,4:15],weights = 1/error.var2.n )
  coefficient.n[,i+2] = as.vector(IRLS.n$coefficient)
  print(coefficient.n)
  if(i >= 5) stop("Exceeded maximum number of iterations")
  error.var2.n = vector()
  for (j in 1:length(IRLS.n$residuals)){
    error.var2.n[j] = var(IRLS.n$residuals[new.nesre$year==new.nesre$year[j]
                                           & new.nesre$month == new.nesre$month[j]])
  }
  i = i+1
}

## Final model 
final.model.n = IRLS.n
final.model.n$coefficient




############################### Data from spot C ########################################
################# Import data 
C = read.table(file.choose(),head=T)    ### only data for year 1997
C[C == -1.23456e+25]  <- NA
head(C)
names(C)
dim(C)

##### Analyzing 
plot(C$VAP.PRES.GRAD,type="n",xlab="time",ylab="vap.press",main="Plot for spot C")
lines(C$VAP.PRES.GRAD,lt=2)
hist(C$VAP.PRES.GRAD)
acf(C$VAP.PRES.GRAD,na.action=na.pass,main="ACF for spot C")
mean(C$VAP.PRES.GRAD,na.rm = T)
var(C$VAP.PRES.GRAD,na.rm = T)


#### converting date values
date4 <-  as.Date(C$DATE,'%Y-%m-%d')
year = as.numeric(format(date4,'%Y'))
month = as.numeric(format(date4,'%m'))
day = as.numeric(format(date4,'%d'))
time4 = as.POSIXct(C$TIME,format='%H:%M')
hour = as.numeric(format(time4,'%H'))
minute = as.numeric(format(time4,"%M"))
head(minute)

fullC = cbind(C,year,month,day,hour,minute)
boxplot(VAP.PRES.GRAD~hour,data=fullC,xlab="hour", ylab="vap.press")
boxplot(VAP.PRES.GRAD~month,data=fullC,xlab="month", ylab="vap.press")
boxplot(VAP.PRES.GRAD~day,data=fullC,xlab="day", ylab="vap.press")

### test of trend ### 
library(Kendall)
trend= MannKendall(fullC$VAP.PRES.GRAD)
trend

C.dayaverage = vector()
C.dayaverage[1]=mean(fullC$VAP.PRES.GRAD[1:86])
for (i in 2:length(unique(fullC$DATE))){
  C.dayaverage[i] = mean(fullC$VAP.PRES.GRAD[fullC$DATE == fullC$DATE[(i-2)*96+87]])
}
plot(C.dayaverage,main="Day average for spot C",xlab="Day",ylab = "Day average")
C.dayvar = vector()
C.dayvar[1]=sqrt(var(fullC$VAP.PRES.GRAD[1:86]))
for (i in 2:length(unique(fullC$DATE))){
  C.dayvar[i] = sqrt(var(fullC$VAP.PRES.GRAD[fullC$DATE == fullC$DATE[(i-2)*96+87]]))
}
plot(C.dayvar,main="Day variance for spot C",xlab="Day",ylab = "Day variance")
lines(lowess(C.dayvar,f=0.3),lw=2,col=4)
C.dayvar.trend = MannKendall(C.dayvar)
C.dayvar.trend

C.houraverage = vector()
for (i in 1:24){
  C.houraverage[i] = mean(fullC$VAP.PRES.GRAD[fullC$hour==(i-1)],na.rm=TRUE)
}
plot(C.houraverage,main="Hour average for spot C",xlab="Hour",ylab = "Hour average")
C.have.trend = MannKendall(C.houraverage)
C.have.trend
lines(lowess(C.houraverage,f=0.3),lw=2,col=4)

acf(C$VAP.PRES.GRAD,na.action = na.pass,main="ACF for vap at ing spot ")

############### linear regression ##############
newC = na.omit(fullC)
dim(newC)
model.c = lm(VAP.PRES.GRAD~.,data=newC[,3:14])
vif(model.c)

## Refined model 
model.c1 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP, data=newC[,3:14])
vif(model.c1)
model.c2 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newC[,3:14])
vif(model.c2)

error.c = model.c2$residuals
plot(error.c,main="Resuduals from linear regression for spot C",xlab="Time Spots",ylab = "Residuals")
boxplot(error.c~newC$month,main="Boxplot of Resuduals from linear regression for C in 1997",xlab="Month",ylab="Residual")


### Check for the error variance in each month
var.c = vector()
for(i in (1:12)){
  var.c[i] = var(error.c[newC$month==i])
}
var.c = as.data.frame(var.c)
var.c   #### variance for month 1,2,3 are larger, variance for the rest of month are smaller

### Check for the error variance for each hour 
houraverage.errorvar.c = vector()
for (i in 1:24){
  houraverage.errorvar.c[i] = var(error.c[newC$hour==(i-1)],na.rm=TRUE)
}
plot(houraverage.errorvar.c,main="Error variance of Hour average for spot C",xlab="Hour",ylab="Error Variance")
lines(houraverage.errorvar.c,lt=2,col=2)


## Calculate variance for each point as weight 
error.var.c = vector()
for (i in 1:length(error.c)){
  error.var.c[i] = var(error.c[newC$month == newC$month[i]])
}
unique(error.var.c)


#### WLS 
IRLS1.c = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newC[,3:14],weights = 1/error.var.c )
IRLS1.c$coefficient
error.var2.c = vector()
for (i in 1:length(IRLS1.c$residuals)){
  error.var2.c[i] = var(IRLS1.c$residuals[newC$month == newC$month[i]])
}

### Iteration WLS  
i=1
coefficient.c = matrix(NA,10,10)
coefficient.c[,i] = as.vector(model.c2$coefficient)
coefficient.c[,i+1] = as.vector(IRLS1.c$coefficient)

while (any(coefficient.c[,i] !=coefficient.c[,i+1])) {
  IRLS.c = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=newC[,3:14],weights = 1/error.var2.c )
  coefficient.c[,i+2] = as.vector(IRLS.c$coefficient)
  print(coefficient.c)
  if(i >= 5) stop("Exceeded maximum number of iterations")
  error.var2.c = vector()
  for (j in 1:length(IRLS.c$residuals)){
    error.var2.c[j] = var(IRLS.c$residuals[newC$month == newC$month[j]])
  }
  i = i+1
}

## Final model 
final.model.c = IRLS.c
final.model.c$coefficient




########### Comparison of differenc regressions ##############
## Chow test. 
final.model.ing
final.model.p
final.model.n
final.model.c
final.model.all

### Use all of the data to fit model ### 
ALL = rbind(newing[,-1],newP[,-1],new.nesre[,-1],newC)
model.all = lm(VAP.PRES.GRAD~., data=ALL[,3:14])
vif(model.all)
model.all1 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP, data=ALL[,3:14])
vif(model.all1)
model.all2 = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=ALL[,3:14])
vif(model.all2)
error.all = model.all2$residuals

error.var.all = vector()
for (i in 1:length(error.all)){
  error.var.all[i] = var(error.all[ALL$year == ALL$year[i]
                                   & ALL$month == ALL$month[i]])
}
unique(error.var.all)


### WLS 
IRLS1.all = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=ALL[,3:14],weights = 1/error.var.all )
IRLS1.all$coefficient
error.var.all2 = vector()
for (i in 1:length(IRLS1.all$residuals)){
  error.var.all2[i] = var(IRLS1.all$residuals[ALL$year == ALL$year[i]
                                              & ALL$month == ALL$month[i]])
}
unique(error.var.all2)
### Iteration WLS  
i=1
coefficient.all = matrix(NA,10,10)
coefficient.all[,i] = as.vector(model.all2$coefficient)
coefficient.all[,i+1] = as.vector(IRLS1.all$coefficient)

while (any(coefficient.all[,i] !=coefficient.all[,i+1])) {
  IRLS.all = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=ALL[,3:14],weights = 1/error.var.all2 )
  coefficient.all[,i+2] = as.vector(IRLS.all$coefficient)
  print(coefficient.all)
  if(i >= 5) stop("Exceeded maximum number of iterations")
  error.var.all2 = vector()
  for (j in 1:length(IRLS.all$residuals)){
    error.var.all2[j] = var(IRLS.all$residuals[ALL$year == ALL$year[j]
                                               & ALL$month == ALL$month[j]])
  }
  i = i+1
}



## Final model 
final.model.all = IRLS.all
final.model.all$coefficient
plot(final.model.all$residuals,xlab="Time Spot",ylab="Residuals")




### Evaluate the different among the four models. Chow test
F_value_all = function(model1,model2,model3,model4,modelfull){
  n1 = length(model1$residuals)
  n2 = length(model2$residuals)
  n3 = length(model3$residuals)
  n4 = length(model4$residuals)
  p = 9
  ss1 = sum(model1$residuals^2)
  ss2 = sum(model2$residuals^2)
  ss3 = sum(model3$residuals^2)
  ss4 = sum(model4$residuals^2)
  ssfull = sum(modelfull$residuals^2)
  f = ((ssfull-ss1-ss2-ss3-ss4)/p)/((ss1+ss2+ss3+ss4)/(n1+n2+n3+n4-4*p))
  p = pf(f,p,(n1+n2+n3+n4-4*p),lower.tail=F)
  out = c(f,p)
}
f.value.all = F_value_all(final.model.ing,final.model.p,final.model.n,final.model.c,final.model.all)
f.value.all    





############################## PCA Regression ########################################
pca_reg = function (data,a,b){
  all_pca=prcomp(data[,(a+1):b],center=T)
  all_x = all_pca$x[,1:2]
  all_y = data[,a]
  pca_ols = lm(all_y~all_x)
  out = pca_ols
}
pca_reg_ing = pca_reg(newing,4,15)
pca_reg_p = pca_reg(newP,4,15)
pca_reg_n = pca_reg(new.nesre,4,15)
pca_reg_c = pca_reg(newC,3,14)
pca_reg_all = pca_reg(ALL,3,14)

pca_cv = function(data,a,b){
  n = dim(data)[1]
  fold = floor(n/10)
  pca_mse = vector()
  for (i in 1:10){
    one = (i-1)*fold+1
    two = i*fold
    train = data[-(one:two),a:b]
    test = data[(one:two),a:b]
    train_pca_x = prcomp(train[,-1],center=T)
    train_pca_y = train[,1]
    test_pca_x = prcomp(test[,-1],center=T)
    test_pca_y = test[,1]
    train_x = train_pca_x$x[,1:2]
    test_x = test_pca_x$x[,1:2]
    train = as.data.frame(cbind(train_pca_y,train_x))
    test = as.data.frame(cbind(test_pca_y,test_x))
    names(test) = c("train_pca_y","train_xPC1","train_xPC2")
    names(train) = c("train_pca_y","train_xPC1","train_xPC2")
    model = lm(train_pca_y~.,data=train)
    pred = predict(model,test)
    error = test_pca_y-pred
    pca_mse[i] = mean(error^2)
  }
  out = mean(pca_mse)
}
site_pca_mse = mean(c(pca_cv(newing,4,15),
                      pca_cv(newP,4,15),
                      pca_cv(new.nesre,4,15),
                      pca_cv(newC,3,14)))
site_pca_mse
reegional_pca_mse = pca_cv(ALL,3,14)
reegional_pca_mse



############################## SELECTING DAYS AND CALCULATE R SQUARE #################
## Evaluate by day  
mix = read.csv(file.choose(),head=T)
mix = na.omit(mix)
DATE = unique(mix$DATE)
DATE
ndate = length(DATE)
ndate
## extract R2 
R2 = function(data,a,b,m){
  set = data[as.character(data$DATE) == as.character(DATE[m]),]
  model = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER,data=set[,a:b])
  R2 = summary(model)$adj.r.squared
  out = R2
}

Rsqr_day.ing = vector()
Rsqr_day.P = vector()
Rsqr_day.N = vector()
Rsqr_day.C = vector()
Rsqr_day.all = vector()
for (i in (1:ndate)){
  Rsqr_day.ing[i] = R2(data=newing,a=4,b=15,m=i)
  Rsqr_day.P[i] = R2(data=newP,a=4,b=15,m=i)
  Rsqr_day.N[i] = R2(data=new.nesre,a=4,b=15,m=i)
  Rsqr_day.C[i] = R2(data=newC,a=3,b=14,m=i)
  Rsqr_day.all[i] = R2(data=ALL,a=3,b=14,m=i)
  rsqr = cbind(Rsqr_day.ing,Rsqr_day.P,Rsqr_day.N,Rsqr_day.C,Rsqr_day.all)
}
spot_mean = colMeans(rsqr)
day_mean = rowMeans(rsqr)
mean_site_rsqr = mean(rsqr[,1:4])
mean_all_rsqr = mean(rsqr[,5])


######## Evaluate by month 
Year = unique(mix$year)
Year   ### So only choose the data in year 1997
MONTH = unique(mix$month)
MONTH
## function calculate month R2 
R2month = function(data,a,b,m){
  DAY = unique(subset(mix,month == m)$day)
  n = length(DAY)
  set = subset(data,year==1997 & month==m & (day>=DAY[1]&day<=DAY[n]))
  model = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER,data=set[,a:b])
  R2 = summary(model)$adj.r.squared
  out = R2
}

Rsqr_month.ing  = vector()
Rsqr_month.P = vector()
Rsqr_month.N = vector()
Rsqr_month.C = vector()
Rsqr_month.all = vector()
for (i in MONTH){
  Rsqr_month.ing[i-5] = R2month(data=newing,a=4,b=15,m=i)
  Rsqr_month.P[i-5] = R2month(data=newP,a=4,b=15,m=i)
  Rsqr_month.N[i-5] = R2month(data=new.nesre,a=4,b=15,m=i)
  Rsqr_month.C[i-5] = R2month(data=newC,a=3,b=14,m=i)
  Rsqr_month.all[i-5] = R2month(data=ALL,a=3,b=14,m=i)
  Rsqr_month = cbind(Rsqr_month.ing,Rsqr_month.P,Rsqr_month.N,Rsqr_month.C,Rsqr_month.all)
}
Rsqr_month
month_mean = colMeans(Rsqr_month)
month_mean



### Evaluate by long term -- OLS 
R2long = function(data,a,b){
  model = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER,data=data[,a:b])
  R2 = summary(model)$adj.r.squared
  out = R2
}
site_ols_R2 = c(R2long(newing,4,15),
                R2long(newP,4,15),
                R2long(new.nesre,4,15),
                R2long(newC,3,14))
mean(site_ols_R2)
reginal_ols_R2 = R2long(ALL,3,14)
reginal_ols_R2

### Eveluate by long term-- wls 
site_rsqr = c((summary(final.model.ing)$adj.r.squared),
              (summary(final.model.p)$adj.r.squared),
              (summary(final.model.n)$adj.r.squared),
              (summary(final.model.c)$adj.r.squared))
mean(site_rsqr)
reginal_rsqr = summary(final.model.all)$adj.r.squared
reginal_rsqr


################################# Cross validation for OLS and extract MSE #################
## extract mse 
mse_cv_day = function(data,a,b,m){
  set = data[as.character(data$DATE) == as.character(DATE[m]),]
  model = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER,data=set[,a:b])
  model.cv = cv.lm(data=set[,a:b],form.lm=model,m=10)
  mse = attr(model.cv, "ms")
  out = mse
}

mse_day_ing = vector()
mse_day_P = vector()
mse_day_N = vector()
mse_day_C = vector()
mse_day_all = vector()
for (i in (1:ndate)){
  mse_day_ing[i] = mse_cv_day(data=newing,a=4,b=15,m=i)
  mse_day_P[i] = mse_cv_day(data=newP,a=4,b=15,m=i)
  mse_day_N[i] = mse_cv_day(data=new.nesre,a=4,b=15,m=i)
  mse_day_C[i] = mse_cv_day(data=newC,a=3,b=14,m=i)
  mse_day_all[i] = mse_cv_day(data=ALL,a=3,b=14,m=i)
  mse_day = cbind(mse_day_ing,mse_day_P,mse_day_N,mse_day_C,mse_day_all)
}
mse_spot_mean = colMeans(mse_day)
mse_spot_mean

mse_mean_site = mean(mse_day[,1:4])
mse_mean_all = mean(mse_day[,5])
mse_mean_site
mse_mean_all

######## Evaluate by month 
Year = unique(mix$year)
Year   ### So only choose the data in year 1997
MONTH = unique(mix$month)
MONTH
## function calculate month R2 
mse_cv_month = function(data,a,b,m){
  DAY = unique(subset(mix,month == m)$day)
  n = length(DAY)
  set = subset(data,year==1997 & month==m & (day>=DAY[1]&day<=DAY[n]))
  model = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER,data=set[,a:b])
  model.cv = cv.lm(data=set[,a:b],form.lm=model,m=10)
  mse = attr(model.cv, "ms")
  out = mse
}

mse_month_ing  = vector()
mse_month_P = vector()
mse_month_N = vector()
mse_month_C = vector()
mse_month_all = vector()
for (i in MONTH){
  mse_month_ing[i-5] = mse_cv_month(data=newing,a=4,b=15,m=i)
  mse_month_P[i-5] = mse_cv_month(data=newP,a=4,b=15,m=i)
  mse_month_N[i-5] = mse_cv_month(data=new.nesre,a=4,b=15,m=i)
  mse_month_C[i-5] = mse_cv_month(data=newC,a=3,b=14,m=i)
  mse_month_all[i-5] = mse_cv_month(data=ALL,a=3,b=14,m=i)
  mse_month = cbind(mse_month_ing,mse_month_P,mse_month_N,mse_month_C,mse_month_all)
}
mse_month_mean = colMeans(Rsqr_month)
mse_month_mean


## Long term OLS 
mse_long = function(data,a,b){
  model = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER,data=data[,a:b])
  model.cv = cv.lm(data=data[,a:b],form.lm=model,m=10)
  mse = attr(model.cv, "ms")
  out = mse
}
site_ols_mse = c(mse_long(newing,4,15),
                 mse_long(newP,4,15),
                 mse_long(new.nesre,4,15),
                 mse_long(newC,3,14))
mean(site_ols_mse)
reginal_ols_mse = mse_long(ALL,3,14)
reginal_ols_mse

### Long term WLS 
# Way one. with implemented functions: cv.lm 
model.cv.ing = cv.lm(data=newing[,4:15],form.lm=final.model.ing,m=10)
mse.ing = attr(model.cv.ing, "ms")

model.cv.p = cv.lm(data=newP[,4:15],form.lm=final.model.p,m=10)
mse.p = attr(model.cv.p, "ms")

model.cv.n = cv.lm(data=new.nesre[,4:15],form.lm=final.model.n,m=10)
mse.n = attr(model.cv.n, "ms")

model.cv.c = cv.lm(data=newC[,3:14],form.lm=final.model.c,m=10)
mse.c = attr(model.cv.c, "ms")
mean(c(mse.ing,mse.p,mse.n,mse.c))

model.cv.all = cv.lm(data=ALL[,3:14],form.lm=final.model.all,m=10)
mse.all = attr(model.cv.all, "ms")

# Way two with implemented functions: cv.glm and glm function  
glm_wls_all = glm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=ALL[,3:14],method = "glm.fit" )
wls.cv.all = cv.glm(data=ALL[,3:14],glmfit=glm_wls_all,K=10)
mse.all = wls.cv.all$delta
mse.all

mse_wls = function(data,a,b){
  wls = glm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER, data=data[,a:b],method = "glm.fit" )
  wls.cv = cv.glm(data=data[,a:b],glmfit=wls,K=10)
  out = wls.cv$delta 
}
site_mse_wls = mean(c(mse_wls(newing,4,15),
                      mse_wls(newP,4,15),
                      mse_wls(new.nesre,4,15),
                      mse_wls(newC,3,14)))
site_mse_wls

# Way 3. calculate the mse by hand 
ten_fold_wls = function(data,a,b,weight){
  n = dim(data)[1]
  fold = floor(n/10)
  error = matrix(NA,fold,10)
  for (i in 1:10){
    one = (i-1)*fold+1
    two = i*fold
    train = data[-(one:two),a:b]
    test = data[(one:two),a:b]
    model = lm(VAP.PRES.GRAD~.-BOT.WAT.TEMP-PYRONOMETER,data=train,weights = 1/weight[-(one:two)])
    pred = predict(model,newdata = test)
    error[,i] = test[,1]-pred
  }
  out = mean(error^2)
}
mse_wls_ing = ten_fold_wls(data=newing,a=4,b=15,weight = error.var2)
mse_wls_p = ten_fold_wls(data=newP,a=4,b=15,weight = error.var2.p)
mse_wls_n = ten_fold_wls(data=new.nesre,a=4,b=15,weight = error.var2.n)
mse_wls_c = ten_fold_wls(data=newC,a=3,b=14,weight = error.var2.c)
mse_wls_all = ten_fold_wls(data=ALL,a=3,b=14,weight = error.var.all2)
mean(c(mse_wls_ing,mse_wls_p,mse_wls_n,mse_wls_c))
mse_wls_all

mse_ols_ing = ten_fold_wls(data=newing,a=4,b=15,weight = rep(1,27655))
mse_ols_p = ten_fold_wls(data=newP,a=4,b=15,weight = rep(1,51407))
mse_ols_n = ten_fold_wls(data=new.nesre,a=4,b=15,weight = rep(1,58578))
mse_ols_c = ten_fold_wls(data=newC,a=3,b=15,weight = rep(1,15454))
mse_ols_all = ten_fold_wls(data=ALL,a=3,b=14,weight = rep(1,153094))
mean(c(mse_ols_ing,mse_ols_p,mse_ols_n,mse_ols_c))
mse_ols_all

################################ Test for resiudal normality ####################
qqPlot(model2$residuals,ylab = "OIH ols Redisuals")
qqPlot(final.model.ing$residuals,ylab = "OIH wls Redisuals")

qqPlot(modelP2$residuals,ylab = "P33 ols Redisuals")
qqPlot(final.model.p$residuals,ylab = "P33 wls Redisuals")

qqPlot(model.nesre2$residuals,ylab = "NES ols Redisuals")
qqPlot(final.model.n$residuals,ylab = "NES wls Redisuals")

qqPlot(model.c2$residuals,ylab = "C111 ols Redisuals")
qqPlot(final.model.c$residuals,ylab = "C111 wls Redisuals")

qqPlot(model.all2$residuals,ylab = "Reginal ols Redisuals")
qqPlot(final.model.all$residuals,ylab = "Reginal wls Redisuals")

qqPlot(pca_reg_ing$residuals,ylab="OIH pca residuals")
qqPlot(pca_reg_p$residuals,ylab="P33 pca residuals") 
qqPlot(pca_reg_n$residuals,ylab="NES pca residuals")
qqPlot(pca_reg_c$residuals,ylab="C111 pca residuals")
qqPlot(pca_reg_all$residuals,ylab="Reginal pca residuals") 




################################ Make up for analyzeing ##########################
### ing 
Monthaverage.ing = vector()
for (i in 1:12){
  Monthaverage.ing[i] = mean(fulling$VAP.PRES.GRAD[fulling$month==i],na.rm=TRUE)
}
plot(Monthaverage.ing,xlab="Month",ylab = "Monthly Average Vaper Pressure Gradient")
### P
Monthaverage.p = vector()
for (i in 1:12){
  Monthaverage.p[i] = mean(fullP$VAP.PRES.GRAD[fullP$month==i],na.rm=TRUE)
}
plot(Monthaverage.p,,xlab="Month",ylab = "Monthly Average Vaper Pressure Gradient")
### N 
Monthaverage.n = vector()
for (i in 1:12){
  Monthaverage.n[i] = mean(fullnesre$VAP.PRES.GRAD[fullnesre$month==i],na.rm=TRUE)
}
plot(Monthaverage.n,,xlab="Month",ylab = "Monthly Average Vaper Pressure Gradient")

### C 
Monthaverage.c = vector()
for (i in 1:12){
  Monthaverage.c[i] = mean(fullC$VAP.PRES.GRAD[fullC$month==i],na.rm=TRUE)
}
plot(Monthaverage.c,,xlab="Month",ylab = "Monthly Average Vaper Pressure Gradient")



############################# plot for other variables ###########################
names(ALL)
## for OIH site 
plot(ALL$AIR.TEMP.GRAD~ALL$DATE,type="l",xlab = "Date",ylab = " AIR.TEMP.GRAD")
plot(ALL$SOIL.HEAT.FLUX~ALL$DATE,type="l",xlab = "Date",ylab = "SOIL.HEAT.FLUX")
plot(ALL$TOP.WAT.TEMP~ALL$DATE,type="l",xlab = "Date",ylab = "TOP.WAT.TEMP")
plot(ALL$BOT.WAT.TEMP~ALL$DATE,type="l",xlab = "Date",ylab = "BOT.WAT.TEMP")
plot(ALL$NET.RAD~ALL$DATE,type="l",xlab = "Date",ylab = "NET.RAD")
plot(ALL$WATER.LEVEL~ALL$DATE,type="l",xlab = "Date",ylab = "WATER.LEVEL")
plot(ALL$AIR.TEMP~ALL$DATE,type="l",xlab = "Date",ylab = "AIR.TEMP")
plot(ALL$WIND.DIR~ALL$DATE,type="l",xlab = "Date",ylab = "WIND.DIR")
plot(ALL$WIND.SPEED~ALL$DATE,type="l",xlab = "Date",ylab = "WIND.SPEED")
plot(ALL$REL.HUMID~ALL$DATE,type="l",xlab = "Date",ylab = "REL.HUMID")
plot(ALL$PYRONOMETER~ALL$DATE,type="l",xlab = "Date",ylab = "PYRONOMETER")



