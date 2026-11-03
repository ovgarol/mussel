# SPDX-FileCopyrightText: 2026 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>

library(readxl)
library(quantreg)
library(scales)
library(lubridate)
library(deSolve)
library(GA)
library(latex2exp)

####

data = NULL
for(i in c(2014,2019,2020,2021,2022,2023)){
  dd = read_xlsx(paste0('./data/mussel_',i,'.xlsx'))
  dd$date = i
  dd = subset(dd,select=c(bank,bank_type,part,area,perimeter,date,date_e))
  dd$date_e = ymd(dd$date_e)
  data = rbind(data,dd)
}

data$bank = substr(data$bank,1,2)

decimal_year = function(date) {
  date = as.Date(date)
  year = as.numeric(format(date, "%Y"))
  start = as.Date(paste0(year, "-01-01"))
  end = as.Date(paste0(year + 1, "-01-01"))
  
  return(year + as.numeric(date - start) / as.numeric(end - start))
}

if(T){
  # Build one annual summary row per bank and bank type by combining the
  # summed perimeter, summed area, and mean partition fraction.
  data0 = merge(
    aggregate(perimeter~date_e+date+bank+bank_type,data=data,sum),
    aggregate(area~date+bank+bank_type,data=data,sum),
    by=c('date','bank','bank_type')
  )
  
  data = merge(
    data0,
    aggregate(part~date+bank+bank_type,data=data,mean),
    by=c('date','bank','bank_type')
  )
  
}

if(interactive()){
  plot(perimeter~area,data=data[data$perimeter<5e5,],log='xy',pch=19)
  abline(s <- lm(log10(perimeter)~log10(area),data=data[data$perimeter<5e5,]),col='red',lwd=2)
  summary(s)
  
  plot(data$area,data$part,log='xy',pch=19)
  abline(s <- lm(log10(part)~log10(area),data=data[data$perimeter<5e5,]),col='red',lwd=2)
  summary(s)
  
  rqfit = rq(log10(part)~log10(area),data=data[data$perimeter<5e5,],tau=0.9)
  abline(rqfit,col='green',lwd=2)
  summary(rqfit)
  
  par(mfrow=c(2,2))
  for(i in unique(data$bank)){
    data0 = data[data$bank == i,]
    
    plot(data0$date,data0$area,col=as.factor(data0$bank_type),pch=19)
    plot(data0$date,data0$perimeter,col=as.factor(data0$bank_type),pch=19)
    mtext(side=4,i,outer=F,line=0)
  }
}

unique(data$bank_type)

####

par(lend=1,las=1)
if(!exists('air.temp')){
  air.temp = read.delim("./data/Cuxhaven_DWD!Lufttemperatur.txt", header=FALSE, comment.char="#")
  air.temp$V1 = as.Date(air.temp$V1)
  air.temp = air.temp[year(air.temp$V1)%in%(2014:2024),]
  air.temp = aggregate(air.temp,by=list(air.temp$V1),function(...)mean(...,na.rm=T))
  air.temp = subset(air.temp, select=c(V1,V2))
  
  salinity = read.delim("./data/LZ_AL!Salzgehalt.txt", header=FALSE, comment.char="#")
  salinity$V1 = as.Date(salinity$V1)
  salinity = salinity[year(salinity$V1)%in%(2014:2024),]
  salinity = aggregate(salinity,by=list(salinity$V1),function(...)mean(...,na.rm=T))
  salinity = subset(salinity, select=c(V1,V2))
  
  # Merge the environmental series and derive a stress index from the
  # combined effect of temperature and salinity thresholds.
  envir = merge(air.temp,salinity, by='V1',all=T,suffixes = c('.air','.sal'))
  envir$stress = 1/(1+exp((10-envir$V2.air)))*1/(1+exp((envir$V2.sal-15))) + 1/(1+exp((20-envir$V2.air))) 
  envir$stress = pmin(envir$stress,1)
  stress.year = aggregate(envir$stress>.5,by=list(year(envir$V1)),function(...)sum(...,na.rm=T))$x
  
  write.csv(data,'./data/mussel_beds.csv',row.names = F)
  colnames(air.temp) = c('date','temperature')
  write.csv(air.temp,'./data/air_temperature.csv',row.names=F)
  colnames(salinity) = c('date','salinity')
  write.csv(salinity,'./data/salinity.csv',row.names=F)
  stress.year.0 = NULL
  stress.year.0$V1 = unique(year(envir$V1))
  stress.year.0$V2 = stress.year
  stress.year.0 = as.data.frame(stress.year.0)
  colnames(stress.year.0) = c('year','days_with_extreme_conditions')
  write.csv(stress.year.0,'./data/stress_year.csv',row.names=F)
}

################################################################################

data = data[data$date>2017,]
#"HE" "HT" "LT" "MB" "NA" "NH" "RL" "SA"
i = 'MB'
print(unique(data$bank))
data = data[data$bank == i & data$bank_type == "Austernbank",]

data$date = data$date_e
data = aggregate(data,by=list(date=data$date),function(...)ifelse(is.numeric(...),sum(...),...))
data = subset(data,select=c(date,area,perimeter))

### new plot

# Create the diagnostic figure that compares the observed bank metrics with
# the environmental stress signal over the study period.
pdf('./fig/data_plot.pdf',width=6,height = 4)
par(mfrow=c(3,1),mgp=c(2,.5,0),tck=0.05,mai=0.1*c(1,4,2,1),oma=0+c(1,1,0,0),las=1)
plot(data$date,data$area/1e6,
     type='o',
     log='',
     ylab=TeX('km$^2$'),
     xlab='',
     xaxt='n',
     col='black',
     pch=19,
     xlim=as.Date(c('2018-01-01','2023-01-01'))
)
axis.Date(side=1,labels = F)
title(main='a. Bank area (A)',line=0.2,adj=0.01)

plot(data$date+0.,data$perimeter/1e3,
     ylab='km',
     xlab='',
     xaxt='n',
     type='o',
     log='',
     col='black',
     pch=19,
     xlim=as.Date(c('2018-01-01','2023-01-01'))
)
axis.Date(side=1,labels = F)
title(main='b. Bank perimeter (P)',line=0.2,adj=0.01)

envir.y=aggregate(envir$stress > 0.25, by = list(year = year(envir$V1)), FUN = function(z) sum(z, na.rm = TRUE))

plot(envir.y$year,envir.y$x/90,type='o',xlim=c(2018,2023),pch=19,ylab='',xlab='')
title(main='c. Summers with stress index > 0.5',line=0.2,adj=0.01)

plot(envir$V1,envir$V2.air,type='l',xlab='',xlim=as.Date(paste0(c(2018,2023),'-01-01')),
     main='',ylab='ºC',xaxt='n')
axis.Date(side=1,labels = F)
title(main='a. Air temperature (T)',line=0.2,adj=0.01)
abline(h=c(10,20),lty=2,col='gray')

plot(envir$V1,envir$V2.sal,type='l', ylab='PSU',xlim=as.Date(paste0(c(2018,2023),'-01-01')),
     xlab='',xaxt='n')
axis.Date(side=1,labels = F)
title(main='b. Salinity (S)',line=0.2,adj=0.01)
abline(h=c(15),lty=2,col='gray')

plot(envir$V1,envir$stress,type='l',main='',ylab='',xlab='',xlim=as.Date(paste0(c(2018,2023),'-01-01')))
envir$stress[is.nan(envir$stress)]=0
title(main='c. Stress index',line=0.2,adj=0.01)
abline(h=c(.5),lty=2,col='gray')

dev.off()

#########

data0 = data
data0$time = round(decimal_date(data0$date),digits=2)
data0$A.real = data0$area/1e6
data0$P.real = data0$perimeter/1e3
data0 = subset(data0,select=c(time,A.real,P.real))

t = round(seq(2018,2023,by=1/100),digits=3)

# Interpolate the environmental time series so the ODE solver can evaluate
# stress and temperature at any requested time point.
stress.approx=approxfun(decimal_date(envir$V1),envir$stress)
temp.approx=approxfun(decimal_date(envir$V1),envir$V2.air)
P0 = 15
system.equations = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    stressor = stress.approx(t)
    temperature = temp.approx(t)
    R = r * 2.^(0.1*(temperature-10))*ifelse(temperature>20,.01,1)
    f =  exp(-P/P0)
    df.dP = -exp(-P/P0)/P0

    dA.dt =R*A*f*(1-4*pi*A/P^2) - c*A*(stressor+.01*k) 
    dP.dt = mu*R*(df.dP*(1-4*pi*A/P^2) + f*2*4*pi*A/P^3) + m0
    list(c(dA.dt, dP.dt))
  })
}

parameters = c(r =58, c= 14, mu = 81, k= 6, m0= 0.61)
state = c(A = aggregate(data0$A.real,by=list(data0$time),median)[1,2], 
          P = aggregate(data0$P.real,by=list(data0$time),median)[1,2])


# The fitness function evaluates how well a candidate parameter set reproduces
# the observed area and perimeter trajectories under the ODE model.
fitness.func =function(x){
  parameters = c(r=x[1],c=x[2],mu=x[3],k=x[4],m0=x[5])
  cc = deSolve::ode(y = state, times = t, func = system.equations, parms = parameters, method = 'euler')
  cc = as.data.frame(cc)
  if(any(is.nan(cc$A))) return(-999)
  if(any(cc$A<=0.1 | cc$A>2.5) ) return(-999)
  
  cc = merge(cc, data0, by='time',all=F)
  n = length(cc$time)

  return((-1*sum((cc$A - cc$A.real)**2) - .1*sum((cc$P - cc$P.real)**2) -12*pi*sum((cc$A/cc$P^2 - cc$A.real/cc$P.real^2)**2))/n ) 
}

x1 = c(25.98483, 17.41540, 63.42554, 49.81319, 33.29349)
x2 = c(25.46246, 17.87520, 65.79566, 47.98589, 34.11554)
x3 = c(25.72086, 17.96333, 65.55289, 47.98679, 34.16258)
x4 = c(25.72086, 17.96333, 65.55289, 47.98679, 34.16035)
suggs = matrix(c(x1,x2,x3,x4),ncol=5,byrow=T) 
#parameters = c(r =58, c= 14, mu = 81, k= 6, m0= 0.61)
low.b = c(0,0,0,1,0)
high.b = 1*c(100,100,100,100,50)

GA1 = ga(type = "real-valued", 
         fitness = fitness.func,
         lower = low.b, 
         upper = high.b,
         maxiter = 10,
         optim=F,
         run=10,
         suggestions = suggs,
         parallel=F,
         pmutation = .95,
         popSize = 100
)

GA1 = ga(type = "real-valued", 
         fitness = fitness.func,
         lower = low.b, 
         upper = high.b,
         maxiter = 100,
         suggestions = GA1@solution,
         optim=F,
         run=100,
         parallel=F,
         popSize = 10
)

if(F)GA1 = ga(type = "real-valued", 
         fitness = fitness.func,
         lower = low.b, 
         upper = high.b,
         maxiter = 100,
         optim=T,
         suggestions = GA1@solution,
         run=100,
         parallel=F,
         popSize = 150
)

plot(GA1)

solution = GA1@solution
x = unname(c(solution[1,]))
parameters = c(r=x[1],c=x[2],mu=x[3],k=x[4],m0=x[5])
cc = ode(y = state, times = t, func = system.equations, parms = parameters,method = 'euler')
cc.center = as.data.frame(cc)

pdf('./fig/model_plot.pdf',width=6,height = 8)
par(mfrow=c(5,1),las=1,mgp=c(2,.5,0),tck=0.05,mai=0.1*c(1,4,2,1),oma=0+c(1,1,0,0))
par(las=1)
col.2 = alpha('gray',.5)

plot(cc.center$time,cc.center$A,main='',type='n',ylab=TeX('km$^2$'),xlab='',xaxt='n',ylim=range(data0$A.real)*c(0.25,1.5))
title(main='a. Area (A)',line=0.2,adj=0.01)
axis.Date(side=1,labels=F)
if(F)for(i in 1:9){
  cc = ode(y = state*(1+0.05*rnorm(2)), times = t, func = system.equations, parms = parameters*(1+0.05*rnorm(5)),method = 'euler')
  cc = as.data.frame(cc)
  lines(cc$time,cc$A,col=col.2,lwd=.5)
}
points(data0$time,data0$A.real,pch=19)
lines(cc.center$time,cc.center$A)
legend('topright',legend=c('observation','model'),pch=c(19,NA),lty=c(NA,1),bty='n')

plot(cc.center$time,cc.center$P,main='',type='n',ylab='km',xaxt='n',xlab='',ylim=range(data0$P.real)*c(0.5,1.5))
title(main='b. Perimeter (P)',line=0.2,adj=0.01)
axis.Date(side=1,labels=F)
if(F)for(i in 1:9){
  cc = ode(y = state*(1+0.05*rnorm(2)), times = t, func = system.equations, parms = parameters*(1+0.05*rnorm(5)),method = 'euler')
  cc = as.data.frame(cc)
  lines(cc$time,cc$P,col=col.2,lwd=.5)
}
points(data0$time,data0$P.real,pch=19)
lines(cc.center$time,cc.center$P)

plot(cc.center$time,4*pi*cc.center$A/cc.center$P^2,main='',type='n',ylab='',xlab='',log='',xaxt='n',
     ylim=c(0,.4))
axis.Date(side=1,labels=F)
title(main=TeX('c. Health index (H = $A/A_{max}$)',bold=T),line=.6,adj=0.01)

if(F)for(i in 1:9){
  cc = ode(y = state*(1+0.05*rnorm(2)), times = t, func = system.equations, parms = parameters*(1+0.05*rnorm(5)),method = 'euler')
  cc = as.data.frame(cc)
  lines(cc$time,4*pi*cc$A/cc$P^2,col=col.2,lwd=.5)
}
points(data0$time,4*pi*data0$A/data0$P^2,pch=19)
lines(cc.center$time,4*pi*cc.center$A/cc.center$P^2)
lines(cc.center$time,cc.center$P/P0/(2+cc.center$P/P0),col='tomato')

plot(cc.center$time,4*pi*cc.center$A/cc.center$P^2,main='',type='n',ylab='',xlab='',log='',
     ylim=c(0,.4))
title(main=TeX('d. Days of summer withs stress index > 0.5',bold=T),line=.6,adj=0.01)
points(stress.year.0$year, stress.year.0$days_with_extreme_conditions/90,type='o',col='black',pch=19)

dev.off()

pdf('./fig/trade_off.pdf',height = 3,width = 3)
par(mfrow=c(1,1),oma=c(1,1,0,0),mai=c(.6,.6,.1,.1),las=1,mgp=c(2,.5,0),tck=0.05)
kk=seq(0,.035,by=0.001)*10
plot(0.05+kk,100*exp(-9*kk),type='l',xlab='perimeter:area ratio, P/A',ylab=TeX('% coverage, A/A$_{max}$'),ylim=c(0,100))
bertolini = read.csv("./data/Bertolini_2020.txt", header=FALSE, comment.char="#")
lines(supsmu(bertolini$V2*10,bertolini$V1),col='gray',lwd=3)
text(.1,80,'Bertolini et al., 2020',col='gray',font=2,adj=c(0,.5))
text(.1,10,'This study',col='black',font=1,adj=c(0,.5))
dev.off()


month.cols = c(
  "#4E79A7", # Jan - coolest (blue)
  "#6B9CC3", # Feb
  "#92C5DE", # Mar
  "#C7E9B4", # Apr
  "#FEE08B", # May
  "#FDAE61", # Jun
  "#F46D43", # Jul - hottest
  "#D73027", # Aug
  "#FC8D59", # Sep
  "#FDD49E", # Oct
  "#B8E186", # Nov
  "#7FBF7B"  # Dec
)


month.cols = c(
  "#6B9CC3", # Jan - coolest (blue)
  "#6B9CC3", # Feb
  "#C7E9B4", # Mar
  "#C7E9B4", # Apr
  "#C7E9B4", # May
  "#D73027", # Jun
  "#D73027", # Jul - hottest
  "#D73027", # Aug
  "#FDD49E", # Sep
  "#FDD49E", # Oct
  "#FDD49E", # Nov
  "#6B9CC3"  # Dec
)


diff0 = function(...)c(diff(...),0)

the.col = 'lightgray'#alpha(month.cols[round(12*(cc.center$time-floor(cc.center$time))+1)],1)

pdf('./fig/stress_health.pdf',height = 6,width = 6)
par(mfrow=c(2,2),oma=c(1,1,0,0),mai=c(.6,.996,.1,.1),las=1,mgp=c(2.5,.25,0),tck=0.05)
plot(stress.approx(cc.center$time),
     diff0(4*pi*cc.center$A/cc.center$P^2)/(4*pi*cc.center$A/cc.center$P^2),
    type='p',
    xlab='stress index',
    ylab='change in Health index, dH/H',
    col=the.col,pch=19
    )
abline(h=0,lty=3)
#lines(seq(0,1,by=0.01),.3-.6/(1+exp(-25*(seq(0,1,by=0.01)-0.45))))
#abline(s <- lm(diff0( 4*pi*cc.center$A/cc.center$P^2)/(4*pi*cc.center$A/cc.center$P^2)~stress.approx(cc.center$time)))
#summary(s)

plot(stress.approx(cc.center$time),
      diff0(4*pi*cc.center$A/cc.center$P^2),
      type='p',
      xlab='stress index',
      ylab='change in Health index, dH',
      col=the.col,pch=19
)
abline(h=0,lty=3)
#lines(seq(0,1,by=0.01),.3-.6/(1+exp(-25*(seq(0,1,by=0.01)-0.45))))
#abline(s <- lm(diff0( 4*pi*cc.center$A/cc.center$P^2)/(4*pi*cc.center$A/cc.center$P^2)~stress.approx(cc.center$time)))
#summary(s)

plot(stress.approx(cc.center$time),
     diff0(cc.center$P),#/(cc.center$P),
     type='p',
     xlab='stress index',
     ylab='change in perimeter, dP',
     col=the.col,pch=19
)
abline(h=0,lty=3)
lines(seq(0,1,by=0.01),-.3+.6/(1+exp(-205*(seq(0,1,by=0.01)-0.45))))


plot(stress.approx(cc.center$time),
     diff0(cc.center$A),#/(cc.center$A),
     type='p',
     xlab='stress index',
     ylab='change in area, dA',
     col=the.col,pch=19
)
abline(h=0,lty=3)

dev.off()



lines(seq(0,1,by=0.01),.3-.6/(1+exp(-25*(seq(0,1,by=0.01)-0.45))))



