# math 189



setwd("C://Users//haw11//Desktop")
data = read.table("gauge-1wb1wa6-2gpel41.txt",header=T)

data.gain = as.vector(data$gain)
data.density = as.vector(data$density)

########   Fitting   #######

fit <- lm(density ~ gain, data = data);fit

summary(fit)

{plot(data$gain,data$density,xlab="gain",ylab="density")
  abline(fit, col="red",lwd=2)
}
#~~~ check for residuals
{plot(residuals(fit))
  abline(0,0, col="red")  # problem is the residuals are not randomly scattered in norma distribution
}
mean(residuals(fit))

hist(residuals(fit),breaks = 20, prob = T,xlim=c(-0.15,0.20))

{qqnorm(residuals(fit))
  qqline(residuals(fit), col="red",lwd=2)
}
#~~~ not exact measure of density blocks     9 ednsity blocks in total



###############
########density on x axis, gain on y axis
#
#
fit.2 <- lm(gain~density , data = data);fit.2
summary(fit.2)
plot(data)
abline(fit.2, col="red",lwd=2)
#~~~ check for residuals
plot(residuals(fit.2))
abline(0,0, col="red")  # problem is the residuals are not randomly scattered in norma distribution

mean(residuals(fit.2))
residuals(fit.2)
hist(residuals(fit.2),breaks = 20, prob = T,xlim=c(-100,150))

qqnorm(residuals(fit.2))
qqline(residuals(fit.2), col="red",lwd=2)


#~~~ not random measure of density blocks
# time might affect the final result of the slope of the regression
# if gain keeps increasing in time, doing measurement in order will bias the 
# measurement of the gain in timely order



#################################################
###################################
###### density on x axis, log(gain) on y axis
###################################
#################################################

newdata<-data
newdata[,2]<-log(newdata[,2]);newdata
fit.3 <- lm(log(gain)~density, data = data);fit.3
summary(fit.3)

plot(data$density,log(data$gain),xlab="density",ylab="log(gain)")
abline(fit.3, col="red",lwd=2)

#~~~ check for residuals
plot(residuals(fit.3))
abline(0,0, col="red",lwd=2) 
mean(residuals(fit.3))

hist(residuals(fit.3),breaks = 20, prob = T,xlab="residuals")
qqnorm(residuals(fit.3))
qqline(residuals(fit.3), col="red",lwd=2) # result is more normal 

#~~~ not exact measure of density blocks     9 ednsity blocks in total
# simulation of write noise

mu1 = 0
sd1 = 0.1

mu2 =0
sd2 =0.03

mu3 = 0
sd3=0.0003
data.density
fit.3
# doing simulation of errors
summary(fit.3)
summary(fit.3)$coefficient[2,1:2]      # Estimate         Std. Error 
# -4.60593690       0.03181727  
CI.2.slope<- c(-4.60593690-1.96*0.03181727 ,-4.60593690+1.96*0.03181727 );CI.2.slope # [-4.668299 -4.543575]

result.2 = rep(0,1000)
lowerbound<-rep(NA,1000)
upperbound<-rep(NA,1000)
for (i in 1 :1000){
  error.2 = c(rnorm(7,mean = mu1, sd = sd1),rnorm(1,mean=mu2,sd=sd2),rnorm(1,mean=mu3,sd=sd3))
  density.2.error = data.density +rep(error.2,each=10)
  
  fit.3.error = lm(log(data.gain)~density.2.error)
  lowerbound[i]<-summary(fit.3.error)$coefficient[2,1]-1.96*summary(fit.3.error)$coefficient[2,2]
  upperbound[i]<-summary(fit.3.error)$coefficient[2,1]+1.96*summary(fit.3.error)$coefficient[2,2]
}
lowerbound;upperbound

for (i in 1:1000)
{
  if(summary(fit.3)$coefficient[2,1]>=lowerbound[i] & summary(fit.3)$coefficient[2,1]<=upperbound[i] )
    result.2[i]<-1
}

mean(result.2)
# 0.28  of the simulation still with in the CI from the original data



summary(fit.3)
fit.3
add3<-predict(fit.3,interval="prediction");add3
add4<-predict(fit.3,interval="confidence");add4

PI = function(a,b,X){
  Xbar = mean(data$density)
  density = as.factor(data$density)
  gain = as.factor(data$gain)
  m = length(levels(density))
  Sisqu = numeric(m)
  j = 1
  for(i in levels(density)){
    Sisqu[j] = var(log(as.numeric(as.character(gain[density == i]))))
    j = j+1
  }
  psigmasqu = mean(Sisqu)
  sum = 0
  for(i in levels(density)){
    sum = sum + (as.numeric(as.character(i))-mean(as.numeric(as.character(density))))^2
  }
  L = (a+b*X)-1.96*sqrt(psigmasqu*(1+1/m+(X-Xbar)^2/sum))
  R = (a+b*X)+1.96*sqrt(psigmasqu*(1+1/m+(X-Xbar)^2/sum))
  return(c(L,R))
}

CI = function(a,b,X){
  Xbar = mean(data$density)
  density = as.factor(data$density)
  gain = as.factor(data$gain)
  m = length(levels(density))
  Sisqu = numeric(m)
  j = 1
  for(i in levels(density)){
    Sisqu[j] = var(log(as.numeric(as.character(gain[density == i]))))
    j = j+1
  }
  psigmasqu = mean(Sisqu)
  sum = 0
  for(i in levels(density)){
    sum = sum + (as.numeric(as.character(i))-mean(as.numeric(as.character(density))))^2
  }
  t = qt(0.975, df = 7)
  L = (a+b*X)-t*sqrt(psigmasqu*(1/m+(X-Xbar)^2/sum))
  R = (a+b*X)+t*sqrt(psigmasqu*(1/m+(X-Xbar)^2/sum))
  return(c(L,R))
}

fit.3                          #predicted: log(gain) = 5.997 -4.606*density
pred1<-matrix(rep(NA,160),ncol=2);pred1
conf1<-matrix(rep(NA,160),ncol=2);conf1
for (i in 1:80)
{pred1[i,]<-PI(5.997,-4.606,-0.1+0.01*i)
conf1[i,]<-CI(5.997,-4.606,-0.1+0.01*i)
}
pred1
conf1

{plot(newdata,xlab="density",ylab="log(gain)",pch=4,col="purple")
abline(fit.3,col="darkgreen")
lines(y=pred1[,1],x=seq(-0.09,0.70,0.01),col="red",lwd=0.5)
lines(y=pred1[,2],x=seq(-0.09,0.70,0.01),col="red",lwd=0.5)
lines(y=conf1[,1],x=seq(-0.09,0.70,0.01),col="orange",lwd=0.5)
lines(y=conf1[,2],x=seq(-0.09,0.70,0.01),col="orange",lwd=0.5)
legend(0.43,5.9, legend=c("fitted","predicted interval", "confidence interval"),
       col=c("darkgreen","red","orange"), lwd =1, cex=0.6)
}


fit.3                          #predicted: log(gain) = 5.997 -4.606*density
PI(5.997,-4.606,0.492)
PI(5.997,-4.606,0.525)
# PI lower [0.492291,0.525586]

# when gain=38.6

 plot(newdata,xlab="density",ylab="log(gain)",pch=4,col="purple",xlim=c(-0.05,0.75),ylim=c(0,6.4))
  abline(fit.3,col="darkgreen")
  lines(y=pred1[,1],x=seq(-0.09,0.70,0.01),col="red",lwd=0.5)
  lines(y=pred1[,2],x=seq(-0.09,0.70,0.01),col="red",lwd=0.5)
  lines(y=conf1[,1],x=seq(-0.09,0.70,0.01),col="orange",lwd=0.5)
  lines(y=conf1[,2],x=seq(-0.09,0.70,0.01),col="orange",lwd=0.5)
  legend(0.43,2.5, legend=c("fitted","predicted interval", "confidence interval"),
         col=c("darkgreen","red","orange"), lwd =1, cex=0.6)
  abline(log(38.6),0)
  abline(log(426.7),0)
  points(x = 0.492291, y = log(38.6), pch = 20, cex = 0.5)
  points(x = 0.525586, y = log(38.6), pch = 20, cex = 0.5)
  
  fPIl = function(x){as.numeric(PI(5.997, -4.606, x)[1])-log(38.6)}
  fPIu = function(x){as.numeric(PI(5.997, -4.606, x)[2])-log(38.6)}
  uniroot(fPIl, c(0,1)) #0.4922937
  uniroot(fPIu, c(0,1)) #0.5255889 PI(38.6) [0.4922937,0.5255889]
  
  fCIl = function(x){as.numeric(CI(5.997, -4.606, x)[1])-log(38.6)}
  uniroot(fCIl,c(0,1))  #0.5011256
  
  fCI2 = function(x){as.numeric(CI(5.997, -4.606, x)[2])-log(38.6)}
  uniroot(fCI2,c(0,1))  #0.5168302  CI(38.6) [0.5011256,0.5168302]
  
  points(x = 0.5011256, y = log(38.6), pch = 17, cex = 0.3,col="brown")
  points(x = 0.5168302, y = log(38.6), pch = 17, cex = 0.3,col="brown")
  points(x = 0.5088467, y = log(38.6), pch = 15, cex = 0.3,col="blue")
  
  fPIl2 = function(x){as.numeric(PI(5.997, -4.606, x)[1])-log(426.7)}
  fPIu2 = function(x){as.numeric(PI(5.997, -4.606, x)[2])-log(426.7)}
  uniroot(fPIl2, c(-0.5,1)) #-0.0309534
  uniroot(fPIu2, c(-0.5,1)) #0.004963428    PI(426.7) [-0.0309534,0.004963428]
  
  fCIl2 = function(x){as.numeric(CI(5.997, -4.606, x)[1])-log(426.7)}
  fCIu2 = function(x){as.numeric(CI(5.997, -4.606, x)[2])-log(426.7)}
  uniroot(fCIl2, c(-0.5,1)) #-0.02434451
  uniroot(fCIu2, c(-0.5,1)) #-0.001797111    CI(426.7) [-0.02434451,-0.001797111]
  
  
  points(x = -0.0309534, y = log(426.7), pch = 20, cex = 0.5)
  points(x = 0.004963428, y = log(426.7), pch = 20, cex = 0.5)
  points(x = -0.02434451, y = log(426.7), pch = 17, cex = 0.3,col="brown")
  points(x = -0.001797111, y = log(426.7), pch = 17, cex = 0.3,col="brown")
  points(x = -0.01282701, y = log(426.7), pch = 15, cex = 0.3,col="blue")
  


######   Cross-Validation    ######

###################
#~~~~ for estimation corresponding to the block of density 0.508
  PInew = function(a,b,X,data){
    Xbar = mean(data$density)
    density = as.factor(data$density)
    gain = as.factor(data$gain)
    m = length(levels(density))
    Sisqu = numeric(m)
    j = 1
    for(i in levels(density)){
      Sisqu[j] = var(log(as.numeric(as.character(gain[density == i]))))
      j = j+1
    }
    psigmasqu = mean(Sisqu)
    sum = 0
    for(i in levels(density)){
      sum = sum + (as.numeric(as.character(i))-mean(as.numeric(as.character(density))))^2
    }
    L = (a+b*X)-1.96*sqrt(psigmasqu*(1+1/m+(X-Xbar)^2/sum))
    R = (a+b*X)+1.96*sqrt(psigmasqu*(1+1/m+(X-Xbar)^2/sum))
    return(c(L,R))
  }
  
  CInew = function(a,b,X,data){
    Xbar = mean(data$density)
    density = as.factor(data$density)
    gain = as.factor(data$gain)
    m = length(levels(density))
    Sisqu = numeric(m)
    j = 1
    for(i in levels(density)){
      Sisqu[j] = var(log(as.numeric(as.character(gain[density == i]))))
      j = j+1
    }
    psigmasqu = mean(Sisqu)
    sum = 0
    for(i in levels(density)){
      sum = sum + (as.numeric(as.character(i))-mean(as.numeric(as.character(density))))^2
    }
    t = qt(0.975, df = 6)
    L = (a+b*X)-t*sqrt(psigmasqu*(1/m+(X-Xbar)^2/sum))
    R = (a+b*X)+t*sqrt(psigmasqu*(1/m+(X-Xbar)^2/sum))
    return(c(L,R))
  }

  
newdata2<-data
newdata2[,2]<-log(newdata2[,2]);dim(newdata2)
ind.1 = which(newdata2$density == 0.508);ind.1
mean(newdata2$gain[ind.1])     # 3.651552

fit.3
newdata2<-newdata2[-ind.1,];dim(newdata2)
dim(newdata)
prediction.1 = lm(gain~density, data = newdata2);prediction.1
summary(prediction.1)

pred2<-matrix(rep(NA,160),ncol=2);pred2
conf2<-matrix(rep(NA,160),ncol=2);conf2
for (i in 1:80)
{pred2[i,]<-PInew(5.997,-4.606,-0.1+0.01*i,newdata2)
conf2[i,]<-CInew(5.997,-4.606,-0.1+0.01*i,newdata2)
}
pred2
conf2
 prediction.1
 fit.3
 plot(newdata2,xlab="density",ylab="log(gain)",pch=4,col="purple",main="without density=0.508")
  abline(prediction.1, col="darkgreen",lwd=1)
  lines(y=pred2[,1],x=seq(-0.09,0.70,0.01),col="red",lwd=0.5)
  lines(y=pred2[,2],x=seq(-0.09,0.70,0.01),col="red",lwd=0.5)
  lines(y=conf2[,1],x=seq(-0.09,0.70,0.01),col="orange",lwd=0.5)
  lines(y=conf2[,2],x=seq(-0.09,0.70,0.01),col="orange",lwd=0.5)
  legend(0.4,5.7, legend=c("fitted","prediction interval", "confidence interval"),
         col=c("darkgreen","red","blue"), lwd =1, cex=0.7)
  abline(log(38.6),0)
  fPIl3 = function(x){as.numeric(PInew(5.997, -4.603, x,newdata2)[1])-log(38.6)}
  fPIu3 = function(x){as.numeric(PInew(5.997, -4.603, x,newdata2)[2])-log(38.6)}
  uniroot(fPIl3, c(0,1)) #0.5041307
  uniroot(fPIu3, c(0,1)) #0.5142443 PI(38.6) [0.5041307,0.5142443]
  
  fCIl3 = function(x){as.numeric(CInew(5.997, -4.606, x,newdata2)[1])-log(38.6)}
  uniroot(fCIl3,c(0,1))  #0.5061798
  
  fCI23 = function(x){as.numeric(CInew(5.997, -4.606, x,newdata2)[2])-log(38.6)}
  uniroot(fCI23,c(0,1))  #0.5115461  CI(38.6) [0.5061798,0.5115461]
  
  points(x = 0.503803, y = log(38.6), pch = 20, cex = 0.5)
  points(x = 0.5139087, y = log(38.6), pch = 20, cex = 0.5)
  points(x = 0.5061798, y = log(38.6), pch = 17, cex = 0.3,col="brown")
  points(x = 0.5115461, y = log(38.6), pch = 17, cex = 0.3,col="brown")
  points(x = 0.5091783, y = log(38.6), pch = 15, cex = 0.3,col="green")
  
  text(x=0.07,y=log(50),labels="given gain=38.6",cex=0.6)

prediction.2                     #predicted: log(gain) = 5.997 -4.603*density

# when gain=38.6, 
## point estimate 0.5091783



##############################
#~~~~for estimation corresponding to the block of density 0.001
newdata3<-newdata
ind.2 = which(newdata3$density == 0.001);ind.2
mean(newdata3[ind.2,2])     #6.056036

newdata3<-newdata3[-ind.2,]

prediction.1
prediction.2 = lm(gain~density, data = newdata3);prediction.2

summary(prediction.2)
fPIl4 = function(x){as.numeric(PInew(5.963, -4.535  , x,newdata3)[1])-log(426.7)}
fPIu4 = function(x){as.numeric(PInew(5.963, -4.535, x,newdata3)[2])-log(426.7)}
uniroot(fPIl4, c(-0.5,1)) #-0.02661253
uniroot(fPIu4, c(-0.5,1)) #-0.01449084 PI(426.7) [-0.02661253,-0.01449084]

fCIl4 = function(x){as.numeric(CInew(5.963, -4.535, x,newdata3)[1])-log(426.7)}
fCIu4 = function(x){as.numeric(CInew(5.963, -4.535, x,newdata3)[2])-log(426.7)}
uniroot(fCIl4, c(-0.5,1)) #-0.02518354
uniroot(fCIu4, c(-0.5,1)) #-0.01596145 CI(426.7) [-0.02518354,-0.01596145]

plot(newdata3$density,newdata3$gain,xlim = c(-0.05,0.75),xlab="density",ylim=c(0,6.3),
      ylab="log(gain)",pch=4,col="blue",main="no density=0.001")
  abline(prediction.2, col="darkgreen",lwd=1)
  
  legend(0.51,6.3, legend=c("fitted","predicted interval", "confidence interval"),
         col=c("darkgreen","red","orange"), lwd =1, cex=0.7)
  abline(log(426.7),0)
  pred3<-matrix(rep(NA,160),ncol=2);pred3
  conf3<-matrix(rep(NA,160),ncol=2);conf3
  for (i in 1:80)
  {pred3[i,]<-PInew(5.963,-4.535,-0.1+0.01*i,newdata3)
  conf3[i,]<-CInew(5.963,-4.535,-0.1+0.01*i,newdata3)
  }
  pred3
  conf3
  lines(y=pred3[,1],x=seq(-0.09,0.70,0.01),col="red",lwd=0.5)
  lines(y=pred3[,2],x=seq(-0.09,0.70,0.01),col="red",lwd=0.5)
  lines(y=conf3[,1],x=seq(-0.09,0.70,0.01),col="orange",lwd=0.5)
  lines(y=conf3[,2],x=seq(-0.09,0.70,0.01),col="orange",lwd=0.5)
  
  points(x = -0.02052507, y = log(426.7), pch = 18, cex = 0.5,col="green")
  points(x = -0.02661253, y = log(426.7), pch = 18, cex = 0.5)
  points(x = -0.02518354, y = log(426.7), pch = 17, cex = 0.3,col="purple")
  points(x = -0.01596145, y = log(426.7), pch = 17, cex = 0.3,col="purple")
  points(x = -0.01449084, y = log(426.7), pch = 18, cex = 0.3)
  text(x=0.12,y=log(400),labels="given gain=426.7",cex=0.3)
  




prediction.2                      #predicted:  log(gain) = 5.963 -4.535*density


