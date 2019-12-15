load("yA")
load("yB")

#1st question
hist(yB,breaks = 100)
hist(yA,breaks = 100, freq = FALSE)

#the _1 means that we are considering the first case K=2

N=100
p = m1_1 = m2_1 = s1_1 = s2_1 = rep(NA,N)

#EM initialisation
p[1] = 0.5

m1_1[1] = 1
m2_1[1] = 2

s1_1[1] = 0.1
s2_1[1] = 0.1


y=yA
tau = rep(NA,n)
for (i in 2:N)
{
  
  tau = p[i-1]/
    ((1-p[i-1])*exp(dnorm(y,m1_1[i-1],s1_1[i-1],log=TRUE)-dnorm(y,m2_1[i-1],s2_1[i-1],log=TRUE))+p[i-1])
  
  p[i] = mean(tau)
  
  m1_1[i] = sum((1-tau)*y)/sum(1-tau)
  m2_1[i] = sum((tau)*y)/sum(tau)
  
  s1_1[i] = sqrt(sum((1-tau)*(y-m1_1[i])^2)/sum(1-tau))
  s2_1[i] = sqrt(sum((tau)*(y-m2_1[i])^2)/sum(tau))
}

thetahat_1 = c(p[i],m1_1[i],m2_1[i],s1_1[i],s2_1[i]); thetahat_1

#we can see the convergence of each parameters.
plot(p[1:i],type="l")
plot(1-p[1:i],type="l")

plot(m1_1[1:i],type="l")  
plot(m2_1[1:i],type="l")  

plot(s1_1[1:i],type="l")
plot(s2_1[1:i],type="l")

x = seq(-2,2,length.out=1000)

plot(x,p[i]*dnorm(x,m1_1[i],s1_1[i])+(1-p[i])*dnorm(x,m2_1[i],s2_1[i]))

#2nd question

hist(y,breaks=100,prob=TRUE)
h=0.3
mean(dnorm(0.5,y,h)) 

s = seq(-2,2,length.out=1000) #evaluate density at these points
l = length(s)
n = length(y)
D = matrix(s,n,l,byrow=TRUE)
D = D-y
C = dnorm(D,0,sd=h)

density = colMeans(C)

plot(x,p[i]*dnorm(x,m1_1[i],s1_1[i])+(1-p[i])*dnorm(x,m2_1[i],s2_1[i]))

lines(s,density,type="l",col="red",lwd=3)


#3rd question

#the _2 means that we are considering the variables of the second case k=3
p1_2 = p2_2 = m1_2 = m2_2 = m3_2 = s1_2 = s2_2 = s3_2 = rep(NA,N)

#EM initialisation
#with the hist we can see that the value of the means of each clusters are around -1, 0 and 1

p1_2[1] = 1/3
p2_2[1] = 1/4
p3_2 = 1-p1_2-p2_2

m1_2[1] = -1
m2_2[1] = 0.1
m3_2[1] = 1

s1_2[1] = 0.1
s2_2[1] = 0.1
s3_2[1] = 0.1


N=100
y=yA
for (i in 2:N)
{
  tau1_2 = p1_2[i-1]*dnorm(y,m1_2[i-1],s1_2[i-1])/
    (p1_2[i-1]*dnorm(y,m1_2[i-1],s1_2[i-1])+p2_2[i-1]*dnorm(y,m2_2[i-1],s2_2[i-1])+
       (1-p1_2[i-1]-p2_2[i-1])*dnorm(y,m3_2[i-1],s3_2[i-1]))
  
  tau2_2 = p2_2[i-1]*dnorm(y,m2_2[i-1],s2_2[i-1])/
    (p1_2[i-1]*dnorm(y,m1_2[i-1],s1_2[i-1])+p2_2[i-1]*dnorm(y,m2_2[i-1],s2_2[i-1])+
       (1-p1_2[i-1]-p2_2[i-1])*dnorm(y,m3_2[i-1],s3_2[i-1]))
  
  tau3_2 = 1-tau1_2-tau2_2
  
  p1_2[i] = mean(tau1_2)
  p2_2[i] = mean(tau2_2)
  p3_2[i] = mean(1 - tau1_2 -tau2_2)
  
  m1_2[i] = sum((tau1_2)*y)/sum(tau1_2)
  m2_2[i] = sum((tau2_2)*y)/sum(tau2_2)
  m3_2[i] = sum((tau3_2)*y)/sum(tau3_2)
  
  s1_2[i] = sqrt(sum((tau1_2)*(y-m1_2[i])^2)/sum(tau1_2))
  s2_2[i] = sqrt(sum((tau2_2)*(y-m2_2[i])^2)/sum(tau2_2))
  s3_2[i] = sqrt(sum((tau3_2)*(y-m3_2[i])^2)/sum(tau3_2))
}

#we can see that all the parameters are converging
plot(p1_2[1:i],type="l")
plot(p2_2[1:i],type="l")
plot(p3_2[1:i],type="l")

plot(m1_2[1:i],type="l")  
plot(m2_2[1:i],type="l")  
plot(m3_2[1:i],type="l")

plot(s3_2[1:i],type="l")
plot(s1_2[1:i],type="l")
plot(s2_2[1:i],type="l")

x = seq(-2,2,length.out=1000)
thetahat_2 = c(p1_2[i], p2_2[i], p3_2[i], m1_2[i],m2_2[i],m3_2[i],s1_2[i],s2_2[i],s3_2[i])
thetahat_2

plot(x,p1_2[i]*dnorm(x,m1_2[i],s1_2[i]) + (p2_2[i])*dnorm(x,m2_2[i],s2_2[i])+
       (1-p2_2[i]-p1_2[i])*dnorm(x,m3_2[i],s3_2[i]))



#4th question 

hist(y,breaks=100,prob=TRUE)
h = 0.1
mean(dnorm(0.5,y,h)) 

s = seq(-2,2,length.out=1000) #evaluate density at these points
l = length(s)
n = length(y)
D = matrix(s,n,l,byrow=TRUE)
D = D-y
C = dnorm(D,0,sd=h)
density = colMeans(C)

hist(yA, breaks=100, freq=F)

lines(x,p1_2[i]*dnorm(x,m1_2[i],s1_2[i]) + (p2_2[i])*dnorm(x,m2_2[i],s2_2[i])+
        (1-p2_2[i]-p1_2[i])*dnorm(x,m3_2[i],s3_2[i]),col='blue')

lines(s,density,type="l",col="black",lwd=3)

#perfect fit!!


#5th question

l_1 = sum(log((1-p[i])*dnorm(y,m1_1[i],s1_1[i])+
                p[i]*dnorm(y,m2_1[i],s2_1[i])))

l_2 = sum(log(p1_2[i]*dnorm(y,m1_2[i],s1_2[i])+
                p2_2[i]*dnorm(y,m2_2[i],s2_2[i])+
                (1-p1_2[i]-p2_2[i])*dnorm(y,m3_2[i],s3_2[i])))

P_1 = 6
P_2 = 9

AIC_1=2*P_1-2*l_1
AIC_2 = 2*P_2 -2*l_2

#since AIC_2<AIC_1, we prefer the second model, with k=3.


#6th question

#BOOTSTRAPPING CONFIDENCE INTERVALS (5)

N = 100 #EM iterations
p1 = p2 = m1 = m2 = m3 = s1 = s2 = s3 = rep(NA,N)
p3 = 1-p1-p2
n=length(y)

#EM initialisation
#with the hist we can see that the value of the means of each clusters are around -1, 0 and 1
p1[1] = 1/3
p2[1] = 1/4
p3 = 1-p1-p2

m1[1] = -1
m2[1] = 0.1
m3[1] = 1

s1[1] = 0.1
s2[1] = 0.1
s3[1] = 0.1

tau1 = tau2 = tau3 = rep(NA,n) #these are updated at each iteration

y = yA
B = 100
pstar1 = pstar2 = pstar3 = rep(NA,B)
for (k in 1:B)
{
  y = sample(yA,replace=TRUE)
  
  for (i in 2:N)
  {
    tau1 = p1[i-1]*dnorm(y,m1[i-1],s1[i-1])/
      (p1[i-1]*dnorm(y,m1[i-1],s1[i-1])+p2[i-1]*dnorm(y,m2[i-1],s2[i-1])+
         (1-p1[i-1]-p2[i-1])*dnorm(y,m3[i-1],s3[i-1]))
    
    tau2 = p2[i-1]*dnorm(y,m2[i-1],s2[i-1])/
      (p1[i-1]*dnorm(y,m1[i-1],s1[i-1])+p2[i-1]*dnorm(y,m2[i-1],s2[i-1])+
         (1-p1[i-1]-p2[i-1])*dnorm(y,m3[i-1],s3[i-1]))
    
    tau3 = 1-tau1-tau2
    
    p1[i] = mean(tau1)
    p2[i] = mean(tau2)
    p3[i] = mean(1 - tau1 -tau2)
   
    m1[i] = sum((tau1)*y)/sum(tau1)
    m2[i] = sum((tau2)*y)/sum(tau2)
    m3[i] = sum((tau3)*y)/sum(tau3)
    
    s1[i] = sqrt(sum((tau1)*(y-m1[i])^2)/sum(tau1))
    s2[i] = sqrt(sum((tau2)*(y-m2[i])^2)/sum(tau2))
    s3[i] = sqrt(sum((tau3)*(y-m3[i])^2)/sum(tau3))
  }
  
  pstar1[k] = p1[i]
  pstar2[k] = p2[i]
  pstar3[k] = p3[i]
}

phat1 = thetahat_2[1]
phat2 = thetahat_2[2]
phat3 = thetahat_2[3]

phat4 = thetahat_2[4]
phat5 = thetahat_2[5]
phat6 = thetahat_2[6]

phat7 = thetahat_2[7]
phat8 = thetahat_2[8]
phat9 = thetahat_2[9]

#normal
phat1+c(-1,1)*1.96*sd(pstar1)
phat2+c(-1,1)*1.96*sd(pstar2)
phat3+c(-1,1)*1.96*sd(pstar3)

phat4+c(-1,1)*1.96*sd(pstar4)
phat5+c(-1,1)*1.96*sd(pstar5)
phat6+c(-1,1)*1.96*sd(pstar6)

phat7+c(-1,1)*1.96*sd(pstar7)
phat8+c(-1,1)*1.96*sd(pstar8)
phat9+c(-1,1)*1.96*sd(pstar9)

#percentiles
quantile(pstar1,c(0.025,0.975))
quantile(pstar2,c(0.025,0.975))
quantile(pstar3,c(0.025,0.975))

quantile(pstar4,c(0.025,0.975))
quantile(pstar5,c(0.025,0.975))
quantile(pstar6,c(0.025,0.975))

quantile(pstar7,c(0.025,0.975))
quantile(pstar8,c(0.025,0.975))
quantile(pstar9,c(0.025,0.975))



#OPTIONAL

Mixture_Gaussian=function(K){
  
  N = 100
  y = yB
  n = length(y)
  p = matrix(NA, nrow = 100, ncol = K)
  m = matrix(NA, nrow = 100, ncol = K) 
  s = matrix (NA , nrow = 100 , ncol = K)
  tau = matrix(NA, nrow = n, ncol = K)
  lk = 0
  
  #Initialization
  for (h in 1:K){
    
    p[1,h] = 1/K
    m[1,h] = min(y)+ h * (max(y)-min(y))/(K+2)  #we initialize every mean such that they cover the whole interval
    s[1,h] = 0.1  
    lk = lk + p[1,i]*dnorm(y,m[1,i],s[1,i])    #the genius comes from the K+2 since its very unlikely that one of the targeted means will be at the upperbound of yA
  }
  
  for (i in  2:N){
    
    tot_sum = rep(0,1000)
    
    for (k in 1:K){
      tot_sum = tot_sum + p[i-1,k]*dnorm(y,m[i-1,k],s[i-1,k])
    }
    
    for (j in  1:K){
      
      tau[,j] = p[i-1,j]*dnorm(y,m[i-1,j],s[i-1,j])/tot_sum
      p[i,j] = mean(tau[,j])
      m[i,j] = sum(tau[,j]*y)/sum(tau[,j])
      s[i,j] = sqrt(sum((tau[,j])*(y-m[i,j])^2)/sum(tau[,j]))
    }
    likelihood = rep(0,1000)
    for (j in seq(K)){
      likelihood = likelihood + p[i,seq(K)[j]]*dnorm(y,m[i,seq(K)[j]],s[i,seq(K)[j]])
    }
    lk = c(lk, sum(log(likelihood)))
  } 

  thetahat = c(p[N,], m[N,], s[N,])
  
  h = c()
  h[["thetahat"]] = thetahat
  h[['likelihood']] = lk
  
  return(h)
}


#Akaike

Akaike=function(values){
  num_par = length(values)-1000
  return(2*length(values[["thetahat"]]-1) - 2*sum(log(values[["likelihood"]])))
}

l=c()

for (k in 1:7){
  values=Mixture_Gaussian(k)
  l = c(l, Akaike(values))
}

plot(l,col='blue',type='l') #we clearly see that 6 is the optimal number of clusters



#plotting

color=c('blue','pink','black','green','brown','yellow','red')

hist(yB, freq=F, breaks=100)
s = seq(-2,2,length.out=1000)

for (k in 1:7){
  values = Mixture_Gaussian(k)[["thetahat"]]
  tot_sum = 0
  
  for (j in 2:k){
    tot_sum = tot_sum + values[j]*dnorm(x,values[j+k],values[j+2*k])
  }
  lines(s,tot_sum, col=color[k])
}

#graphically k=6

Mixture_Gaussian(k)[["thetahat"]]



