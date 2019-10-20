#### Stats Assignment No. 2
# Team Members:
#   Leonie Intat
#   Hugo Paolini
#   Konrad Eilers

#EXERCISE 1
#Question1
# a)
## (i)
c=4/3
# pdf
f1<-function(x){
  ind.zero_to_one = (0 < x & x <= 1)
  ind.other = (x > 1 & x < 0)
  u = c()
  u[ind.zero_to_one] = 4/3 * (x+x**3)
  u[ind.other] = 0
  return(u)
}

## (ii)
# cdf
cdf<-function(x){
  ind.zero_to_one = (0 < x & x <= 1)
  ind.other = (x > 1 & x < 0)
  u=c()
  u[ind.zero_to_one]=4/3*(1/2*x^2 + 1/4*x^4)
  u[ind.other]=0
  return (u)
  }

## (iii) / b (i)
invcdf<-function(x){
  ind.zero_to_one = (0 < x & x <= 1)
  ind.other = (x > 1 & x < 0)
  u = c()
  u[ind.zero_to_one] = sqrt(-1+sqrt(1+3*x))
  u[ind.other] = 0
  return(u)
}

## b (ii)
u=runif(10000)
#using transformation method
x=invcdf(u)

## b (iii)
hist(x,breaks=100,freq=F)
xx=seq(0,1,0.01)
yy=f1(xx)
lines(xx,yy, col='blue')  #looks good

## b (iv)
var(x)
mean(x)

#Question2
# (a) (i) (ii)
f2<-function(x){
  return((sin(pi*x)**2)*(x+x**3))
}
yyy=f2(xx)

plot(xx,yyy,col='red',type='l')

# A normal distribution fits well as a rejection sampler
lines(xx, dnorm(xx, mean = 0.7, sd = 0.2), type = "l")

# b) Rejection Sampling (i-ii)
# Proposal function: Normal Distribution with mean = 0.7 and sd = 0.2
max(yyy)
# M = 1
RSSample <- function(n){
  M = 1
  x = rnorm(n,mean = 0.7, sd = 0.2)
  u = runif(n)
  x = x[u<f2(x)/(dnorm(x,0.7,0.2)*M)]
  return(x)
}

sample = RSSample(50000)
hist(sample[sample<1],breaks=100,freq=F)
lines(xx,yyy,col='red')  #looks good as well

# (iii)
mean(sample)
var(sample)

# (iv)
#lets check the value of g through integration
integrate(f2,0,1)
g=0.337

f2_normalised<-function(x){
  return(1/g *(sin(pi*x)**2)*(x+x**3))
}

yyyy=f2_normalised(xx)
hist(sample[sample<1],breaks=100,freq=F)
lines(xx,yyyy,col='red')  

# Estimation of c
accepted = length(sample) # Accepted
rejected = 50000 - length(sample) # Rejected
p = accepted / (accepted + rejected)
# P(accepted) = 1/M
# M = c*M'
# M' = 1
# c = 1/P(accepted)
c = 1/p
1/g # We observe c to be close to the "true" value

#EXERCISE 2
#a)
dualf<-function(x,y){
  return((sin(x+y)**2 * cos(x**2)**2) * exp(-1/2 * (abs(x)+ abs(y))))
}

xx1<-seq(-2,2,0.04)
xx2<-seq(-2,2,0.04)
z=outer(xx1,xx2,dualf)

#b)
persp(xx1,xx2,z,theta=10,phi=20, expand=0.5,col='lightblue')
filled.contour(xx1,xx2,z)

#c)
#we sample from two laplace (0,2) because we notice thath the function exp(...) is the product
#of two laplace(0,2) so the two variables are independant.

#since we obtain a 1/4 coefficient by sampling from a laplace(0,2) we set M=4
#we can now proceed with the rejection sampling

#install.packages('rmutil')
plot.new()
library('rmutil')
RSSsample_Laplace <- function(n) {
  samplex_f = c()
  sampley_f = c()
  count = 0
  for (i in 1:10000){
    u = 10000
    xp = 0
    yp = 0
    while (u>(sin(xp+yp)**2 * cos(xp)**2))
    {
      u = runif(1)
      count = count + 1
      xp=rlaplace(1,0,2)
      yp=rlaplace(1,0,2)
    }
    samplex_f[i] = xp
    sampley_f[i] = yp
  }
  data = data.frame(samplex_f,sampley_f)
  return(data)
}

data = RSSsample_Laplace(n=1000)

contour(xx1,xx2,z)
points(data[,1],data[,2],col='red') #looks good

# Expected Value
mean(data[,1]) 
mean(data[,2])
# Variance
var(data[,1])
var(data[,2])

#d) Importance Sampling
ISsample_Laplace <- function(n)
{
  xp = rlaplace(n,m = 0, s = 2)
  yp = rlaplace(n,m = 0, s = 2)
  val = list(xp, yp)
  w = dualf(xp,yp) / (dlaplace(xp,0,2)*dlaplace(yp,0,2))
  return(data.frame(val,w))
}
sampleIS = ISsample_Laplace(10000)
# mean
e_x = sum(sampleIS[,1]*sampleIS[,3]) / sum(sampleIS[,3])
e_y = sum(sampleIS[,2]*sampleIS[,3]) / sum(sampleIS[,3])
# var
var_x = mean((sampleIS[,1]-e_x)^2)
var_y = mean((sampleIS[,2]-e_y)^2)

# e)
# Rejection Sampling
g_proposed <-function(x_1, x_2){
  (1/(2*pi))*exp(-1/2*(x_1 + x_2))
}
x = seq(-10,2,0.1)
y = seq(-10,2,0.1)
z_f =dualf(x,y)
z_g_prop = g_proposed(x,y)
sum(z_g_prop < z_f)
# The proposed function is below  f. 
# This violates the conditions of RS. Rejection sampling does not work.
# While this could theoretically be fixed through M, we can clearly observe
# g function does not behave well, and converges towards infinity where x,y -> infinity
# This is not useful for g.
z=outer(x,y,g_proposed)
filled.contour(x,y,z)

# Importance Sampling
ISsample_Laplace_Gauss <- function(n)
{
  xp = rnorm(n,mean = 0, sd = 1)
  yp = rnorm(n,mean = 0, sd = 1)
  val = list(xp, yp)
  w = dualf(xp,yp) / g_proposed(xp,yp)
  return(data.frame(val,w))
}
sampleIS_gauss = ISsample_Laplace_Gauss(10000)
# mean
e_x = sum(sampleIS_gauss[,1]*sampleIS_gauss[,3]) / sum(sampleIS_gauss[,3])
e_y = sum(sampleIS_gauss[,2]*sampleIS_gauss[,3]) / sum(sampleIS_gauss[,3])
# var
var_x = mean((sampleIS_gauss[,1]-e_x)^2)
var_y = mean((sampleIS_gauss[,2]-e_y)^2)

#f) 
j = 100
n = 20
meanRS_x = c() 
meanRS_y = c()
RS_var_x = c()
RS_var_y = c()

meanIS_LA_x = c()
meanIS_LA_y = c()
IS_var_LA_x = c()
IS_var_LA_y = c()

meanIS_GA_x = c()
meanIS_GA_y = c()
IS_var_GA_x = c()
IS_var_GA_y = c()

for (k in 1:n)
{
  sampleRS = RSSsample_Laplace(j)
  sampleIS = ISsample_Laplace(j)
  sampleIS_gauss = ISsample_Laplace_Gauss(j)
  
  meanRS_x[k] = mean(sampleRS[,1])
  meanRS_y[k] = mean(sampleRS[,2])
  meanIS_LA_x[k] = sum(sampleIS[,1]*sampleIS[,3]) / sum(sampleIS[,3])
  meanIS_LA_y[k] = sum(sampleIS[,2]*sampleIS[,3]) / sum(sampleIS[,3])
  meanIS_GA_x[k] = sum(sampleIS_gauss[,1]*sampleIS_gauss[,3]) / sum(sampleIS_gauss[,3])
  meanIS_GA_y[k] = sum(sampleIS_gauss[,2]*sampleIS_gauss[,3]) / sum(sampleIS_gauss[,3])
  
  RS_var_x[k] = var(sampleRS[,1])
  RS_var_y[k] = var(sampleRS[,2])
  IS_var_LA_x[k] = mean((sampleIS[,1]-meanIS_LA_x[k])^2)
  IS_var_LA_y[k] = mean((sampleIS[,2]-meanIS_LA_y[k])^2)
  IS_var_GA_x[k] = mean((sampleIS_gauss[,1]-meanIS_GA_x[k])^2)
  IS_var_GA_y[k] = mean((sampleIS_gauss[,2]-meanIS_GA_y[k])^2)

}

boxplot(meanRS_x,meanIS_LA_x, meanIS_GA_x, names=c("RS IIc","IS E. IId","IS (Gauss) E. IIe"))
boxplot(meanRS_y,meanIS_LA_y, meanIS_GA_y, names=c("RS IIc","IS E. IId","IS (Gauss) E. IIe"))
boxplot(RS_var_x,IS_var_LA_x, IS_var_GA_x, names=c("RS IIc","IS E. IId","IS (Gauss) E. IIe"))
boxplot(RS_var_y,IS_var_LA_y, IS_var_GA_y, names=c("RS IIc","IS E. IId","IS (Gauss) E. IIe"))

# We use boxplots to compare the relative performance. We can clearly see that the 
# rejection sampling method is more efficient at estimating the mean & variance

# Further, we can observe that the Importance Sampling using Gauss leads to biased estimates for the 
# variance and expected value. 



