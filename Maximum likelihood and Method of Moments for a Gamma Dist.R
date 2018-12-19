# setwd("C:/Users/dli73/Desktop")




# Maximum Likelihood
data <- read.delim("Neyman.dat", header=FALSE)[,1]
n <- length(data)
xbar <- mean(data)
sumlog <- sum(log(data))

# alpha is 0.4408
alpha <- optimize(function(a){n*a*log(a/xbar)-n*log(gamma(a))+(a-1)*sumlog-n*a},c(0,10),maximum=TRUE)$maximum

# lambda is 1.9643
lambda <- alpha/xbar

# Parametric Bootstrap
alphaStar <- numeric(0)
lambdaStar <- numeric(0)
for( i in 1:1000) {
  sample <- rgamma(227, shape=alpha, rate = lambda)

  n <- length(sample)
  xbar <- mean(sample)
  sumlog <- sum(log(sample))

  alphaB <- optimize(function(a){n*a*log(a/xbar)-n*log(gamma(a))+(a-1)*sumlog-n*a},c(0,10),maximum=TRUE)$maximum
  lambdaB <- alpha/xbar

  alphaStar <- c(alphaStar,alphaB)
  lambdaStar <- c(lambdaStar,lambdaB)
  rm(sample,n,xbar,sumlog,alphaB,lambdaB)
}

SEalpha <- sqrt((1/1000)*sum((alphaStar-mean(alphaStar))^2))
SElambda <- sqrt((1/1000)*sum((lambdaStar-mean(lambdaStar))^2))
alphaCI <- c(sort(alphaStar)[25], sort(alphaStar)[975])
lambdaCI <- c(sort(lambdaStar)[25], sort(lambdaStar)[975])















# Method of Moments
data <- read.delim("Neyman.dat", header=FALSE)[,1]
firstMoment <- (1/length(data))*sum(data)
secondMoment <- (1/length(data))*sum(data^2)

# lambda is 1.6842
lambda <- firstMoment/(secondMoment-firstMoment^2)

# alpha is 0.3779
alpha <- (firstMoment^2)/(secondMoment-firstMoment^2)

# Parametric Bootstrap
alphaStar <- numeric(0)
lambdaStar <- numeric(0)
for( i in 1:1000) {
  sample <- rgamma(227, shape=alpha, rate = lambda)
  firstMomentB <- (1/length(sample))*sum(sample)
  secondMomentB <- (1/length(sample))*sum(sample^2)
  lambdaB <- firstMomentB/(secondMomentB-firstMomentB^2)
  alphaB <- (firstMomentB^2)/(secondMomentB-firstMomentB^2)

  alphaStar <- c(alphaStar,alphaB)
  lambdaStar <- c(lambdaStar,lambdaB)
}

SEalpha <- sqrt((1/1000)*sum((alphaStar-mean(alphaStar))^2))
SElambda <- sqrt((1/1000)*sum((lambdaStar-mean(lambdaStar))^2))
alphaCI1 <- c(sort(alphaStar)[25], sort(alphaStar)[975])
lambdaCI <- c(sort(lambdaStar)[25],sort(lambdaStar)[975])






#Minimum chi square
data <- read.delim("Neyman.dat", header=FALSE)[,1]
n <- length(data)
partition <- c(-999,.1,.5,1,1.5,2,999)
# number of bins
k <- length(partition)-1
# probability of being in a bin
pi <- list()
count <- numeric(k)
#this holds the chi-sq term (obs-exp)^2/exp for each bin
temp <- list()
for (i in 1:k) {
  count[i] <- sum(partition[i]<data & data <= partition[i+1])
  pi[i] <-eval(substitute(
    list(function(a,l){pgamma(upper,shape=a,rate=l)-pgamma(lower,shape=a,rate=l)}),
    list(lower=partition[i],upper=partition[i+1])))

  temp[i] <- eval(substitute(
    list(function(a,l){(ct-n*pi(a,l))^2/(pi(a,l))}),
    list(ct=count[i],pi=pi[[i]],n=n)))
}
chisq <- function(a,l) {
  sum(sapply(temp,function(x){x(a,l)}))
}
estimate<-optim(par=c(.4,2),fn=function(x){chisq(x[1],x[2])})
#alpha is 0.3439
alpha<-estimate$par[1]
#lambda is 1.5068
lambda<-estimate$par[2]
rm(pi,count,temp,i,chisq,estimate)

#Parametric Bootstrap
alphaStar <- numeric(0)
lambdaStar <- numeric(0)
for( i in 1:1000) {
  sample <- rgamma(227, shape=alpha, rate = lambda)
  n <- length(sample)

  # probability of being in a bin
  pi <- list()
  count <- numeric(k)
  #this holds the chi-sq term (obs-exp)^2/exp for each bin
  temp <- list()
  for (i in 1:k) {
    count[i] <- sum(partition[i]<sample & sample <= partition[i+1])
    pi[i] <-eval(substitute(
      list(function(a,l){pgamma(upper,shape=a,rate=l)-pgamma(lower,shape=a,rate=l)}),
      list(lower=partition[i],upper=partition[i+1])))

    temp[i] <- eval(substitute(
      list(function(a,l){(ct-n*pi(a,l))^2/(pi(a,l))}),
      list(ct=count[i],pi=pi[[i]],n=n)))
  }
  chisq <- function(a,l) {
    sum(sapply(temp,function(x){x(a,l)}))
  }
  estimate<-optim(par=c(.4,2),fn=function(x){chisq(x[1],x[2])})
  alphaB<-estimate$par[1]
  lambdaB<-estimate$par[2]
  rm(pi,count,temp,i,chisq,estimate)

  alphaStar <- c(alphaStar,alphaB)
  lambdaStar <- c(lambdaStar,lambdaB)
}

SEalpha <- sqrt((1/1000)*sum((alphaStar-mean(alphaStar))^2))
SElambda <- sqrt((1/1000)*sum((lambdaStar-mean(lambdaStar))^2))
alphaCI <- c(sort(alphaStar)[25], sort(alphaStar)[975])
lambdaCI <- c(sort(lambdaStar)[25], sort(lambdaStar)[975])
