cap_charge <- function(x){
  K0=1
  lamb=1
  x0=1
  K<-K0-lamb*(x-x0)**2
  #return(K)
}

aa <- function(x1,x2){
  sigma=1
  a<- exp(1/2*(x1-x2)^2/(sigma^2))
  #return (a)
}

ss <-function(x1,x2){
  r=1
  s <- r*(1-aa(x1,x2)*(cap_charge(x1))/(cap_charge(x2)))
  return(s)
}

ss(x1,x2)
library(rgl)

x1=seq(0.1,1.9,0.01)
x2=seq(0.1,1.9,0.01)
length(x1)

recup= rep(NA,19*19)
n=1
X1_long <-rep(0,length(x2)**2)
X2_long <-rep(0,length(x1)**2)

for(k in x1){
  for (i in x2){
    X1_long[n] <- k
    X2_long[n] <- i
    
    ss(k,i)
    recup[n]=ss(k,i)
    n=n+1
    
  }
}
X2_long
length(recup)
length(X1_long)
plot3d(X1_long,X2_long,recup)



###PACKAGE RGL POUR LES PLOT 3D
plot3d(ss,xlim=c(0.1,1.9), ylim=c(0.1,1.9))
