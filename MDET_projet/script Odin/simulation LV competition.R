

library(ggplot2)
library(tidyverse)


N = 10

sigma_a_ = 2
lambda_ = 1
K0_ = 1
r = 1

xmin = 0
xmax = 1

x0_ = xmin + (xmin + xmax)/2

X = seq(xmin, xmax, length.out = N)

a = function(x1,x2, sigma_a = sigma_a_){ 
  exp( (-1/2*(x2-x1)**2)/sigma_a**2 )
}

fonction_K = function(x, K0 = K0_, lambda=lambda_, x0 = x0_){
  max((K0 - lambda*(x-x0)**2), 0)
}

# A

A = matrix(0, nrow = N, ncol = N)

for (i in 1:N){
  for (j in 1:N){
    A[i,j] = a(X[i],X[j])
  }
}

# K

K = sapply(X, FUN = fonction_K)

B = matrix(rep(1/K, times = N),byrow = T, nrow = N, ncol = N)

# M
M = A*B

#



# Equation generale

# N_0

n = rep(K0_/N, times = N)

UN = matrix(1, ncol = N, nrow = 1)

dn_dt = r * n * ( UN - n %*% M )

# Evolution

duree = 300

n_histo = matrix(n, byrow = T,ncol = 1, nrow = N )

for (t in 1:duree){
  
  n_ = matrix(n_histo[, ncol(n_histo)], byrow = T,ncol = 1, nrow = N ) 
  
  dn_dt = r * t(n_) * ( UN - t(n_) %*% M )
  
  n_ = n_ + t(dn_dt)
  
  n_histo = cbind(n_histo, n_)
    
    
}



p = ggplot()
for (i in 1:N){
  data_i = tibble(temp = 1:(duree+1), 
                  n_i = n_histo[i,])
  p = p + geom_line(data = data_i, aes(x =  temp, y = n_i))
}
p

p2 = ggplot()
for (i in 1:N){
  data_i = tibble(temp = 1:(duree+1), 
                  x_i = rep(X[i], times = (duree+1)))
  p2 = p2 + geom_line(data = data_i, aes(x = x_i , y = temp))
}
p2



