

library(ggplot2)
library(tidyverse)


N = 70

sigma_a_ = 1
lambda_ = 1
K0_ = 10
r = 1

xmin = 0
xmax = 1

# X0
x0_ = xmin + (xmin + xmax)/2

# Les traits
X = seq(xmin, xmax, length.out = N)

X_X = matrix(X, byrow = T, nrow = N, ncol = N )
delta = X_X - t(X_X)



# les fonction de competitions et d environnment

function_a = function(delta, sigma_a = sigma_a_){ 
  exp( (-1/2*(delta)**2)/sigma_a**2 )
}

function_K = function(x, K0 = K0_, lambda=lambda_, x0 = x0_){
  output = max((K0 - lambda*(x-x0)**2), 0) + 1E-9
}

# A

A = sapply(delta, FUN = function_a)

# K

K = sapply(X, FUN = function_K)

# B
B = matrix(rep(1/K, times = N),byrow = T, nrow = N, ncol = N)

# M
M = A*B


# Equation general

UN = matrix(1, ncol = N, nrow = 1)

# Evolution

duree = 4000

n = matrix(rep(0, times = N), nrow= 1, ncol = N) #n0
n[20] = 1


GLV = function(t, n, params){
  with(as.list(params), {
    
    dn = r * n * ( UN - n %*% M )
    
    pheno_maj = which.max(n)
    
    if (runif(1)<0.5){
      dn[pheno_maj+1] = dn[pheno_maj+1] + 0.1 * dn[pheno_maj]
      
    }
    else{
      dn[pheno_maj-1] = dn[pheno_maj+1] + 0.1 * dn[pheno_maj]
      
    }
    
    
    return(list(dn))
  })
}


params = list(r = r, M = M)

dyna_glv = ode(y = n, times = 1:duree, func = GLV, params , method = "euler")



#
dyna_glv_long = as.data.frame(dyna_glv)
str(dyna_glv_long)
dyna_glv_long = dyna_glv_long%>%
  pivot_longer(-time, values_to = "dens", names_to = "sp_ID")

dyna_glv_long$trait = X[ as.numeric(dyna_glv_long$sp_ID)]

ggplot(dyna_glv_long)+
  geom_raster(aes(trait, time,  fill=dens))+
  scale_fill_gradient2(low = "white" ,
                       high = "red")+
  main_theme
#


