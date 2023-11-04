

library(ggplot2)
library(tidyverse)


N = 30

sigma_a_ = 0.5
lambda_ = 1
K0_ = 1
r = 1

xmin = 0
xmax = 1

x0_ = xmin + (xmin + xmax)/2

X = seq(xmin, xmax, length.out = N)

X_X = matrix(X, byrow = T, nrow = N, ncol = N )

delta = X_X - t(X_X)

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

duree = 300

n = rep(K0_/N, times = N) #n0
n_histo = matrix(n, byrow = T,ncol = N, nrow = 1 )

for (t in 1:duree){
  
  dn_dt = r * n * ( UN - n %*% M )
  
  n = n + dn_dt
  
  n_histo = rbind(n_histo, n)
  
}


n_histo_ =  data.frame(n_histo)
colnames(n_histo_) = paste0("Species",1:N)

n_histo_2 <- n_histo_ %>%
  tibble() %>% 
  gather(key = "variable", value = "value") %>% 
  mutate(temps = rep(1:(duree+1), times = N))

# Visualization
ggplot(n_histo_2, aes(x = temps, y = value)) + 
  geom_line(aes(color = variable))


#faire avec filter

n_histo_3 = n_histo_2 %>% 
  mutate(x0_value = rep(X, each = duree+1)) %>% 
  subset(value > 1/N)

ggplot(n_histo_3, aes(x = x0_value, y =  temps)) + 
  geom_line(aes(color = variable))




# sigma et sa relation à l'étalement de la competition

# plus sigma faible, plpus un une faible différence de trait limite la compétition, et limite la compétition 
#et permet théotiquement d'acceuillir plus d'espèces dans ce modèle

library("plot3D")


n_histo_ =  data.frame(n_histo)
colnames(n_histo_) = 1:N 

n_histo_3 = n_histo_ %>%
  tibble() %>% 
  gather(key = "variable", value = "value") %>% 
  mutate(temps = rep(1:(duree+1), times = N)) %>% 
  mutate(x0_value = rep(X, each = duree+1)) %>% 
  subset(value > 1/N) %>% 
  mutate(species = as.numeric(variable)) %>% 
  mutate(temps    = as.numeric(temps))
  
x = n_histo_3$species
y = n_histo_3$temps
z = n_histo_3$value

scatter3D(x, y, z)

# Evolution ---------------------------------------------------------------




