

library(ggplot2)
library(tidyverse)
library(rgl)


sigma_a_ = 1
lambda_ = 1
K0_ = 1
x0_ = 1

r_ = 1


a = function(x1,x2, sigma_a = sigma_a_){ 
  exp( (-1/2*(x2-x1)**2)/sigma_a**2 )
}

a(x1 =  4, x2 = 5)

K = function(x, K0 = K0_, lambda_=1, x0 = x0_){
   (K0 - lambda_*(x-x0)**2)
}

K(100)

s = function(x1, x2, r = r_){
  r*(1-a(x1, x2)*(K(x1)/K(x2)))
}


s(x1=1, x2=5, r = 1)


X = 1:100

data_ = tibble(X=X)

data_ =data_ %>% 
  mutate(Y = s(X, x2 = 5, r = 1))

ggplot()+
  geom_line(data = data_, aes(x = X, y = Y))




## Set range and domain of plot
x1_interval  <- seq(0.1, 1.9, length.out = 25);
x2_interval  <- seq(0.1, 1.9, length.out = 25);

## Interpolate surface

z  <- outer(x1_interval,x2_interval,
            FUN = s)

p  <- persp(x1_interval,x2_interval,z, theta = 30, phi = 20,
            col = "lightblue", shade = 0.4, ticktype = "detailed")



library(plotly)


x <- x1_interval
y <- x2_interval
plot_ly() %>% add_surface(x = ~x, y = ~y, z = ~z)
