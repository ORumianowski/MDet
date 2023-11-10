library(rgl)
library(deSolve)
library(tidyverse)

# clear
#Penser à commenter

main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=28))

#Parametres
i = 20 # Nombre d'especes
Xmin = 0
Xmax = 2
X0 = (Xmax - Xmin)/2
K0 = 1
lambda = 1
sigma = 0.2 #etalement de la competition
r = 1
############ Fonctions #############

calcul_K <- function (X){
  K = max(0, K0 - lambda*(X - X0)^2) + 10^(-9) # On rajoute 10^-9 pour ne pas avoir de pb de division par 0
  return (K)
}
a <- function (Xi,Xj){
  a = exp(-0.5*(Xi-Xj)^2/sigma^2)
  return (a)
}
fitness <- function (X1,X2){
  s = r*(1-a(X1,X2)*calcul_K(X1)/calcul_K(X2))
  return (s)
}

############## Plot fitness ####################

plot3d(fitness, xlim = c(0.1,1.9), ylim = c(0.1,1.9), n=100, col = colorRampPalette(c('green','yellow','red')) )

#plot du signe de la fonction (permet de visualiser les eq ???) Pairwise invasibility plot PIP
varx = seq(0.1,1.9,length.out = 100)
vary = seq(0.1,1.9,length.out = 100)

############### Simuler LV pour un nombre arbitraire d'espece ############



#Conditions initiales
x = seq(from = 0, to = 2, length.out = i) # Valeurs des phenotypes
k = c()
for (i in 1:length(x)){
  k = c(c(k),c(calcul_K(x[i])))
}
 # K de chaque espece
n = rep(1/i,i) #densite de pop initiale

# On met sous forme de matrice
N0 = t(as.matrix(n))
K = matrix(rep(k,i), ncol = i, nrow = i, byrow = FALSE)
X = matrix(rep(x,i), ncol = i, nrow = i, byrow = TRUE)

#Calcul de la matrice M
D = t(X) - X
A = exp(-0.5*D^2/(sigma^2))
M = A/t(K)

#Resolution du systeme d'ODE
model <- function(t,N,sigma){
  dN <- r * N * (1 - N %*% M)
  return(list(dN))
}

ode <-ode(N0,0:100,model,parms = sigma)

# Representation
#barplot(ode[100,1:i+1])

ode <- as.data.frame(ode)
#colnames(ode) <- c("t",'x1','x2','x3')
ode_plot = ode%>%
  pivot_longer(-time, values_to = "dens", names_to = "sp_ID")

ggplot(ode_plot)+
  geom_line(aes(time, dens, col = sp_ID))

ode_plot$sp_ID <- as.numeric(ode_plot$sp_ID)
ode_plot$trait <- x[as.numeric(ode_plot$sp_ID)]

ggplot(ode_plot)+
  geom_raster(aes(trait, time,  fill=dens))+
  scale_fill_gradient2(low = "white" ,
                       high = "red")
######## Nombre d'espece qui persiste en fonction de sigma (pour l'instant ne marche pas..) ##############
seuil <- 0.05
nb_especes <- c()
for (s in seq(0.3,5,by = 0.1)){
  ode <- ode(N0,0:100,model,parms = s)
  nb = 0
  for (m in as.vector(tail(ode,1)[,1:i+1])){
    if (m>seuil){
      nb = nb + 1
    }
  }
  nb_especes <- c(c(nb_especes),c(nb))
  
}
plot(seq(0.3,5,by = 0.1),nb_especes)
######## Implementation de mutation ##########

#parametres
p = 0.1 #taux de la pop qui mute
Temps = 10000 #temps de la simulation
t = 200 #pas de temps entre chaque mutation
n = rep(0,i) 
n[i/6] <- 1 # densites de pop initiale (une seule espece)
N0 <- t(as.matrix(n))

#Initialisation des valeures de K
k = c()
for (i in 1:length(x)){
  k = c(c(k),c(calcul_K(x[i])))
}
K = matrix(rep(k,i), ncol = i, nrow = i, byrow = FALSE)

M = A/t(K)

# Evolution sans mutation sur le premier pas de temps t
ode <-ode(N0,1:t,model,parms = sigma, method = 'lsoda')
ode_mutation <- ode

#trouver un moyen de faire muter tout le monde... (peut etre metre un seuil de densité min..)

seuil = 0.05 #densité min on considère qu un trait en dessous a une densité de 0


for (z in 2:(Temps/t)){
  
  #Elimination des traits sous le seuil
  tail <- as.matrix(tail(ode,1)[,1:i+1])
  for (j in 1:length(tail)){
    if(tail[j] < seuil){
      tail[j] <- 0
    }
  }
  N0 <- t(tail)
  
  # Mutation des traits restants
  for (index in which(tail != 0)){
    a = runif(1)
    if (a<0.5){ #mutation vers la gauche
      N0[max(0,index-1)] <- p*N0[index]
      N0[index] <- (1-p)*N0[index]
    } else { #mutation vers la droite
      N0[min(i,index+1)] <- p*N0[index]
      N0[index] <- (1-p)*N0[index]
    }
   
  }
  
  
  ode <-ode(N0,(z*t):(z*t+t),model,parms = sigma, method='lsoda')
  
  ode_mutation <- rbind(ode_mutation,ode) 

}

ode_mutation_plot = as.data.frame(ode_mutation)%>%
  pivot_longer(-time, values_to = "dens", names_to = "sp_ID")

ode_mutation_plot$sp_ID <- as.numeric(ode_mutation_plot$sp_ID)
ode_mutation_plot$trait <- x[as.numeric(ode_mutation_plot$sp_ID)]

ggplot(ode_mutation_plot)+
  geom_raster(aes(trait, time,  fill=dens))+
  scale_fill_gradient2(low = "white" ,
                       high = "red") +
  main_theme

# ggplot(ode_mutation_plot)+
#   geom_tile(aes(trait, time,  fill=dens))+
#   scale_fill_gradient2(low = "white" ,
#                        high = "red") +
#   main_theme
