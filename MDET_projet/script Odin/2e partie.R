
```{r}
############### Simuler LV pour un nombre arbitraire d'espece ############

#Conditions initiales
x = seq(from = 0, to = 2, length.out = i) # Valeurs des phenotypes
k = c()
for (i in 1:length(x)){
  k = c(c(k),c(K(x[i])))
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
```


```{r}
ode_plot$sp_ID <- as.numeric(ode_plot$sp_ID)
ode_plot$trait <- x[as.numeric(ode_plot$sp_ID)]

ggplot(ode_plot)+
  geom_raster(aes(trait, time,  fill=dens))+
  scale_fill_gradient2(low = "white" ,
                       high = "red")
```


```{r}
######## Nombre d'espece qui persiste en fonction de sigma (pour l'instant ne marche pas..) ##############
seuil <- 0.05
nb_especes <- c()
for (sigma in seq(0.3,5,by = 0.1)){
  ode <-ode(N0,0:100,model,parms = sigma)
  nb = 0
  for (m in tail(ode,1)[,1:i+1]){
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
T = 100000
t = 200
n = rep(0,i)
n[i/6] <- 1 # densites de pop initiale (une seule espece)
N0 <- t(as.matrix(n))

#Initialisation des valeures de K
k = c()
for (i in 1:length(x)){
  k = c(c(k),c(K(x[i])))
}
K = matrix(rep(k,i), ncol = i, nrow = i, byrow = FALSE)

M = A/t(K)

ode <-ode(N0,1:t,model,parms = sigma)
ode_mutation <- ode

#trouver un moyen de faire muter tout le monde... (peut etre metre un seuil de densité min..)

for (z in 2:(T/t)){
  N0 <- t(as.matrix(tail(ode,1)[,1:i+1]))
  a = runif(1)
  index = which.max(N0)
  if (a<0.5){ #mutation vers la gauche
    N0[max(0,index-1)] <- p*N0[index]
    N0[index] <- (1-p)*N0[index]
  } else { #mutation vers la droite
    N0[min(i,index+1)] <- p*N0[index]
    N0[index] <- (1-p)*N0[index]
  }
  ode <-ode(N0,(z*t):(z*t+t),model,parms = sigma)
  
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

```