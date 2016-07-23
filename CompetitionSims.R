library(deSolve)
library(RColorBrewer)

#Simulate human-monster dynamics
source("R/Comp_func.R")

state <- c(
  H = 100,        #Humans
  M = 1          #Monsters
  )               

parameters <- c(
  Rh = 0.5,
  Rm = 0.5,
  a_hm = 0.8, #Effect of monsters on humans (higher means monsters outcompete humans more)
  a_mh = 0.8,   #Effect of humans on monsters
  Kh = 200,
  Km = 200
)

EqH <- (Kh + a_hm*Km)/(1-a_mh*a_hm)
EqM <- (Km + a_mh*Kh)/(1-a_mh*a_hm)


time <- seq(0,5000, by = 0.1)

out <- ode(y = state, times = time, func = Comp_func, parms = parameters)

plot(out[,1],out[,2],type="l",ylim=c(0,200))
lines(out[,1],out[,3],lty=2)



