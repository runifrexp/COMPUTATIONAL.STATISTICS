#LABORATORIO 3


#-------------------------------------ES 1--------------------------------------

#creazione del campione (faccio finta che sia un campione qualunque, che io non
#sappia la sua distribuzione)
dati <- rnorm(50, 0, 1)

#media campionaria 
mean(dati)
n <- length(dati)

#creo 1000 campioni bootstrap 
B <- 1000
bts <- function(dati, B){
  cmp.bts <- rep(NA, B)
  n <- length(dati)
  for (i in 1:B){
    cmp.bts[i] <- mean(dati[floor(1+n*runif(n))])
  }
  return(cmp.bts)
}

#creazione campioni bootstrap
ybar.bts1 <- bts(dati,B)

#oppure piu' velocemente
ybar.bts2 <- replicate(B,mean(dati[floor(1 + n * runif(n))])) 

#stimo la distorsione
mean(ybar.bts1) - mean(dati)
mean(ybar.bts2) - mean(dati)


#invece la varianza e'
var(ybar.bts1)
var(ybar.bts2)

par(mfrow = c(1,2))
hist(ybar.bts1, nclass = 50)
hist(ybar.bts2, nclass = 50)


# <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 
#la distorsione vale molto contanto che e' nello stesso ordine di grandezza 
#della distibuzione da cui si campiona bootstrap

#verifichiamola meglio con MC sopra Bootstrap
MC <- function(Y, f, B, n = 40){
  v <- rep(NA, B)
  for (i in 1:B){
    v[i] <- mean(camp(Y, n))
  }
  hist(v, nclass = 50)
  return(f(v))
}

B <- 1000
#distorsione
MC(Y, mean, B) - mean(Y)
#la media campionaria e' un stimatore poco distorto

#mentre la varianza e'
MC(Y, var, B)
#la media campionaria e' uno stimatore tutto sommato efficiente
# <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 



#punto b

#aumento la varianza per vedere cosa succede
dati2 <- rnorm(50, 0, 10)

ybar.bts3 <- bts(dati2, B)
#oppure
ybar.bts4 <- replicate(B,mean(dati2[floor(1 + n * runif(n))])) 

#correttezza
mean(ybar.bts3) - mean(dati2)
mean(ybar.bts4) - mean(dati2)
#la distorsione aumenta da prima, anche se non di molto

var(ybar.bts3)
var(ybar.bts4)
#la varianza aumenta di molto 

hist(ybar.bts3, nclass = 50)
hist(ybar.bts4, nclass = 50)
#si nota varianza piu' alta


#punto c

#lascio invariata media e varianza e aumento la numerosita' del campione da cui
#faccio la stima bootstrap 
B2 <- 10000

ybar.bts5 <- bts(dati, B2)

mean(ybar.bts5) - mean(dati)
#la distorsione diminusice (aumentare B nel bts e' come aumentare n nel 
#mondo reale)

var(ybar.bts5)
#la varianza non diminuisce di molto

hist(ybar.bts5, nclass = 50)
#la distribuzione e' molto incentrata vicino alla media


#<- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 
MC(Y, mean, B2) - mean(Y)

MC(Y, sd, B2)
#si nota che lo stimatore media campionaria bootstrap per B grande aumenta la 
#sua correttezza
#vale lo stesso per l'efficienza
#<- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 





#-------------------------------------ES 2--------------------------------------

dati <- rnorm(50, 0, 1)
dati

med <- median(dati)
med
library(sm)
sm.density(dati)

z <- sm.density(dati, eval.points = med)$estimate
z

#std error asintotico
se.med <- 1/(2 * sqrt(n) * z) 
se.med



#parte b
B <- 1000

med.bts1 <- bts(dati, B)

#oppure piu' semplice
med.bts2 <- replicate(B,median(dati[floor(1 + n * runif(n))]))

mean(med.bts1) - med # Distorsione bootstrap.
mean(med.bts2) - med
sd(med.bts1) #
sd(med.bts2) #








#-------------------------------------ES 3--------------------------------------

#parte a

data(faithful)
dati <- faithful$eruptions
head(dati)

#stimo theta.cap: mediana campionaria 
theta.cap <- median(dati)
theta.cap

#ripeto B volte la stima di theta.cap, ergo ne trovo la distribuzione

f.faith <- function(rf, B){
  sapply(rep(1,B), function(d) median(rf()))
}

#funzione per generare casualmente dai dati e poterne dedurre la distribuzione 
#mettendola in input alla funzione f.faith
n <- length(dati)
r.faith <- function(){
  index <- floor(1 + n * runif(n))
  dati[index]
}

B <- 1000
s.np <- f.faith(r.faith, B)

ic <- 2 * theta.cap - quantile(s.np, c(0.975, 0.025))
ic


#parte b

B1 <- 500
B2 <- 1000
B3 <- 5000
B4 <- 10000

#funzione che fa MC sopra bootstrap
MC.bs <- function(N, B, dati){
  theta.cap <- median(dati)
  ics.inf <- rep(NA, N)
  ics.sup <- rep(NA, N)
  for (i in 1:N){
    ics.inf[i] <- 2 * theta.cap - quantile(f.faith(r.faith, B), 0.975)
    ics.sup[i] <- 2 * theta.cap - quantile(f.faith(r.faith, B), 0.025)
  }
  cat("Deviazione standard dell'estremo inferiore", sd(ics.inf),
      "\nDeviazione standard dell'estremo superiore", sd(ics.sup))
  if (sd(ics.sup) < 0.005 & sd(ics.inf) < 0.005){
    cat('\n',B, "e' un valore sufficientemente grande per avere una deviazione 
standard degli estremi dell'intervallo di confidenza inferiore a 0.005")
  }
  else{
    cat("\nUna delle condizioni non è soddisfatta")
  }
}

MC.bs(20, B1, dati)
MC.bs(20, B2, dati)
MC.bs(20, B3, dati)
MC.bs(20, B4, dati)




#-------------------------------------ES 4--------------------------------------

data(faithful)
dati <- faithful$eruptions
n <- length(dati)

#theta.cap (vero valore del parametro nel mondo bootstrap)
media <- mean(dati)

#funzione che genera un campione bootstrap 
r.bs <- function(n, dati){
  c.bts <- rep(NA, n) #genera un vettore vuoto di 272 elementi per riempirlo con
  #il campione bootstrap
  for (i in 1:n){
    index <- floor(1 + n * runif(1)) #crea l'indice (1 indice) che mi indica che 
    #elemento estrarre dal vettore di dati
    c.bts[i] <- dati[index] #estrae l'elemento con indice generato da index
  }
  return(c.bts)
}

#ripeto B volte la generazione del campione e calcolo le medie dei campioni, 
#quindi l'output e' un vettore di B = 1000 medie di campioni bootstrap
means <- function(r, dati, n, B = 1000){
  mns <- rep(NA, B)
  for (i in 1:B){
    mns[i] <- mean(r(n, dati))
  }
  return(mns)
}

#distribuzione bootstrap degli errori di stima assoluti
err.stm.boot <- means(r.bs(n, dati), dati, n) - media
err.stm.boot
hist(err.stm.boot, prob = T, nclass = 100)

#trovo la probabilita' richiesta
mean(abs(err.stm.boot) > 0.2)



# <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 
#molto piu' facile
data(faithful)
dati <- faithful$eruptions

n = length(dati)

m <- mean(dati) # Media campionaria.

B <- 10000
m.bts <- replicate(B,mean(dati[floor(1 + n * runif(n))]))

mean(abs(m.bts - m) > 0.2) # Stima della probabilit� richiesta.
# La media campionaria, nel mondo
# bootstrap, gioca il ruolo di vera media.
# <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 


#-------------------------------------ES 5--------------------------------------

impianto <- read.table(file.choose(), header = T) #impianto.dat
head(impianto)
y <- impianto$Perdita.Prod
y
x <- impianto$Temperatura
x
c(x,y)
plot(x,y)
#si notano gli outliers
n <- length(x)
n
B <- 1000


#stima bootstrap NP

library(MASS)

b1.cap <- as.numeric(coef(ltsreg(y ~ x))[1])
b1.cap


b2.cap <- as.numeric(ltsreg(y~x)$coefficients[2])
b2.cap

stm.np <- function(x, y, n = length(x), B = 1000){
  params <- matrix(NA, B, 2)
  for (i in 1:B){
    idx <- floor(1 + n * runif(n))
    params[i,] <- coef(ltsreg(x[idx], y[idx]))
  }
  params
}

dim(stm.np(x,y))

#distorsione
#beta1
mean(stm.np(x,y)[,1]) - b1.cap

#beta2
mean(stm.np(x,y)[,2]) - b2.cap






#stima bootstrap SP: assumo la distribuzione della variabile risp. 

m <- ltsreg(y ~ x)

res <- residuals(m)
res

mods.s <- function(b1.cap, b2.cap, x, res){
  index <- floor(1 + n * runif(1))
  return(b1.cap + b2.cap * x[index] + res[index])
}


stm.sp <- function(x, y){
  params <- matrix(NA, B, 2) #vettore contenente le B stime bts di b1 e b2
  m.bts <- matrix(NA, B, n)
  for (i in 1:B){
    for (j in 1:n){
      m.bts[j] <- ltsreg(mods.s(b1.cap, b2.cap, x, res) ~ x)
    }
    params[i,] <- coef(m.bts)
  }
  params
}

stm.sp(x,y)




# <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 
mods.stm <- function(){
  y.cap <- matrix(NA, B, n)
  for (i in 1:n){
    index <- floor(1 + n * runif(1))
    y.cap[i] <- b1.cap + b2.cap * x[index] + res[index]
  }
  return(y.cap)
}
# <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 



res <- residuals(m)




























#-------------------------------------ES 6--------------------------------------
dati <- read.table(file.choose(), header = T) #esponenziale.dat
y <- dati$y
y
n <- length(y)
n
#ic con bts parametrico

media.cap <- mean(y)
media.cap
lambda.cap <- 1/media
lambda.cap

alpha <- 0.05

#stima bootstrap di lambda
B <- 1000
lambda.bts <- replicate(B, 1/mean(rexp(n,lambda.cap)))
lambda.bts

#stima bootstrap degli errori di stima 
err.stm <- lambda.bts - lambda.cap
err.stm

#ic costruito sulla distribuzione degli errori di stima
lambda.cap - c(quantile(err.stm, 1-alpha/2), quantile(err.stm, alpha/2))
lambda.cap - quantile(err.stm, c(1-alpha/2, alpha/2))

#ic costruito sulla distribuzione delle stime bootstrap di lambda
2 * lambda.cap - quantile(lambda.bts,c(1-alpha/2,alpha/2))











