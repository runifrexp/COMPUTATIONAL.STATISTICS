##############################################################################
######## PARTE 1: Simulazione e generazione di numeri pseudo-casuali  ########
##############################################################################

##############################################################################
### 1.1 generazione mediante Generatori Congruenti Lineari
##############################################################################

options(scipen=999) # Per i numeri superiori a scipen si preferisce la notazione
# esponenziale, per gli altri la notazione fissa.

### Numeri Pseudo-Casuali Uniformi in [0,1]
### Algoritmo Wichmann-Hill (1982)

.current.seed <- c(123,456,789) # Seme iniziale

runif.wh <- function(n){
  a <-c(171,172,170)
  b <- c(30269,30307,30323) # Sono tre numeri primi.
  s <- .current.seed
  u <- rep(0, n)
  for(i in 1:n){
    s <- (a*s) %% b # Aggiornamento sequenza (n volte)
    u[i] <- sum(s/b) %% 1 # Parte frazionari della somma
    # della sequenza dei numeri in s 
  }
  .current.seed <<- s
  u
} 

### Esercizio n. 1
###
### a. Utilizzare il metodo della simulazione per inversione
###    per simulare da una v.c. esponenziale.
### b. Verificare l'adattamento su n=10000 valori simulati.

# La funzione di ripartizione F di Y ~ Exp(lambda) è F(y)=1-e^{-lambda*y}

# Funzione densità v.c. esponenziale.

dExp <- function(y, lambda=1){
  lambda*exp(-lambda*y)
}

curve(dExp(x),0,10,lty="solid",ylab="f(x)")
curve(dExp(x,0.5),lty="dashed",add=TRUE)
curve(dExp(x,1.5),lty="dotted",add=TRUE)

# Funzione di ripartizione v.c. esponenziale

pExp <- function(y, lambda=1){
  1-exp(-lambda*y)
}

curve(pExp(x),0,10,lty="solid",ylab="F(x)")
curve(pExp(x,0.5), lty="dashed",add=TRUE)
curve(pExp(x,1.5), lty="dotted",add=TRUE)

# La funzione quantile e' F^−1(u) = -log(1-u)/lambda, 0 < u < 1

qExp <- function(u, lambda=1){
  -log(1-u)/lambda
}

curve(qExp(x),0,1,lty="solid",xlab="u",ylab=expression(paste(F^{-1},(u)))) # Uso "expression" per le
# formule nelle etichette
curve(qExp(x,0.5), lty="dashed",add=TRUE)
curve(qExp(x,1.5), lty="dotted",add=TRUE)

rExp <- function(n, lambda=1){
  y <- qExp(runif.wh(n),lambda)
  #  y <- -log(1 - runif.wh(n))/lambda # Equivalente
  
  return(y)
}

y <- rExp(5, 2)
y

# Verifichiamo la bontà di adattamento con n=10000

y <- rExp(10000,2) # 10000 numeri simulati
plot(ecdf(y),do.points=FALSE,xlim=c(0,10)) # FdR empirica
curve(pExp(x,2),0,10,add=TRUE,col="red") # FdR teorica

# Per l'incorrelazione:

acf(y,100,ylim=c(-1,1))
Box.test(y,type="L") # Test per il ritardo 1 (lag=1)
# ?Box.test per informazioni


##############################################################################
### 1.2 generazione mediante Accettazione/rifiuto
##############################################################################

### Esercizio n. 2: Simulare mediante accettazione/rifuto da una Beta(2,2)

## Si parta da una Uniforme(0,1)

## Scriviamo l'algoritmo A/R

r.accetto.rifiuto <- function(n, f, g, rg, k, report=TRUE){
  y <- double(n)
  ntry <- 0
  for (i in 1:n){
    done <- FALSE
    while (!done){
      ntry <- ntry + 1
      z <- rg(1)
      u <- k * g(z) * runif.wh(1)
      if (u <= f(z)) done <- TRUE
    }
    y[i] <- z
  }
  if (report)
    cat("r.accetto.rifiuto: ho fatto",ntry,"tentativi\n")
  y
}

## Si trovi dapprima il max della Beta(2,2) in questo intervallo
## A tale scopo scriviamo la distribuzione di probabilità della beta definendo:

Beta <- function(alfa,beta){
  (factorial(alfa-1)*factorial(beta-1))/factorial(alfa+beta-1)
}
dBeta <- function(y,alfa=2,beta=2){
  (y^(alfa-1)*(1-y)^(beta-1))/Beta(alfa,beta)
}

# Si noti che:
# - factorial(x) coincide con la funzione gamma(x+1) e, per x intero,
#   e' x fattoriale.
# - la funzione Beta definita da noi coincide con la funzione beta
#   presente in R.

x <- seq(0,1,0.001)
k = max(dBeta(x,2,2))

# Usando optimize:

optimize(dBeta,c(0,1),maximum=T) # Uso: optimize(funzione, range di x, devo trovare il massimo?)
# $maximum e' la x a cui corrisponde il massimo
# $objective e' il massimo valore della funzione
# Quindi, potremmo scrivere k = optimize(dBeta,c(0,1),maximum=T)$objective

# Si noti che in questo caso k puo' anche essere trovato analiticamente,
# derivando la densita' beta.

r.beta.ar <- r.accetto.rifiuto(10000, dBeta, dunif, runif, k)

## Adesso verifichiamo che il nostro algoritmo funzioni bene:

# Adattamento:

hist(r.beta.ar, prob = T) # "freq = F" e "prob = T" sono equivalenti
curve(dBeta(x,2,2), add = T)

plot(ecdf(r.beta.ar),do.points=FALSE) # FdR empirica
curve(pbeta(x,2,2),add=TRUE,col="red") # FdR teorica

# Incorrelazione:

acf(r.beta.ar,100,ylim=c(-1,1))
Box.test(r.beta.ar,type="L")

### Proviamo a cambiare il parametro k

r.beta.ar.k2 <- r.accetto.rifiuto(10000, dBeta, dunif, runif, k=2)
r.beta.ar.k3 <- r.accetto.rifiuto(10000, dBeta, dunif, runif, k=3)
r.beta.ar.k1 <- r.accetto.rifiuto(10000, dBeta, dunif, runif, k=1) # k troppo piccolo!

# Il numero medio di tentativi richiesti per generare n
# valori e' n*k. Dunque, k deve essere il piu' piccolo
# possibile, ma non troppo piccolo: la densita' da cui generiamo (che
# nei lucidi chiamiamo g), moltiplicata per k, deve dominare la
# densita' obiettivo (che nei lucidi chiamiamo f).
# Nel caso presente, in cui generiamo da una densita' g uniforme, cio'
# significa che k non puo essere piu' piccolo del massimo della
# densita' f (che e' una beta(2,2)).

plot(ecdf(r.beta.ar.k2),do.points=FALSE)
curve(pbeta(x,2,2),add=TRUE,col="red") # OK

plot(ecdf(r.beta.ar.k3),do.points=FALSE)
curve(pbeta(x,2,2),add=TRUE,col="red") # OK

plot(ecdf(r.beta.ar.k1),do.points=FALSE)
curve(pbeta(x,2,2),add=TRUE,col="red") # NO!


### Esercizio n. 3: Simulare una mezza normale standard (normale troncata) con
###                 un algoritmo A/R basato su una distribuzione esponenziale.

## Si consideri la seguente distribuzione di probabilità:

dHalfNorm <- function(y){ # Densità di una normale troncata (y > 0).
  (2/sqrt(2*pi))*exp(-y^2/2)
}

curve(dHalfNorm(x),0,10)

# Consideriamo la distribuzione esponenziale di parametro lambda = 1:
#
# g(x) = exp(-x).
#
# Dobbiamo trovare il valore di k tale che
#
# k * g(x) >= f(x) per ogni x, dove f(x) è la densità della normale
# troncata ai soli valori positivi.
#
# Deve essere dunque:
#
# k * exp(-x) >= 2/sqrt(2*pi) * exp(-x^2/2) per ogni x.
#
# Equivalentemente:
#
# k >= sqrt(2/pi) * exp(x-x^2/2)
#
# E' facile vedere che x-x^2/2 raggiunge un massimo per x = 1.
# Dunque, possiamo definire k = sqrt(2/pi) * exp(1/2) = sqrt(2*exp(1)/pi).

k <- sqrt(2*exp(1)/pi)

curve(dExp(x)*k,0,10,lty="dashed",add=TRUE)

# Usiamo la funzione rExp, creata nell'esercizio n. 1, all'interno della
# funzione per la generazione A/R

r.halfNorm.ar <- r.accetto.rifiuto(10000, dHalfNorm, dExp, rExp, k)

## Controllo simulazione

hist(r.halfNorm.ar, prob= T, nclass=50)
curve(dHalfNorm(x), add = T)


##############################################################################
### 1.3 generazione mediante Trasformazione (Box-Muller)
##############################################################################

## Esercizio n. 4

# Funzione per generare una Normale Standard mediante transformazione
# (algoritmo di Box-Muller).

sim.norm.bm <- function(n){
  u1 <- runif(n)
  u2 <- runif(n)
  R <- sqrt(-2*log(u1))
  theta <- 2*pi*u2
  x <- R*cos(theta)
  y <- R*sin(theta)
  x # Naturalmente, scritto cosi' l'algoritmo e' inefficiente
  # perche' buttiamo via meta' dei valori generati.
}

out1 <- sim.norm.bm(2000)

hist(out1,br=25,xlab="",ylab="",main="2,000 simulazioni da una Normale",freq=FALSE)
curve(dnorm(x),from=-4,to=4,col="blue",add=TRUE)

# Esercizio n. 4.1: si scriva sim.norm.bm in maniera piu' efficiente.

# Esempio di soluzione.

sim.norm.bm.1 <- function(n){
  n1 <- floor((n-1)/2) + 1 # Per n pari n1 = n/2;
  # per n dispari n1 = (n + 1)/2
  #                  = (parte intera di (n-1)/2) + 1
  u1 <- runif(n1)
  u2 <- runif(n1)
  R <- sqrt(-2*log(u1))
  theta <- 2*pi*u2
  x <- R*cos(theta)
  y <- R*sin(theta)
  c(x,y)[1:n] # Considero solo i primi n valori (se n e' dispari ne genero n+1).
}

out1 <- sim.norm.bm.1(2000)

hist(out1,br=25,xlab="",ylab="",main="2,000 simulazioni da una Normale",freq=FALSE)
curve(dnorm(x),from=-4,to=4,col="blue",add=TRUE)

# Confrontiamo l'efficienza di sim.norm.bm e sim.norm.bm.1.

system.time(sim.norm.bm(1E6))
system.time(sim.norm.bm.1(1E6))

# Esercizio n. 4.2: si confronti l'efficienza di Box-Muller
# (implementazione efficiente) con l'algoritmo per generare da una
# normale basato sul metodo accettazione/rifiuto, a partire da una
# distribuzione di Laplace, visto a lezione.