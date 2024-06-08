#LAB 5


#-------------------------------------ES 1--------------------------------------

cmp <-  c(1,1,1,1,0,0,0,0,0,0)
n <- length(cmp)
w <- sum(cmp)


#scrittura 'analitica' della veros. e della log veros.
likelihood <- function(p, w, n) p^w*(1-p)^(n-w) #(non necessaria)

llik.1 <- function(p, w, n) w*log(p)+n*log(1-p)-w*log(1-p)

loglik.1 <- function(p, w, n) { #per poter prendere input vettoriali (e quindi 
  sapply(p, function(p) llik.1(p, w, n)) #poter plottare la funzione)
}

loglik.1(c(0.1,0.2,0.5,0.8), w, n)


########################################trucchetto input vettoriali
g <- function(x){
  if (x<3)
    y = 0
  else
    y = 1
  y
}
g(3)
g(c(4,2,5,1,1,5))
#cosi potra' prendere input vettoriali
gv <- function(x) sapply(x, g)
gv(c(3,4,2,5,1))
##################################################################


#scrittura della logverosimiglianza con il 'trucchetto'
dens <- function(p, w, n) (p^w) * (1-p)^(n-w)

llik.2 <- function(p, w, n) sum(log(dens(p, w, n)))

loglik.2 <- function(p, w, n){ #per poter prendere input vettoriali (e quindi 
  sapply(p, function(p) llik.2(p, w, n)) #poter plottare la funzione)
}

loglik.2(c(0.1,0.2,0.5,0.8), w, n)


#le due definizioni di logverosimiglianza sono equivalenti 
curve(loglik.1(x, w, n), 0,1)
curve(loglik.2(x, w, n), 0,1)

#dato che sappiamo che dobbiamo stimare una probabilita', possiamo dare in input 
#un'intervallo preconfezionato alla funzone optimize, senza dover 
#riparametrizzare
optimize(llik.1, w , n, interval = c(0,1), maximum = T)
optimize(llik.2, w , n, interval = c(0,1), maximum = T)

optimize(llik.1, w , n, interval = c(0,1), maximum = T)$maximum
optimize(llik.2, w , n, interval = c(0,1), maximum = T)$maximum

#confronto con la stima anaitica
p.cap <- w/n
p.cap


####parte b

alpha = 0.05

#@@WALD@@: ic basato sulla distribuzione asintotica dello SMV 
#metodo 'analitico'
p.cap <- w/n
p.cap

v.p.cap <- p.cap*(1-p.cap)/n

ic.1 <- p.cap + c(-1,1) * qnorm(1-alpha/2) * sqrt(v.p.cap)
ic.1


#metodo 'alternativo' usando l'informazione osservata
library(nlme) #per poter usare fdHess
v = -1/fdHess(p.cap, function(x) llik.1(x,w,n))$Hessian 
s = sqrt(v)
s
s = as.double(s) #al psoto di vedere l'output come matrice 1x1, lo teniamo come
#scalare

ic.2 <- p.cap + c(-1,1) * s * qnorm(1-alpha/2)
ic.2


 


#@@WILKS@@ (doppia caduta di log-veros.) (preferibile perche' ha un approssimazione
#in meno), (distribuzione asintotica del rapporto di verosimiglianza)

curve(llik.1(x, w, n)- llik.1(p.cap, w, n), from = 0, to = 1)
abline(h = -qchisq(1-alfa,1)/2, lty=2)

g <- function(l, w, n){
  abs(llik.1(l, w, n) - llik.1(p.cap, w, n) + (qchisq(0.975/2, 1)))
}

lower.bound <- optimize(g, w, n, interval = c(0,p.cap))$minimum
upper.bound <- optimize(g, w, n, interval = c(p.cap,1))$minimum

c2 <- c(lower.bound, upper.bound)
c2





#-------------------------------------ES 2--------------------------------------
rm(list=ls())

#Y rappresenta il numero di insuccessi che precedono il k-mo successo, in una 
#successione di prove indipendenti con probabilita' di successo p
f <-  c(6, 6, 6, 8, 6, 4, 5, 4, 1, 2, 1, 1, 1, 1)
y <- rep(0:13, f)



###punto a

func <- function(r,k){ #definizione del fattoriale 
  gamma(r+k)/(gamma(k)*factorial(r))
}

llik <- function(par,f){ #log veros. (input parametri e frequenze)
  k <- par[1]
  p <- par[2]
  r <- 0:(length(f)-1) #valori osservati 
  dens <- func(r,k) * p^k * (1-p)^r 
  sum(f*log(dens))
}

llik(c(2, 0.4), f)




###punto b


#disegno le curve di livello (in quanto la veros. e' bidimensionale)

#per capire gli intervalli in cui trovo le SMV, faccio prima il ###punto c,
#oppure vado a tentativi.
#p sappiamo essere in [0,1], mentre k sappiamo essere > 1

cat('p=',p.tilde,'\nk=',k)

k.seq1 <- seq(1, 100, by = 1)
p.seq1 <- seq(0.01, 0.99, by = 0.01)
 
z1 <- matrix(NA, length(k.seq1), length(p.seq1))
for (i in 1:length(k.seq1)){
  for (j in 1:length(p.seq1)){
    z1[i, j] <- llik(c(k.seq1[i], p.seq1[j]), f)
  }
}

range(z1)

livelli = quantile(z1, prob = c(0.5, 0.7, 0.8,
                               seq(0.9,1,by=0.01)))
contour(k.seq1, p.seq1, z1, levels = livelli, xlab="k", ylab="p")

#faccio un grafico piu' vicino alla SMV

#ridefinisco le sequenze
k.seq2 <- seq(1, 15, by = 1)
p.seq2 <- seq(0.01, 0.5, by = 0.001)


z2 <- matrix(NA, length(k.seq2), length(p.seq2))
for(i in 1:length(k.seq2)){
  for (j in 1:length(p.seq2)){
    z2[i, j] <- llik(c(k.seq2[i], p.seq2[j]), f)
  }
}

dim(z2)

livelli = quantile(z2, prob = c(0.5, 0.7, 0.8,
                               seq(0.9,1,by=0.01)))
contour(k.seq2, p.seq2, z2, levels = livelli, xlab="k", ylab="p")



###punto c

#metodo dei momenti: sostituisco momento teorico a momento empirico 

m <- mean(y)
v <- var(y)

p.tilde <- m/v
p.tilde 

#trovo k dato che conosco p.tilde

k <- (m*p.tilde)/(1-p.tilde)
k
cat('p =',p.tilde,'\nk =',k) #stime preliminari di p e k



###punto d


#MAPPE E MAPPE INVERSE

#@p
#mappa da [0,1] a R: logit 
theta <- log(p/(1-p))
#mappa inversa 
p <- e^theta/(1+e^theta)

#@k
#mappa da (a, +inf) a R: traslazione e log
#traslazione 
phi <- k - a
#log
psi <- log(phi)
psi <- log(k-a)
#mappa inversa
k <- a + e^psi


psi <- log(k-1)

#minimizzo llik cambiata di segno con input le mappe inverse
#definisco la nuova veros. 
llik.rip.neg <- function(par, f){
  k <- 1 + exp(par[1])
  p <- exp(par[2])/(1+exp(par[2]))
  - llik(c(k,p), f)
}

llik.rip.neg(c(2, 0.4), f)
llik(c(2, 0.4), f)



#valori di partenza passati per la mappa
vp <- c(log(p.tilde/(1-p.tilde)), log(k-1))


#ottimizzo rispetto ai valori passati per la mappa inversa e con valori iniziali
#i valori stimati passati per la mappa
opt <- optim(fn = llik.rip.neg, par = vp, f = f)
opt

#ritorno ai parametri iniziali con la mappa inversa
psi <- opt$par[1]
psi
theta <- opt$par[2]
theta
p <- exp(theta)/(1+exp(theta))
p
k <- 1 + exp(psi)
k


#-------------------------------------ES 4--------------------------------------










