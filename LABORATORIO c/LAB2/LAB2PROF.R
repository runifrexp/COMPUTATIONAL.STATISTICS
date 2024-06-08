#LABORATORIO 2 PROF

##################################################
######## PARTE 2: Simulazione Monte Carlo ########
##################################################

### Esercizio n. 1

h=function(x){(cos(50*x)+sin(20*x))^2}

B <- 10^4

I <- h(runif(B))
mean(I)

s2 = var(I)
se = sqrt(s2/B)
mean(I) + c(-1,1) * qnorm(0.975) * se  # Intervallo di confidenza con grado di
# fiducia 0.95. Forniamo cosi' una
# valutazione della precisione.

integrate(h,0,1) # Per un'integrale unidimensionale, anche la funzione
# integrate risolve il problema.


# Valutazione per B = 2,...,10^4.
# Si noti che B deve essere > 1 per poter valutare l'incertezza.

estint = icl = icu = rep(NA,10^4-1) # estint conterra' le stime dell'integrale.
# icl conterra' gli estremi inferiori degli i.c.
# icu conterra' gli estremi superiori degli i.c.
for (B in 2:10^4){
  print(B)
  I  <- h(runif(B))
  estint[B-1] = mean(I)
  s2 = var(I)
  se = sqrt(s2/B)
  icl[B-1] = mean(I) - qnorm(0.975) * se
  icu[B-1] = mean(I) + qnorm(0.975) * se
}

par(mfrow=c(2,1))
plot(estint,xlab="B",ty="l",ylim=c(min(icl),max(icu)))
lines(icl,col="red")
lines(icu,col="green")

# Metodo piu' furbo, senza ciclo for.

I <- h(runif(10^4))
estint= cumsum(I)/(1:10^4) # estint contiene le medie di B valori generati per h(x),
# con B = 1,...,10^4.
# A differenza del ciclo for, in cui generavo prima B=2, poi B=3,
# poi B=4, ... valori, qui genero B=10^4 valori una volta sola.
csI2 = cumsum((I-estint)^2) # Somme degli scarti quadratici dalla media, per B=1,...,10^4.
s2 = csI2[2:10^4]/(1:(10^4-1)) # Stime della varianza, basate su B=2,3,...,10^4
# simulazioni (al denominatore ho B-1).
se=sqrt(s2/(2:10^4)) # Valori degli standard error, per B=2,...,10^4.

estint = estint[2:10^4] # Considero solo B=2,...,10^4 anche per estint.

plot(estint,xlab="B",ty="l",lwd=2,ylim=c(min(estint+2*se,estint-2*se), # qnorm(0.975) � circa 2
                                         max(estint+2*se,estint-2*se)))
lines(estint+2*se,col="red",lwd=2)
lines(estint-2*se,col="green",lwd=2)

# Diamo un'occhiata alla funzione h(x) in (0,1).

h=function(x){(cos(50*x)+sin(20*x))^2}

par(mfrow=c(2,1))
curve(h,from=0,to=1,xlab="x",ylab="h(x)",lwd="2") # Di default, la curva viene calcolata in
# n=101 punti, che in questo caso sono
# troppo pochi per disegnarla bene.
curve(h,from=0,to=1,xlab="x",ylab="h(x)",lwd="2",n=1E3) # ?curve per dettagli.

par(mar=c(2,2,2,1)) # Con altri margini.
curve(h,from=0,to=1,xlab="x",ylab="h(x)",lwd="2")
curve(h,from=0,to=1,xlab="x",ylab="h(x)",lwd="2",n=1E3)

par(mar=c(5, 4, 4, 2) + 0.1,mfrow=c(1,1)) # Torniamo ai valori di default.


###############################################################################

### Esercizio n. 2

B <- 1000
y1 <- runif(B) # Trasformiamo le variabili con yi = xi/3,
y2 <- runif(B) # ottenendo Y_i con distr. U(0,1), per i = 1,2,3.
y3 <- runif(B)
z <- 27 * exp(-abs(3 * (y1 - y2))^(2 + sqrt(3 * y3)))
cat("Stima:", mean(z), "\nIntervallo:", mean(z)+c(-1,1)*qnorm(0.975)*sd(z)/sqrt(B), "\n")


###############################################################################

## Esercizio n. 3

B <- 1000
x1 <- runif(B, 0, 3)
x2 <- runif(B, 0, 3)
x3 <- runif(B, 0, 3)
z <- 27*exp(-abs(x1 - x2)^(2 + sqrt(x3))) # Devo moltiplicare per 3^3, visto che la
# densita' di una U(0,3) e' 1/3.
cat("Stima:", mean(z), "\nIntervallo:", mean(z)+c(-1,1)*qnorm(0.975)*sd(z)/sqrt(B), "\n")

# Si noti che con questo approccio stiamo applicando, in pratica, il campionamento per
# importanza, con phi_i(x) = 1 e psi_i(x) = 1/3 per i = 1,2,3.
# In questo caso, tuttavia, il nome "per importanza" non ha molto senso, in quanto stiamo usando
# delle densita' psi_i(x) costanti, che non tengono dunque in alcun conto l'andamento della
# funzione da integrare.


###############################################################################

### Esercizio n. 4

f <- function(x){exp(-sqrt(abs(x)))*sin(x^2/50)}

curve(f,-60,60,n=1000) # La funzione f(x) non ha un massimo in 0,
# come la densita' N(0,400). Tuttavia, e'
# chiaro che f(x) e' simmetrica attorno a 0 e
# che in generale assume valori sempre piu'
# piccoli allontanandosi da 0. Inoltre, e' chiaro
# che il contributo all'integrale dei valori f(x),
# con |x| > 40 = 2 * dev. standard, e' piccolo.
# Queste considerazioni rendono la densita' N(0,400)
# adeguata per il campionamento per importanza.

B <- 10^4
x <- rnorm(B,0,20)
z <- f(x)/dnorm(x,0,20) # Campionamento per importanza, con phi(x) = 1 e psi(x)
# densita' di una N(0,400).
cat("Stima:", mean(z), "\nIntervallo:", mean(z)+c(-1,1)*qnorm(0.975)*sd(z)/sqrt(B), "\n")

integrate(f,-Inf,Inf) # Per confronto.

# Proviamo a risolvere il problema senza campionamento per importanza,
# riportandoci all'intervallo (0,1) con una trasformazione.
# In particolare, usiamo y = 1/(1+exp(-x)),
# che assume valori in (0,1) quando x varia sulla retta reale.

# Otteniamo x = -log((1-y)/y) e dx/dy = 1/(y*(1-y)).
# Quindi integreremo sull'intervallo (0,1) la funzione

fy <- function(y){
  x = -log((1-y)/y)
  1/(y*(1-y)) * exp(-sqrt(abs(x)))*sin(x^2/50)
}

w <- fy(runif(B))
cat("Stima:", mean(w), "\nIntervallo:", mean(w)+c(-1,1)*qnorm(0.975)*sd(w)/sqrt(B), "\n")

# La precisione e' molto inferiore a quella ottenuta con il
# campionamento per importanza. Inoltre, possiamo avere problemi
# di stabilita' per valori generati da runif troppo vicini a 0 o a 1
# (per y=0 e y=1, si ha che x diverge).

integrate(fy,1E-16,1-1E-16) # Dobbiamo escludere 0 e 1 dall'intervallo di
# integrazione, altrimenti x diverge.

###############################################################################

### Esercizio n. 5

B <- 10000

lambda <- 1

# Procedura di stima "diretta".

x <- rexp(B, lambda)
y <- sapply(x, function(x) runif(1, 0, x)) # Generiamo B valori per Y, che ha distrib. U(0,X).
cat("Stima di E(Y):", mean(y), "\nIntervallo:", mean(y)+c(-1,1)*qnorm(0.975)*sd(y)/sqrt(B), "\n")

# Possiamo usare Rao-Blackwell, dato che sappiamo che la media di una
# U(0,X), condizionata a X, e' X/2.

cat("Stima di E(Y):", mean(x/2), "\nIntervallo:", mean(x/2)+c(-1,1)*qnorm(0.975)*sd(x/2)/sqrt(B), "\n")

# Confrontiamo le precisioni delle due procedure via Monte Carlo,
# invece che usando intervalli di confidenza basati sul teorema del
# limite centrale come abbiamo fatto sopra.

M <- 500
stima.diretta <- stima.rb <- rep(NA,M)
for (i in 1:M){
  print(i)
  x <- rexp(B, lambda)
  y <- sapply(x, function(x) runif(1, 0, x))
  
  stima.diretta[i] <- mean(y)
  stima.rb[i] <- mean(x/2)
}

mean(stima.diretta) # In media il risultato e' lo stesso.
mean(stima.rb)

sd(stima.diretta) # Con Rao-Blackwell la precisione aumenta.
sd(stima.rb)

# Se, oltre a Rao-Blackwell, usiamo anche le variabili antitetiche,
# cambia solo il modo di generare x.

stima.rbant <- rep(NA,M)
for (i in 1:M){
  u <- runif(B/2)
  x1 <- qexp(u,lambda)
  x2 <- qexp(1-u,lambda)
  x <- c(x1,x2)
  stima.rbant[i] <- mean(x/2)
}

mean(stima.rbant) # In media il risultato e' lo stesso.
sd(stima.rbant) # La precisione e' superiore a quella ottenuta usando solo
# Rao-Blackwell.

boxplot(cbind(stima.diretta,stima.rb,stima.rbant))

###

# Stima di P(Y > 1).

# Procedura "diretta".

x <- rexp(B, lambda)
y <- sapply(x, function(x) runif(1, 0, x))
y.1 <- y > 1 # y.1 vale TRUE quando y > 1 e vale 0 altrimenti.
p.hat = mean(y.1) # Stima di P(Y > 1).
cat("Stima di P(Y > 1):", p.hat,
    "\nIntervallo:", p.hat+c(-1,1)*qnorm(0.975)*sqrt(p.hat*(1-p.hat)/B), "\n")

# Sappiamo che P(Y > 1 | X) e' (X-1)/X se X>1 e 0 altrimenti. Quindi,
# possiamo applicare Rao-Blackwell come segue.

f <- function(x) (x-1)/x * (x>1)
p.hat <- mean(f(x))
cat("Stima di P(Y > 1):", p.hat,
    "\nIntervallo:", p.hat+c(-1,1)*qnorm(0.975)*sqrt(p.hat*(1-p.hat)/B), "\n")

# Valutiamo la precisione via Monte Carlo.

M <- 500
stima.diretta <- stima.rb <- rep(NA,M)
for (i in 1:M){
  print(i)
  x <- rexp(B, lambda)
  y <- sapply(x, function(x) runif(1, 0, x))
  
  stima.diretta[i] <- mean(y>1)
  stima.rb[i] <- mean(f(x))
}

mean(stima.diretta) # In media il risultato e' lo stesso.
mean(stima.rb)

sd(stima.diretta) # Con Rao-Blackwell la precisione aumenta.
sd(stima.rb)

# Come sopra, usando le variabili antitetiche cambia solo il modo di generare x.

stima.rbant <- rep(NA,M)
for (i in 1:M){
  u <- runif(B/2)
  x1 <- qexp(u, lambda)
  x2 <- qexp(1-u, lambda)
  x <- c(x1,x2)
  
  stima.rbant[i] <- mean(f(x))
}

mean(stima.rbant) # Stessa media
sd(stima.rbant) # Maggiore precisione.

boxplot(cbind(stima.diretta,stima.rb,stima.rbant))

###############################################################################

### Esercizio n. 6 (Birthday problem)

# Come funziona "sample".

sample(1:10,3) # Campionamento di 3 valori, nell'insieme {1,2,...,10},
# senza reinserimento.
sample(1:10,3,replace=TRUE) # Campionamento di 3 valori, nell'insieme {1,2,...,10},
# con reinserimento.

# Come funziona "duplicated".

duplicated(1:10) # Nessun valore ripetuto.
duplicated(c(1:4,4,5,6,6)) # TRUE per i valori ripetuti.

# Come funziona "any".

any(rep(F,8)) # "any" restituisce TRUE se almeno uno degli elementi e' TRUE.
any(c(rep(F,8),T))

# Gia' che ci siamo: come funziona "all".

all(rep(F,8)) # "all" restituisce TRUE se tutti gli elementi sono TRUE.
all(c(rep(F,8),T))
all(c(rep(T,8)))

bpSim <- function(n, B){
  # is.rep restituisce un vettore di B codici TRUE e FALSE. Ottengo TRUE se, tra gli n valori
  # generati con sample, c'e' almeno un valore ripetuto.
  is.rep <- function(n,B) sapply(rep(n,B),function(n) any(duplicated(sample(1:365, n, replace=TRUE))))
  
  sapply(n, function(n) mean(is.rep(n,B))) # In questo modo, n puo' essere vettoriale.
}

B <- 10000
n <- 2:50
prob <- bpSim(n,B) # Probabilita' stimata di compleanni ripetuti, quando si
# prendono n persone a caso.

# Disegniamo i risultati.

plot(n, prob, pch=20)
abline(h=0.5,lty=2)
min(n[prob>0.5]) # Con n almeno pari a 23, la probabilita' che tra n persone
# prese a caso ce ne siano almeno due nate lo stesso giorno
# e' maggiore della probabilita' complementare (ossia la
# probabilita' che siano nate tutte in giorni diversi).


# Se non conosco any e duplicated, devo lavorare in maniera inefficiente:

any.duplicated <- function(x){
  n = length(x)
  ris <- FALSE
  
  for(i in 1:(n-1)){ # I cicli for sono inefficienti.
    cat("i =",i,"\n") # Togliere il commento per il debugging.
    for(j in (i+1):n){
      cat("j =",j,"\n") # Togliere il commento per il debugging.
      if(x[i] == x[j]){
        ris <- TRUE
        break # Fatto: usciamo dal ciclo for per j.
      }
    }
    if(ris) break # Fatto: usciamo dal ciclo for per i.
  }
  
  ris
}

# Proviamo se funziona.

any.duplicated(c(1:100,100))
any.duplicated(rep(7,10))
any.duplicated(3:107)
any.duplicated(c(3,3:107))

# Funziona. Scriviamo dunque una funzione analoga a bpSim, ma
# basata su any.duplicated, e quindi inefficiente.

bpSim1 <- function(n, B){
  is.rep <- function(n,B) sapply(rep(n,B),function(n) any.duplicated(sample(1:365, n, replace=TRUE)))
  sapply(n, function(n) mean(is.rep(n,B)))
}

# Otteniamo gli stessi risultati, ma il tempo impiegato e'
# sensibilmente superiore.

prob1 <- bpSim1(n,B)
plot(n, prob1, pch=20)
abline(h=0.5,lty=2)
min(n[prob1>0.5])

system.time(bpSim(n,B))
system.time(bpSim1(n,B))