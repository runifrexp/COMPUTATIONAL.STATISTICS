#LABORATORIO 3 PROF 

#### Bootstrap

### Esercizio n. 1

## [a]

# Generiamo i dati.

n <- 50
mu <- 0
sigma <- 1
dati <- rnorm(n,mu,sigma)

# Lo stimatore che usiamo per la media della popolazione
# � la media campionaria

mean(dati) # Media campionaria osservata.

# Stimiamo la varianza della media campionaria con il bootstrap.
# Usiamo 1000 pseudo-campioni bootstrap.

B <- 1000
y.bar.b <- replicate(B,mean(dati[floor(1 + n * runif(n))])) # Bootstrap non
# parametrico.

hist(y.bar.b,nclass=40)

mean(y.bar.b) - mean(dati) # Distorsione bootstrap (deve essere
# prossima a 0).
var(y.bar.b) # Varianza bootstrap.
sigma^2/n # La varianza bootstrap stima questo valore.

## [b]

sigma <- 10 # Aumentiamo sigma di 10 volte.
dati <- rnorm(n,mu,sigma) # Generiamo i dati.

mean(dati) # La media campionaria � affetta da un'incertezza molto
# maggiore e quindi � molto pi� facile avere una media
# campionaria osservata lontana dalla vera media.

y.bar.b <- replicate(B,mean(dati[floor(1 + n * runif(n))]))

hist(y.bar.b,nclass=40) # La distribuzione bootstrap � pi� dispersa.

mean(y.bar.b) - mean(dati) # Distorsione bootstrap. La vera distorsione
# bootstrap � sempre zero, anche se abbiamo
# aumentato la varianza. Tuttavia, la
# distribuzione di y.bar.b � pi� dispersa,
# e questo si riflette, a parit� di altre
# condizioni (tipicamente, a parit� di B),
# nel fatto che la distorsione calcolata
# approssima peggio il vero valore zero.
var(y.bar.b) # Varianza bootstrap (deve aumentare, visto che sta
# stimando un valore pi� grande).
sigma^2/n # La varianza bootstrap stima questo valore.

## [c]

B <- 10000 # Aumentiamo B di 10 volte, lasciando
# invariato tutto il resto.

y.bar.b <- replicate(B,mean(dati[floor(1 + n * runif(n))]))

mean(y.bar.b) - mean(dati) # Distorsione bootstrap (tendenzialmente,
# aumentando B, la distorsione bootstrap
# sar� pi� prossima al suo vero valore, che
# � zero). Il punto � che, come mean(dati)
# tende a mu per n che tende ad infinito,
# allo stesso modo mean(y.bar.b) tende a
# mean(dati) per B che tende ad infinito.
var(y.bar.b) # Aumentare B non implica una riduzione della varianza
# bootstrap. Implica solo che questa varianza
# approssimer� meglio la "vera" varianza della media
# campionaria nel mondo bootstrap, che � sigma2.hat/n,
# dove sigma2.hat � la varianza della funzione di
# ripartizione empirica, calcolata dando peso 1/n a tutte
# le osservazioni (in altri termini, sigma2.hat � la
# stima di massima verosimiglianza della vera varianza
# sigma^2).

# Per riassumere, l'analisi precedente ci dice che, qualsiasi sia il
# parametro stimato, il comportamento, all'interno del mondo
# bootstrap, di theta.hat.star rispetto a theta.hat, � rappresentato
# tanto meglio 1) quanto pi� piccola � l'incertezza sigma nel mondo
# reale e 2) quanto pi� alto � il numero B di replicazioni bootstrap.


### Esercizio n. 2

## [a]

# I dati sono gli stessi dell'esercizio n. 1.

# La teoria ci dice che la mediana (Med) ha distribuzione
# asintoticamente normale con media uguale alla mediana della
# popolazione e varianza pari a:

# 1/[4*n*f(Med)^2]

# La distorsione, dunque, � asintoticamente nulla. Per stabilire
# l'accuratezza dello stimatore, dobbiamo tuttavia valutare anche
# la varianza: abbiamo il problema che essa dipende dalla funzione
# di densit� f, che � incognita. 
# Stimiamo dunque f con il metodo del nucleo, calcolandone poi il
# valore in corrispondenza della mediana campionaria.

# Generiamo i dati (come per l'esercizio 1).

n <- 50
mu <- 0
sigma <- 1
dati <- rnorm(n,mu,sigma)

med = median(dati) # Mediana campionaria.

library(sm) # Biblioteca per la stima della densit�.

sm.density(dati) # Grafico della densit� stimata.

# Con l'argomento "eval.points" decidiamo dove la funzione
# "sm.density" calcola la stima della densit�. Il valore della stima �
# contenuto nella componente "estimate" restituita dalla funzione.
#
# display="none" serve per evitare che la funzione dia un risultato
# grafico, che qui non ci serve.

f.med <- sm.density(dati, eval.points=med, display="none")$estimate
f.med # Valore che la densit� stimata assume in corrispondenza della
# mediana campionaria.
se.med <- 1/(2 * sqrt(n) * f.med) # Standard error asintotico.
se.med

## [b]

# Soluzione via bootstrap.

B <- 1000
med.bts <- replicate(B,median(dati[floor(1 + n * runif(n))]))

mean(med.bts) - med # Distorsione bootstrap.
sd(med.bts) # Deviazione standard bootstrap.

# Si noti che la distorsione della mediana campionaria e' nulla
# (almeno asintoticamente) indipendentemente dalla densita' f da cui
# provengono i dati (la mediana campionaria e' asintoticamente
# corretta). Al contrario, la deviazione standard della mediana
# campionaria dipende da f. Questo significa che la valutazione della
# distorsione fatta nel mondo bootstrap assomigliera' alla distorsione
# reale indipendentemente dal fatto che F_n assomigli a F oppure no
# (dove indichiamo con F_n e F la funzione di ripartizione empirica e
# la vera funzione di ripartizione, rispettivamente). Non possiamo
# dire la stessa cosa per la valutazione bootstrap della deviazione
# standard, che sara' prossima alla vera deviazione standard della
# mediana campionaria solo se F_n e' prossima a F.

# Questo, naturalmente, vale ogniqualvolta le proprieta' di interesse
# di uno stimatore non dipendono (come la distorsione della media
# campionaria, che e' sempre zero) o dipendono (come la varianza della
# media campionaria, che e' funzione della vera varianza) dalla
# distribuzione da cui provengono i dati.


### Esercizio n. 3

# Carichiamo i dati.

data(faithful)

summary(faithful)
?faithful

## [a]

# A noi interessa studiare le durate delle eruzioni:

dati <- faithful$eruptions
n = length(dati)

med <- median(dati) # Mediana campionaria.

B <- 1000
med.bts <- replicate(B,median(dati[floor(1 + n * runif(n))]))

mean(med.bts) - med # Distorsione bootstrap.
sd(med.bts) # Deviazione standard bootstrap.

hist(med.bts)

# Intervallo di confidenza.

alpha = 0.05
2 * med - quantile(med.bts,c(1-alpha/2,alpha/2))

## [b]

# Procedura per scegliere B.

n.rep <- 30 # Ripetiamo la procedura n.rep volte.
# Aumentando n.rep, dobbiamo aspettare di pi�!
B.try <- c(500,1000,5000,10000) # Valori di B tra cui scegliere.
crit <- rep(NA,length(B.try)) # Criterio di valutazione: deviazione
# standard massima degli estremi degli
# intervalli di confidenza.
for (i in 1:length(B.try)){
  B <- B.try[i] # Assegno B.
  ic <- matrix(NA,n.rep,2) # Calcoler� n.rep intervalli di
  # confidenza.
  for (j in 1:n.rep){
    cat("i = ",i,", j = ",j,", B = ",B,"\n",sep="") # A che punto siamo?
    med.bts <- replicate(B,median(dati[floor(1 + n * runif(n))]))
    ic[j,] <- 2 * med - quantile(med.bts,c(1-alpha/2,alpha/2))
  }
  crit[i] = max(apply(ic,2,sd)) # Prendo il massimo tra le due
  # deviazioni standard dei due estremi
  # degli intervalli di confidenza.
}

min(B.try[crit < 0.005]) # Valore di B considerato sufficiente.


### Esercizio n. 4

data(faithful)
dati <- faithful$eruptions

n = length(dati)

m <- mean(dati) # Media campionaria.

B <- 10000
m.bts <- replicate(B,mean(dati[floor(1 + n * runif(n))]))

mean(abs(m.bts - m) > 0.2) # Stima della probabilit� richiesta.
# La media campionaria, nel mondo
# bootstrap, gioca il ruolo di vera media.


### Esercizio n. 5

impianto <- read.table("impianto.dat",header=TRUE)

x <- impianto[,1] # Temperatura.
y <- impianto[,2] # Perdita di produzione.

plot(x,y) # Sono evidenti dei valori molto lontani da un'ipotetica
# retta di regressione.

library(MASS)
m <- ltsreg(y ~ x) # Stima del modello di regressione robusta.

beta <- coef(m) # Coefficienti di regressione stimati.

abline(beta) # La retta di regressione � evidentemente stata stimata
# ignorando i valori anomali.

n <- length(x)
B <- 1000

# Bootstrap non parametrico.

blts.np <- function(n){ # Funzione che calcola i coefficienti stimati.
  idx <- floor(1 + n * runif(n)) # Bootstrap non parametrico.
  m.bts <- ltsreg(y[idx] ~ x[idx])
  coef(m.bts)
}

beta.bts.np <- sapply(rep(n,B),blts.np) # Bootstrap non parametrico.

dist.np <- apply(beta.bts.np,1,mean) - beta # Calcolo della distorsione.

beta.corr.np <- beta - dist.np # Stime corrette per la distorsione.

abline(beta.corr.np,lty=2) # Retta di regressione corretta.

# Bootstrap semi-parametrico.

res <- residuals(m)

blts.sp <- function(n){ # Funzione che calcola i coefficienti stimati.
  idx <- floor(1 + n * runif(n))
  y.bts <- beta[1] + beta[2] * x + res[idx] # Bootstrap semi-parametrico.
  m.bts <- ltsreg(y.bts ~ x)
  coef(m.bts)
}

beta.bts.sp <- sapply(rep(n,B),blts.sp) # Bootstrap semi-parametrico.

dist.sp <- apply(beta.bts.sp,1,mean) - beta # Calcolo della distorsione.

beta.corr.sp <- beta - dist.sp # Stime corrette per la distorsione.

abline(beta.corr.sp,lty=3) # Retta di regressione corretta.

# Bootstrap parametrico.

# Per il bootstrap parametrico non ha senso assumere la normalit�.
# Una distribuzione spesso usata per rappresentare la presenza di
# valori anomali � la mistura di due normali. In questo caso,
# tuttavia, i valori anomali sono causati da un arresto parziale
# dell'impianto e vanno quindi solo in una direzione: non ha senso
# pensare che provengano da una normale, ossia che possano essere sia
# positivi che negativi. Genereremo i valori anomali da una
# distribuzione esponenziale.

gamma <- 0.14 # Dal grafico, vedo 3 valori anomali su 21. Stimiamo la
# probabilit� di avere un valore anomalo con 3/21, che
# � circa 0.14.

ind <- order(abs(res)) # Indici che servono per ordinare i residui assoluti.
abs(res)[ind] # Residui assoluti ordinati.
res[ind] # Residui semplici, ordinati in base al valore assoluto.
q <- 3 # Numero di valori anomali.
s2.ord <- var(res[ind][1:(n-q)]) # Stimo la varianza "ordinaria", della
# distribuzione normale che genera i
# valori non anomali, usando gli n-q
# residui con valori assoluti pi� piccoli.
m.an <- abs(mean(res[ind][(n-q+1):n])) # Stimo la media "anomala", della
# distribuzione esponenziale che genera i
# valori anomali, usando i q residui
# con valori assoluti pi� grandi.

# Equivalentemente, visto che i 3 residui anomali sono gli unici
# inferiori a -7:

s2.ord <- var(res[res > -7])
m.an <- abs(mean(res[res < -7]))

rres <- function(n){ # Generazione parametrica dei residui.
  v <- rbinom(n, 1, gamma) # v vale 1 con probabilit� gamma e vale 0
  # con probabilit� 1-gamma.
  (1-v) * rnorm(n,0,sqrt(s2.ord)) - v * rexp(n,1/m.an) # I residui
  # anomali sono negativi
  # (causati dall'arresto
  # parziale dell'impianto).
}

y.sim <- beta[1] + beta[2] * x + rres(n) # Verifico graficamente se i
plot(x,y.sim)                            # dati simulati assomigliano
# a quelli reali.

# Si noti che, in questa implementazione, i valori simulati per y
# possono essere anche negativi, mentre non � cos� per i valori di y
# reali. In un'implementazione migliore dovremmo preoccuparcene, ad
# esempio generando gli errori anomali da un'esponenziale troncata in
# beta[1] + beta[2] * x (� facile per inversione).

blts.par <- function(n){ # Funzione che calcola i coefficienti stimati.
  res.bts <- rres(n) # Bootstrap parametrico.
  y.bts <- beta[1] + beta[2] * x + res.bts
  m.bts <- ltsreg(y.bts ~ x)
  coef(m.bts)
}

beta.bts.par <- sapply(rep(n,B),blts.par) # Bootstrap parametrico.

dist.par <- apply(beta.bts.par,1,mean) - beta # Calcolo della distorsione.

beta.corr.par <- beta - dist.par # Stime corrette per la distorsione.

plot(x,y)
abline(beta)
abline(beta.corr.par,lty=4) # Retta di regressione corretta.

dist.np  # I tre metodi danno correzioni diverse. In particolare, il
dist.sp  # metodo non parametrico fornisce correzioni di entit�
dist.par # maggiore. Questo suggerisce che l'assunzione di linearit� �
# da valutare attentamente.


### Esercizio n. 6

dati <- read.table("esponenziale.dat", header=TRUE)
y <- dati$y

lambda.hat = 1/mean(y) # Stima di massima verosimiglianza di lambda.

n = length(y)
B = 1000

lambda.bts <- replicate(B,1/mean(rexp(n,lambda.hat))) # Bootstrap parametrico.

alpha = 0.05

# [a.] Usando la funzione "quantile".

2 * lambda.hat - quantile(lambda.bts,c(1-alpha/2,alpha/2))

# [b.] Senza usare la funzione "quantile".

quantile.mio <- function(x,p){ # x � il vettore di osservazioni.
  # p � il vettore di probabilit� che
  #   caratterizza  i quantili.
  n <- length(x)
  idx <- round(n*p) # La cosa pi� semplice � arrotondare all'intero pi� vicino.
  # La funzione "quantile" lavora in maniera pi� complessa,
  # usando diverse forme di interpolazione.
  sort(x)[idx]
}

2 * lambda.hat - quantile.mio(lambda.bts,c(1-alpha/2,alpha/2))

# I risultati sono leggermente diversi da quelli forniti dalla
# funzione "quantile" (pi� affidabile).