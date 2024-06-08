#LAB 4 PROF



########################################################################################
#################### Esercizio n. 1 ####################################################

d <- read.table("asma.dat", header=TRUE)

d[1:5,]

# Per semplicita' creiamo due vettori x e y.

x <- d[d[,"trattamento"]==1,"risposta"]
y <- d[d[,"trattamento"]==2,"risposta"]

#### [a.]

# Guardiamo i dati.

boxplot(x,y,main="Dati osservati")

summary(x)
summary(y)

# Ambedue le terapie sembrano funzionare, visto che la differenza tra i
# picchi di flusso e' sempre positiva. Dal grafico, il primo
# trattamento sembra dare una migliore risposta in media, ma presenta
# anche una variabilita' superiore. Per di piu', alcuni pazienti
# sembrano reagire pochissimo al primo trattamento. L'effetto del
# secondo trattamento e' piu' stabile.

# Vediamo ora se esistono differenze significative in media.
# Si noti che stiamo usando il test di Welch, senza assumere
# uguaglianza delle varianze.

t.test(x,y)

# La differenza tra le medie delle risposte sembrerebbe non essere
# significativa.

# Il risultato puo' pero' dipendere dalla assunzione di normalita',
# che forse non e' del tutto lecita, dato che

shapiro.test(x)
shapiro.test(y)

#### [b.]

# Proviamo quindi a fare un test bootstrap non parametrico sulla
# differenze tra le due medie e basato sulla stessa statistica test.

# Dobbiamo ragionare sotto H0, quindi le due funzioni di ripartizione
# empiriche da cui genereremo i dati nel mondo bootstrap devono avere
# la stessa media. Senza perdita di generalita' (la funzione test non
# cambia al variare del valore assunto dalla media comune), fissiamo
# la media comune pari a zero.

x0 <- x-mean(x) # Useremo le FdR empiriche di
y0 <- y-mean(y) # x0 e y0.
mean(x0)
mean(y0)

B <- 10000 # Numero di replicazioni bootstrap.

nx <- length(x) # Numerosita' campionarie.
ny <- length(y)

# Funzione test (e' quella del test di Welch):

tstat <- function(x,y) (mean(x)-mean(y))/sqrt(var(x)/nx+var(y)/ny)


# Applichiamo ora la procedura bootstrap.

t.boot <- replicate(B,tstat(x0[floor(runif(nx,1,nx+1))],y0[floor(runif(ny,1,ny+1))]))

# Comando equivalente, basato su "sample" invece che su "runif".

#t.boot <- replicate(B,tstat(sample(x0, replace=TRUE),sample(y0, replace=TRUE)))

hist(t.boot,nclass=50,prob=TRUE) # Rappresentazione grafica della
# distribuzione (bootstrap) della
# funzione test.

t.dati <- tstat(x,y) # Valore osservato della funzione test.
points(t.dati,0,cex=3,pch="*") # Il valore osservato della funzione
# test e' "strano" quando H0 e' vera?

# Confrontiamo la distribuzione bootstrap della funzione test con la
# distribuzione sotto assunzione di normalita' (con la correzione di
# Welch per le varianze possibilmente diverse tra le due popolazioni).

curve(dt(x,31.354),-4,10,col="red",lw=2,add=TRUE) # I gradi di
# liberta' (31.354) li abbiamo
# recuperati dall'output di
# "t.test(x,y)".

# La densita' bootstrap sembra, almeno in parte, differente da quella
# basata sull'assunzione di normalita'.

# Calcoliamo il livello di significativita' osservato basato sulla
# distribuzione bootstrap:

cat("p=",mean(abs(t.boot)>=abs(t.dati)),"\n")

# Piu' correttamente, senza assumere la simmetria attorno a zero della
# distribuzione, il livello di significativita' osservato e':

cat("p=",2*min(mean(t.boot<=t.dati),mean(t.boot>=t.dati)),"\n")

# In ogni caso, i risultati non cambiano: la differenza tra
# le due medie non e' significativa.

##### [c.]

# La conclusione sulla non significativita' della differenza tra le
# medie potrebbe dipendere dalla presenza di alcuni pazienti che
# rispondono male al primo trattamento, ovvero dalla presenza di due
# outlier verso il  basso nel primo gruppo. Possiamo quindi provare a
# confrontare una misura piu' robusta della posizione di una
# distribuzione; ad esempio, la mediana.

mstat <- function(x,y) median(x)-median(y) # Funzione test.

# Mettiamoci sotto H0 (che in questo caso ipotizza uguaglianza delle
# mediane). Senza perdita di generalita', fissiamo le mediane delle
# FdR empiriche pari a zero (un diverso valore non cambia la funzione
# test).

x0 <- x-median(x)
y0 <- y-median(y)
median(x0)
median(y0)

# Applichiamo ora la procedura bootstrap.

m.boot <- replicate(B,mstat(x0[floor(runif(nx,1,nx+1))],y0[floor(runif(ny,1,ny+1))]))

# Comando equivalente, basato su "sample" invece che su "runif".

#m.boot <- replicate(B,mstat(sample(x0, replace=TRUE),sample(y0, replace=TRUE)))

m.dati <- mstat(x,y) # Valore osservato della funzione test.

hist(m.boot,nclass=50,prob=TRUE,xlim=range(m.boot,m.dati))
# Rappresentazione grafica della
# distribuzione (bootstrap)
# della funzione
# test. Specifichiamo "xlim"
# in modo di essere certi di
# poter tracciare m.dati sul grafico.

points(m.dati,0,cex=3,pch="*") # Il valore osservato della funzione
# test e' "strano" quando H0 e' vera?
# In questo caso si'.

# Calcoliamo il livello di significativita' osservato basato sulla
# distribuzione bootstrap:

cat("p=",mean(abs(m.boot)>=abs(m.dati)),"\n")

# Piu' correttamente, senza assumere la simmetria attorno a zero della
# distribuzione, il livello di significativita' osservato e':

cat("p=",2*min(mean(m.boot<=m.dati),mean(m.boot>=m.dati)),"\n")

# Dunque, la differenza tra le due mediane sembra essere altamente
# significativa. In altri termini, i due trattamenti danno una
# risposta significativamente diversa in posizione, se nella
# valutazione si usa un indicatore come la mediana, robusto rispetto
# ai valori anomali.

# Esempio di conclusione con suggerimenti per i medici.
# Il primo farmaco sembra garantire una migliore risposta ad almeno
# il 50% dei pazienti.
# In alcuni rari casi pero' la risposta dei pazienti e'
# ridotta.
# Supponendo che i due farmaci abbiano costi comparabili, i risultati
# suggeriscono di
# 1. usare il primo farmaco come farmaco di prima scelta, ovvero, come
# farmaco elettivo per i nuovi pazienti;
# 2. monitorare la risposta. Se l'aumento del picco di flusso e'
# inferiore a 55, proporre al paziente l'uso del secondo farmaco, che
# diventa quindi il farmaco di seconda scelta.


########################################################################################
#################### Esercizio n. 2 ####################################################

rm(list=ls()) # Cancello tutto cio' che ho in memoria, per non fare
# confusione con variabili precedenti.

# Definiamo i parametri delle due distribuzioni gamma da cui
# genereremo, nelle diverse iterazioni, i dati.

nx = 15
ny = 15

alfa = 0.05

mu.x = 3
mu.y = 3
#mu.y = 5 # Togliere il commento per valutare la potenza.

alpha.x = 2
alpha.y = 2

M = 400 # Numero di iterazioni Monte Carlo.
B = 1000 # Numero di replicazioni bootstrap.

t.pval <- b.pval <- rep(NA,M) # Inizializziamo i vettori che
# conterranno i livelli di
# significativita' osservati (lso).
for (i in 1:M){ # Iterazioni Monte Carlo.
  
  if(i %% 10 == 0) print(i) # Ogni 10 iterazioni stampo il punto a cui
  # sono arrivato.
  
  x = rgamma(nx,alpha.x,alpha.x/mu.x) # Dati x.
  y = rgamma(ny,alpha.y,alpha.y/mu.y) # Dati y.
  
  t.pval[i] <- t.test(x,y)$p.value # lso: test di Welch.
  
  tstat <- function(x,y) (mean(x)-mean(y))/sqrt(var(x)/nx+var(y)/ny) # Funzione test.
  
  x0 <- x-mean(x) # Useremo le FdR empiriche di
  y0 <- y-mean(y) # x0 e y0.
  
  # Bootstrap:
  
  t.boot <- replicate(B,tstat(x0[floor(runif(nx,1,nx+1))],y0[floor(runif(ny,1,ny+1))]))
  t.dati <- tstat(x,y) # Valore osservato della funzione test.
  
  b.pval[i] <- 2*min(mean(t.boot<=t.dati),mean(t.boot>=t.dati)) # lso: bootstrap.
}

mean(t.pval<alfa) # Potenze del test di Welch
mean(b.pval<alfa) # e del test bootstrap.

# I livelli di significativita' dei due test sono, sostanzialmente,
# corretti. Il test bootstrap e' piu' potente.


########################################################################################
#################### Esercizio n. 3 ####################################################

rm(list=ls()) # Cancello tutto cio' che ho in memoria, per non fare
# confusione con variabili precedenti.

# Dati di Newcomb.

x <- c(28, -44, 29, 30, 26, 27, 22, 23, 33, 16, 24, 29, 24, 40 , 21, 31, 34, -2, 25, 19)
hist(x)
abline(v=33.02,lty=2) # Valore attualmente noto per la velocita' della luce.

shapiro.test(x) # I dati non sembrano proprio avere distribuzione normale.


t.test(x,mu=33.02) # Il lso e' un po' piu' grande di 0.01, ma non ci
# possiamo fidare, a causa della non normalita' dei
# dati.

# Creiamo la FdR empirica, sotto H0: mu=33.02.

x0 <- x - mean(x) + 33.02
mean(x0)

B <- 100000
n <- length(x) # Numerosita' campionaria.

# Funzione test:

tstat <- function(x) (mean(x)-33.02)/sqrt(var(x)/n)

# Applichiamo ora la procedura bootstrap.

t.boot <- replicate(B,tstat(x0[floor(1+n*runif(n))]))

hist(t.boot,nclass=50,prob=TRUE) # La distribuzione della funzione
# test e' molto diversa da quella
# sotto assunzione di normalita'.
curve(dt(x,n-1),-4,8,col="red",lwd=2,add=TRUE)

t.dati <- tstat(x) # Valore osservato della funzione test.

points(t.dati,0,cex=3,pch="*") # Il valore osservato della funzione
# test e' poco comune, se H0 e' vera.

2*min(mean(t.boot<t.dati),mean(t.boot>t.dati)) # Rifiutiamo H0. I dati
# di Newcomb non sembrano
# essere, in media,
# compatibili con la velocita'
# della luce attualmente nota.


# Costruiamo ora un test bootstrap per verificare che la mediana della
# popolazione da cui provengono i dati di Newcomb e' 33.02.

boxplot(x)
abline(h=33.02,lty=2) # Valore attualmente noto per la velocita' della luce.

# FdR sotto H0: mediana=33.02.

x0 <- x - median(x) + 33.02
median(x0)

med.stat <- function(x) median(x)-33.02 # Funzione test.

med.boot <- replicate(B,med.stat(x0[floor(1+n*runif(n))])) # Bootstrap.

med.stat.dati <- med.stat(x) # Funzione test osservata.

hist(med.boot,nclass=50,prob=TRUE)
points(med.stat.dati,0,cex=3,pch="*")

2*min(mean(med.boot<med.stat.dati),mean(med.boot>med.stat.dati))
# Il lso, pur essendo piu' grande
# di quello ottenuto per la
# media, porta comunque a
# rifiutare H0.


########################################################################################
#################### Esercizio n. 4 ####################################################

rm(list=ls()) # Cancello tutto cio' che ho in memoria, per non fare
# confusione con variabili precedenti.

x <- c(28, -44, 29, 30, 26, 27, 22, 23, 33, 16, 24, 29, 24, 40 , 21, 31, 34, -2, 25, 19)

n <- length(x) # Numerosita' campionaria.
alfa = 0.01

# Test t di Student.

mean(x) + c(-1,1) * sd(x)/sqrt(n) * qt(1-alfa/2,n-1) # 33.02 e' nell'i.c.
t.test(x,mu=33.02,conf.level=1-alfa) # Equivalente.

# Intervallo di confidenza bootstrap per la media.

B <- 100000
m.boot <- replicate(B,mean(x[floor(1+n*runif(n))]))
hist(m.boot)

2*mean(x) - quantile(m.boot,c(1-alfa/2,alfa/2)) # 33.02 e' nell'intervallo!

# Quindi, basandoci su un intervallo di confidenza, anche il bootstrap
# suggerisce, come il test t di Student, di accettare H0: mu=33.02,
# almeno al livello alfa=0.01.
# Come mai si verifica questa discrepanza tra il risultato ottenuto
# con il test bootstrap sulla media ed il risultato ottenuto con
# l'intervallo di confidenza bootstrap per la media?
# Il fatto e' che la funzione test equivalente a questo intervallo non
# e' la funzione tipo t di Student che abbiamo usato nell'esercizio
# precedente, ma e' semplicemente la differenza tra la media
# campionaria ed il valore 33.02 ipotizzato (senza alcuna
# standardizzazione).
# Implementiamo un test bootstrap basato su questa funzione test:

tstat1 <- function(x) mean(x)-33.02 # Funzione test.
x0 <- x-mean(x)+33.02
t.boot1 <- replicate(B,tstat1(x0[floor(1+n*runif(n))])) # Bootstrap.
2*min(mean(t.boot1<tstat1(x)),mean(t.boot1>tstat1(x)))

# Come si vede, in questo caso il test porta ad accettare H0, al
# livello alfa=0.01.
# Quest'ultimo test e' meno affidabile, in quanto la standardizzazione
# applicata nella statistica t di Student permette, in qualche modo,
# di "accorgersi" della presenza di valori anomali e di valutare
# meglio gli allontanamenti dal valore 33.02.
# Gli intervalli di confidenza bootstrap studentizzati, analoghi al
# test t di Student, in questo caso sono preferibili: provare ad applicarli.

# Intervallo di confidenza bootstrap per la mediana.

med.boot <- replicate(B,median(x[floor(1+n*runif(n))]))
2*median(x) - quantile(med.boot,c(1-alfa/2,alfa/2)) # 33.02 non e' nell'intervallo.

# L'insegnamento che traiamo dall'esercizio e' che gli intervalli di
# confidenza bootstrap non studentizzati vanno valutati con cautela,
# specialmente se si basano su stimatori non robusti e siamo in
# presenza di valori anomali. Il problema e' meno grave se si
# considerano stimatori robusti, come la mediana.


########################################################################################
#################### Esercizio n. 5 ####################################################

rm(list=ls()) # Cancello tutto cio' che ho in memoria, per non fare
# confusione con variabili precedenti.

houseprices <- read.table("houseprices.dat",header=TRUE)

x <- houseprices$area
y <- houseprices$sale.price

plot(x,y)
mod <- lm(y~x)
abline(coef(mod))
mod.s <- summary(mod)

mod.s

n <- length(x)

# Calcoliamo un intervallo di confidenza per beta2, assumendo la normalita'.

alfa <- 0.05
coef(mod.s)[2,1] + c(-1,1) * coef(mod.s)[2,2] * qt(1-alfa/2,n-2)

# L'intervallo include 0.4, ma non ci possiamo fidare. Infatti:

shapiro.test(rstandard(mod)) # Rifiutiamo l'ipotesi di normalita'.

B <- 20000

# Bootstrap non parametrico.

b.np <- function(n){ # Funzione che calcola il coefficiente stimato.
  idx <- floor(1 + n * runif(n)) # Bootstrap non parametrico.
  m.bts <- lm(y[idx] ~ x[idx])
  coef(m.bts)[2]
}

b.bts.np <- sapply(rep(n,B),b.np)

2*coef(mod)[2] - quantile(b.bts.np,c(1-alfa/2,alfa/2)) # L'intervallo
# non include 0.4.

# Bootstrap parametrico.

b.par <- function(n){ # Funzione che calcola il coefficiente stimato.
  yy <- mod.s$coef[1] + mod.s$coef[2] * x + mod.s$sigma * rnorm(n) # Bootstrap parametrico.
  m.bts <- lm(yy ~ x)
  coef(m.bts)[2]
}

b.bts.par <- sapply(rep(n,B),b.par)

2*coef(mod)[2] - quantile(b.bts.par,c(1-alfa/2,alfa/2)) # L'intervallo
# include 0.4.

# Bootstrap semiparametrico.

b.sp <- function(n){ # Funzione che calcola il coefficiente stimato.
  idx <- floor(1 + n * runif(n)) # Bootstrap semi-parametrico.
  yy <- mod.s$coef[1] + mod.s$coef[2] * x + mod.s$res[idx]
  m.bts <- lm(yy ~ x)
  coef(m.bts)[2]
}

b.bts.sp <- sapply(rep(n,B),b.sp)

2*coef(mod)[2] - quantile(b.bts.sp,c(1-alfa/2,alfa/2)) # L'intervallo
# include 0.4.

# Conclusioni. La distribuzione dei disturbi nel modello di
# regressione lineare semplice non sembra proprio gaussiana. La
# numerosita' campionaria non e' cosi' elevata da poter ricorrere in
# tutta tranquillita' a risultati asintotici. Abbiamo dunque applicato
# tecniche bootstrap. L'ipotesi nulla di interesse viene accettata se
# si assume che la relazione tra le variabili sia lineare. Altrimenti,
# viene rifiutata.


########################################################################################
#################### Esercizio n. 6 ####################################################

rm(list=ls()) # Cancello tutto cio' che ho in memoria, per non fare
# confusione con variabili precedenti.

n <- 100   # Numerosita' campionaria.
B <- 1000 # Numero di campioni bootstrap.
M <- 200 # Numero di replicazioni Monte Carlo.
alfa <- 0.3 # Probabilita' che definisce il quantile da studiare.

# Bootstrap non parametrico con FdR empirica.

q.boot.MC.emp <- rep(NA,M)
for (i in 1:M){ # Iterazioni Monte Carlo.
  
  if(i %% 10 == 0) print(i) # A che punto siamo?
  
  y <- rnorm(n) # Generiamo i dati via Monte Carlo.
  q.stima <- quantile(y,alfa) # Stima del quantile: gioca il ruolo di
  # vero quantile nel mondo bootstrap in
  # cui la vera distribuzione e' la FdR empirica.
  q.boot <- sapply(rep(n,B), function(n) quantile(y[1+n*runif(n)], alfa))
  dist <- mean(q.boot) - q.stima
  q.boot.MC.emp[i] <- q.stima - dist # Stima del quantile,
  # corretta con la stima
  # bootstrap della distorsione.
}

eqm1=sqrt(mean((q.boot.MC.emp-qnorm(alfa))^2)) # Radice del'errore
# quadratico medio. Usiamo
# il fatto, solo in questa
# istruzione, che sappiamo
# che il vero quantile e'
# "qnorm(alfa)".
eqm1

rnucleo <- function(n, y, h) { # Generazione da una densita' stimata con
  # il metodo del nucleo, usando un
  # nucleo N(0,1).
  ny <- length(y)
  idx <- floor(1 + ny * runif(n))
  y[idx] + h * rnorm(n)
}

library(sm) # Ci serve la biblioteca "sm" per "h.select".

q.boot.MC.nuc <- rep(NA,M)
for (i in 1:M){
  
  if(i %% 10 == 0) print(i) # A che punto siamo?
  
  y <- rnorm(n) # Generiamo i dati via Monte Carlo.
  q.stima <- quantile(y,alfa) # Questo e' il vero quantile nel mondo bootstrap?
  h <- h.select(y,method="sj")
  yb <- matrix(rnucleo(n * B, y, h), n, B)
  q.boot <- apply(yb, 2, function(x) quantile(x, alfa))
  dist <- mean(q.boot) - q.stima
  q.boot.MC.nuc[i] <- q.stima - dist
}

eqm2=sqrt(mean((q.boot.MC.nuc-qnorm(alfa))^2))
eqm2
# L'EQM, almeno nelle condizioni studiate, peggiora.
# In queste condizioni quindi il metodo del nucleo non e' utile,
# almeno dal punto di vista della stima puntuale.
# Il metodo del nucleo richiede di stimare il parametro di lisciamento
# e questo puo' creare incertezza.
# Quanto estremo e' il quantile che stiamo cercando di stimare (i
# quantili per alfa piu' piccoli sono piu' difficili da stimare),
# assieme alla numerosita' campionaria, giocano sicuramente un ruolo
# importante. Anche la forma della vera distribuzione dei dati e'
# importante.

# Notiamo, tuttavia, quanto segue.

# "quantile(y,alfa)" e', sostanzialmente (viene applicata una qualche
# interpolazione), basato sulla FdR empirica.
# Quindi, se i dati bootstrap provengono dalla densita' stimata
# con il metodo del nucleo, "quantile(y,alfa)" non e' il vero quantile
# alfa del mondo bootstrap. Invece di usare "quantile(y,alfa)" come
# riferimento per stimare la distorsione, come abbiamo appena fatto,
# proviamo a calcolare il quantile alfa basato sulla densita' stimata
# con il metodo del nucleo.
# Purtroppo, la FdR basata sul metodo del nucleo non e'
# invertibile. Dobbiamo applicare una procedura di ottimizzazione.

pnucleo <- function(x, y, h) { # FdR basata sul metodo del nucleo.
  a <- outer(x, y, "-")/h
  apply(pnorm(a), 1, mean)
}

# La funzione qnucleo minimizza, rispetto ad x, |pnucleo(x, y, h) - alfa|.

qnucleo <- function(alfa, y, h){
  optimize(function(x,y,h,alfa) abs(pnucleo(x,y,h)-alfa),
           y=y,h=h,alfa=alfa,
           interval=c(-4,4),tol=1E-5)$minimum
}
pnucleo(qnucleo(alfa,y,h),y,h) # Proviamo se funziona: deve restituire alfa.

q.boot.MC.nuc1 <- rep(NA,M)
for (i in 1:M){
  
  if(i %% 10 == 0) print(i) # A che punto siamo?
  
  y <- rnorm(n) # Generiamo i dati via Monte Carlo.
  q.stima <- qnucleo(alfa,y,h) # Cambiamo la definizione di vero
  # quantile nel mondo bootstrap.
  # Il resto della procedura non cambia.
  h <- h.select(y,method="sj")
  yb <- matrix(rnucleo(n * B, y, h), n, B)
  q.boot <- apply(yb, 2, function(x) quantile(x, alfa))
  dist <- mean(q.boot) - q.stima
  q.boot.MC.nuc1[i] <- q.stima - dist
}

eqm3=sqrt(mean((q.boot.MC.nuc1-qnorm(alfa))^2)) # L'EQM diminuisce.
eqm3

c(eqm1,eqm2,eqm3) # I tre EQM sotto radice.

# N.B.: per diversi valori di n e alfa, o per un'altra distribuzione
# dei dati, le cose potrebbero andare diversamente.