#LABORATORIO 4 



#-------------------------------------ES 1--------------------------------------

dati <- read.table(file.choose(), header = T) #asma.dat
head(dati)
attach(dati)



#### [a.]
x <- dati[dati[, 'trattamento'] == 1, 'risposta']
y <- dati[dati[, 'trattamento'] == 2, 'risposta']
boxplot(x,y,main="Dati osservati")
#da un'analisi grafica preliminare sembra che le medie siano diverse.
#sembra che il trattamento 1 sia migliore

t.test(x,y)
#noto che il pvalue e' grande, quindi accetto H_0 ipotesi di uguaglianza delle 
#medie di efficacia dei due trattamenti diversi (anche se dal boxplot ho il 
#sospetto che non sia aeffettivamente cosi')

#verifico la normalita', su cui si basa il risultato del t test
shapiro.test(x)
shapiro.test(y)
#noto che i pvalue sono molto bassi quindi devo rifiutare H_0, percio' non 
#posso assumere normalita' nei dati
#Quindi anche il test t a due campioni non sara' molto affidabile... 



#### [b.]
#faccio test bootstrap NP per verificare la stessa ipotesi, facendo parlare 
#solamente i dati

#funzione per generare campioni bootsrap NP
r.trat <- function(x, y, B = 10000){
  nx <- length(x)
  ny <- length(y)
  T <- rep(NA, B)
  for (i in 1:B){
    xb <- x[floor(1 + nx * runif(nx))]
    yb <- y[floor(1 + ny * runif(ny))]
    T[i] <- (mean(xb) - mean(yb))/sqrt(var(xb)/nx+var(yb)/ny)
  }
  T
}


t.boot <- r.trat(x - mean(x), y - mean(y)) #mettendo questi dati in input non 
#faccio l'errore di assumere indipenenza tra le due variabili, in questo modo 
#trovo una soluzione all'indicizzazione diversa portata dalla numerosita' 
#diversa di x e y
#inoltre (cosa piu' importante), con questo input mi metto sotto H_0

hist(t.boot, nclass = 50, prob = T)
curve(dt(x,31.354),-4,10,col="red",lw=2,add=TRUE) #df trovati da t.test(x,y)
t.dati <- tstat(x,y) # Valore osservato della funzione test.
points(t.dati,0,cex=3,pch="*") 

nx <- length(x)
ny <- length(y)
t.dati <- (mean(x) - mean(y)) / sqrt(var(x)/nx+var(y)/ny) #stat. test dei dati

#alpha oss:
mean(abs(t.boot) >= abs(t.dati))
#accetto H_0: le due medie non hanno una differenza significativa, quindi il 
#test t a due campioni, nonostante la non normalita' dei dati, dava un 
#risultato corretto


#### [c.]
#la validita' dei test bootstrap con ipotesi sulla media puo' risultare scarsa 
#per via di possibili outliers (vedi boxplot), puo' essere quindi utile fare dei
#test sulla mediana 
mstat <- function(x,y) median(x)-median(y) # Funzione test.

x0 <- x-median(x)
y0 <- y-median(y)
median(x0)
median(y0)
#ok, siamo sotto H_0, ricorda che non si perde generalita' se si pone la mediana
#uguale a zero, si poteva eguagliare a qualsiasi valore, pero' zero e' il piu'
#pratico 

#creo il campione bootstrap
B <- 10000
m.boot <- replicate(B, mstat(sample(x0, replace=TRUE),sample(y0, replace=TRUE)))

hist(m.boot, nclass = 50, prob = T)

m.dati <- mstat(x,y)

2*min(mean(m.boot<=m.dati, mean(m.boot>=m.dati)))

cat("p=",2*min(mean(m.boot<=m.dati),mean(m.boot>=m.dati)),"\n")
#e' zero, quindi rifiuto H_0.
#cio' significa che la differenza tra le due mediane e' significativa



#-------------------------------------ES 2--------------------------------------

#### [a.]
rm(list = ls())

x <- rgamma(15, 2, 2/3)
y <- rgamma(15, 2, 2/3)

nx <- length(x)
ny <- length(y)

mu.x = mu.y = 3

alpha.x = alpha.y = 2

t.test(x,y, var.equal = F)
#non assumo omoschedasticita', test di Welch

#pvalue grande quindi si pensa di poter accettare.
shapiro.test(x)
shapiro.test(y)
#la normalita' non viene rispettata ne' da x ne' da y, quindi il test t di 
#di student a due campioni NON ha un esito affidabile

welch.MC <- function(M = 400){
  v <- rep(NA, M)
  for (i in 1:M){
    x <- rgamma(15, alpha.x, alpha.x/mu.x)
    y <- rgamma(15, alpha.y, alpha.y/mu.y)
    v[i] <- t.test(x,y)$p.value #alphaoss del test di Welch
    cat(v[i], '\n')
  }
  v
}
MC.pvalue <- welch.MC()

#alphaoss MC
mean(MC.pvalue)
#no, l' alphaoss e' distante dall' alpha dichiarato (=0.05)



#### [b.]

#monte carlo sul bootstrap


#stima bootstrap NP delle stat. test, ne creo un
welch.bts <- function(x,y, B = 1000) {
  nx <- length(x)
  ny <- length(y)
  test <- rep(NA, B)
  for (i in 1:B){
    xb <- x[floor(1+nx*runif(nx))]
    yb <- y[floor(1+ny*runif(ny))]
    test[i] <- (mean(xb)-mean(yb))/sqrt((var(xb)/nx)+(var(yb)/ny))
  }
  test
} 

welch.bts(x,y)

hist(welch.bts(x,y), nclass = 100)
#Mi accorgo che la media non risiede in zero, come postulerebbe l'ipotesi nulla
#È imperativo intervenire con solerzia, modificare le circostanze per poter 
#continuare a vagare nel vasto mondo del bootstrap. Altrimenti, se non rispetto 
#le loro tre leggi auree, sarò gettato nell'oblio. 
#Presto! Devo agire prima che mi smascherino...
welch.bts(x-mean(x),y-mean(y))
hist(welch.bts(x-mean(x),y-mean(y)), nclass = 100)
#fiu, per fortuna sono stato lesto e nessuno ha notato che non eravamo sotto
#l'egida dell'ipotesi nulla. Ora, per avere piena fiducia nelle mie 
#affermazioni, mi concedo una sontuosa escursione a Monte Carlo, affinché possa 
#ripetere l'esperimento innumerevoli (non proprio innumerevoli in quanto non
#dispongo di Leonardo) volte, ottenendo una certezza granitica che un singolo 
#risultato non sia semplicemente frutto del capriccio del caso.
welch.MC.bts <- function(M = 400){
  v <- rep(NA, M)
  for (i in 1:M){
    xb <- rgamma(15, alpha.x, alpha.x/mu.x)
    yb <- rgamma(15, alpha.y, alpha.y/mu.y)
    st.test <- (mean(xb)-mean(yb))/sqrt((var(xb)/nx)+(var(yb)/ny))
    t.boot <- welch.bts(x-mean(x),y-mean(y))
    v[i] <- mean(abs(t.boot)>abs(st.test))
  }
  v
}

b.pvalue <- welch.MC.bts() 
mean(b.pvalue) #alphaoss con MC su bootstrap



#### [c.]

alpha <- 0.05
mu.x <- 3
mu.y <- 5
mean(MC.pvalue < alpha) # Potenze del test di Welch
mean(b.pvalue < alpha)




#-------------------------------------ES 3--------------------------------------

#### [a.]

rm(list=ls())
x <- c(28, -44, 29, 30, 26, 27, 22, 23, 33, 16, 24, 29, 24, 40 , 21, 31, 34, 
       -2, 25, 19)
n <- length(x)

#verifico la normalita'
shapiro.test(x)
#alpha oss piccolo: rifiuto H_0: il campione non sembra avere 
#distribuzione normale

#test media singolo campione
t.test(x,mu = 33.02) 
#dato che i dati non sono normali, non mi posso fidare dell'esito del test.
#oltretutto alphaoss non e' un valore che ci permette di accettare o rifiutare 
#H_0 con sicurezza



#### [b.]

#provo con il bootstrap: creo il mondo

#stessa FdR della teorica sotto H_0: mu = 33.02
x0 <- x - mean(x) + 33.02
mean(x0) #ok, siamo sotto H_0

tstat <- function(x) (mean(x)-33.02)/sqrt(var(x)/n) #statistica test

#statistica test del bootstrap
B <- 100000
t.boot <- replicate(B, tstat(x0[floor(1+n*runif(n))]))

#statistica test del mondo reale
t.dati <- tstat(x)

#alphaoss trovato via bootstrap
2*min(mean(t.boot<t.dati),mean(t.boot>t.dati))
#non posso accettare H_0: la vera media non e' 33.02





#### [c.]

#dato che ci sono certi valori anomali si prova a usare la mediana in quanto 
#stimatore piu' robusto

#creo il mondo bootrsap: STAI SOTTO H_0! Devi starci...

x0 <- x - median(x) + 33.02
median(x0)

#definisco funzione mediana utile per il bootstrap
med <- function(x) median(x)-33.02

med.rw <- med(x) #rw = real world
med.rw

#mediana bootstrap
med.bts <- replicate(B, med(x[floor(1+n*runif(n))]))
med.bts

hist(med.bts, nclass = 50, prob = T)
points(med.rw,0,cex=3,pch="*")

2*min(mean(med.bts<med.rw),mean(med.bts>med.rw))
#accetto H_0 con mooooolta forza, la mediana e' praticamente 33.02



#-------------------------------------ES 4--------------------------------------


rm(list=ls())

x <- c(28, -44, 29, 30, 26, 27, 22, 23, 33, 16, 
       24, 29, 24, 40 , 21, 31, 34, -2, 25, 19)

n <- length(x) # Numerosita' campionaria.
alpha = 0.01




#### [a.]

shapiro.test(x)
#alpha oss piccolo: rifiuto la normalita' dei dati

#test media singolo campione
t.test(x,mu = 33.02) 
#dato che i dati non sono normali non posso fidarmi della risposta di questo
#test, inoltre se mi dovessi fidare 
t.test(x,mu = 33.02)$p.value
#e' un valore piccolo, quindi non accetterei H_0

#costruisco degli ic per verificare se 33.02 rientra negli intevalli

#ic 'classico'

mu.cap <- mean(x)

s2 <- var(x)/n

ic <- mu.cap +  sqrt(s2) * qt(c(alpha/2, 1-alpha/2), n-1)
ic

t.test(x,mu=33.02,conf.level=1-alpha) # Equivalente.

#mu0 rientra nell'ic, quindi accettiamo H_0



#via bootstrap

B <- 10000

#ic con theta.cap.star
x0 <- x - mean(x) +33.02 #ora siamo sotto H_0

mu.cap.star <- replicate(B, mean(x[floor(1+n*runif(n))]))

hist(mu.cap.star, nclass = 100)

points(mu.cap, pch = '*', cex = 3)

ic <- 2 * mean(x) - quantile(mu.cap.star, c(1-alpha/2, alpha/2))
ic
#mu0 rientra nell'ic anche nel bootstrap, quindi accettiamo H_0




#### [b.]

#ic basato su errore di stima

tstat1 <- function(x) mean(x)-33.02 # Funzione test.

x0 <- x - mean(x) + 33.02 #campione che rispetta le regole del mondo
#bootatrap

T.bts <- replicate(B, tstat1(x0[floor(1+ n *runif(n))])) #distribuzione
#bootstrap dell'errore di stima

hist(T.bts, nclass = 50)

points(mean(x0), pch = '*', cex = 3)

2*min(mean(T.bts<tstat1(x)),mean(T.bts>tstat1(x)))

#per la mediana 
med.boot <- replicate(B,median(x[floor(1+n*runif(n))]))
2*median(x) - quantile(med.boot,c(1-alfa/2,alfa/2))
#33.02 non e' nell'intervallo: rifiuto H_0








#-------------------------------------ES 5--------------------------------------

rm(list = ls())

prices <- read.table(file.choose(), header = T)
head(prices)

x <- prices[,1] #area, v. exp
x

y <- prices[,2]
y

n <- length(x)

m <- lm(y~x)

plot(y~x)
abline(coef(m))

summary(m)


shapiro.test(res)
#non si puo' accettare l'ipotesi di normalita' nei residui (errori)

xbar <- mean(x)
ybar <- mean(y)

b2.cap <- sum((x-xbar)*(y-ybar))/sum((x-xbar)^2)

res <- residuals(m)

s2 <- sum(res^2)/(n-2)

stima.var.b2 <- s2/sum((x-mean(x))^2)

alpha <- 0.05

ic <- b2.cap + c(-1,1) * sqrt(stima.var.b2) * qt(1-alpha/2, n-2)
ic


#### [b.]

#@@@ic NP

#ic su distribuzione della stima di theta.star

B <- 20000

b.np <- function(n){ # Funzione che calcola il coefficiente stimato.
  idx <- floor(1 + n * runif(n)) # Bootstrap non parametrico.
  m.bts <- lm(y[idx] ~ x[idx])
  coef(m.bts)[2]
}

b.bts.np <- sapply(rep(n,B),b.np)

2*coef(mod)[2] - quantile(b.bts.np,c(1-alfa/2,alfa/2)) 
#ic NON include 0.4


#@@@ic P

#ricampiono solo dalle y, le x le tengo fissate
b.par <- function(n){ # Funzione che calcola il coefficiente stimato.
  yy <- mod.s$coef[1] + mod.s$coef[2] * x + mod.s$sigma * rnorm(n) # Bootstrap parametrico.
  m.bts <- lm(yy ~ x)
  coef(m.bts)[2] #estraggo b2
}

b.bts.par <- replicate(B, b.par(n))
#oppure
b.bts.par <- sapply(rep(n,B),b.par)

2*coef(mod)[2] - quantile(b.bts.par,c(1-alfa/2,alfa/2)) 
#ic include 0.4


#@@@ ic SP

b.sp <- function(n){ # Funzione che calcola il coefficiente stimato.
  idx <- floor(1 + n * runif(n)) # Bootstrap semi-parametrico.
  yy <- mod.s$coef[1] + mod.s$coef[2] * x + mod.s$res[idx]
  m.bts <- lm(yy ~ x)
  coef(m.bts)[2]
}

b.bts.sp <- sapply(rep(n,B),b.sp)

#oppure
b.bts.sp <- replicate(B, b.sp(n))

2*coef(mod)[2] - quantile(b.bts.sp,c(1-alfa/2,alfa/2))
#ic contiene 0.4


