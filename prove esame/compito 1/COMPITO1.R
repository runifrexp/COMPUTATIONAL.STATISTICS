#COMPITO 1 

#-------------------------- ESERCIZIO 1

rm(list = ls())
###parte a

f <- function(x){
  10 * exp(-x*10 *(0.5 +abs(sin(x*10))))
}

B <- 5000

mc <- rep(NA, B)
K <- 10
u <- runif(B)
u
f(u)
I <- mean(f(u))
I

se <- sqrt(var(f(u))/length(u))

ic <- mean(f(u)) + c(-1, 1) * qnorm(0.975) * se
  
ic



###parte b

dens <- function(x) exp(-x)

f <- function(x){
  exp(-x *(0.5 +abs(sin(x))))
}


x <- rexp(B)
z <- f(x)/dexp(x)
z
I <- mean(z)
I



###parte c

dens <- function(x) (1/I) * f(x)

r.acc.rif <- function(n, f, g, rg, k){
  y <- rep(NA, n)
  for(i in 1:n){
    done <- F
    while(!done){
      z <- rg(1)
      u <- k * g(z) *runif(1)
      if (u <= f(z)){
        done <- T
      }
    }
    y[i] <- z
  }
  y
}


g <- function(x) 0.5 *exp(-0.5*x)

h <- function(x) f(x)/g(x)

curve(h, from = -1, to = 10)
curve(k*g(x), from = -1, to = 10, col = 'darkgreen')
curve(f(x), from = -1, to = 10, add = T,  col = 'magenta')
#cosi' capisco da dove posso far partire l'ottimizzatore per ottenere k, sembra 
#essere un valore prossimo allo zero 

k <- optimize(f = h, lower = 0, upper = 1, maximum = T)$objective
k

rg <- function(n) (rexp(n, 0.5))

n <- 1000
vals <- r.acc.rif(n, f, g, rg, k)
hist(vals, nclass = 100, prob = T)
curve(f(x), from = -1, to = 10, add = T,  col = 'magenta')




#-------------------------- ESERCIZIO 2

rm(list = ls())

dati <- read.table(file.choose(), header = T)
dati

x <- dati$x
x

y <- dati$y
y

###parte a

#repliche bootstrap non parametrico: non assumo ne' normalita' ne' linearita' 
#dei parametri 

mod.vero <- lm(y~x)
summary(mod.vero)
plot(x,y)
abline(mod.vero)
n <- length(x)
n
B <- 10000

b0.bts <- function(n){
  idx <- (floor(1+n*runif(n)))
  m.b <- lm(y[idx]~x[idx])
  coef(m.b)[1]
}

b0.b <- sapply(rep(n, B), b0.bts)
b0.bts <- mean(b0.b) #beta0 bootstrap
b0.bts

b1.bts <- function(n){
  idx <- (floor(1+n*runif(n)))
  m.b <- lm(y[idx]~x[idx])
  coef(m.b)[2]
}

b1.b <- sapply(rep(n, B), b1.bts)
b1.bts <- mean(b1.b) #beta1 bootstrap
b1.bts


plot(x,y)
abline(mod.vero)
abline(b0.bts - (b0.bts - coef(mod.vero)[1]), b1.bts -(b1.bts - 
                                          coef(mod.vero)[2]) , col = 'magenta')


###parte b

#distorsione 
#b0
b0.bts - coef(mod.vero)[1]

#b1
b1.bts - coef(mod.vero)[2]

#varianza
#b0
var(b0.b)
#b1
var(b1.b)

hist(b0.b, nclass = 100, prob = T) #distribuzione di b0
curve(dnorm(x, b0.bts, sqrt(var(b0.b))), add = T, col = 'magenta', lwd = 2)
hist(b1.b, nclass = 100, prob = T) #distribuzione di b1
curve(dnorm(x, b1.bts, sqrt(var(b1.b))), add = T, col = 'gold', lwd = 2)



###punto c

alpha <- 0.01

vero.b1 <- coef(mod.vero)[2]
ic <- 2 * vero.b1 - quantile(b1.b, c(1-alpha/2, alpha/2))
ic

alpha2 <- 0.1
ic2 <- 2 * vero.b1 - quantile(b.bts[,2], c(1-alpha2/2, alpha2/2))
ic2

#si nota che l'intervallo si stringe in quanto vengono presi valori piu' piccoli
#da quantile in quanto alpha e' piu' piccolo, ed essendo che questi quantili 
#vengono sottratti a due volte la stima nel 'mondo reale', ossia il vero 
#parametro per il mondo bootstrap, l'intervallo risulta essere piu' stretto
#



#-------------------------- ESERCIZIO 3

rm(list = ls())
dati <- read.table(file.choose(), header = T)
dati


###parte a
dens <- function(par, x){
  a <- par[1]
  b <- par[2]
  (gamma(a+b)/(gamma(a)*gamma(b)))*(x^(a-1))*((1-x)^(b-1))
}

llik <- function(par, x){
  sum(log(dens(par, x)))
}


###parte b
  
a.seq <- seq(0.1,10, by = 0.1) 
la <- length(a.seq)

b.seq <- seq(0.1,10, by = 0.1) 
lb <- length(b.seq)
  
zz <- matrix(NA, la, lb)
dim(zz)
for(i in 1:la){
  for(j in 1:lb){
    zz[i,j] <- llik(c(a.seq[i], b.seq[j]), dati)
  }
}
zz
contour(a.seq, b.seq, zz, levels = quantile(zz, seq(0.8, 1, by = 0.02)))
#il grafico della verosimiglianza evidenzia che le SMV per a e b sembrano essere 
#rispettivamente per a 2, mentre per b 4



###parte c

start <- log(c(2,4))
start

nllik.rip <- function(par, x){
  - llik(exp(par), x)
}

smv.da.riparaetrizzare <- nlm(f = nllik.rip, x = dati, p = start)$estimate
smv.da.riparaetrizzare

SMV <- exp(smv.da.riparaetrizzare)
SMV



###parte d

library(nlme)
H <- fdHess(SMV, function(par) llik(par, dati))
H
vcov <- solve(-H$Hessian)
stde <- sqrt(diag(vcov))

alpha <- 0.05
cbind(SMV - qnorm(1-alpha/2) *stde, SMV + qnorm(1-alpha/2) *stde)

