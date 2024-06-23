#COMPITO 3

rm(list = ls())
# -------------------------- ESERCIZIO 1

###parte a
amianto <- read.table(file.choose(), header = T, sep = ',')
amianto
dim(amianto)

anfibolo <- amianto[amianto$X == 1, 1]
anfibolo
length(anfibolo)

crisotilo <- amianto[amianto$X == 0, 1]
crisotilo
length(crisotilo)

t.test(anfibolo, crisotilo)
#risultato incerto, non sono sicuro se accettare o non accettare

#verifico la normalita' con shapiro.test
shapiro.test(anfibolo)
shapiro.test(crisotilo)
#alphaoss praticamente zero quindi rifiuto la normalita' dei dati quindi non
#mi posso fidare dei risultati ottenuti con il test t a due campioni



###parte b 

#trovo la 'vera' mediana del mondo bootstrap, ossia la stima della mediana dai dati
ma <- median(anfibolo)
ma
mc <- median(crisotilo)
mc

B <- 5000
na <- length(anfibolo)
nc <- length(crisotilo)
a.bts <- rep(NA, B*na)
c.bts <- rep(NA, B*nc)

m.t.bts <- rep(NA, B)

a0 <- anfibolo - median(anfibolo)
c0 <- crisotilo - median(crisotilo)

m.test <- function(x, y){
  median(x) - median(y)
}

median.dati <- m.test(anfibolo, crisotilo)



for (i in 1:B){
  idxa <- floor(1+na*runif(na))
  idxc <- floor(1+nc*runif(nc))
  m.t.bts[i] <- m.test(a0[idxa], c0[idxc])
}

2 * min(mean(m.t.bts<= median.dati), mean(m.t.bts>= median.dati))



hist(m.t.bts, nclass = 100, prob = T)



###parte c

ic <- 






# -------------------------- ESERCIZIO 1
rm(list=ls())

cardio <- read.table(file.choose(), header = T)
cardio

x <- cardio$age
y <- cardio$y
n <- cardio$n

pi <- function(beta, x){
  b0 <- beta[1]
  b1 <- beta[2]
  (exp(b0+b1*x))/(1+exp(b0+b1*x))
}

bin <- function(n, y){
  factorial(n)/(factorial(y)*factorial(n-y))
}

d <- function(x, n, y, beta){
  bin(n,y) * (pi(beta, x)^y) * (1-pi(beta, x))^(n-y)
}

llik <- function(beta, x, n, y){
  sum(log(d(x, n, y, beta)))
}



#punto b

b0.seq <- seq(-10, 0, by = 0.1)
b1.seq <- seq(0, 0.2, by = 0.002)
n0 <- length(b0.seq)
n0
n1 <- length(b1.seq)
n1

zz <- matrix(NA, n0, n1)
dim(zz)

for (i in 1:n0){
  for (j in 1:n1){
    zz[i,j] <- llik(c(b0.seq[i], b1.seq[j]), x, n, y )
  }
}
zz


contour(b0.seq, b1.seq, zz, levels = quantile(zz, seq(0.2, 0.99, by =0.02)))



###punto c

start <- c()

nllik <- function(beta, x, n, y){
  - llik(beta, x, n, y)
}

start <- c(-5, 0.13)
start
nlm(f = nllik, p = start, n = n, y = y, x = x)

SMV <- nlm(f = nllik, p = start, n = n, y = y, x = x)$estimate
SMV

glm(cbind(y, n - y) ~ x, family = binomial(logit))



###punto d

library(nlme)

H <- fdHess(SMV, function(beta) llik(beta, x, n, y))$Hessian
H

vcov <- solve(-H)

stde <- sqrt(diag(vcov))

cbind(SMV - qnorm(0.975) *stde, SMV + qnorm(0.975) *stde)


