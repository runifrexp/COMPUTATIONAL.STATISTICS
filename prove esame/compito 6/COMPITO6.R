#COMPITO 6


#-------------------------- ESERCIZIO 1

#funzione di probabilita' (densita' per le discrete)
fdistr <- function(x, p){
  p *(1-p)^(x-1) * (x>0)
}

#FdR
Fdistr <- function(x,p){
  v <- 0
  for (i in 1:x){
    v <- v + fdistr(i, p)
  }
  return(v)
}

#per generare da una discreta devo avere la sua funzione di ripartizione,
#la ottengo sommando 

###parte a
gen.distr <- function(p){
  u <- runif(1)
  val <- 0
  while (Fdistr(val, p) < u){
    val <- val + 1
  }
  val
}

rdistr <- function(n, p){
  sapply(rep(p,n), function(p) gen.distr(p))
}




###parte b

p <- 0.6
n <- 100
y <- rdistr(n, p)
y
hist(y, nclass = 10, prob = T)
curve(fdistr(x, p), col = 'magenta', lwd = 2, from = 0, to = 100)

#la funzione generatrice sembra essere corretta in quanto l'istogramma di 1000 
#valori generati assomiglia al grafico della funzione di probabilita' descritta
#nell'esercizio




###parte c
x <- rpois(n, )

curve(fdistr(x, p), col = 'magenta', lwd = 2, from = 20, to = 30)
curve(dpois(x, 1/20), add = T, col = 'lightskyblue', lwd = 2)

lambda <- 20
x <- rpois(n, lambda)
x

h <- (fdistr(x, p)/dpois(x, 30)) *(x>20)
h
mean(h)

alpha <- 0.01




###parte d

x <- rdistr(100, p)
x



B <- 100
y <- replicate(B, function(n, x) rpois(n, x))
y <- sapply(rep(n, B), function(n) rpois(n,x))
y

r <- function(){
  x <- rep(NA, B)
  for (i in 1:B){
    x[i] <- rdistr(1, p)
  }
  x
}
r()

n <- 100
y <- sapply(rep(n, B), function(n) rpois(n, r()))


y
mean(y)






#-------------------------- ESERCIZIO 2


rm(list =ls())
dati <- read.table(file.choose(), header = T) #dati del compito 4 esercizio 2
dati
x <- dati$x
x
y <- dati$y
y
n <- length(x)
n

###parte a

#H0 : med(x) = med(y)


#bootstrap np

B <- 10000
v.diff.med <- median(x) - median(y)
v.diff.med
  
x0 <- x - median(x)  
median(x0)

y0 <- y -median(y)
median(y0)

med <- function(x,y) median(x) - median(y)

med.bts <- rep(NA, B)
for (i in 1:B){
  idx <- floor(1+n*runif(n))
  xb <- x0[idx]
  yb <- y0[idx]
  med.bts[i] <- med(xb,yb)
}

hist(med.bts, prob = T, nclass = 40)
2*min(mean(med.bts <v.diff.med), mean(med.bts >v.diff.med))



###parte b

v.diff.var <- var(x) / var(y)
v.diff.var

x0v <- x / sd(x)
var(x0v)

y0v <- y/sd(y)
var(y0v)
#con x0v e y0v sono sotto H0 in quanto 1-1 = 0 (diff varianze)

f.var.b <- function(x,y) var(x) / var(y)

v.b <- rep(NA, B)
for (i in 1:B){
  idx <- floor(1+n*runif(n))
  xb <- x0v[idx]
  yb <- y0v[idx]
  v.b[i] <- f.var.b(xb, yb)
}

hist(v.b, nclass=100, freq = F)
curve(df(x, 30, 10), add = T, col = 'magenta')


2*min(mean(v.b>=v.diff.var), mean(v.b<=v.diff.var))



###punto c

alpha <- 0.05

#mediana 

med.bts <- rep(NA, B)
for (i in 1:B){
  idx <- floor(1+n*runif(n))
  xb <- x[idx]
  yb <- y[idx]
  med.bts[i] <- med(xb,yb)
}


ic <- 2* v.diff.med - quantile(med.bts, c(1-alpha/2, alpha/2))
ic
#rifiuto l'ipotesi nulla



#varianza

v.b <- rep(NA, B)
for (i in 1:B){
  idx <- floor(1+n*runif(n))
  xb <- x[idx]
  yb <- y[idx]
  v.b[i] <- f.var.b(xb, yb)
}

ic <- 2* v.diff.var - quantile(v.b, c(1-alpha/2, alpha/2))
ic
#accetto l'ipotesi nulla


###punto d
alpha <- 0.01

vera.var.x <- var(x)
vera.var.x

v.xb <- rep(0,B)
for (i in 1:B){
  idx <- floor(1+n*runif(n))
  v.xb[i] <- var(x[idx])
}

v.xbb <- replicate(B, var(x[floor(1+n*runif(n))]))

hist(v.xbb, nclass = 50)
a <- 2* vera.var.x - quantile(v.xbb,0.01)
c(-Inf, a)


#-------------------------- ESERCIZIO 3

rm(list = ls())
x <- c(1.1670504, 2.2437352, 4.3123664, 1.9894887, 2.7824318, 2.0157447, 
          1.6567842, 4.0169583, 4.4108051, 3.9200230,
                    1.8994196, 1.9990551, 4.9023921, 0.7810407, 1.9128419, 
          1.5853165, 0.6539528, 4.7494736, 2.2244362, 0.8042446,
                    3.6359439, 1.4133516, 3.9277407, 10.4341082, 17.7799186, 
          2.5969365, 3.1621217, 6.7595466, 4.5241052, 2.4502098,
                    2.9652919, 3.6072830, 13.2651019, 4.4439039, 3.8883708, 
          4.9557474, 2.1419572, 1.8494534, 0.7431344, 4.9798833,
                    4.5669072, 3.8204000, 2.7582161, 0.4754341, 7.3623801,
          0.5583623, 0.4732613, 2.3972103, 2.6671306, 0.5089649,
                    8.1234074, 6.4645910, 1.1770131, 1.1734521, 1.0066231, 
          4.0890214, 11.2765140, 10.1184573, 1.9521476, 3.1410062,
                    5.6810402, 0.6775010, 4.8434583, 1.3026003, 1.4806974, 
          1.2361622, 4.2062373, 6.9679000, 3.7586663, 3.1912299,
                    4.5404899, 3.0939591, 6.9711306, 4.1833636, 1.3181758, 
          1.6938775, 0.2154538, 5.9633446, 3.2677088, 1.0557894,
                    4.7932988, 6.3460110, 14.8074798, 6.9406986, 4.1221336,
          3.3714387, 1.9493682, 5.5225306, 5.0483030, 6.0667487,
                    5.4012607, 1.1069708, 2.1543786, 3.2537756, 0.3849553, 
          4.3934407, 6.7692799, 0.8035184, 4.2874418, 2.9766753)


###parte a

fdens <- function(par, x){
  mu <- par[1]
  sigma <- par[2]
  ex <- exp((x-mu)/sigma)
  (1/sigma) * ex/((1+ex)^2)
}

llik <- function(par, x){
  sum(log(fdens(par,x)))
}

#verifico se funziona con un esempio
llik(c(0-3,4), x)



###parte b

#metodo dei momenti: eguaglio i momenti teorici con i momenti empirici: 
#dai dati dell'esercizio sappiamo che E(X) = mu e Var(X) = ((pi^)2/3)*(sigma^2)

mu <- mean(x)
mu

#si trova analiticamente 

sigma <- sqrt((var(x)/pi^2)*3)
sigma




#parte c

#stima di massima verosimiglianza

#riparametrizzo sigma in modo tale da avere un nuovo parametro in R, non solo in 
#R+
start = c(mu, log(sigma))
start

nllik.rip <- function(par, x){
  - llik(c(par[1], exp(par[2])), x)
}

SMV.rip <- nlm(f = nllik.rip, p = start, x = x)

#devo ritornare alle stime iniziali 
SMV <- list(SMV.rip$estimate[1], exp(SMV.rip$estimate[2]))
SMV



###parte d

alpha <- 0.05

m.seq <- seq(0.1, 7, by = 0.5)
lm <- length(m.seq)
lm

s.seq <- seq(0.1, 7, by = 0.5)
ls <- length(s.seq)
ls


zz <- matrix(NA, lm, ls)

for (i in 1:lm){
  for (j in 1:ls){
    zz[i,j] <- llik(c(m.seq[i], s.seq[j]), x)
  }
}

contour(m.seq, s.seq, zz, levels = quantile(zz, seq(0.8, 0.99, length=lm)))
contour(m.seq, s.seq, zz, levels = (llik(SMV, x) - (qchisq(0.99, 2)/2)), col = 'magenta', add = T)
points()

contour(a.seq, t.seq, zz, levels = (llik(SMV, x) - 
                                      (qchisq(0.99, 2)/2)), col = 'magenta', add = T)










