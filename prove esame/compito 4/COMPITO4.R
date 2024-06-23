#COMPITO 4



#-------------------------- ESERCIZIO 1

g <- function(y){
  exp(-y)
}

f <- function(y){
  sqrt(2/pi)*exp((-y^2)/2)
}

h <- function(y){
  f(y)/g(y)
}

k <- optimize(f = h, lower = -10, upper = 10, maximum = T)$objective
k


r.acc.rif <- function(n, f, g, rg, k){
  y <- rep(NA, n)
  ntry <- 0
  for (i in 1:n){
    done <- F
    while(!done){
      ntry <- ntry +1
      z <- rg(1)
      u <- k * g(z) * runif(1)
      if (u <= f(z)){
        done <- T
      }
    }
    y[i] <- z
  }
  y
}

rg <- function(n){
  rexp(n)
}

n <- 100
r.acc.rif(n , f, g, rg, k)




###punto b

n <- 1000
y <- r.acc.rif(n , f, g, rg, k)
y
hist(y, nclass = 100, prob = T)
curve(f, col = 'magenta', add = T)
mean(y)


###punto c

B <- 10000
n <- 1000
xbars <- function(n, y){
  t <- rexp(n, 1/y)
  mean(sqrt(t))
}

mc <- rep(NA, B)
for (i in 1:B){
  mc[i] <- xbars(n, y)
}
mc
mean(mc)

ic <- mean(mc) + c(-1, 1) *sd(mc) *qnorm(0.975)
ic



###punto d

raob <- function(n, y){
  t <- rexp(n, (1/2) * sqrt(pi*y))
  mean(t)
}

rb <- rep(NA, B)
for (i in 1:B){
  rb[i] <- raob(n, y)
}
rb
mean(rb)








#-------------------------- ESERCIZIO 2

###punto a 

dati <- read.table(file.choose(), header = T)
dati

x <- dati$x
x

y <- dati$y
y

t.test(x,y)

shapiro.test(x)
shapiro.test(y)
#notiamo che la normalita' non viene accettata ne' in x ne' in y, percio' non ci 
#si puo' fidare del risultato ottenuto dal test t (oltre al fatto che 
#alphaoss = 0.02 e' un risultato che non mi fa ne' accettare ne' rifiutare 
#l'ipotesi di uguaglianza delle medie)




###punto b
#dato che il test t non e' affidabile data la mancanza di assunzione di 
#normalita' nei dati, posso effettuare il test con approccio bootstrap non
#parametrico

#genero delle statistica test sotto H_0 per valutare l'alphaosservato
B <- 5000

tboot <- function(x, y){
  nx <- length(x)
  ny <- length(y)
  (mean(x)-mean(y))/sqrt((var(x)/nx)+(var(y)/ny))
}

bts <- replicate(B, tboot((x-mean(x))[floor(1 + nx *runif(nx))], (y - mean(y))
                          [floor(1 + ny *runif(ny))]))
bts

hist(bts, prob = T, nclass = 100)
#si nota che siamo sotto H_0 perche' la differnza delle medie e' centrata in 0

ic <- 2*min(mean(bts <= tboot(x, y)), mean(cmp.boot>= t(x,y)))
ic


###parte c

t <- function(x, y) var(x)/var(y) 

x0 <- x/sd(x)
x0

y0 <- y/sd(y)
y0

cmp.boot <- replicate(B, t(x0[floor(1+nx*runif(nx))], y0[floor(1+ny*runif(ny))]))
cmp.boot
hist(cmp.boot, prob = T,  nclass = 100)
curve(df(x, nx-1, ny-1), add = T, col = 'magenta')

ic <- mean(cmp.boot) + c(-1,1) *sd(cmp.boot) *qt(0.975, c( nx-1, ny-1))
ic
mean(cmp.boot)

alphaoss <- 2*min(mean(cmp.boot<= t(x,y)), mean(cmp.boot>= t(x,y)))
alphaoss



###parte d

alpha <- 0.05

theta <- mean(x) - mean(y)
theta

diff <- function(x, y){
  mean(x) - mean(y)
}


theta.bts <- replicate(B, diff(x[floor(1+nx*runif(nx))], y[floor(1+ny*runif(ny))]))
theta.bts
mean(theta.bts)

ic <- 2* mean(theta.bts) - quantile(theta.bts, c(0.975, 0.025))
ic
