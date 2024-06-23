#COMPITO 2


#-------------------- ESERCIZIO 1

### punto a
pEE <- function(x, lambda, alpha){
  (1-e^(-lambda *x))^alpha
}

qEE <- function(p, lambda, alpha) {
  (log(1-p^(1/alpha))/(-lambda))
}
#queste due funzioni sono state trovate analticamente, integrando la funzione di
#denstià da zero a x, e la funzione quantile e' stata trovata invertendo la 
#funzione di ripartizione eguagliandola a p ed esplicitandola in funzione di p


### punto b
rEE <- function(n, lambda, alpha){
  qEE(runif(n), lambda, alpha)
}


###punto c
alpha <- 2
lambda <- 3
x <- rEE(1000, alpha, lambda)
hist(x, prob = T, nclass = 100)

dEE <- function(x, alpha, lambda) {
  alpha * lambda *((1-exp(-lambda*x))^(alpha-1))*exp(-lambda*x)
}
curve(dEE(x, alpha, lambda), col = 'magenta', add = T)


###parte d

B <- 10000
z <- rEE(B, alpha, lambda)
h <- cos(z)/dEE(z, alpha, lambda)
mean(h)






#-------------------- ESERCIZIO 2

d <- read.table(file.choose(), header = T)
d
y <- d$y
y
x <- d$x
x
n <- length(y)
n

mod0 <- lm(y ~ x + I(x^2))
summary(mod)
res <- residuals(mod0)
res
b0 <- coefficients(mod0)[1]
b0

b1 <- coefficients(mod0)[2]
b1

b2 <- coefficients(mod0)[3]
b2

B <- 5000

m <- matrix(NA, B, 3)
mods <- rep(0, B)

for (i in 1:B){
    index <- floor(1 + n *runif(n))
    y.b <- b0 + b1 *x + b2 *x^2 + res[index]
    mod <- lm(y.b ~ x + I(x^2))
    m[i,] <- coefficients(mod)
}

m

b0.bts <- apply(m, 2, mean) #cosi' trovo i beta bootstrap
b0.bts


###parte b

#distorsione b0
distr.distorsione.b0 <- m[,1] - coefficients(mod0)[1]
hist(distr.distorsione.b0, nclass = 50)
distorsione.b0 <- mean(m[,1]) - coefficients(mod0)[1]
distorsione.b0

distr.distorsione.b1 <- m[,2] - coefficients(mod0)[2]
hist(distr.distorsione.b1, nclass = 50)
distorsione.b1 <- mean(m[,2]) - coefficients(mod0)[2]
distorsione.b1


distr.distorsione.b2 <- m[,3] - coefficients(mod0)[3]
hist(distr.distorsione.b2, nclass = 50)
distorsione.b2 <- mean(m[,3]) - coefficients(mod0)[3]
distorsione.b2

#deviazione standard 
#b0
sd(m[1,])

#b1 
sd(m[2,])

#b2
sd(m[3,])



###parte c

b0.bts.vector <- rep(NA, B)
for (i in 1:B){
  index <- floor(1 + n *runif(n))
  y.b <- b0 + b1 *x + b2 *x^2 + res[index]
  mod <- lm(y.b ~ x + I(x^2))
  b0.bts.vector[i] <- coefficients(mod)
}

b0.bts.vector
hist(b0.bts.vector, nclass = 50)
b0.bts <- mean(b0.bts.vector)
b0.bts

ic <- 2 * b0.bts - quantile(b0.bts.vector, c(0.975, 0.025))
ic





#-------------------- ESERCIZIO 3


dati <- read.table(file.choose(), header = T)
dati
x <- x$x
x

a <- -2
b <- 4

#punto a

dNorm <- function(par, x, a = -2, b = 4){
  mu <- par[1]
  sigma <- par[2]
  y1 <- dnorm((x-mu)/sigma)
  y2 <- pnorm((b-mu)/sigma)
  y3 <- pnorm((a-mu)/sigma)
  (y1)/(y2-y3)*(1/sigma) *(x>=a)*(x<=b)
}

llik <- function(par, x){
  sum(log(dNorm(par, x)))
}
llik(c(4, 2.4), x)


###parte b
seq.alpha <- seq(0.1, 5, by = 0.1)
seq.sigma <- seq(0.1, 5, by = 0.1)
la <- length(seq.alpha)
ls <- length(seq.sigma)

zz <- matrix(NA, la, ls)
dim(zz)
for (i in 1:la){
  for (j in 1:ls){
    zz[i,j] <- llik(c(seq.alpha[i], seq.sigma[j]), x)
  }
}
zz

contour(seq.alpha, seq.sigma, zz, levels = quantile(zz, seq(0.2, 0.999, length= ls)))



#punto c

negllik <- function(ripar, x){
  mu = ripar[1]
  sigma = exp(ripar[2])
  - llik (c(mu, sigma), x)
}
negllik(c(1,3), x)

start = c(4, log(2.4)) #visto dal contour

nlm(f = negllik, p = start, x = x)
ris.nlm = nlm(f = negllik, p = c(1, 1), x = x, hessian = TRUE)

mu.cap <- nlm(f = negllik, p = c(1, 1), x = x, hessian = TRUE)$estimate[1]
mu.cap
sigma.map <- nlm(f = negllik, p = c(1, 1), x = x, hessian = TRUE)$estimate[2]
sigma.map
sigma.cap <- exp(sigma.map)
sigma.cap
SMV <- c(mu.cap, sigma.cap)
SMV

#punto d
library(nlme)
hess <- fdHess(SMV, function(par) llik(par, x))$Hessian
hess
vcov = solve(-hess)
vcov

stde <- sqrt(diag(vcov))
stde

cbind(SMV - qnorm(0.975) *stde, SMV + qnorm(0.975) *stde)





