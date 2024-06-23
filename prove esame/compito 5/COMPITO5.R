#COMPITO 5 (20-07-2023)


#-------------------------- ESERCIZIO 1


###punto a

FdR <- function(x, k) {
  1 - exp(-(x/3)^k)
}

f.quant <- function(u, k) {
  3 * (-log(1-u))^(1/k)
}

r <- function(n, k = 5) {
  u <- runif(n)
  f.quant(u, k)
}

###punto b

n <- 1000
k <- 5
y <- r(n)
y

plot(ecdf(y))
curve(FdR(x, k), add = T, col = 'magenta', from = -10, to = 10)



###punto c









###punto d

#genero monte carlo E(Y*Z)

B <- 5000

x <- r(n)
x
#E(Y|X=x)
rpois(n, x)

#E(Z|X=x)
rexp(n, x)

b <- replicate(B, mean(rpois(n, x) * rexp(n, x)))
b

I <- mean(b)
I

se <- sqrt(var(b)/n)

ic <- (I +c(-1,1) * qnorm(0.975) * se)
ic[2]-ic[1] #ampiezza intervallo di confidenza


#-------------------------- ESERCIZIO 2

dati <- read.table(file.choose(), header = T) #esercizio 2.dat compito 1
dati
y <- dati$y
y
x <- dati$x
x
###parte a

#bootstrap sp per b0, b1, b2

#modello dai dati 'reali'
mod <- lm(y ~ x + I(log(x)))
summary(mod)
plot(x,y)
abline(mod)

b1 <- lm(y ~ x + I(log(x)))$coefficients[1]
b1

b2 <- lm(y ~ x + I(log(x)))$coefficients[2]
b2

b3 <- lm(y ~ x + I(log(x)))$coefficients[3]
b3

n <- length(x)
n



#calcolo i residui 
res <- residuals(mod)
res

B <- 5000
m <- matrix(NA, B, 3)
yb <- matrix(NA, B, n)
dim(yb)
for (i in 1:B){
  index <- floor(1+n*runif(n))
  yb[i,] <- b1 + b2 * x + b3 *log(x) + res[index]
}

beta <- matrix(NA, B, 3)
for (i in 1:B){
  beta[i,] <- coef(lm(yb[i,]~x + I(log(x))))
}
coeff.bts <- apply(beta, 2, mean)
coeff.bts

#distorsione 
coeff.mondoreale <- coef(lm(y ~ x + I(log(x))))
coeff.bts - coeff.mondoreale

#deviazione standard
sd.coeff.bts <- apply(beta, 2, sd)
sd.coeff.bts

#valori delle stime corrette
#tolgo il bias ai valori bootstrap
bias <- coeff.bts - coeff.mondoreale
bias
coeff.mondoreale - bias 



###punto c

beta2.cap <- coeff.mondoreale[2]
beta2.cap
#distr. b2 
db2 <- beta[,3]

#ic e' del tipo (-inf, a] (unilaterale sinistro)
a <- 2 * beta2.cap - quantile(db2 +1, 0.01)

ic <- c(-Inf, a)
ic








#-------------------------- ESERCIZIO 3
rm(list = ls())
x <- c(1.1670504, 2.2437352, 4.3123664, 1.9894887, 2.7824318, 2.0157447, 
            1.6567842, 4.0169583, 4.4108051, 3.9200230,
            1.8994196, 1.9990551, 4.9023921, 0.7810407, 1.9128419, 1.5853165, 
            0.6539528, 4.7494736, 2.2244362, 0.8042446,
            3.6359439, 1.4133516, 3.9277407, 10.4341082, 17.7799186, 2.5969365,
            3.1621217, 6.7595466, 4.5241052, 2.4502098,
            2.9652919, 3.6072830, 13.2651019, 4.4439039, 3.8883708, 4.9557474, 
            2.1419572, 1.8494534, 0.7431344, 4.9798833,
            4.5669072, 3.8204000, 2.7582161, 0.4754341, 7.3623801, 0.5583623, 
            0.4732613, 2.3972103, 2.6671306, 0.5089649,
            8.1234074, 6.4645910, 1.1770131, 1.1734521, 1.0066231, 4.0890214, 
            11.2765140, 10.1184573, 1.9521476, 3.1410062,
            5.6810402, 0.6775010, 4.8434583, 1.3026003, 1.4806974, 1.2361622, 
            4.2062373, 6.9679000, 3.7586663, 3.1912299,
            4.5404899, 3.0939591, 6.9711306, 4.1833636, 1.3181758, 1.6938775, 
            0.2154538, 5.9633446, 3.2677088, 1.0557894,
            4.7932988, 6.3460110, 14.8074798, 6.9406986, 4.1221336, 3.3714387, 
            1.9493682, 5.5225306, 5.0483030, 6.0667487,
            5.4012607, 1.1069708, 2.1543786, 3.2537756, 0.3849553, 4.3934407, 
            6.7692799, 0.8035184, 4.2874418, 2.9766753)
x


###punto a

dens <- function(par, x){
  alpha <- par[1]
  theta <- par[2]
  (1/(gamma(alpha)*theta^alpha))*x^(alpha-1)*exp(-x/theta)
}

llik <- function(par, x){
  alpha <- par[1]
  theta <- par[2]
  sum(log(dens(c(alpha, theta), x)))
}



###punto b
#metodo dei momenti: eguaglio momenti teorici ed empirici

sigma2.emp <- var(x)
sigma2.emp
mu.emp <- mean(x)
mu.emp

#analiticamente si e' trovato che alpha e theta valgono:
alpha <- mu.emp^2/sigma2.emp
alpha

theta <- mu.emp /alpha
theta


###punto c

#riparametrizzazione: dato che alpha e theta sono in [0, +inf) adotto la 
#riparametrizzazione logaritmica 
phi <- log(alpha)
phi

psi <- log(theta)
psi

c(phi, psi)
#phi e psi saranno i punti di partenza dell'ottimizzatore

nllik.rip <-  function(par, x){
  alpha <- par[1]
  theta <- par[2]
  - llik(c(exp(alpha), exp(theta)), x)
}

nlm(f = nllik.rip, p = c(phi, psi), x = x)

#SMV 
SMV.ripar <- nlm(f = nllik.rip, p = c(phi, psi), x = x)$estimate
SMV <- exp(SMV.ripar) 


###punto d

a.seq <- seq(0.01, 4, by = 0.01)
la <- length(a.seq)
t.seq <- seq(0.01, 4, by = 0.01)
lt <- length(t.seq)

zz <- matrix(NA, la, lt)

for(i in 1:la){
  for(j in 1:lt){
    zz[i,j] <- llik(c(a.seq[i], t.seq[j]), x)
  }
}
zz


contour(a.seq, t.seq, zz, levels = quantile(zz, seq(0.8, 1, by =0.01)))

contour(a.seq, t.seq, zz, levels = (llik(SMV, x) - 
                              (qchisq(0.99, 2)/2)), col = 'magenta', add = T)


