#LABORATORIO 1



#-------------------------------------ES. 1-------------------------------------

.current.seed <- c(123,456,789) # Seme iniziale
runif.wh <- function(n){
  a <- c(171,172,170)
  b <- c(30269,30307,30323)
  s <- .current.seed
  u <- rep(0, n)
  for(i in 1:n){
    s <- (a*s) %% b # Aggiornamento sequenza (n volte)
    u[i] <- sum(s/b) %% 1 # Parte frazionaria della somma
    # della sequenza dei numeri in s
  }
  .current.seed <<- s # Seme finale
  u
}




#densità exp 
dExp <- function(y, lambda=1){
  lambda*exp(-lambda*y)
}

#plots di alcune densità
curve(dExp(x),0,10,lty="solid",ylab="f(x)")
curve(dExp(x,0.5),lty="dashed",add=TRUE)
curve(dExp(x,1.5),lty="dotted",add=TRUE)



#FdR
pExp <- function(y, lambda=1){
  1-exp(-lambda*y)
}

#plots di alcune FdR
curve(pExp(x),0,10,lty="solid",ylab="F(x)")
curve(pExp(x,0.5), lty="dashed",add=TRUE)
curve(pExp(x,1.5), lty="dotted",add=TRUE)



#funzione quantile exp
qExp <- function(u, lambda){
  (-log(1-u))/lambda
}

#plots di alcune f. quantile
curve(qExp(x,1),0,1,lty="solid",xlab="u",ylab=expression(paste(F^{-1},(u))))
curve(qExp(x,0.5), lty="dashed",add=TRUE)
curve(qExp(x,1.5), lty="dotted",add=TRUE)



#funzione generatrice exp
rExp <- function(n, lambda = 1) { 
  qExp(runif(n), lambda)
}

#genero 10000 valori
ex <- rexp(10000)

#li printo
ex



#verifica di casualità

#FdR 
plot(ecdf(ex), do.points = FALSE) #FdR empirica
curve(pexp(x), col = "red", add = TRUE) #FdR teorica

#densità
hist(ex, prob = TRUE, nclass =100) #istogramma dei valori generati, nclass = 100
#indica il numero di classi in cui viene suddiviso l'istogramma
curve(dexp(x), add = TRUE) #densità teorica

#diagramma di dispersione
plot(ex[-length(ex)], ex[-1], pch = '.')

#ACF (funzione di autocorrelazione) e Box_Pierce test
acf(ex, ylim = c(-1, 1))
Box.test(u,type="L")













#-------------------------------------ES. 2-------------------------------------



#definisco la funzione di densità della Beta(2,2)
f <- function(x){return(6*x-6*x^2)}
#oppure
curve(dbeta(x, 2, 2), add = TRUE, col = 'red')


#per massimizzare la densità per trovare il k che cerchiamo uso:
massimo <- optimize(f, interval = c(0, 1), maximum = TRUE)
k = optimize(f, interval = c(0, 1), maximum = TRUE)$objective

#$maximum
#[1] 0.5

#$objective
#[1] 1.5

r.acc.rif <- function(n, f ,g ,rg ,k ,report = TRUE){
  y <- double(n) 
  ntry <- 0
  for (i in 1:n){  
    done <- FALSE 
    while (!done) { 
      ntry <- ntry +1. 
      z <- rg(1)       
      u <- k*g(z) *runif(1) 
      if (u <= f(z)) 
        done <- TRUE 
      maxtry = 1e+9
      if (ntry >= maxtry)
        stop('troppi tentativi! controlla il codice')
    }
    y[i] <- z 
  }
  if (report)
    cat('tentativi', ntry) 
  return(y)
}




#con laplace
rbeta.ar.lap <- r.acc.rif(10000, dbeta(x,2,2), dlaplace, rlaplace, k)

hist(rbeta.ar.lap, prob = T, nclass = 100) 
curve(dbeta(x,2,2), add = T)
#l'istogramma NON si adatta bene alla densità della beta 


# prove con laplace
rbeta.ar.lap1 <- r.acc.rif(10000, dbeta(x,2,2), dlaplace, rlaplace, k)
#tentativi 35608
rbeta.ar.lap2 <- r.acc.rif(10000, dbeta(x,2,2), dlaplace, rlaplace, k)
#tentativi 35748
rbeta.ar.lap3 <- r.acc.rif(10000, dbeta(x,2,2), dlaplace, rlaplace, k)
#tentativi 35408
rbeta.ar.lap4 <- r.acc.rif(10000, dbeta(x,2,2), dlaplace, rlaplace, k)
#tentativi 35631
rbeta.ar.lap5 <- r.acc.rif(10000, dbeta(x,2,2), dlaplace, rlaplace, k)
#tentativi 36309



#con uniforme 
rbeta.ar.unif <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, k)

hist(rbeta.ar.unif, prob = T, nclass = 100) 
curve(dbeta(x,2,2), add = T)
#l'istogramma si adatta bene alla densità della beta


# prove con uniforme
rbeta.ar.unif1 <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, k)
#tentativi 35608
rbeta.ar.unif2 <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, k)
#tentativi 35748
rbeta.ar.unif3 <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, k)
#tentativi 35408
rbeta.ar.unif4 <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, k)
#tentativi 35631
rbeta.ar.unif5 <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, k)
#tentativi 36309




#efficienza e tentativi
system.time(r.acc.rif(1e+6, dbeta(x,2,2), dlaplace, rlaplace, k))
#user  system elapsed 
#17.317   0.046  17.353
system.time(r.acc.rif(1e+6, dbeta(x,2,2), dunif, runif, k))
#user  system elapsed 
#5.552   0.017   5.564


#prove con k = 1,2,3


#k=1
x1 <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, 1)
plot(ecdf(x1), do.points = FALSE, col = 'red') #FdR empirica
curve(pbeta(x, 2, 2), add = TRUE) #FdR teorica
#le FdR non si sovrappongono


#k=2
x2 <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, 2)
plot(ecdf(x2), do.points = FALSE, add = FALSE, col = 'green')
curve(pbeta(x, 2, 2), add = TRUE)
#le FdR si sovrappongono


#k=3
x3 <- r.acc.rif(10000, dbeta(x,2,2), dunif, runif, 3)
plot(ecdf(x3), do.points = FALSE, add = FALSE, col = 'violet')
curve(pbeta(x, 2, 2), add = TRUE)
#le FdR si sovrappongono





#-------------------------------------ES. 3-------------------------------------


#densità della normale troncata ai valori positivi
dNorm <- function(y){    #f
  (2/sqrt(2*pi)) * exp(-y^2/2)
}

#densità dell'esponenziale (lambda = 1)
dExp2 <- function(y) { #g
  exp(-y) 
}

#funzione h: trovo k(y) t.c. k>=h(y) dove h(y) = f(y)/g(y)
h <- function(y) { #h=f/g
  dNorm(y)/dExp2(y) 
}

curve(dNorm(x), 0, 10)
#trovo il k 
k <- optimize(h,c(0,100),maximum=T)$objective


curve(k*dExp2(x), lty = 2, add = T)

#genero 10000 valori
r.Norm.ar <- r.acc.rif(10000, dNorm, dexp, rExp, k)


#verifico che la densità empirica assomigli alla densità teorica 
hist(r.Norm.ar, prob = T, nclass = 50)
curve(dNorm, add = T)





#-------------------------------------ES. 4-------------------------------------



box.muller <- function(n){
  v <- matrix(NA, n, 2) #preallocazione
  for (i in 1:n){
    u1 <- runif(1) #genero u1: un valore generato casualmente da un'uniforme
    u2 <- runif(1) #genero u2: un valore generato casualmente da un'uniforme
    R <- sqrt(-2*log(u1)) #R
    theta <- 2*pi*u2      #theta
    X <- R*cos(theta)     #X
    Y <- R*sin(theta)     #Y
    v[i,] <- c(X,Y) #riempio il vettore vuoto con le coppie (X,Y) generate 
  }  
  return (v)
}

#verifcare se i valori generati assomigliano alla densità di una normale std
z <- box.muller(10000)
z
hist(z, prob = T, nclass = 1000)
curve(dnorm(x), add = T, col = 'red', lwd = 4)
#si nota che l'istogramma assomiglia alla densità teorica della normale std.
