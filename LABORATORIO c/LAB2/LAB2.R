#LABORATORIO 2



#-------------------------------------ES. 1-------------------------------------

#per ic uso l'ic definito dalla distribuzione asintotica (se TLC valido) della 
#stima dell'integrale
monte.carlo1 <- function(n, f, B = 10000, alpha = 0.01){
  y <- sapply(rep(n,B), function(n) f(runif(n)))
  ybar <- mean(y)
  ic <- qnorm(1-alpha/2)*sd(y)/sqrt(B)
  cat ("La stima di I(f) e'", ybar,
       "\nIl vero I(f) e' compreso tra ", ybar - ic, "e", ybar + ic,
       "\n con probabilita'", 1 - alpha, "\n")
  plot(y, pch = '.')
}

h <- function(x) (cos(50*x)+sin(20*x))^2


monte.carlo1(1000,h)


#metodo più veloce ed efficiente 

h <- function(x) (cos(50*x)+sin(20*x))^2

B <- 10000

I <- h(runif(B))

mean(I)
s2 = var(I)
se = sqrt(s2/B)
mean(I) + c(-1,1) * qnorm(0.975) * se #ic 

#oppure direttamente
mean(I) + c(-1,1) * qnorm(0.975) *sd(I) /sqrt(B)

#verifico con la funzione nativa di R dato che l'integrale è univariato
integrate(h, 0, 1)
#la stima fatta via MC è molto vicina al valore restituito dalla funzioen nativa

#parte b

ris <- matrix(NA, 10000, 4)
colnames(ris) <- c('B', 'mean', 'inf', 'sup')


for (B in 2:10000){
  x <- runif(B)
  ris[B-1, 1] <- B
  ris[B-1, 2] <- mean(h(x))
  ris[B-1, 3] <- ris[B-1, 2]-qnorm(0.975)*sd(h(x))/sqrt(B)
  ris[B-1, 4] <- ris[B-1, 2]+qnorm(0.975)*sd(h(x))/sqrt(B)
}
plot(ris[,2], pch = '.')
points(ris[,3], pch = '.', col = 'red')
points(ris[,4], pch = '.', col = 'green')
points(ris[,2], pch = '.') #per mettere in evidenza i valori generati 

#-------------------------------------ES. 2-------------------------------------

monte.carlo2 <- function(n, f, B = 10000, alpha = 0.01) {
  y <- sapply(rep(n,B), function(n) f(runif(n), runif(n), runif(n)))
  ybar <- mean(y)
  ic <- qnorm(1-alpha/2)*sd(y)/sqrt(B)
  cat ("La stima di I(f) e'", ybar,
       "\nil vero valore di I(f) e' compreso tra ", ybar - ic ,
       "e", ybar + ic,
       "\n con probabilita'", 1 - alpha)
}

#per ovviare agli estemi di integrazione 'scomodi' e quindi dover generare da 
#un'unif(0,3), effettuo un cambio di variabili: yi = xi/3, passando 
#dalla funzione g1
g1 <- function(x1,x2,x3) exp(-abs(x1-x2)^(2+sqrt(x3)))
#alla funzione g2
g2 <- function(y1,y2,y3) 27 * exp(-abs(3*y1-3*y2)^(2+sqrt(3*y3)))

monte.carlo2(1000, g2)


#-------------------------------------ES. 3-------------------------------------

monte.carlo3 <- function(n, f, B = 10000, alpha = 0.01){
  y <- sapply(rep(n,B), function(n) f(runif(n,0,3),runif(n,0,3),runif(n,0,3)))
  #genero da uniformi (0,3)
  ybar <- mean(y)
  ic <- qnorm(1-alpha/2)*sd(y)/sqrt(B)
  cat ("Il valore previsto per I(f) e'", ybar,
       "\n il vero valore di I(f) e' tra", ybar-ic, "e", ybar+ic,
       "\n con probabilità", 1-alpha)
}

g1 <- function(x1,x2,x3) exp(-abs(x1-x2)^(2+sqrt(x3))) * 27 #*27 perchè la 
#densità dell'uniforme (0,3) è 1/3

monte.carlo3(1000, g1)



#-------------------------------------ES. 4-------------------------------------

t <- function (x) (exp(-sqrt(abs(x))) * sin(x^2/50))

B <- 10000

x <- rnorm(B, 0, 20)

z <- t(x)/dnorm(x,0,20)

g <- function(x) dnorm(x, 0,20)
integrate(g, -Inf, +Inf)
#qui si nota che (giustamente) la densità della normale integra ad 1, quindi 
#nella stima dell'integrale non è influente.

#intervallo di confidenza
cat(mean(z), mean(z)+c(-1,1)*qnorm(0.975)*sd(z)/sqrt(B), "\n")

#dato che l'integrale è univariato lo posso stimare con la funzione integrate()
integrate(t, -Inf, +Inf)


#senza campionamento per importanza: faccio un cambio di variabile per poter
#generare da un'unif(0,1):
# y = 1/(1 + e −x )
fy <- function(y){
  x = -log((1-y)/y)
  1/(y*(1-y)) * exp(-sqrt(abs(x)))*sin(x^2/50)
}

w <- fy(runif(B))
cat("Stima:", mean(w), "\nIntervallo:", 
    mean(w)+c(-1,1)*qnorm(0.975)*sd(w)/sqrt(B), "\n")
#la stima ha una bassissima precisione, questo sinpuo' capire dall'ampiezza 
#dell' ic

integrate(fy,1E-16,1-1E-16) # Dobbiamo escludere 0 e 1 dall'intervallo di
# integrazione, altrimenti si possono toccare valori 'scomodi' per la funzione,
#ossia 0 e 1



#-------------------------------------ES. 5-------------------------------------

#media di Y

# - procdeura diretta
B <- 10000

x <- rexp(B)

diretta <- sapply(x, function(x) runif(1, 0, x))

cat("Stima diretta di E(Y):", mean(diretta), 
    "\nIntervallo:", mean(diretta)+c(-1,1)*qnorm(0.975)*sd(diretta)/sqrt(B), "\n")

sd(diretta) #deviazione standard molto alta



# - Rao-Blackwell
#uso la media condizionata di Y|X, essendo X un'uniforme ha media X/2, quindi 
#calcolo la media di Y tramite la media condizionata di Y|X, ossia X/2
#la varianza usando rb è piu' bassa (vedi pag. 176 2 PARTE)
xrb <- rexp(B)
rb <- mean(xrb/2)

sd(xrb/2) #deviazione standard molto alta
cat("Stima di E(Y) con RB:", mean(xrb/2), "\nIntervallo:", 
    mean(xrb/2)+c(-1,1)*qnorm(0.975)*sd(xrb/2)/sqrt(B), "\n")


# - RB+anti
#genero solo al 1° livello (x~exp(n)) (RB)
#creo correlazione negativa alla generazione fatta al 1° livello (anti)
M <- 500
rbant <- rep(NA,M)
for (i in 1:M){
  u <- runif(B/2)
  x1 <- qexp(u)
  x2 <- qexp(1-u) #variabile antitetica 1-u per ottenere correlazione 
  x <- c(x1,x2)                                                 #negativa
  rbant[i] <- mean(x/2)
}

mean(rbant) 
sd(rbant) #ottima precisione

boxplot(cbind(diretta, rb, rbant))

cat('Deviazioni standard: \n',
    '\ndiretta', sd(diretta),
    '\nRB', sd(xrb/2),
    '\nrbant', sd(rbant))

mean(y>1)



#-------------------------------------ES. 5-------------------------------------

# - diretta 
B <- 1000
x <- rexp(B)
y <- sapply(x, function(x) runif(1, 0, x))
I1 <- mean(y)
P1 <- mean(y>1)
cat("La stima diretta di E(Y) e'", I1,
    "\nE' contenuta tra", I1 - qnorm(0.975)*sd(y)/sqrt(B),
    "e", I1 + qnorm(0.975)*sd(y)/sqrt(B),
    "con una probabilita' del 95%",
    "\nLa stima di P(Y>1) è", P1)

# - RB
#devo fare le medie delle varie exp, poi faccio la media delle medie
#condizionate
I2 <- mean(x/2)
P2 <- mean((x/2)>1)
cat("La stima RB di E(Y) e'", I2,
    "\nE' contenuta tra", I2 - qnorm(0.975)*sd(x/2)/sqrt(B),
    "e", I2 + qnorm(0.975)*sd(x/2)/sqrt(B),
    "con una probabilita' del 95%",
    "\nLa stima di P(Y>1) è", P2)


# - RB+anti
#RB perchè genero solo dalla variabile di primo livello (X)
#anti perche' introduco correlazione negativa tra la vechia variabile e una 
#nuova variabile. 
u <- runif(B/2) #genero da uniforme per poi applicare la funzione quantile di X
x1 <- qexp(u) #u e 1-u per apportare correlazione negativa che fa abbassare
x2 <- qexp(1-u) #la varianza della stima
x <- c(x1,x2) #media marginale delle variabili antitetiche 
I3 <- mean(x/2) #media delle medie marginali delle v. anti.
P3 <- mean((x/2)>1)
cv <- cov(x1, x2)
cat("La stima RB+anti di E(Y) e'", I3,
    "\nE' contenuta tra", I3 - qnorm(0.975)*sd(x/2)/sqrt(B),
    "e", I3 + qnorm(0.975)*sd(x/2)/sqrt(B),
    "con una probabilita' del 99%",
    "\nHa senso usare le variabili antitetiche in quanto cov(x1,x2)=", cv, 
    "\nLa stima di P(Y>1) è", P3)


#confronto dell'efficienza delle stime: faccio MC su MC
stima.diretta <- function(){
  B <- 10000
  x <- rexp(B)
  y <- sapply(x, function(x) runif(1, 0, x))
  mean(y)
}
stima.RB <- function(){
  B <- 10000
  x <- rexp(B)
  mean(x/2)
}
stima.RBanti <- function(){
  B <- 10000
  u <- runif(B/2) 
  x1 <- qexp(u)  
  x2 <- qexp(1-u)  
  x <- c(x1,x2) 
  mean(x/2)
}

M <- 500
ris.diretta <- ris.RB <- ris.RBanti <- rep(NA, M)
for (i in 1:M){
  cat('\n',i)
  ris.diretta[i] <- stima.diretta()
  ris.RB[i] <- stima.RB()
  ris.RBanti[i] <- stima.RBanti()
}
mean(ris.diretta) #stima MC su MC ottenuta con metodo diretto
mean(ris.RB) #stima MC su MC ottenuta con metodo Rao-Blackwell
mean(ris.RBanti) #stima MC su MC ottenuta con v. antitetiche e RB
sd(ris.diretta) #stima MC su MC ottenuta con metodo diretto
sd(ris.RB) #stima MC su MC ottenuta con metodo Rao-Blackwell
sd(ris.RBanti) #stima MC su MC ottenuta con v. antitetiche e RB
boxplot(cbind(ris.diretta, ris.RB, ris.RBanti))


#stima di P(Y>1)

B <- 10000
x <- rexp(B)
y <- sapply(x, function(x) runif(1, 0, x))

# - diretto
p1 <- mean(y>1)
cat("La stima diretta di P(Y>1) e'", p1,
    "\nSta nell'intervallo [", p1+c(-1,1)*qnorm(0.975)*sqrt(p1*(1-p1)/B), "]",
    "con probabilita' 0.95")

# - RB
f <- function(x) (x-1)/x * (x>1)
p2 <- mean(f(x))
cat("La stima di P(Y>1) con RB e'", p2,
    "\nSta nell'intervallo [", p2+c(-1,1)*qnorm(0.975)*sqrt(p2*(1-p2)/B), "]",
    "con probabilita' 0.95")

# - RB+anti
u <- runif(B/2)
x1 <- qexp(u)
x2 <- qexp(1-u)
cov(x1,x2)
x <- c(x1,x2)
p3 <- mean(f(x))
cat("La stima di P(Y>1) con RB e'", p3,
    "\nSta nell'intervallo [", p3+c(-1,1)*qnorm(0.975)*sqrt(p3*(1-p3)/B), "]",
    "con probabilita' 0.95")

#MC su MC 
M <- 5000
p.diretta <- p.RB <- p.RBanti <- rep(NA, M)

for (i in 1:M){
  cat('\n', i)
  x <- rexp(B)
  y <- sapply(x, function(x) runif(1, 0, x))
  p.diretta[i] <- mean(y>1)
  p.RB[i] <- mean(f(x))
}
for (i in 1:M){
  cat('\n', i)
  u <- runif(B)
  x1 <- qexp(u)
  x2 <- qexp(1-u)
  x <- c(x1,x2)
  p.RBanti[i] <- mean(f(x))
}

mean(p.diretta)
mean(p.RB)
mean(p.RBanti)

sd(p.diretta)
sd(p.RB)
sd(p.RBanti)

boxplot(cbind(p.diretta, p.RB, p.RBanti))



#-------------------------------------ES. 6-------------------------------------

#funzione che genera da una uniforme discretizzata
f <- function(){
  x <- runif(365)
  y <- sapply(x, function(x) rbinom(1, 364, x))
  1+y
}

#verifica di casualità
b <- seq(1,365)
n <- length(b)

#densità
hist(y/n, nclass = 100, prob = T)
curve(dunif(x), add = T, col = 'blue')

#grafico di dispersione
plot(b,y)

#confronto FdR empirica e teorica
plot(ecdf(y/n))
curve(punif(x), add = T, col = 'green')


#oppure piu facilmente
x <- sample(1:365, 1)



sameday <- function(n){
  x <- sample(1:365, 1)
  y <- rep(NA,n)
  for (i in 1:n){  
    y[i] <- sample(1:365, 1)
  }
  v <- x == y
  sum(v)
}

#MC su MC per stimare la frazione di persone che compie gli anni lo stesso 
#giorno tramite simulazione
mc.sd <- function(n){
  B <- 10000
  sm <- 0
  for (i in 1:B){
    sm <- sm + sameday(n)
  }
  sm/B
}





