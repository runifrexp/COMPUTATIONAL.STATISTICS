#LAB 5 PROF


#### Ottimizzazione

###############################################################################
### Esercizio n. 1 ############################################################
###############################################################################

rm(list=ls())

# Dati.

n = 10 # Numero di prove.
w = 4 # Numero di successi.

#### a.

# La distribuzione del campione, sulla base dell'informazione
# disponibile, e' descritta da una distribuzione binomiale per W.

# Funzione di verosimiglianza.
L = function(p,w,n) p^w*(1-p)^(n-w)

# Funzione di log-verosimiglianza.
l = function(p,w,n) w*log(p) + (n-w)*log(1-p)

# Disegniamo la distribuzione di probabilita' dei dati come funzione di p.
curve(dbinom(w,n,x), 0, 1)

# Disegniamo la funzione di verosimiglianza. A meno di una costante di
# proporzionalita', coincide con la probabilita' binomiale che abbiamo
# appena tracciato.
curve(L(x,w,n), 0, 1)

# Disegniamo la funzione di log-verosimiglianza.
curve(l(x,w,n), 0, 1)

# Usiamo "optimize" (NB: possiamo farlo solo perche' il parametro e'
#                    unidimensionale)

optimize(f=l, interval = c(0,1), w=w, n=n, maximum=TRUE)
optimize(f=L, interval = c(0,1), w=w, n=n, maximum=TRUE) # Equivalente, o quasi.

# Ovviamente conosciamo la soluzione analitica.
p.hat = w/n # SMV
p.hat

#### b.

# Derivando la log-verosimiglianza, con qualche passaggio si ottiene
# che il reciproco dell'informazione osservata, calcolata in
# corrispondenza della SMV, e'

v = p.hat*(1-p.hat)/n # Stima della varianza asintotica.
s = sqrt(v) # Errore standard asintotico.

# Equivalentemente:

library(nlme)
v = -1/fdHess(p.hat, function(x) l(x,w,n))$Hessian
s = sqrt(v)
s
s = as.double(s) # Invece che come matrice (1 x 1), vediamo s come scalare

# Intervallo di confidenza basato sulla distribuzione asintotica dello SMV.

alfa = 0.05
p.hat + c(-1,1) * s * qnorm(1-alfa/2)

# Usiamo ora la distribuzione asintotica del rapporto di
# verosimiglianza.

curve(l(x,w,n) - l(p.hat,w,n), 0, 1)   # Grafico.
abline(h = -qchisq(1-alfa,1)/2, lty=2)

# Valori degli estremi, con ottimizzazione.

g = function(p,w,n) abs(l(p,w,n) - l(p.hat,w,n) + qchisq(1-alfa,1)/2)

optimize(g,interval=c(0,p.hat), w=w, n=n) # Estremo sinistro.
optimize(g,interval=c(p.hat,1), w=w, n=n) # Estremo destro.


###############################################################################
### Esercizio n. 2 ############################################################
###############################################################################

rm(list=ls())

# Dati.

f <- c(6,6,6,8,6,4,5,4,1,2,1,1,1,1) # Frequenze osservate.

# Invece che ragionare su modalita' e frequenze, nel seguito sarebbe
# anche possibile usare il vettore di osservazioni generato con
# y = rep(0:13,f)

#### [a.]

choose.cont <- function(r,k) gamma(r+k)/(gamma(k)*factorial(r))

l.bn <- function(par,f){ # Log-verosimiglianza.
  k <- par[1]
  p <- par[2]
  r <- 0:(length(f)-1) # Valori osservati 0:13.
  
  pr <- choose.cont(r,k) * p^k * (1-p)^r # "choose.cont" calcola il coefficiente binomiale,
  # ammettendo valori di k non interi.
  # Siccome r e' un vettore, "choose.cont" qui mi
  # restituisce un insieme di coefficienti
  # binomiali, calcolati per 0, 1, 2,..., 13
  # insuccessi prima del k-mo successo.
  sum(f*log(pr))
}

#### [b.]

# Definisco i vettori kvec e pvec, contenenti sequenze di valori dei
# parametri k e p, rispettivamente. Calcolo poi i valori assunti dalla
# log-verosimiglianza in corrispondenza di tutte le coppie di punti
# (kvec[i],pvec[j]), per i = 1,...,length(kvec) e j = 1,...,length(pvec).

kvec <- seq(1, 100, by=1)
pvec <- seq(0.01, 0.99, by=0.01)

mat.l.bn <- matrix(NA,length(kvec),length(pvec)) # Matrice di valori
# della log-ver.
for(i in 1:length(kvec)){
  for(j in 1:length(pvec)){
    mat.l.bn[i,j] <- l.bn(c(kvec[i], pvec[j]), f)
  }
}

range(mat.l.bn)

contour(kvec,pvec,mat.l.bn,
        levels=quantile(mat.l.bn,seq(0.8,1,by=0.01)),xlab="k", ylab="p")

# Non si vede molto bene, ma sembra che il massimo sia attorno a
# (k,p)=(5,0.35). Concentriamo il grafico su quella regione.

kvec <- seq(1, 15, by=1)
pvec <- seq(0.01, 0.5, by=0.001)

mat.l.bn <- matrix(NA,length(kvec),length(pvec)) # Matrice di valori
# della log-ver.
for(i in 1:length(kvec)){
  for(j in 1:length(pvec)){
    mat.l.bn[i,j] <- l.bn(c(kvec[i], pvec[j]), f)
  }
}

range(mat.l.bn)

contour(kvec,pvec,mat.l.bn,
        levels=quantile(mat.l.bn,seq(0.8,1,by=0.01)),xlab="k", ylab="p")

# La forma della funzione non sembra molto semplice. Capiamo che i
# punti di partenza dell'algoritmo di ottimizzazione sono importanti.

#### [c.]

# Dalle espressioni di media e varianza, vediamo che possiamo scrivere
# V(Y) = E(Y)/p. Ne segue che p = E(Y)/V(Y).
# Sostituendo questa espressione nella formula per E(Y), e risolvendo
# rispetto a k, troviamo k = [E(Y)]^2/[V(Y)-E(Y)].

# Troviamo la media e la varianza campionarie.

r <- seq(0,length(f)-1) # Modalita' osservate.
m <- sum(r*f)/sum(f) # Media campionaria osservata.
v <- sum((r-m)^2*f)/(sum(f)-1) # Varianza campionaria osservata.

# Equivalentemente

y <- rep(0:13,f)
m <- mean(y)
v <- var(y)

# Stime mediante il metodo dei momenti.

km <- m^2/(v-m)
km
pm <- m/v
pm

#### [d.]

# Definiamo la funzione di log-verosimiglianza negativa, adottando una
# riparametrizzazione che ci assicuri di rimanere dentro lo spazio
# parametrico durante la ricerca dell'ottimo.

nll.rip <- function(ripar, f) { # ripar[1] = psi; ripar[2] = phi.
  k = 1 + exp(ripar[1])
  p = exp(ripar[2])/(1+exp(ripar[2]))
  
  -l.bn(c(k,p),f)
}

init = c(log(km-1),log(pm/(1-pm))) # Valori di partenza per psi e phi.
opt <- optim(fn = nll.rip, par=init, f=f)

psi.hat <- opt$par[1] # SMV di psi.
phi.hat <- opt$par[2] # SMV di phi.
k.hat <- 1+exp(psi.hat) # SMV di k.
p.hat <- exp(phi.hat)/(1+exp(phi.hat)) # SMV di p. 
theta.hat <- c(k.hat,p.hat) # SMV bidimensionale. La usiamo nel seguito.
theta.hat

#### [e.]

# Qui non usiamo riparametrizzazioni, sfruttando la capacita' di
# nlminb di cercare il massimo dentro una "scatola".

nll <- function(par, f) { # Log-verosimiglianza negativa, senza
  # riparametrizzazioni.
  -l.bn(par,f)
}

nlminb(objective = nll, start = c(km,pm), lower=c(1,0), upper=c(Inf,1), f=f)

# Le SMV sono molto simili a quelle trovate con optim.

#### [f.]

# Regioni di confidenza.

# Usiamo la distribuzione asintotica dello SMV.

library(nlme)
H <- fdHess(par=theta.hat, function(par) l.bn(par,f))
H # Matrice hessiana.

vcov <- solve(-H$Hessian) # Matrice di varianze e covarianze asintotica.
vcov

alpha = 0.1 # Livello di significativita' prescelto.

k.seq = seq(0.5, 6, length = 50) # Sequenza di valori per k e p.
p.seq = seq(0.1, 0.6, length = 50) # Usiamo valori attorno alle SMV.
zn = matrix(NA, 50, 50)
for (i in 1:50) {
  for (j in 1:50) {
    theta = c(k.seq[i],p.seq[j])
    zn[i, j] = t(theta.hat-theta) %*% solve(vcov) %*% (theta.hat-theta)
  }
}
contour(k.seq, p.seq, zn, levels = qchisq(1-alpha,2))

# Usiamo ora la distribuzione asintotica del rapporto di verosimiglianza.

z = matrix(NA, 50, 50)
for (i in 1:50) {
  for (j in 1:50) {
    z[i, j] = l.bn(c(k.seq[i], p.seq[j]), f)
  }
}
contour(k.seq, p.seq, z, levels = l.bn(theta.hat, f) -
          qchisq(1-alpha, 2)/2, add = TRUE, col = "red")

points(k.hat,p.hat,pch="*",cex=2) # Stima di massima verosimiglianza.

#### [g.]

# Intervalli di confidenza.

# Distribuzione asintotica degli SMV.

se <- sqrt(diag(vcov)) # Standard error.
se

cbind(theta.hat - qnorm(1-alpha/2) * se, theta.hat + qnorm(1-alpha/2) * se)

# Verosimiglianze profilo.

# Intervallo di confidenza per k.

# Scriviamo la verosimiglianza profilo per k.
# Dobbiamo dunque, per ogni valore di k, massimizzare rispetto a p. Lo
# facciamo usando la riparametrizzazione che ci assicura che l'ottimo
# verra' cercato in (0,1). Potremmo anche usare la capacita' di nlminb
# di cercare l'ottimo in un intervallo, senza riparametrizzazioni.

nll.rip.k <- function(phi, k, f) { # Verosimiglianza negativa
  # riparametrizzata, vista come
  # funzione solo di phi (e quindi di p).
  p = exp(phi)/(1+exp(phi))
  
  -l.bn(c(k,p),f)
}

sup.lik.k <- function(k, p.start, f) { # Verosimiglianza profilo per k.
  m <- nlminb(objective = nll.rip.k, start = log(p.start/(1-p.start)),
              k=k, f=f)
  -m$objective
}

# Equivalentemente, senza riparametrizzazioni:

#nll.k <- function(p, k, f) { # Verosimiglianza negativa
#                             # vista come funzione solo di p).
#    -l.bn(c(k,p),f)
#}

#sup.lik.k <- function(k, p.start, f) { # Verosimiglianza profilo.
#  m <- nlminb(objective = nll.k, start = p.start, lower=0, upper=1, k=k, f=f)
#  -m$objective
#}

# Rappresentazione grafica dell'intervallo cercato.

k.v <- seq(0.5, 10, length = 50)
lik.prof.k <- rep(NA, length(k.v))
for (i in 1:length(k.v))
  lik.prof.k[i] <- sup.lik.k(k.v[i], p.start=pm, f)

plot(k.v, lik.prof.k, type = "l")
target <- l.bn(theta.hat, f) - qchisq(1-alpha, 1)/2
abline(h = target, lty = 2)
grid()

# Troviamo numericamente gli estremi dell'intervallo.

g <- function(k){
  abs(sup.lik.k(k, p.start=pm, f) - target)
}

optimize(g, interval=c(0.5, k.hat), tol=1e-5)
optimize(g, interval=c(k.hat, 6), tol=1e-5)

# Intervallo di confidenza per p.

nll.rip.p <- function(psi, p, f) {  # Verosimiglianza negativa
  # riparametrizzata, vista come
  # funzione solo di psi (e quindi di k).
  k = 1 + exp(psi)
  
  -l.bn(c(k,p),f)
}

sup.lik.p <- function(p, k.start, f) { # Verosimiglianza profilo per p.
  m <- nlminb(objective = nll.rip.p, start = log(k.start - 1),
              p=p, f=f)
  -m$objective
}

# Equivalentemente, senza riparametrizzazioni:

#nll.p <- function(k, p, f) {  # Verosimiglianza negativa vista come
#                                # funzione solo di k.
#    -l.bn(c(k,p),f)
#}

#sup.lik.p <- function(p, k.start, f) { # Verosimiglianza profilo per p.
#  m <- nlminb(objective = nll.p, start = k.start, lower = 1, upper = Inf,
#              p=p, f=f)
#  -m$objective
#}

# Rappresentazione grafica dell'intervallo cercato.

p.v <- seq(0.01, 0.9, length = 50)
lik.prof.p <- rep(NA, length(p.v))
for (i in 1:length(p.v))
  lik.prof.p[i] <- sup.lik.p(p.v[i], k.start=km, f)

plot(p.v, lik.prof.p, type = "l")
target <- l.bn(theta.hat, f) - qchisq(1-alpha, 1)/2
abline(h = target, lty = 2)
grid()

# Troviamo numericamente gli estremi dell'intervallo.

g <- function(p){
  abs(sup.lik.p(p, k.start=km, f) - target)
}

optimize(g, interval=c(0.01, p.hat), tol=1e-5)
optimize(g, interval=c(p.hat, 0.9), tol=1e-5)


###############################################################################
### Esercizio n. 3 ############################################################
###############################################################################

rm(list=ls())

# Dati.

y <- c(3, 5, 5, 13, 14, 15, 22, 22, 23, 30, 36, 39, 44, 46, 50, 72,
       79, 88, 97, 102, 139, 188, 197, 210)

#### a.

wei.loglik <- function(y, alpha, q) {
  
  beta = log(-log(0.9))/(log(q)-log(alpha)) # Dalla funzione quantile
  # discende che beta e' questa
  # funzione di alpha e del
  # quantile 10% q.
  
  sum(log(beta/alpha) + (beta-1) * log(y/alpha)
      - (y/alpha)^beta) # Log-verosimiglianza.
}

wei.nll <- function(par) -wei.loglik(y, par[1], par[2]) # Log-ver. negativa.

alpha.0 = 63 # Valore iniziale per alpha.
q.0 = 63 * ((-log(0.9))^(1/1.1)) # Dal valore iniziale per beta,
# ricaviamo il valore iniziale per q.

# Ottimizziamo, senza riparametrizzare, ma vincolando lo spazio
# parametrico con nlminb.

m <- nlminb(objective = wei.nll, start = c(alpha.0,q.0), lower = c(0, 0), upper = c(Inf, Inf))
m
theta.hat <- m$par
theta.hat

# Naturalmente, almeno fino a qui, sarebbe stato possibile trovare le
# stime di massima verosimiglianza di alfa e beta e poi applicare l'equivarianza
# per trovare la stima di q.

#### b.

# Intervallo di confidenza per alpha.

wei.nll.a <- function(q, alpha, y)
  -wei.loglik(y, alpha, q)

sup.lik.a <- function(alpha, q.start, y) { # Log-ver. profilo per alpha.
  m <- nlminb(objective = wei.nll.a, start = q.start, lower = 0, upper = Inf,
              alpha=alpha, y=y)
  -m$objective
}

# Soluzione numerica.

target <- -m$objective - qchisq(0.95, 1)/2
g <- function(alpha){
  abs(sup.lik.a(alpha, q.start=q.0, y) - target)
}

optimize(g, interval=c(40, theta.hat[1]), tol=1e-5) # Intervallo di
optimize(g, interval=c(theta.hat[1],120), tol=1e-5) # confidenza per alpha.

# Intervallo di confidenza per q.

wei.nll.q <- function(alpha, q, y)
  -wei.loglik(y, alpha, q)

sup.lik.q <- function(q, alpha.start, y) { # Log-ver. profilo per q.
  m <- nlminb(objective = wei.nll.q, start = alpha.start, lower = 0,
              upper = Inf, q=q, y=y)
  -m$objective
}

# Rappresentazione grafica.

q.v <- seq(1, 30, length = 50)
lik.prof.q <- rep(NA, length(q.v))
for (i in 1:length(q.v))
  lik.prof.q[i] <- sup.lik.q(q.v[i], alpha.start=alpha.0, y)

plot(q.v, lik.prof.q, type = "l")
abline(h = target, lty = 2)
grid()

# Soluzione numerica.

g <- function(q){
  abs(sup.lik.q(q, alpha.start=alpha.0, y) - target)
}

optimize(g, interval=c(1, theta.hat[2]), tol=1e-5) # Intervallo di confidenza
optimize(g, interval=c(theta.hat[2],20), tol=1e-5) # per q.


###############################################################################
### Esercizio n. 4 ############################################################
###############################################################################

rm(list=ls())

# Dati.

y = c(0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1)
x = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)

#### a.

nll <- function(par, y, x) { # Log-verosimiglianza negativa.
  alpha <- par[1]
  beta <- par[2]
  p <- exp(alpha+beta*x)/(1+exp(alpha+beta*x))
  -sum(y * log(p) + (1-y) * log(1-p))
}

m.o <- optim(fn = nll, par = c(0,0), y=y, x=x)
m.o
m.n <- nlminb(objective = nll, start = c(0,0), y=y, x=x) # Risultato analogo.
m.n

# Determiniamo le probabilita' di successo stimate per x=0 e x=1.
# Il predittore lineare alfa + beta * x, stimato in x=0 e x=1, e'

pl = m.n$par[1] + m.n$par[2] * c(0,1)

# Quindi, le probabilita' di successo per x=0 e x=1 sono

exp(pl)/(1+exp(pl))

# Le probabilita' di successo sono 4/6 e 3/6, ossia le frazioni
# di successi nei due gruppi.

#### b.

# Usiamo ora la distribuzione asintotica del rapporto di
# verosimiglianza per tracciare una regione di confidenza.

alpha.v <- seq(-4,4,length=50)
beta.v <- seq(-4,4,length=50)
z = matrix(NA, 50, 50)
for (i in 1:50) {
  for (j in 1:50) {
    z[i, j] = -nll(c(alpha.v[i], beta.v[j]), y, x)
  }
}
contour(alpha.v, beta.v, z, levels = -m.n$objective - qchisq(0.95, 2)/2)

points(m.n$par[1],m.n$par[2],pch="*",cex=2) # Stima di massima verosimiglianza.


###############################################################################
### Esercizio n. 5 ############################################################
###############################################################################

rm(list=ls())

y <- scan("cauchy.dat") # Usare "file.choose()".

l <- function(alpha,y) { # Log-verosimiglianza con beta=1.
  beta <- 1
  -length(y)*log(beta) - sum(log(1+((y-alpha)/beta)^2))
}

nll <- function(alpha,y) -l(alpha,y) # Log-verosimiglianza negativa,
# con beta=1.

init <- c(-11,-1,0,1.5,4,4.7,7,8,38) # Valori di partenza.
smv <- conv <- obj <- rep(NA,length(init))
for (i in 1:length(init)){
  m <- nlminb(objective = nll, start = init[i], y=y)
  smv[i] <- m$par
  conv[i] <- m$conv
  obj[i] <- m$obj
}

smv # Ottimi individuati.
obj # Minimi raggiunti.
conv # E' sempre 0: la convergenza viene sempre raggiunta, qualsiasi
# sia il valore iniziale.

# In tutti i casi viene raggiunta la convergenza. Tuttavia, vediamo
# che per valori iniziali di grandezza intermedia, l'algoritmo
# converge ad un minimo locale. Valori iniziali sufficientemente
# grandi o piccoli portano ad individuare il massimo corretto.

init[abs(obj-min(obj))<1e-8] # Valori iniziali che portano al minimo corretto.
smv[abs(obj-min(obj))<1e-8] # Ottimo globale.
smv[abs(obj-min(obj))>1e-8] # Ottimo locale.

# Si noti che non sono necessariamente i valori iniziali piu'
# vicini all'ottimo globale che portano ad individuarlo correttamente.

# Per completezza, tracciamo l'andamento della log-verosimiglianza per
# alpha (con beta fissato pari a 1).

# Si noti che "curve" passa un vettore alla funzione che si desidera
# disegnare, aspettandosi come risposta un vettore di valori della
# funzione. La funzione l(alpha,y) che abbiamo scritto non funziona
# con alpha vettoriale. Riscriviamo questa funzione in modo che
# accetti alpha vettoriali.

l.1 <- function(alpha,y) { # Log-verosimiglianza con beta=1.
  # Accetta anche alpha vettoriali.
  beta <- 1
  -length(y)*log(beta) - sapply(alpha,FUN = function(a) sum(log(1+((y-a)/beta)^2)))
}

l(0:2,y) # Errore!
l.1(0:2,y) # OK.

curve(l.1(x,y),-4,8) # Vediamo che la funzione di log-verosimiglianza
# e' multimodale, con i massimi locale e globale
# che avevamo individuato.