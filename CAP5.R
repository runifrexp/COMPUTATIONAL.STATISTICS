#OTTIMIZZAZIONE 





#UN PRIMO ESEMPIO 



campione = c(0.25, 0.01, 0.45, 0.92, 0.55, 0.03, 1.51, 0.9, 0.2, 0.04)
n = length(campione)
campione

#funzione di verosimiglianza per la stima di lambda.cap
#il primo argomento deve sempre essere
lik <- function(lambda, y) { 
  lambda^length(y) * exp(-lambda * sum(y))
}
lik(1, campione)

#log-verosimiglainza
llik <- function(lambda, y) {
  length(y) * log(lambda) - lambda * sum(y)
}

par(mfrow = c(1, 2), mar = c(5, 4, 0.2, 0.2))
curve(lik(x, campione), from = 0, to = 5, lwd = 3)
#stampo la funzione di veros. per verificare che non ci siano dei massimi locali
#che possono disturbare le funzioni di ottimizzazione per la stima dello SMV
curve(llik(x, campione), from = 0, to = 5, lwd = 3)

#dopo aver trovato analiticamwente (vedi appunti) la SMV tramite lo score 
#eguagliato a zero abbiamo ottenuto 
smv.es = 1/mean(campione)
smv.es

#la stima della varianza invece vale 
var.es = smv.es^2/n
var.es



#ora usiamo la funzione optimize: devo specificare maximize = T perche' il 
#default e' maximize = F, quindi di default optimize trova un minimo 
optimize(llik,y=campione,interval=c(0,5),
         maximum=T,tol=1E-5)
#argomenti obbligatori: funzione da ottimizzare, intervallo
#La funzione obiettivo deve avere come primo argomento il vettore dei parametri 
#rispetto ai quali ottimizzare, eventuali altri argomenti (tipicamente il 
#campione) vengono passati con il nome con cui compaiono nella funzione obiettivo

#sintassi alternativa 
optimize(llik,y=campione,lower=0,upper=5,maximum=TRUE,tol=1E-5)




#INTERVALLI DI CONFIDENZA

#IC basato sulla distribuzione asintotica dello SMV: approccio WALD
a = 0.05
ic.norm = smv.es + c(-1, 1) * qnorm(1 - a/2) * sqrt(var.es)
ic.norm

#IC basato sulla distribuzione asintotica del rapporto di verosimiglianza: 
#approccio WILKS
curve(llik(x, campione) - llik(smv.es, campione), from = 0,
      to = 5)
abline(h = -qchisq(1 - a, 1)/2, lty = 2)

g = function(l, campione) {
  abs(llik(l, campione) - llik(smv.es, campione) + qchisq(1 - a, 1)/2)
}
ris.low = optimize(g, lower = 0, upper = smv.es, campione = campione)
ris.up = optimize(g, lower = smv.es, upper = 5, campione = campione)
c(ris.low$minimum, ris.up$minimum)


#confronto graifco tra WALD e WILKS
curve(llik(x, campione) - llik(smv.es, campione), from = 0,
      to = 5) #WILKS 
llikquad = function(lambda, campione){ 
  llik(smv.es, campione) - 0.5 * (n/smv.es^2) * (lambda - smv.es)^2
}
curve(llikquad(x, campione) - llik(smv.es, campione),
      add = TRUE, col = "red") #WALD (si nota la forma parabolica)
abline(h = -qchisq(1 - a, 1)/2, lty = 2)
abline(v = ris.low$minimum, lty = 2)
abline(v = ris.up$minimum, lty = 2)
abline(v = ic.norm[1], col = "red", lty = 2)
abline(v = ic.norm[2], col = "red", lty = 2)





#SMV PER UN MODELLO GAMMA CON DUE PARAMETRI



n= 100
alfa = 2
mu = 1.5
y = rgamma(n, alfa, alfa/mu)

#funzione di log-veros. per la gamma
ll=function(t,y){
  -t[1]*sum(y)/t[2]+(t[1]-1)*sum(log(y))+
    length(y)*t[1]*log(t[1]/t[2])-
    length(y)*log(gamma(t[1]))
}

#griglia di valori di (alfa, mu) in corrispondenza dei quali la 
#vogliamo rappresentare.
t1.seq = seq(0.01, 5, length = 50)
t2.seq = seq(0.01, 5, length = 50)


z = matrix(NA, 50, 50) 
for (i in 1:50) {
  for (j in 1:50) {
    z[i, j] = ll(c(t1.seq[i], t2.seq[j]), y)
  } 
}

range(z)

contour(t1.seq, t2.seq, z, x,
        levels = seq(-3000, -1000, length = 100))
#si nota che la rappresentazione non e' chiara

#per disegnare soltanto le curve vicino al massimo:
livelli = quantile(z, prob = c(0.5, 0.7, 0.8,
                               seq(0.9,1,by=0.01)))
contour(t1.seq, t2.seq, z, levels = livelli)
#In questo caso i livelli non sono equispaziati ed il grafico non puo' essere 
#usato per studiare la velocita' di crescita della funzione.



#OPTIM
#optim e' simile a optimize solo che vuole input multidimensionali e non un 
#intervallo, bensi' un punto di partenza per cercare l'ottimo
#optim puo' soltanto minimizzare ua funzione obiettivo, quindi bisogna definire
#la funzione di log verosimiglianza cambiata di segno, in quanto 
#il minimo di una funzione cambiata di segno corrisponde al suo massimo
#e' utile ricordare che le funzioni di ottimizzazione vogliono come primo 
#input il vettore di parametri da ottimizzare
nll = function(t, y) -ll(t, y)
ris = optim(fn = nll, par = c(1, 1), y = y)
ris

ris$convergence
#dato che e' zero vuol dire che l'ottimizzazione dovrebbe essere andata a 
#buon fine
#se fosse diverso da zero significherebbe che l'ottimizzatore non e' sicuro di 
#essere arrivato all'ottimo

optim(fn=nll, par=c(1,1), y=y, method="Nelder-Mead",
      control=list(maxit=1000, reltol=0.01))

#funzione alternativa a optim():
ris.nlm = nlm(f = nll, p = c(1, 1), y = y)
ris.nlm





#STIMA DELLA DISTRIBUZIONE ASINTOTICA 

#stima degli ic WALD
#mi serve l'informazione osservata per ottenere la matrice di 
#varianza e covarianza (sigma), che mi serve per poter determinare l'ic con 
#l'approccio alla WALD
#
ris = optim(fn = nll, par = c(1, 1), y = y, hessian = TRUE) 
round(ris$hessian, 4)

#ecco la matrice di varianze e covarianze ottenuta con l'inversione della
#matrice di informazione osservata
vcov = solve(ris$hessian) 
round(vcov, 4)

#comando analogo con nlm
ris.nlm = nlm(f = nll, p = c(1, 1), y = y, hessian = TRUE)
ris.nlm


#calcolo analitico della matrice di varianze e covarianze (R come calcolatrice)
alfa.hat = ris$par[1] #alfa
mu.hat = ris$par[2] #mu
n/alfa.hat - n * trigamma(alfa.hat) #sigma11
-2 * alfa.hat * sum(y)/mu.hat^3 + n * alfa.hat/mu.hat^2# sigma21 = sigma12
sum(y)/mu.hat^2 - n/mu.hat #sigma22



#rappresentazione grafica della regione di confidenza per theta

#WALD
liv.conf = 0.9
t1.seq = seq(1, 3.5, length = 50) 
t2.seq = seq(1, 1.8, length = 50) 
zn = matrix(NA, 50, 50)
theta.hat = c(alfa.hat,mu.hat) 
for (i in 1:50) {
  for (j in 1:50) {
    theta = c(t1.seq[i],t2.seq[j])
    zn[i, j] = t(theta.hat-theta) %*% solve(vcov) %*% (theta.hat-theta)
  } 
}
contour(t1.seq, t2.seq, zn, levels = qchisq(liv.conf,2), col = "red")
#si nota la simmetria della regione di confidenza


z = matrix(NA, 50, 50) 
for (i in 1:50) {
  for (j in 1:50) {
    z[i, j] = ll(c(t1.seq[i], t2.seq[j]), y)
  } 
}
contour(t1.seq, t2.seq, z, levels = ll(ris$par, y) - 
          qchisq(liv.conf, 2)/2, add = TRUE)
#preferisco WILKS a WALD in quanto non ha l'ulteriore vincolo 
#dell'approssimazione parabolica che porta con se' simmetria


#nota sui tempi di calcolo
f <- function(x) sum((x-1:length(x))^2)
a <- rep(1,5)
system.time(replicate(500, optim(fn=f, par=a)))
system.time(replicate(500, nlm(f=f, p=a)))
#guardare solo elapsed: si nota che nlm e' piu' veloce di optimize





#FREQUENZA DI FIGLI ALBINI


f = c(25, 23, 10, 1, 1) 
y = 1:5
m=5
n = sum(f)
y.bar = sum(y * f)/n

#definisco la log verosimiglianza
l <- function(p) {
  n * (y.bar * log(p) + (m-y.bar) * log(1-p) - log(1-(1-p)^m))
}

#grafico della log verosimiglianza
curve(l(x), from = 0, to = 1, lwd = 3)

#soluzione con optimize (dato che stimo un solo parametro)
ris = optimize(l, lower = 0, upper = 1, maximum = TRUE)
ris

#soluzione alternativa (vedi appunti pg 64)
g = function(p) (y.bar/m) * (1-(1 - p)^m) 
p0 = 0.5
p1 = g(p0)
while (abs(p1 - p0) > 1e-05) {
  p0 = p1
  p1 = g(p1)
  print(p1)
}

ris$maximum - p1 
l(ris$maximum) - l(p1)
#si evince che optimize e' piu' preciso dell'algoritmo definito a mano









#TEMPI DI GUASTO


y <- c(3, 5, 5, 13, 14, 15, 22, 22, 23, 30, 36, 39, 44,
       46, 50, 72, 79, 88, 97, 102, 139, 188, 197, 210)

#definizione della densita' Weibull, alcuni esempi 
dweib <- function(y, alpha, beta)
  beta/alpha * (y/alpha)^(beta - 1) * exp(-(y/alpha)^beta)
#@@@@@@@  e' fondamentale scrivere correttamente la densita', poi basta usare il 
#@@@@@@@  'trucchetto'

alpha <- 1
beta <- 0.5
curve(dweib(x,alpha,beta), 0, 10, ylab = "")

alpha <- 1
beta <- 1
curve(dweib(x,alpha,beta), 0, 10, add = TRUE, lty = "dotted")

alpha <- 1
beta <- 2
curve(dweib(x,alpha,beta), 0, 10, add = TRUE, lty = "dashed")

alpha <- 3
beta <- 2
curve(dweib(x,alpha,beta), 0, 10, add = TRUE, lty = "dotdash")

#si evince che alpha e' il parametro di scala, beta e' il parametro di forma.

legend(6, 1, lty = c("solid", "dotted", "dashed", "dotdash"),
       legend = c("a=1, b=0.5", "a=1, b=1", "a=1, b=2", "a=3, b=2"),
       lwd = 2)


#log-verosimiglianza
wei.loglik <- function(y, alpha, beta) { 
  sum(log(beta/alpha) + (beta-1) * log(y/alpha) - (y/alpha)^beta)
}

#usando il 'trucchetto:
llik <- function(y, par){
  alpha = par[1]
  beta = par[2]
  sum(log(dweib(y, alpha, beta)))
}

llik(y, c(4,0.3))
wei.loglik(y, 4, 0.3)
#si nota che le scritture sono equivalenti, ma e' piu' comodo procedere con il 
#il secondo metodo





#le procedure di ottimizzazione funzionano meglio se dispongono di buoni 
#punti di partenza per la ricerca dell’ottimo

#cerchiamo quindi delle stime preliminari per alpha e beta 
#al posto di p metto alcuni quantili
p <- c(0.25, 0.5, 0.75)
w <- log(quantile(y, p))
x <- log(-log(1 - p))
v <- coef(lm(w ~ x)) #coefficienti stimati
start <- c(exp(v[1]), 1/v[2]) #stime preliminari
start

#nel nostro caso, notiamo che alpha e beta sono numeri positivi 

#se non imponiamo dei vincoli, niente vieta all'algoritmo di ottimizzazione
#di estendere la sua ricerca al di fuori dello spazio parametrico
#e' bene fare in modo che all'algoritmo sia impossibile uscire dallo spazio 
#parametrico

#un modo per imporre i vincoli richiesti consiste nel “riparametrizzare” 
#opportunamente il problema

#dopo l'opportuna riparametrizzazione, va minimizzata la log-verosimiglianza
#cambiata di segno con argomenti l'inverso della mappa che abbiamo definito
#come riparametrizzazione
wei.nll.rip <- function(lp){
  -wei.loglik(y,exp(lp[1]),exp(lp[2]))
}

#la funzione nlminb (in box), non necessita di riparametrizzazioni (vedremo
#in seguito)
m <- nlminb(objective = wei.nll.rip, start = log(start))
m

par <- exp(m$par)
par

max.loglik <- -m$objective 
max.loglik

#un'alternativa consiste nel "dichiarare" a nlminb che il minimo deve 
#essere cercato in un sottoinsieme del piano cartesiano
wei.nll <- function(p) -wei.loglik(y, p[1], p[2])

m <- nlminb(objective = wei.nll, start = start, 
            lower = c(0, 0), upper = c(Inf, Inf))
m            #ecco il sottoinsieme del piano


#verifichiamo che la stima sia sensata confrontando la FdR teorica con quella 
#empirica
plot(ecdf(y)) #FdR teorica
alpha <- par[1]
beta <- par[2]
curve(1 - exp(-(x/alpha)^beta), 0, 220, add = TRUE) #FdR empirica





###INTERVALLI E REGIONI DI CONFIDENZA

##WALD

#calcolo della matrice Hessiana per poi trovare la matrice di varianza e 
#covarianza invertendo l'Hessiana per costruire gli ic
library(nlme)
H <- fdHess(par, function(par) wei.loglik(y, par[1], par[2]))
H

#inverto l'Hessiana cambiata di segno per trovare la matrice di varianza e 
#covarianza
vcov <- solve(-H$Hessian) 
vcov

#ecco lo standard error per costruire gli ic (radice delle varianze, quindi 
#degli elementi sulla diagonale)
se <- sqrt(diag(vcov))
se

#ic per alpha (intercept) e per beta (x)
cbind(par - qnorm(0.975) * se, par + qnorm(0.975) * se)


#PASSAGGIO DA WALD A WILKS

#una regione di confidenza non e', generalmente, il rettangolo (cubo, per tre 
#parametri, o ipercubo, per piu' di tre parametri) dato dal prodotto cartesiano 
#degli intervalli di confidenza unidimensionali

#costruiamo allora una regione di confidenza, ad usando la distribuzione 
#asintotica del rapporto di verosimiglianza

#e' importante notare pero' che l'implementazione della log verosimiglianza
#profilo e' molto piu' dispensiosa in quanto va fissato alpha e calcolato
#SMV per ogni alpha fissato, e lo stesso va fatto per beta

a <- seq(20, 150, length = 50)
b <- seq(0.1, 2, length = 50)
llik <- matrix(NA, length(a), length(b))
for (ia in 1:length(a))
  for (ib in 1:length(b)) {
    llik[ia, ib] <- wei.loglik(y, a[ia], b[ib])
  }
contour(a, b, llik, levels = max.loglik - qchisq(c(0.95, 0.99), 2)/2)
points(par[1], par[2], pch = "*", cex = 2) #SMV


#@@@ottimizzazione per alpha
wei.nll.rip.a <- function(lb, alpha, y) - wei.loglik(y, alpha, exp(lb))

sup.lik.a <- function(alpha, beta.start, y) {
  m <- nlminb(objective = wei.nll.rip.a, start = log(beta.start),
              alpha=alpha, y=y)
  -m$objective 
}

alpha.v <- seq(30, 140, length = 50)
lik.prof.a <- rep(NA, length(alpha.v))
for (i in 1:length(alpha.v))
  lik.prof.a[i] <- sup.lik.a(alpha.v[i], beta.start=1, y)

plot(alpha.v, lik.prof.a, type = "l")
target <- max.loglik - qchisq(0.95, 1)/2
abline(h = target, lty = 2)
grid()


g <- function(alpha){
  abs(sup.lik.a(alpha, beta.start=1, y) - target)
}
optimize(g, interval=c(40, par[1]), tol=1e-5)

optimize(g, interval=c(par[1],120), tol=1e-5)


#@@@Ripetiamo la stessa procedura per beta

wei.nll.rip.b <- function(la, beta, y)
  -wei.loglik(y, exp(la), beta)

sup.lik.b <- function(beta, alpha.start, y) {
  m <- nlminb(objective = wei.nll.rip.b, start = log(alpha.start),
              beta=beta, y=y)
  -m$objective 
}

beta.v <- seq(0.5, 1.5, length = 50)
lik.prof.b <- rep(NA, length(beta.v))
for (i in 1:length(beta.v))
  lik.prof.b[i] <- sup.lik.b(beta.v[i], alpha.start=mean(y), y)

plot(beta.v, lik.prof.b, type = "l")
target <- max.loglik - qchisq(0.95, 1)/2
abline(h = target, lty = 2)
grid()

g <- function(beta){
  abs(sup.lik.b(beta, alpha.start=mean(y), y) - target)
}
optimize(g, interval=c(0.5, par[2]), tol=1e-5)

optimize(g, interval=c(par[2],1.5), tol=1e-5)






