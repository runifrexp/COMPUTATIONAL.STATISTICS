#GENERAZIONE NUMERI CASUALI



# ALGORITMO DI GENERAZIONE CONGRUENZIALE IN R
# P. 7-8 'Generazione di numeri (pseudo-)casuali

# La seguente funzione genera n interi da una v.a. uniforme.
# La funzione prende in input :
# - Il valore iniziale x0
# - I parametri a, b, m;

generan = function(n,x0,a,b,m){
  ris=rep(NA,n)
  ris[1]=x0
  for (i in 2:n){
    ris[i]=(a*ris[i-1]+b)%%m
    }
  return(list(int=ris, unif=ris/m))
  }

#esempi di plots e numeri casuali generati 
p = generan(200,40,1,8,81)
plot(p$int, type ="l")

p = generan(200,40,1,8,2^32)
plot(p$int, type ="l")

p = generan(200,2,24325343,25453523,2^32)
plot(p$int, type ="l")



#VERIFICA DI CASUALITÀ 

#funzione di ripartizione
u = generan(200, 2, 1664525, 1013904223, 2^32)$unif
plot(ecdf(p), do.points = FALSE) #funzione di ripartizione empirica uniforme
curve(punif(x), col = "red", add = TRUE) #funzione di rip. teorica uniforme

#quantili
plot(sort(u), ((1:length(u)) - 0.5)/length(u)) 
abline(0,1)

#test indipendenza e chi quadrato 

#creo la suddivisione 
b <- seq(0,1,by=0.2)
#numero di intervalli 
q <- length(b)-1 
#creo un nuovo grafico 
plot(b,b,type="n") 
#linee orizz. tratteggiate 
for ( i in 1:length(b) ) abline(h=b[i],lty=2) 
#linee vert. tratteggiate
for ( i in 1:length(b) ) abline(v=b[i],lty=2) 
#tabella con le appartenenze rispettivamente a righe e colonne definite 
#dagli intervalli dei cut 
freq <- table(cut(u[-1], breaks = b), cut(u[-length(u)], breaks = b))
#numero di coppie di valori successivi (n-1)
nm1 <- sum(freq) 
#arrotondamento grafico freq: ora ho le freq. relative (4 cifre decimali)
round(freq/nm1, 4)
#frequenze assolute attese
attese <- nm1/q^2 
#statistica test 
X2 <- sum((freq - attese)^2/attese) 
#vettore con statistica test e alfa osservato in base alle ipotesi assunte in 
#partenza rifiuto o meno l'ipotesi di indipendenza 
c(X2, 1 - pchisq(X2, q^2 - 1)) 

#funzione di autocorrelazione (ACF function)
acf(u, ylim = c(-1, 1))

#diagramma di dispersione
plot(u[-length(u)], u[-1], pch = ".")





#ALGORITMO DI WHICHMANN-HILL
#generazione di vettori aleatori sfruttando l'algoritmo di Whichmann-Hill

.current.seed <- c(123, 456, 789) runif.wh <- function(n) {
  a <- c(171, 172, 170)
  b <- c(30269, 30307, 30323) s <- .current.seed
  u <- rep(0, n)
  for (i in 1:n) {
    s <- (a * s)%%b
    u[i] <- sum(s/b)%%1 }
  .current.seed <<- s
  u
}


#la funzione nativa di R che implementa la generazione di numeri casuali è 
runif(3)


#il corrispondente di .current.seed è .Random.seed, ma non conviene modificarlo
#per modificare il seme conviene usare 
set.seed(1)
#in questo modo verranno generati gli stessi numeri se verrà chiamata due volte
#runif()


#VERIFICA CHE EFFETTIVAMENTE I NUMERI GENERATI SIANO (PSEUDO)CASUALI

#1.1 verificare che la funzione di ripartizione empirica corrisponda alla 
#funzione di ripartizione teorica di una distribuzione uniforme:
u <- runif.wh(1e+05)
plot(ecdf(u), do.points = FALSE) 
curve(punif(x), col = "red", add = TRUE)

#1.2 verificare che i quantili empirici siano vicini a quelli teorici di una
#variabile uniforme
plot(sort(u), ((1:length(u)) - 0.5)/length(u)) 
abline(0, 1) #disegna una retta del tipo y=x

#2 test di indipendenza e test di bontà di adattmento(contemporaneamente)
b <- seq(0,1,by=0.2) #creo una suddivisione di passo 0.2
b
q <- length(b)-1
q
plot(b,b,type="n") #creo un grafico vuoto
for ( i in 1:length(b) ) abline(h=b[i],lty=2) #linee orizz. tratteggiate 
for ( i in 1:length(b) ) abline(v=b[i],lty=2) #linee vert. tratteggiate
#per confrontare le frequenze osservate e le attese con test chi-quadrato
freq <- table(cut(u[-1], breaks = b), cut(u[-length(u)], breaks = b))
#questo input crea una tabella con asse x gli elementi che appartengono 
#alla prima classe definita nel primo cut, e sulla y gli elementi che 
#appartengono alla classe definita nel secondo cut 
nm1 <- sum(freq) # Questo e' il numero di coppie di valori successivi 
# (u_i, u_{i-1}) ossia n-1
round(freq/nm1, 4) # Dovrebbero essere valori vicini a 1/q^2. (0.04)
attese <- nm1/q^2
X2 <- sum((freq - attese)^2/attese)
c(X2, 1 - pchisq(X2, q^2 - 1)) #livello di significatività osservato 

#3 funzione di autocorrelazione
acf(u, ylim = c(-1, 1))

#4 diagramma di dispersione 
plot(u[-length(u)], u[-1], pch = '.')



#funzione cut. (Per convertire una variabile continua in una categoriale)
x <- 1:10
x
a <- cut(x, breaks = c(0, 5, 10))
a
#restituisce se il numero è contenuto nell'intervallo (0,5] o in (5,10]
is.factor(a) #verifica che a sia una variabile di tipo factor (return: bool)




#GENERAZIONE PER INVERSIONE

#basta mettere come input della funzione quantile un vettore di numeri ( in 
#[0,1]generati casualmente da un algoritmo qualsiasi. Il codominio della 
#funzone di ripartizione, ossia [0,1], sarà il dominio della funzione quantile,
#mentre il dominio, ossia {ℝ} della funzione di ripartizione diventerà il 
#codominio della funzione quantile.
rlaplace <- function(n, a = 0, b = 1) { 
  qlaplace(runif.wh(n), a, b)
}

#verifico la casualità:
#FdR empirica assomiglia a FdR teorica?
y <- rlaplace(1e+05)   #1e+05 = 100'000
plot(ecdf(y), do.points = FALSE, xlim = c(-6, 6)) #FdR empirica 
curve(plaplace(x), -10, 10, add = TRUE, col = "red") #FdR teorica 




#DISTRIBUZIONE DI LAPLACE
#fondamentale per generare numeri con il metodo dell'inversione: ha ottime 
#caratteristiche in quanto assomiglia a una normale, ma ha una FdR esplicita 
#e facilmente invertibile.

dlaplace <- function(y, a=0, b=1){
  exp(-abs(y-a)/b)/(2*b)
}

curve(dlaplace(x, 0, 1), -10, 10, lty = 'solid')
curve(dlaplace(x, 0, 4), -10, 10, lty = "dashed", add = TRUE)
curve(dlaplace(x, 4, 1), -10, 10, lty = "dotted", add = TRUE)


#Laplace FdR (trovata algebricamente, vedi appunti p. 33)
plaplace <- function(y, a=0, b=1){
  (1+sign(y-a)*(1-exp(-abs(y-a)/b)))/2
}
curve(plaplace(x), -10, 10)

#funzione quantile Laplace (trovata algebricamente, vedi appunti p. 34)
#non si può usare if() per argomenti vettoriali, per questo si usa sign()
qlaplace <- function(p, a = 0, b=1){
  a - b * sign(p - 0.5) * log(1 - 2 * abs(p - 0.5))
}

#alcuni grafici
#densità
curve(plaplace(x,0,1),-10,10,lty=1)
curve(plaplace(x,0,4),lty=2,add=TRUE)
curve(plaplace(x,4,1),lty=2,add=TRUE)
#funzione quantile
curve(qlaplace(x,0,1),0,1,lty=1,ylim=c(-10,10))
curve(qlaplace(x,0,4),lty=2,add=TRUE)
curve(qlaplace(x,4,1),lty=2,add=TRUE)



#GENERAZIONE DA UNA DISTRIBUZIONE DISCRETA

x = c(1,2,4,7)

pr = c(0.2,0.4,0.1,0.3)

xx=cut(runif(1000),breaks=c(0,cumsum(pr)),
       labels=c(1,2,4,7))

table(xx)


#generazione da una bernouolliana 

rbern1 = function(n,p) {
  ris = rep(NA, n)
  for (i in 1:n) {
    u = unif(1)
    if (u<p)
      ris[i] = 1
    else 
      ris[i] = 0
  }
  return (ris)
}


#per rendere l'algoritmo (molto) più efficiente:
rbernv2 = function (n,p){
  ris = rep(NA, n)
  u = runif(n)
  ris[u<p] = 1
  ris[u>=p] = 0
  return(ris)
}


#facciamo finta di dimenticarci dell'esistenza di cut() e as.double
rx = function (){
  u= runif(1)
  if (u<0.2)
    x = 1
  else if (u<0.6)
    x = 2
  else if (u<0.7)
    x = 2
  else 
    x = 7
  return (x)
}


#GENERAZIONE DA UNA DISCRETA CON SUPPORTO NON FINITO 

#FdR Poisson
F = function(x) ppois(x, 3) 

#funzione generatrice 
generaF = function() {
  u = runif(1)
  valore = 0
  while (F(valore) < u) {
    valore = valore + 1 
  }
  return(valore) 
}

#il programma funziona ma è molto inefficiente, infatti se dovessimo generare
#più osservazioni e verificassimo il tempo di esecuzione
x = rep(NA, 10000) #creo il vettore per la preallocazione
#riempio il vettore usando generaF
system.time(for (i in 1:10000) x[i] = generaF())

#Ecco una versione più efficiente:
#Con questa funzione genero contemporaneamente più numeri, con quella precedente
#ne generavo uno e la ripetevo con un ciclo for per ottenere un vettore.
#L'indicizzazione logica aiuta molto l'efficienza di questo algoritmo.
F = function(x) ppois(x, 3) 
generaFn = function(n) {
  u = runif(n) 
  valore = rep(0, n) 
  repeat {
    ind = F(valore) < u
    valore[ind] = valore[ind] + 1
    if (all(!ind))
      break
  }
  return(valore) 
}


#è evidente controllando i tempi di esecuzione degli algoritmi che il secondo 
#ci mette molto meno

#tempo 1° algoritmo
x = rep(NA, 10000)
system.time(for (i in 1:10000) x[i] = generaF())

#tempo 2° algoritmo
system.time(generaFn(10000))



#verifica di casualità#confronto la FdR empirica (a.1) con la FdR teorica(a.2)plot(ecdf(y),do.points=FALSE,xlim=c(-4,4)) #(a.1)curve(pnorm(x),-4,4,col="red",add=TRUE) #(a.2)#confronto la densità teorica (b.1) con la densità empirica (b.2) sotto forma di #istogramma delle frequenze relaativehist(y,prob=TRUE,nclass=100) #(b.1)curve(dnorm(x),-5,5,col='green',add=TRUE) #(b.2)acf(y, 100, ylim = c(-1, 1), demean = FALSE) #verifico l'ndipendenza con ACF#ALGORITMO A FETTE (SLICE SAMPLING)f=function(x) dbeta(x,2,3) #distribuzione beta x=y=rep(NA,100)x[1]=0.5y[1]=0.1 simsufy=function(y) {  u=runif(1,0,1)  while (f(u)<y) u=runif(1,0,1)  u}for (i in 2:10000){  y[i]=runif(1,0,f(x[i-1]))  x[i]=simsufy(y[i]) }y #stampo i valori generati#per distribuzione beta vedi: https://it.wikipedia.org/wiki/Distribuzione_Beta


