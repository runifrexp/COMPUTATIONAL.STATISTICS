#SIMULAZIONI E APPLICAZIONI



#INTEGRAZIONE MONTE CARLO

monte.carlo <- function(n, f, B = 50000, print.details = TRUE, alpha = 0.01) {
  y <- sapply(rep(n, B), function(n) f(runif(n))) #sapply è un alternativa a un 
  #ciclo for 
  ybar <- mean(y) #è già l'approssimazione dell'integrale 
  if (print.details) { 
    ic <- qnorm(1 - alpha/2) * sd(y)/sqrt(B)
    cat("La stima di I(f) e'", ybar,
        "\nIl vero I(f) e' compreso tra ", ybar - ic, "e", ybar + ic,
        "\n con probabilita'",
        1 - alpha, "\n")
  }
  invisible(ybar) #resstituisce un valore senza stamparlo
}



#esempio 1 
monte.carlo(100, function(x) 3 * sum(x^2))
#con 1000 elementi si ha più precisione (ci mette molto)
monte.carlo(1000, function(x) 3 * sum(x^2))



#esempio 2
f <- function(x) { 
  x <- x - 0.5 
  sum(outer(x, x)) - sum(x^2) 
}
monte.carlo(100, f)



#esempio 3 (calcolare pi tramite simulazioni MC)
f.pi <- function(x) {
  d <- (x[1] - 0.5)^2 + (x[2] - 0.5)^2 
  if (4 * d < 1)
    4 else 0
}

#disegno del grafico 
xx <- seq(0, 1, length = 200) #crea una sequenza di 200 numeri da 0 a 1
yy <- seq(0, 1, length = 200) #crea una sequenza di 200 numeri da 0 a 1
zz = matrix(NA, 200, 200) #matrice che verrà riempita con le coppie dei 200 
for (i in 1:200) {        #valori di xx e i 200 valori di yy 
  for (j in 1:200) {
    zz[i, j] = f.pi(c(xx[i], yy[j])) 
  } 
}

par(mfrow = c(1, 2), mar = c(2, 2, 0.3, 0.3)) #mfrow prepara la finestra grafica
#I quattro valori corrispondono ai margini inferiore,sinistro,superiore e destro
image(xx, yy, zz) #stampa la matrice (a sx)   #mar prepara i margini
contour(xx, yy, zz) #curva di livello (a dx)


monte.carlo(10000, f.pi)

monte.carlo(2, f.pi)

monte.carlo(2, f.pi)










#esempio 3
f.pi <- function(x) {
  d <- (x[1] - 0.5)^2 + (x[2] - 0.5)^2 
  if (4 * d < 1)
    4 else 0
}

xx <- seq(0, 1, length = 200) 
yy <- seq(0, 1, length = 200) 
zz = matrix(NA, 200, 200)
for (i in 1:200) {
  for (j in 1:200) {
    zz[i, j] = f.pi(c(xx[i], yy[j]))
  } 
}

par(mfrow = c(1, 2), mar = c(2, 2, 0.3, 0.3)) #mfrow prepara la finestra grafica
image(xx, yy, zz)                             #mar prepara i margini
contour(xx, yy, zz) #curva di livello 


monte.carlo(10000, f.pi)

monte.carlo(2, f.pi)

monte.carlo(2, f.pi)





#TEOREMA DEL LIMITE CENTRALE

#dato che vengono usati risultati asintotici per procedimenti inferenzaili
#con n finito, bisogna verificare che l'alfa nominale sia effettivamente 
#l'alfa reale. 
#per fare ciò definiamo i vari psi, ossia le vere probabilità di rifiutare H_0
#per n finito
#in questo caso il calcolo è esatto in quanto consociamo la distribuzione per n 
#finito
psi.basso.normale=function(n,alfa) pt(qnorm(alfa),n-1)
psi.alto.normale=function(n,alfa) 1-pt(qnorm(1-alfa),n-1)
psi.bilat.normale=function(n,alfa) 2*pt(qnorm(alfa/2),n-1)

#verifica grafica per capire per che n gli psi si avvicinano agli alfa 
#l'unica differenza tra le prime tre righe di codice è l'alfa dichiarato
plot(2:30,psi.basso.normale(2:30,0.1), ylim=c(0,0.25),pch=20, col = 'red') 
points(2:30,psi.basso.normale(2:30,0.05),pch=20, col = 'yellowgreen')
points(2:30,psi.basso.normale(2:30,0.01),pch=20, col = 'sienna')
abline(h = c(0.1, 0.05, 0.01), lty = 2)



#ora affrontiamo una situazione in cui il calcolo non è esatto:
#voglio conoscere la distribuzione della statistica test. 
#procediamo generando MC n T_n: la frazione dei casi in cui rifiutiamo H_0 è psi

#partiamo dalla generazione di B T_n
molte.t <- function(n, rf, mu, B = 10000) { 
  y <- matrix(rf(n * B), n, B) #genero n*B valori da una Fisher
  ybar <- apply(y, 2, mean) #applica media alle colonne di y
  s <- apply(y, 2, sd) #applica standard deviation alle colonne di y
  sqrt(n) * (ybar - mu)/s #formula di T_n
}

#funzionamento di apply()
a <- matrix(1:6, 2, 3)
a
apply(a, 2, mean) 
apply(a, 1, mean)
#il primo apply restituisce la media delle colonne (2) della matrice
#il secondo apply restituisce la media delle righe (1) della matrice


#le righe importanti sono 3, il resto è una ripetizione per altri casi.
#verranno indicate con #***
conta.t <- function(n, rf, mu, B = 10000, alpha = c(0.01, 0.05, 0.1)) {
  tn <- molte.t(n, rf, mu, B) #*** generazione delle statistica test
  zbasso <- qnorm(alpha) #*** soglia di rifiuto (f. quantile normale std)
  zalto <- -zbasso
  zbil <- qnorm(1 - alpha/2)
  basso <- alto <- bil <- rep(0, length(alpha)) 
  for (i in 1:length(alpha)) {
    basso[i] <- mean(tn <= zbasso[i]) #*** frazione dei casi in cui rifiuto sui
    alto[i] <- mean(tn >= zalto[i])   #casi totali
    bil[i] <- mean((tn <= -zbil[i]) | (tn >= zbil[i]))
  }  
  ans <- list()
  ans$basso <- cbind(alpha, basso, sqrt(basso * (1 - basso)/B))                                                
  colnames(ans$basso) <- c("alpha", "psi", "se")
  ans$alto <- cbind(alpha, alto, sqrt(alto * (1 - alto)/B))
  colnames(ans$alto) <- c("alpha", "psi", "se")
  ans$bilaterale <- cbind(alpha, bil, sqrt(bil * (1 - bil)/B))
  colnames(ans$bilaterale) <- c("alpha", "psi", "se")
  ans
}

#alcuni esempi
conta.t(5, runif, 0.5)
conta.t(30, runif, 0.5)
conta.t(100, runif, 0.5)

conta.t(5, rexp, 1)
conta.t(30, rexp, 1)
conta.t(10000, rexp, 1)

#verifica se la distribuzione delle T_n generate assomiglia a una Z
hist(molte.t(10, runif, 0.5), nclass = 100, prob = TRUE) #densità empirica
curve(dnorm(x), -4, 4, col = "red", add = TRUE) #densità teorica
qqnorm(molte.t(10, runif, 0.5))#quantili teorici rispetto a quantili empirici:
abline(0,1, col = 'purple', lwd = 2)    #sono allineati?

#aumentando n l'approssimazione è migliore
hist(molte.t(30, runif, 0.5), nclass = 100, prob = TRUE)
curve(dnorm(x), -4, 4, col = "red", add = TRUE)
qqnorm(molte.t(30, runif, 0.5))
abline(0,1, col = 'purple',lwd = 2)

#verifica con esponenziale
hist(molte.t(10, rexp, 1), nclass = 100, prob = TRUE)
curve(dnorm(x), -4, 4, col = "red", add = TRUE)
qqnorm(molte.t(10, rexp, 1))
abline(0,1, col = 'purple',lwd = 2)





#MISTURA DI NORMALI

#FdR di una mistura di normali
pmist2n <- function(y, mu = 0, sigma = 1, pi = 0, eta = 0, lambda = 1) {
  (1 - pi) * pnorm(y, mu, sigma) + pi * pnorm(y, eta, lambda)
}

#densità di una mistura di normali
dmist2n <- function(y, mu = 0, sigma = 1, pi = 0, eta = 0, lambda = 1) {
  (1 - pi) * dnorm(y, mu, sigma) + pi * dnorm(y, eta, lambda)
}

#generazione da una mistura di normali
rmist2n <- function(n, mu = 0, sigma = 1, pi = 0, eta = 0, lambda = 1) {
  v <- rbinom(n, 1, pi) + 1 #genero 1 o 2 da una bernoulliana
  c(mu, eta)[v] + c(sigma, lambda)[v] * rnorm(n) 
}

#confronto tra la densità teorica della nromale e quella della mistura
curve(dnorm(x), -4, 14)
curve(dmist2n(x, pi = 0.03, eta = 10, lambda = 2), add = TRUE, col = "red")

#c'è però da stare attenti ai classici stimatori media campionaria e varianza
#campionaria corretta perchè sono influenzati molto alla presenza di pochi 
#valori anomali
#per verificarlo generiamo da una normale std e da una distribuzione contaminata
#generazione di 10000 campioni di dimensione 50 
B <- 10000
n <- 50
y.buoni <- matrix(rnorm(n * B), n, B) #campione da nomrale std
y.cont <- matrix(rmist2n(n * B, pi = 0.03, eta = 10, lambda = 2), n, B)
#campione da distribuzione contaminata

#applico media e d std ai due campioni 
ybar.buoni <- apply(y.buoni, 2, mean)
ybar.cont <- apply(y.cont, 2, mean)
sd.buoni <- apply(y.buoni, 2, sd)
sd.cont <- apply(y.cont, 2, sd)


#confronto tra densità della media campionaria non contaminata e contaminata
par(mfrow = c(1, 2))
hist(ybar.buoni, nclass = 30, xlim = c(-1, 2))
hist(ybar.cont, nclass = 30, xlim = c(-1, 2))


#confronto tra densità della deviazione standard campionaria non contaminata 
#e contaminata
par(mfrow = c(1, 2))
hist(sd.buoni, nclass = 30, xlim = c(0, 5))
hist(sd.cont, nclass = 30, xlim = c(0, 5))

#confronto con qqplot
qqplot(ybar.buoni,ybar.cont); abline(0,1)
qqplot(sd.buoni,sd.cont); abline(0,1)




#NUOVI STIMATORI 

B <- 10000
n <- 50
y.buoni <- matrix(rnorm(n * B), n, B) #campione da nomrale std
y.cont <- matrix(rmist2n(n * B, pi = 0.03, eta = 10, lambda = 2), n, B)
#campione da distribuzione contsaminata

#applico mediana e mad(median absolute deviation)
med.buoni <- apply(y.buoni, 2, median)
med.cont <- apply(y.cont, 2, median)
mad.buoni <- apply(y.buoni, 2, mad)
mad.cont <- apply(y.cont, 2, mad)


par(mfrow = c(1, 2))
hist(med.buoni, nclass = 30, xlim = c(-1, 2))
hist(med.cont, nclass = 30, xlim = c(-1, 2))
#si nota che il coportamento non cambia tanto anche se sono state modificate le 
#condizioni iniziali: questo è sintomo di robustezza di mediana e mad

par(mfrow = c(1, 2))
hist(mad.buoni, nclass = 30, xlim = c(0, 5))
hist(mad.cont, nclass = 30, xlim = c(0, 5))
#si nota lo stesso anche per mad

qqplot(med.buoni,med.cont); abline(0,1)
qqplot(mad.buoni,mad.cont); abline(0,1)


#scarsa efficienza di med e mad (ma ottima robustezza)

sqrt(c(mean(ybar.buoni^2), mean(med.buoni^2), mean(ybar.cont^2),
       mean(med.cont^2)))
#si nota che in condizioni ottimali si comportano meglio gli stimatori classici
#per med e mad non cambia RMSE con o senza contaminazione
#per gli stimatori classici aumenta RMSE se contaminati

sqrt(c(mean((sd.buoni - 1)^2), mean((mad.buoni - 1)^2),
       mean((sd.cont - 1)^2), mean((mad.cont - 1)^2)))


#SETACCIO
#la scarsa efficienza rende pcoo attraenti med e mad, si possono però combinare 
#con gli stimatori classici per ottenere stimatori che uniscono robustezza
#ed efficienza

#per la media
mu.cappello <- function(y, k = 2.5) {
  idx.buone <- abs(y - median(y)) <= k * mad(y) 
  mean(y[idx.buone])
}

#per la varianza
sigma.cappello <- function(y, k = 2.5) {
  idx.buone <- abs(y - median(y)) <= k * mad(y)
  sd(y[idx.buone]) 
}
#il setaccio opera così:
#calcolare mediana e mad → determinare “osservazioni buone”; 
#calcolare media e s.q.m. delle osservazioni “dichiarate” buone.

#verifico se la promessa della combinazione di efficienza e robustezza viene 
#mantenuta
mc.buoni <- apply(y.buoni, 2, mu.cappello)
mc.cont <- apply(y.cont, 2, mu.cappello)
sqrt(c(mean(ybar.buoni^2), mean(med.buoni^2), mean(mc.buoni^2),
       mean(ybar.cont^2), mean(med.cont^2), mean(mc.cont^2)))

sc.buoni <- apply(y.buoni, 2, sigma.cappello)
sc.cont <- apply(y.cont, 2, sigma.cappello)
sqrt(c(mean((sd.buoni - 1)^2), mean((mad.buoni - 1)^2),
       mean((sc.buoni - 1)^2), mean((sd.cont - 1)^2),
       mean((mad.cont - 1)^2), mean((sc.cont - 1)^2)))

#sembra di si: si alza di poco il RMSE dal primo all'ultimo elemento del vettore
#l'inefficienza del setaccio è dovuta al fatto che anche esso può sbagliare
#valutando come non contaminate delle variabili che invece lo sono




#TEST MULTIPLI, ESPONENZIALE 


#verifica d'ipotesi: 
#H_0: mu_1 = mu_2
#statistica test W: doppia caduta di log-verosimiglianza: misura la distanza di
#r = mu_1 / mu_2 da 1, dove 1 è l'ipotesi nulla mu_1 = mu_2
curve(4 * log((sqrt(x) + 1/sqrt(x))/2), 0, 6)


#genero molte W
rv.sim <- function(n) {
  r <- sqrt(mean(rexp(n))/mean(rexp(n))) 
  4 * n * log((r + 1/r)/2)
}
W <- sapply(rep(5, 1e+05), rv.sim)


#verifica se i risultati asintotici valgono anche per n=5, si.
hist(W, nclass = 100, prob = TRUE, xlim = c(0, 10))
curve(dchisq(x, 1), 0, 10, add = TRUE, col = "red")



table(W <= qchisq(0.95, 1))/1e+05
#si nota che l'alfa dichiarato è molto vicino alla vera frazione di valori di W
#inferiori al quantile 1-alfa (=0.95)


#inizializzo la matrice y
y <- rbind(c(1,2,10,18,19,20,26,27,35,82,136,369),
           c(6,13,13,39,51,56,59,69,77,97,132,164),
           c(1,2,37,38,41,46,80,96,118,195,230,324),
           c(18,40,45,68,79,151,152,179,195,205,276,303),
           c(0,28,50,101,108,120,123,210,223,260,266,310),
           c(9,40,44,53,61,106,115,139,147,149,443,515),
           c(2,18,43,48,67,75,83,161,229,352,588,615))
rownames(y) <- c("A","B","C","D","E","F","G") #assegno i nomi alle colonne


apply(y, 1, mean) #applica la media alle righe della matrice


#la seguente funzione calcola tutti i test rapporto di verosimiglianza a coppie
#ossia tutte le statistica test di Wilks per tutti i confronti a coppie
rv2exp <- function(y) {
  ybar <- apply(y, 1, mean)
  r <- outer(ybar, ybar, function(yi, yj) sqrt(yi/yj)) 
  4 * NCOL(y) * log((r + 1/r)/2)
}

round(rv2exp(y), 2) #arrotonda la matrice 7x7 a due cifre decimali dopo la
#virgola dopo aver applicato rv2exp
#la matrice è ridondante: la triangolare superiore è uguale alla triangolare 
#inferiore 

c95 <- qchisq(0.95, 1) #FdR della chi-quadrato 1 grado di libertà, quantile 0.95
c95

rv2exp(y) > c95 #restituisce una matrice 7x7 di T e F in base al fatto se i 
#valori della matrice sono maggiori di c95

#si nota che 7 valori superano il valore critico, quindi dato che bisogna 
#rifiutare H_0 se un solo valore supera la soglia del quantile 0.95, rifiuto H_0



#con questa simulazione si viene a sapere la vera probabilità di errore di primo
#tipo 

#pwmax è una FdR
#x è il valore in corrispondenza al quale calcolo la FdR
#k è il numero di gruppi = 7
#n è il numero di colonne, ossia il numero di osservazioni per gruppo, quindi 12
#B è il numero di volte che genero il valore 
pwmax <- function(x, k, n, B = 10000) {
  f <- function(n) {
    y <- matrix(rexp(n * k), k, n)
    max(rv2exp(y)) #calcola la statistica test per ogni coppia di numeri e ne 
    #trova il max
  }
  maxW <- sapply(rep(n, B), f) #per B volte calcolo W_max
  sapply(x, function(x) mean(maxW <= x))#mean(maxW <= x) calcola la frazione di 
  #valori a sinistra di x, ossia l'approssimazione di avere valori a sx di x sotto
  #H_0, mean vuole input scalari, quindi bisogna prima passare la funzione sapply,
  #che perette di prendere input vettoriali
}

#in questo modo possiamo possiamo ottenere stime dell’errore di I tipo: per 2 
#gruppi con 12 osservazioni per gruppo, la prob. di errore del I tipo con soglia
#c95 è
1 - pwmax(c95, 2, 12)

c1 <- 1 - pwmax(c95, 3, 12)

c2 <- 1 - pwmax(c95, 5, 12)

c3 <- 1 - pwmax(c95, 7, 12) #vale la pena tirare una moneta...

cbind(c1,c2,c3)
#è facile norare che la stima della P di errore di primo tipo aumenta con 
#l'aumentare del numero di gruppi, non dal numero di osservazioni dentro 
#ogni singolo gurppo



#il punto è che stiamo usando la soglia sbagliata,
#la soglia corretta con cui confrontare le varie statistiche test W_ij è il
#quantile della distribuzione di W_max.
#la funzione qwmax stima i quantili di W_max

qwmax <- function(p, k, n, B = 10000) { 
  f <- function(n) {
    y <- matrix(rexp(n * k), k, n)
    max(rv2exp(y)) #mi da un valore di W_max
  }
  quantile(sapply(rep(n, B), f), p) #ordina i valori e restituisci quello che 
  #lascia a sx una frazione p di valori 
}

qwmax(0.95, 7, 12)
#Utilizzando questo valore nessuna delle differenze nei dati di prima può 
#essere considerata significativa. Non abbiamo elementi, in altre parole, per 
#dubitare che l’affidabilità delle componenti prodotte dalle 7 marche sia 
#differente.


wmax <- max(rv2exp(y))
lso <- 1 - pwmax(wmax, 7, 12, B = 10000)
c(wmax, lso, sqrt(lso * (1 - lso)/10000))
#wmax è stat. test 
#lso è liv. signigficatività osservato
#ultimo elemento è l'incertezza (P di avere valori più grandi di W_max_oss |H_0)



#PROVA CLINICA

#esempio su tossicità in una prova clinica

#funzione di monitoraggio sequenziale basata su test alla Wald su v.a. binomiali
#questa funzione serve per verificare che man mano che arrivano nuove 
#osservazioni va fermato il test o meno (test sequenziale)
monitor.binom <- function(y, N, theta0, k) { 
  tn <- seq(1, N) #scorrere del tempo (indice i)
  media <- theta0 * tn #cambia all'aumentare di tn (sequenziale)
  sqm <- sqrt(tn * theta0 * (1 - theta0)) #evolversi dello SQM
  limiti <- media + k * sqm #evolversi della soglia g_n
  sn <- cumsum(y) #somma cumulata dei pazienti intossicati
  #da qui in poi ci sono istruzioni grafiche
  plot(sn, xlim = c(0, N), ylim = range(0, max(sn, limiti)))
  lines(limiti, lty = "dashed", lwd = 2)
  grid()
}

#prova sotto H_0, si nota che l'andamento è quello desiderato
y <- rbinom(200, 1, 0.1)
monitor.binom(y, 200, 0.1, 4)


#in questo caso non siamo sotto H_0 (alfa = 0.5 ≠ alfa_dichiarato = 0.1 )
y <- rbinom(200, 1, 0.5)
monitor.binom(y, 200, 0.1, 4)
#l'andamento, come ci si poteva aspettare è che la tossicità della prova
#è troppo alta




#La soglia può portarci a sbagliare, ossia di interrompere la prova anche quando
#non andrebbe interrotta. 
#non vogliamo ciò, quindi va scelto un k t.c. 
#P_theta_0(bloccare erroneamente lo studio) = alpha

#SCELTA DELLA SOGLIA 

#La seguente funzione di R calcola, via Monte Carlo, una approssimazione di 
#F^−1(alpha), in quanto ci serve che k sia >= quantile (1-alpha) della
#distribuzione max{T_n}

qmonitor.binom <- function(p, N, theta0, B = 20000) { 
  tn <- seq(1, N)
  media <- tn * theta0
  sqm <- sqrt(tn * theta0 * (1 - theta0))
  y <- matrix(rbinom(N * B, 1, theta0), N, B) #matrice NxB di campioni generati
  #da variabili uniformi sotto H_0
  sn <- apply(y, 2, cumsum) #statistiche sufficienti, n° di pazienti cumulati 
  #ad ogni passo
  a <- apply(sn, 2, function(col) max((col - media)/sqm)) #ecco le statistica 
  #test max{T_n}
  quantile(a, p)
}

alpha <- 0.05
qmonitor.binom(1 - alpha, 200, 0.1)

#per verificare l'accuratezza del risultato, ripetiamo.
#riptetendo il comando si nota che il valore restituito è simile a quello 
#precedente. Da ciò capiamo che B è abbastanza grande 
qmonitor.binom(1 - alpha, 200, 0.1)


#possiamo inoltre effettuare una verifica grafica che ci fa capire qunata 
#densità c'è nei vari valori plausibili, e quanto grande è il range nel quale 
#vengono generati questi valori:
#questo controllo serve per avere sotot contorllo l'incertezza Monte Carlo, 
#si chiama Monte Carlo su Monte Carlo. 
#è molto lento poichè genera 20'000 valori per 1000 volte 
par(mar = c(2, 2, 0.2, 0.2), cex = 1.3)
rip = rep(NA, 1000)
h <- for (i in 1:1000) rip[i] = qmonitor.binom(1 - alpha, 200, 0.1) 
hist(rip, cex = 1.3)




#ecco un esempio in cui B è troppo basso
qmonitor.binom( 1-alpha , 200, 0.1, B=100)

#ripeto il comando 
qmonitor.binom( 1-alpha , 200, 0.1, B=100)
#noto che il risultato è molto diverso dal primo ottenuto


#con il controllo grafico la variabilità è ancora più evidente  
par(mar=c(2,2,0.2,0.2),cex=1.3) rip=rep(NA,1000)
for(i in 1:1000)
  rip[i]=qmonitor.binom( 1-alpha , 200, 0.1,B=100) 
hist(rip,cex=1.3)



#questa funzione da un'idea visiva dell'andamento della tossicità nei pazienti
#basata sul test alla Wilks (confronto tra verosimiglianze, no SMV o 
#distribuzioni asintotiche)
monitor.rv <- function(y, N, theta0, theta1, C) { 
  V <- log((1 - theta1)/(1 - theta0))
  U <- log(theta1/theta0) - V
  limiti <- (C - seq(1, N) * V)/U
  sn <- cumsum(y)
  plot(sn, xlim = c(0, N), ylim = range(0, max(sn,
                                               limiti)), xlab = "n", ylab = "sn")
  lines(limiti, lty = "dashed", lwd = 2)
  grid()
}



#esempi come quelli effettuati sui test alla Wald

#sotto H_0
y <- rbinom(200, 1, 0.1)
monitor.rv(y, 200, 0.1, 0.2, 2)
#sta quasi sempre sotto la linea di tolleranza


#non sotto H_0
y <- rbinom(200, 1, 0.5)
monitor.rv(y, 200, 0.1, 0.2, 2)
#non sta mai sotto la linea di tolleranza




#pure in questo caso è necessario verificare la scelta della soglia come nel
#test alla Wald, fa più o meno la stessa cosa della funzione qmonitor.binom

qmonitor.rv <- function(p, N, theta0, theta1, B = 20000) { 
  V <- log((1 - theta1)/(1 - theta0))
  U <- log(theta1/theta0) - V
  Vn <- V * seq(1, N)
  y <- matrix(rbinom(N * B, 1, theta0), N, B)
  sn <- apply(y, 2, cumsum)
  a <- apply(sn, 2, function(col) max(col * U + Vn))
  quantile(a, p)
}


#noto che 20'000 è sufficiente in quanto i valori generati non si discostano 
#troppo tra loro.
qmonitor.rv(1 - alpha, 200, 0.1, 0.2)
qmonitor.rv(1 - alpha, 200, 0.1, 0.2)
qmonitor.rv(1 - alpha, 200, 0.1, 0.2)



#ritorna rl, run lenght, ossia il punto fin cui arriva lo studio prima di 
#fermarsi perchè rileva un'anomalia

#per il test alla Wald
rlmonitor.binom <- function(N, theta0, k, vero.theta, B = 20000) {
  tn <- seq(1, N)
  media <- tn * theta0
  sqm <- sqrt(tn * theta0 * (1 - theta0)) 
  limiti <- media + k * sqm
  rl <- function(y) {
    sn <- cumsum(y)
    a <- which(sn > limiti)
    if (length(a))
      min(a) else N + 1
  }
  y <- matrix(rbinom(N * B, 1, vero.theta), N, B)
  apply(y, 2, rl) }



#per il test alla Wilks

rlmonitor.rv <- function(N, theta0, theta1, C, vero.theta, B = 20000) {
  V <- log((1 - theta1)/(1 - theta0)) 
  U <- log(theta1/theta0) - V
  limiti <- (C - seq(1, N) * V)/U
  rl <- function(y) {
    sn <- cumsum(y)
    a <- which(sn > limiti)
    if (length(a))
      min(a) else N + 1
  }
  y <- matrix(rbinom(N * B, 1, vero.theta), N, B)
  apply(y, 2, rl) 
}


#Le soglie che fanno sì che la prob. di allarme sia 0.05 se 
#theta = theta_0 = 0.1 sono

k <- 3.27 #per Wald
C <- 2.74 #per Wilks

#
a <- rlmonitor.binom(200, 0.1, k, 0.1)
b <- rlmonitor.rv(200, 0.1, 0.2, C, 0.1)
c(mean(a < 201), mean(b < 201))
#ritorno la stima di alpha fatta nei due modi


#ripeto l'esperimento per osservare l'incertezza MC
a <- rlmonitor.binom(200, 0.1, k, 0.1)
b <- rlmonitor.rv(200, 0.1, 0.2, C, 0.1)
c(mean(a < 201), mean(b < 201))



#da qui in poi verranno trattati casi con theta != theta_0, non sotto H_0
#theta è P(tossico)

#theta = 0.12
a <- rlmonitor.binom(200, 0.1, k, 0.12)
b <- rlmonitor.rv(200, 0.1, 0.2, C, 0.12)
#stima di gamma (0.12)
gamma <- c(mean(a < 201), mean(b < 201))
gamma
#rv migliore
#stima di rho(0.12)
rho <- c(mean(a[a < 201]), mean(b[b < 201]))
rho
#binom migliore 

#risultato complessivo 
(1 - gamma) * 200 + gamma * rho
#è meglio rv



#theta = 0.15
a <- rlmonitor.binom(200, 0.1, k, 0.15)
b <- rlmonitor.rv(200, 0.1, 0.2, C, 0.15)
gamma <- c(mean(a < 201), mean(b < 201))
gamma
rho <- c(mean(a[a < 201]), mean(b[b < 201]))
rho
#meglio rv in entrambi i casi
(1 - gamma) * 200 + gamma * rho



#theta = 0.2
a <- rlmonitor.binom(200, 0.1, k, 0.2)
b <- rlmonitor.rv(200, 0.1, 0.2, C, 0.2)
gamma <- c(mean(a < 201), mean(b < 201))
gamma
rho <- c(mean(a[a < 201]), mean(b[b < 201]))
rho
(1 - gamma) * 200 + gamma * rho



#theta = 0.25
a <- rlmonitor.binom(200, 0.1, k, 0.25)
b <- rlmonitor.rv(200, 0.1, 0.2, C, 0.25)
gamma <- c(mean(a < 201), mean(b < 201))
gamma
rho <- c(mean(a[a < 201]), mean(b[b < 201]))
rho
(1 - gamma) * 200 + gamma * rho



#theta = 0.3
a <- rlmonitor.binom(200, 0.1, k, 0.3)
b <- rlmonitor.rv(200, 0.1, 0.2, C, 0.3)
gamma <- c(mean(a < 201), mean(b < 201))
gamma
rho <- c(mean(a[a < 201]), mean(b[b < 201]))
rho
(1 - gamma) * 200 + gamma * rho



#theta = 0.5
a <- rlmonitor.binom(200, 0.1, k, 0.5)
b <- rlmonitor.rv(200, 0.1, 0.2, C, 0.5)
gamma <- c(mean(a < 201), mean(b < 201))
gamma
rho <- c(mean(a[a < 201]), mean(b[b < 201]))
rho
(1 - gamma) * 200 + gamma * rho


theta1 = c(seq(0.1, 0.25, by = 0.01), 0.3, 0.35, 0.4, 0.45, 0.5)
gamma = rho = pazienti = matrix(NA, length(theta1), 2)
for (i in 1:length(theta1)) {
  a <- rlmonitor.binom(200, 0.1, k, theta1[i])
  b <- rlmonitor.rv(200, 0.1, 0.2, C, theta1[i]) 
  gamma[i, ] <- c(mean(a < 201), mean(b < 201)) 
  rho[i, ] <- c(mean(a[a < 201]), mean(b[b < 201])) 
  pazienti[i, ] = (1 - gamma[i, ]) * 200 + gamma[i, ] * rho[i, ]
}


#confronto grafico 
par(mar = c(2, 2, 0.2, 0.2), mfrow = c(1, 3))
matplot(theta1, gamma[, c(1, 1, 2, 2)], type = "pl",
        col = c("black", "red")[c(1, 1, 2, 2)], pch = 20,
        lty = 1, ylim = c(0, 1), yaxs = "i")
matplot(theta1, rho[, c(1, 1, 2, 2)], type = "pl",
        col = c("black", "red")[c(1, 1, 2, 2)], pch = 20,
        lty = 1)
matplot(theta1, pazienti[, c(1, 1, 2, 2)], type = "pl",
        col = c("black", "red")[c(1, 1, 2, 2)], pch = 20,
        lty = 1)





#CAMPIONAMENTO PER IMPORTANZA


#Per aumentare l'efficienza di una procedura Monte Carlo posso sfruttare il 
#campionamento per importanza

#genero non più uniformemente ma da una distribuzione che segue la forma 
#dell'integrale che voglio stimare, in modo da dare più importanza alle parti
#della funzione che ha più peso nel calcolo dell'integrale 


#faccio MC su MC per abbassare l'incertezza
f = function(x) (5 * (1 - x)^4 * (1 + 0.3 * abs(sin(12 * pi * x)))) * (x < 1)
M = 3000
ris = matrix(NA, 3, M) 
B = 1000
for (i in 1:M) {
  ris[1, i] = mean(f(runif(B)))           #unif(0,1)
  ris[2, i] = 2 * mean(f(runif(B, 0, 2))) #unif(0,a)
  u = rbeta(B, 1, 4)
  ris[3, i] = mean(f(u)/dbeta(u, 1, 4))   #beta(1,4)
}


a = apply(ris, 1, mean)


b = apply(ris, 1, sd)

par(mfrow = c(3, 1))  # Imposta i margini destro e sinistro
# Traccia gli istogrammi dei tre metodi (usando distribuzioni diverse)
hist(ris[1,], xlim = c(min(ris), max(ris)), main = "Histogram of unif(0,1)")
hist(ris[2,], xlim = c(min(ris), max(ris)), main = "Histogram of unif(0,a)")
hist(ris[3,], xlim = c(min(ris), max(ris)), main = "Histogram of beta(1,4)")
#si nota che il terzo grafico è quello con una deviazione std minore 



#ESEMPIO SULLA MINI RETE DI COMUNICAZIONE

#X_i è una variabile bernoulliana che indica se l'i-esimo arco funziona o meno
#la funzione 'funziona' indica se dato in input un vettore binario X di 
#lunghezza 8 (ossia il numero di archi) la rete può funzionare o meno

funziona <- function(x) {
  a <- x[1] * (x[7] + x[4] * x[6]) + x[2] * x[6] + x[3] *(x[5] * x[6] + x[8]) 
  a>0
}
#i * vanno interpretati come 'e', mentre i + come 'oppure'
#ritorna T se funziona e F se non funziona

#esempi
funziona(c(1, 0, 0, 1, 0, 1, 0, 0))
funziona(c(1, 0, 0, 1, 0, 0, 0, 0))


#supponendo che:
# - il funzionamento di un arco sia indipendente dal funzionamento degli altri
# - la probabilità con cui un arco funziona in un determinato istante di tempo è
p.arco <- 0.99

#Il problema di cui ci occupiamo è di calcolare la probabilità con cui, in un
#certo istante di tempo, la rete non funzioni.

#questa funzione enumera tutti i 2^8 possibilis tati della rete  
stati <- function(M) { if (M == 1) {
  c(0, 1) } else {
    altri <- stati(M - 1)
    rbind(cbind(0, altri), cbind(1, altri)) }
}

#esempi
stati(1)
stati(2)
stati(3)

#quindi la probabilità che la rete non funzioni è 
a <- stati(8)
s <- apply(a, 1, sum)
f <- apply(a, 1, funziona)
p.stato <- (p.arco^s) * (1 - p.arco)^(8 - s)
p.rete <- sum((1 - f) * p.stato)
p.rete


#se avessimo invece di 8 archi un numero grande M, che rende il problema molto 
#più realistico, non si potrebbe più fare un calcolo esatto per via della 
#pesantezza computazionale del calcolo da fare.
#per questo si opera simulando via Monte Carlo.


#APPROCCIO NAIVE (INGENUO)
#la funziona 'prob.guasta' funziona così:
# - generiamo casualmente B possibili stati per la rete;
# - poniamo la stima della probabilità che ci interessa uguale alla 
#   proporzione di casi in cui abbiamo osservato, nel “campione”, uno stato 
#   corrispondente a “rete guasta”.

prob.guasta <- function(p.arco, M, funziona, B = 10000) { 
  x <- matrix(rbinom(M * B, 1, p.arco), B, M)
  f <- apply(x, 1, funziona)
  mean(1 - f)
}

#esempi (per vedere che l'approccio naive non è molto utile)
prob.guasta(p.arco, 8, funziona)
prob.guasta(p.arco, 8, funziona)
prob.guasta(p.arco, 8, funziona)
#gli esiti sono quasi tutti zero in quanto si stima una probsabilità molto bassa
#(a meno di non usare un B enorme)



#APPROCCIO CON CAMPIONAMENTO PER IMPORTANZA
#Generiamo possibili stati della rete facendo finta che ogni singolo arco si 
#guasti con una probabilità più elevata di quella vera.

#la probabilità con cui un arco funzioni sia più piccola di quella vera
q.arco <- 0.75

#la seguente funzione implementa una simulazione MC con campionamento per imp.
guasta.ma.importante <- function(p.arco, q.arco, M, funziona, B = 10000) {
  x <- matrix(rbinom(M * B, 1, q.arco), B, M)
  s <- apply(x, 1, sum)
  w <- (p.arco^s) * (1 - p.arco)^(M - s)/((q.arco^s) * (1 - q.arco)^(M - s))
  f <- apply(x, 1, funziona)
  mean((1 - f) * w)
}

#esempi (funzionati)
guasta.ma.importante(p.arco, q.arco, 8, funziona)
guasta.ma.importante(p.arco, q.arco, 8, funziona)
guasta.ma.importante(p.arco, q.arco, 8, funziona)
#si nota che avendo cambiato la densitò da cui si genera la precisione della 
#stima aumenta di molto 



#CONFRONTO DEI DUE METODI
#Per confrontare le due stime usiamo la media dell’errore relativo percentuale
N <- 100
a <- matrix(0, N, 2)
colnames(a) <- c("naive", "importante") 
for (i in 1:N) {
  a[i, 1] <- prob.guasta(p.arco, 8, funziona)
  a[i, 2] <- guasta.ma.importante(p.arco, q.arco, 8, funziona) 
}
a <- 100 * abs(a - p.rete)/p.rete
apply(a, 2, mean)
#si nota che l'errore commesso dalla procedura naive è estremamente più grande 
#dell'errore fatto 



#VARIABILI ANTITETICHE

#altro metodo per aumentare l'efficienza della simulazione MC

#Prendiamo X_1, X_2 ∼ IID(f)
#Assumendo B pari per comodità, generiamo poi B/2 valori x_1,1, . . . , x_1,B/2 
#e x_2,1, . . . , x_2,B/2 dalle due variabili casuali.

#si fa in questo modo perchè la varianza della media campionaria 'classica' è 
#maggiore della varianza della somma delle due medie campionarie se le due 
#variabili hanno una correlazione negativa (vedi teoria pag. 165)



g = function(x) x
Finv = function(x) qchisq(x, 3)

#Effettueremo M prove, con B simulazioni MC ad ogni prova. I risultati saranno 
#memorizzati in ris.
B = 1000
M = 1000
ris = array(NA, dim = c(2, M))

#otteniamo tramite generazione per inversione B repliche da una chi-quadro 
for (i in 1:M) { 
  u = runif(B)
  ris[1, i] = mean(g(Finv(u))) 
}
mean(ris[1, ])
sd(ris[1, ])

#ora precediamo con il emtodo di dividere la generazione in B/2 proveniente da 
#U ‰
correlaz = rep(NA, M) 
for (i in 1:M) {
  u = runif(B/2)
  ris[2, i] = mean(c(g(Finv(u)), g(Finv(1 - u))))
  correlaz[i] = cor(g(Finv(u)), g(Finv(1 - u)))
}
mean(ris[2, ])
sd(ris[2, ])
#in confronto al risultato precedente si abbassa per visa della correlazione 
#negativa

mean(correlaz)
#se vedo che ènegativa vuol dire che sto aumentando l'efficienza della stima




#questa funzione effettua alcune prove 
antitetiche = function(g = function(x) x, Finv = function(x) qchisq(x, 3), 
                       B = 1000, M = 1000) {
  ris = array(NA, dim = c(3, M))
  for (i in 1:M) {
    u = runif(B)
    ris[1, i] = mean(g(Finv(u)))
    ris[2, i] = mean(c(g(Finv(u[1:(B/2)])), g(Finv(1 - u[1:(B/2)]))))
    ris[3, i] = cor(g(Finv(u[1:(B/2)])), g(Finv(1 - u[1:(B/2)])))
    return(ris)
  } 
}

r = antitetiche(g=function(x) sin(0.8*pi*x),
                Finv=function(x) qbeta(x,2,2),
                B=1000,
                M=1000)

round(apply(r,1,mean),4)
#la media non cambia 
#output:
# - media calcolata con il metodo naive 
# - media calcolata con le variabili antitetiche 
# - correlazione
round(apply(r,1,sd),4)
#la differenza tra le variazne è molto bassa, infatti la correlazione è quasi 
#nulla 


r=antitetiche(g=function(x) sin(0.5*pi*x),
              Finv=function(x) qbeta(x,2,2),
              B=1000,
              M=1000)
round(apply(r,1,mean),4)

round(apply(r,1,sd),4)

#in questo caso invece, data la correlazione fortemente negativa, l'efficienza 
#aumenta di moto in quanto si abbassa sensibilmente la varianza passando dal
#procedimento naive a quello per variabili antitetiche 





#RAO-BALCKWELL E ANTITETICHE


#verrà affrontata la stima della media con 5 modi diversi:
# - MC diretta
# - variabili antitetiche
# - Rao-Blackwell
# - RB + anti.
# - funzione nativa di R


#stima MC diretta
diretta <- function(beta0 = 4, beta1 = -0.1, B = 10000) { 
  s <- rnorm(B)
  t <- exp(beta0 + beta1 * s) * rexp(B)
  mean(t) 
}



#variabili antitetiche

#verifica di correlazione negativa, quindi se ha senso oppure no usare le 
#variabili antitetiche 
s <- rnorm(10000)
mu <- exp(4 - 0.1 * s)
u <- runif(10000)
cor(-mu * log(u), -mu * log(1 - u))
#si, ha senso

antitetiche <- function(beta0 = 4, beta1 = -0.1, B = 10000){
  s <- rnorm(B/4)
  u <- runif(B/2)
  mu <- c(exp(beta0 + beta1 * s), exp(beta0 - beta1 * s))
  t <- -c(mu, mu) * log(c(u, 1 - u))
  mean(t)
}
antitetiche()



#Rao-blackwell
rao.blackwell <- function(beta0 = 4, beta1 = -0.1, B = 10000) { 
  s <- rnorm(B)
  mean(exp(beta0 + beta1 * s)) 
}


#RB+anti
rb.antitetica <- function(beta0 = 4, beta1 = -0.1, B = 10000) { 
  s <- rnorm(B/2)
  (sum(exp(beta0 + beta1 * s)) + sum(exp(beta0 - beta1 * s)))/B
}



#MC su MC per verificare che le stime siano simili
M <- 500
a <- matrix(0, M, 4)
colnames(a) <- c("diretta", "anti", "rb", "rb.anti") 
for (i in 1:M) {
  a[i, "diretta"] <- diretta()
  a[i, "anti"] <- antitetiche()
  a[i, "rb"] <- rao.blackwell()
  a[i, "rb.anti"] <- rb.antitetica()
}

apply(a, 2, mean)

#verifico quale dei metodi si porta dietro la minore variabilità
apply(a, 2, sd)
#è la RB+anti

#con i boxplot rendo ancora più evidente la variabilità dei vari metodi
boxplot(list(a[, "diretta"], a[, "anti"], a[, "rb"], a[,"rb.anti"]))



#dato che l'integrale è univariato, non è necessario utilizzare simulazioni MC,
#ecco la funzione nativa di R che esegue il calcolo integrale 
integrate(function(x) exp(4 - 0.1 * x) * dnorm(x), -Inf, Inf)


#questo integrale si può fare anche 'a mano' in modo analitico, possiamo farlo 
#grazie al fatto che la densità di una variabile integra a 1 tra -inf e +inf:
#l'integrale si riduce a 
exp(4 + 0.1^2/2)














 



































