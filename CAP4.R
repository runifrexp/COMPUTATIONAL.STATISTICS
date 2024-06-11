#BOOTSTRAP - APPLICAZIONI 



#UN TEST BASATO SUL BOOTSTRAP

d <- read.table(file.choose(), header = TRUE)
head(d)

x <- d[d[, 1] == 1, 2]
y <- d[d[, 1] == 2, 2]


boxplot(x,y, horizontal = T)
#si nota che le medie sono diverse: si tratta di una casualita' oppure e' un 
#risultato costante

#test t-stud con assunzione di omoschedasticita'
t.test(x, y, var.equal = T)

#test di Welch: test t-stud senza assunzione di omoschedasticita'
t.test(x,y)


#questi due test si basano sull'assunzione di normalita': la verifico
shapiro.test(x)
#accetto normalita' (H_0)
shapiro.test(y)
#sono incerto: pvalue basso, non so se accettare


#quindi: quanto dipende la conclusione del test da quest'assunzione?
nx <- length(x)
ny <- length(y)
#statistica test per il test di Shapiro_Wilk
t.dati <- (mean(x) - mean(y))/sqrt(var(x)/nx + var(y)/ny)
t.dati


rtduecampioni <- function( x , y , B=10000) { 
  nx <- length(x)
  ny <- length(y)
  T <- rep( 0 , B )
  for ( i in 1:B ) {
    xb <- x[floor(1+nx*runif(nx))] 
    yb <- y[floor(1+ny*runif(ny))] 
    T[i] <- (mean(xb)-mean(yb))/sqrt(var(xb)/nx+var(yb)/ny)
  }
  T
}

#stima lso sbagliata:
t.boot <- rtduecampioni(x, y)
mean(abs(t.boot) >= abs(t.dati))
#dovrebbe essere centrata in zero la distribuzione delle statistica test in 
# in quanto le statistica test sono calcolate sotto H_0, quindi dovrebbero 
#essere 0 in quanto sotto H_0 le medie sono uguali quindi si annullano


hist(t.boot, nclass = 50)
points(t.dati, 0, pch = "*", cex = 5) #cex rende l'asterisco piu' grande
#lo si nota anche dall'istogramma delle statistica test genereate 


#verifichiamo se effettivamente H_0 e' vera o meno:
mean(x) - mean(y)
#no, non e' vera


#dato che non siamo sotto H_0 creo un mondo bootstrap nel quale è valida H_0:
t.boot <- rtduecampioni(x - mean(x), y - mean(y))
mean(abs(t.boot) >= abs(t.dati))
#si, ora siamo (circa) sotto H_0



hist(t.boot, prob = TRUE, nclass = 50)
curve(dt(x, 21), -5, 5, lty = "dashed", lwd = 2, add = TRUE)
#ora notiamo che la distribuzione e' simmetrica e centrata in 0, e l'istogramma
#assomiglia alla forma della densita' di una tstudent con 21 gradi di liberta'
#come richiesto nel test di Welch (diverso da t-student in quanto nel test di 
#Welch non viene assunta omoschedasticita', e con GdL diversi (21 vs 22))
points(t.dati, 0, pch = "*", cex = 5)# questo asterisco e' t^oss





#BOOTSTRAP LISCIATO


#formulazione del problema:
#la durata del 90% delle componenti e', con alta probabilita', superiore a a(y)

y.fake <-read.table(file.choose(), header = T) #durate.dat

head(y.fake)
y <- y.fake$Durate

library(sm)

#due metodi equivalenti per calcolare l'h migliore per il trade off tra 
#efficienza e correttezza
h.select(y, method = "cv")
h.select(y, method = "sj")


h <- h.select(y, method = "sj")

#plot della densita' stiamta 
sm.density(y, h = h)


pnucleo <- function(x, y, h) { 
  a <- outer(x, y, "-")/h 
  apply(pnorm(a), 1, mean)
}


plot(ecdf(y), do.points = FALSE)

curve(pnucleo(x, y, h), xlim = range(y) +c (-10, 10),
      lty = "dashed", lwd = 2, add = TRUE)



rnucleo <- function(n, y, h) { 
  ny <- length(y)
  idx <- floor(1 + ny * runif(n))
  y[idx] + h * rnorm(n) 
}

z <- rnucleo(1e+05, y, h)
sm.density(y, h = h)
hist(z, nclass = 100, prob = TRUE, add = TRUE)

##############scelta ottimale degli intervalli (nclass) per hist################
library(KernSmooth)
bin_width <- dpih(z$V1)
nbins <- seq(min(z$V1) - bin_width,
             max(z$V1) + bin_width,
             by = bin_width)
hist(z$V1, breaks = nbins, add = T, prob = T)
#altro metodo
hist(z$V1, breaks = "FD", add = T, prob = T) 
################################################################################


#stima del quantile 10% (decimo percentile)
q10.stima <- quantile(y, 0.1)
q10.stima
#quindi circa 23 ore 



#ora usiamo la funzione di densita' stimata per generare 10'000 valori di q10
#e per poterne tracciare quindi la distribuzione 

B <- 10000
n <- length(y)
yb <- matrix(rnucleo(n * B, y, h), n, B)
#ecco la distribuzione del decimo percentile via bootstrap
q10.boot <- apply(yb, 2, function(x) quantile(x, 0.1))


#una rappresentazione degli errori che si possono fare e'
hist(q10.boot - quantile(y, 0.1), nclass = 100, prob = TRUE)
#Gli errori di stima non sono perfettamente centrati intorno allo zero: lo 
#stimatore ha una seppur piccola distorsione, eccone la stima
dist <- mean(q10.boot) - q10.stima
dist


#Una stima del quantile di interesse "corretta" via bootstrap e' quindi
quantile(y, 0.1) - dist

#ecco il risultato finale: quello che cercavamo era stimare via bootstrap 
#la durata del 90% delle componenti e', con alta probabilita', superiore a a(y)

a <- 2 * quantile(y, 0.1) - quantile(q10.boot, 0.95)
cat("la durata del 90% delle componenti e', con alta probabilita', 
    superiore a", a)












