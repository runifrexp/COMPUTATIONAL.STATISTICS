#BOOTSTRAP




#GENERALITÀ



# - 1° caso: conosciamo la distribuzione del campione: la distribuzione, e tutto della stima 
#            della stima della media campionaria: la sua distribuzione e la 
#            distribuzione della sua quantità pivotale
#l'obiettivo e' la stima di mu, lo stimatore com'e' noto e' la media 
#campionaria, della quale conosciamo la distribuzione campionaria
n = 10
mu = 0
y = rnorm(n, mu, sqrt(sigma2))
y.bar = mean(y)
y.bar

s2 = var(y) 
s2



# - 2° caso: se NON consocessimo la distribuzione della media campionaria o 
#            della Q (quantità pivotale)
#genero via MC da F, dato che è nota, stimo le medie campionarie e in seguito 
#le quantità pivotali Q*. Per trovare la distribuzione di Q possiamo sostituire 
#il p-quantile della vera distribuzione di Q con Q*
B = 10000
stim.mc = replicate(B, mean(rnorm(n, mu, sqrt(sigma2)))) #stima della media 
#verifica casualità hist, qqplot
hist(stim.mc, freq = FALSE, main = "", nclass = 50)
curve(dnorm(x, mu, sqrt(sigma2/n)), add = TRUE, col = 'green')
qqplot(stim.mc, rnorm(B, mu, sqrt(sigma2/n)))
abline(0, 1, col = 'red')

#stima delle q. pivotali,genera un campione di 10 osservazioni per 10000 volte 
#e ne calcola le Q*
Q.mc = sapply(rep(n, B), FUN = function(n) { 
  x = rnorm(n, mu, sqrt(sigma2)) 
  (mean(x) - mu)/sqrt(var(x)/n)
})
#verifica casualità
hist(Q.mc, freq = FALSE, main = "", nclass = 50) 
curve(dt(x, df = n - 1), add = TRUE, col = 'green') 
qqplot(Q.mc, rt(B, df = n - 1)) 
abline(0, 1, col = 'red')



# - 3° caso: non conosco nulla del campione, ne' distribuzione ne' parametri:
#            bootstrap non parametrico
#stimo F via bootstrap non parametrico generando un campione casuale per avere 
#F.cap, ossia la stima di F per poi proseguire via MC come fatto nel primo caso
#per stimare F.cap si faranno delle estrazioni con reinserimento
#trunc restituisce un numero intero (arrotondamento)
ybar.bts = replicate(B, mean(y[trunc(runif(n, 1, n + 1))])) #come sapply
hist(ybar.bts, freq = FALSE, main = "", nclass = 50)
curve(dnorm(x, mu, sqrt(sigma2/n)), add = TRUE, col = 'green')
qqplot(ybar.bts, rnorm(B, mu, sqrt(sigma2/n)))
abline(0, 1, col = 'red')
#si nota facilmente che l'approssimazione non e' buona, questo perche' i
#campioni hanno numerosita' troppo basse (n=10)


#per studiare il comportamento in media (correttezza) di ybar.bts usiamo come 
#vero valore della media la media campionaria y.bar
mean(ybar.bts)
y.bar
mean(ybar.bts) - y.bar
#si nota che lo stimatore è corretto, com'e' la media campionaria per la 
#(vera) media di una distribuzione, questo perche' il comportamento in media 
#di una distribuzione viene traslato al mondo bootstrap 

#Analogamente a quanto succedeva per la media campionaria, il comportamento 
#relativo in media di sigma2*.hat rispetto a sigma2.hat (osservato) coincide con
#il comportamento relativo in media di sigma2.hat rispetto a sigma2 (vero valore
#di sigma2)
shat.bts = sapply(rep(n, B), FUN = function(n, y) { 
  x = y[trunc(runif(n, 1, n + 1))]
  (mean(x^2) - mean(x)^2) #stima varianza
}, y = y) 
mean(shat.bts) #media delle stime della varianza via bootstrap

s.hat = mean(y^2) - mean(y)^2
s.hat

mean(shat.bts) - s.hat

#stima della distribuzione di Q tramite Q*
Q.bts = sapply(rep(n, B), FUN = function(n, y) { x = y[trunc(runif(n, 1, n + 1))]
(mean(x) - mean(y))/sqrt(var(x)/n)
}, y = y)
#verifica di casualità
hist(Q.bts, freq = FALSE, main = "", nclass = 50) 
curve(dt(x, df = n - 1), add = TRUE, col = 'green') 
qqplot(Q.bts, rt(B, df = n - 1)) 
abline(0, 1, col = 'red')
#e' facile notare che Fhat stima male F in quanto n e' piccolo (n=10)



# - 4° caso: assumo che il campione si distribuisca come una normale di 
#            parametri ignoti: bootstrap parametrico
#per il resto procedo come nel MC (2° caso)

#stima cella media
ybar.parbts = replicate(B, mean(rnorm(n, y.bar, sqrt(s2))))
#verifica casualità
hist(ybar.parbts, freq = FALSE, main = "")
curve(dnorm(x, mu, sqrt(sigma2/n)), add = TRUE)
qqplot(ybar.parbts, rnorm(B, mu, sqrt(sigma2/n)))
abline(0, 1)





#BOOTSTRAP E SALMONI


d <- read.table(file.choose(), header = TRUE)
d[1:10,]

plot(d[, "adulti"], d[, "giovani"])
abline(a = 0, b = 1, lty = "dashed")
#si nota che la relazion tra givoani ed adulti non e' proporzionale 
#sa popolazione sembra quindi avere una qualche forma di “autoregolazione”


#verifica di corrlelazione, in quanto si tratta dell'evoluzione di un 
#comportamento nel tempo 
par(mar = c(2, 2, 0.2, 0.2), mfrow = c(1, 3))
acf(d[, "adulti"], ylim = c(-1, 1))
acf(d[, "giovani"], ylim = c(-1, 1))
#cross correlation: legame tra giovsani e adulti nello stesso anno
ccf(d[, "giovani"], d[, "adulti"], ylim = c(-1, 1))


#modello di Beverton-Holt per studiare il comportamento del numero di giovani 
#al variare del numero di adulti
plot(1/d[, "adulti"], 1/d[, "giovani"])
#ora la relazione sembra lienare


#costruico il modello lineare semplice di questa popolazione 
m <- lm(I(1/d[, "giovani"]) ~ I(1/d[, "adulti"])) 
summary(m)


plot(residuals(m))
#non si nota struttura

plot(fitted(m), residuals(m))
#ora sembra che la varianza delle prime unità è molto superiore a quella delle 
#unita' alla fine, per verificare se effettivamente è cosi devo fare un test F


beta <- coef(m)
plot(d[,"adulti"],d[,"giovani"])
curve(1/(beta[1]+beta[2]/x),
      min(d[,"adulti"]), 
      max(d[,"adulti"]), add=TRUE)
abline(0,1, lty = 2)
#l'intersezione tra la retta tratteggiata e la curva è E cappello


equilibrio <- function(d) {
  beta <- coef(lm(I(1/d[, "giovani"]) ~ 
                    I(1/d[, "adulti"]))) 
  E <- (1 - beta[2])/beta[1]
  names(E) <- "E"
  E
}


equilibrio(d)
#la stima della popolazione di equilibrio è di 150'000 salmoni





#tre modi per stimare E:
# - non parametrico
# - semi parametrico
# - parametrico



#non parametrico
rEconilcappello <- function(rsalmoni, B=1000) { 
  sapply(rep(1,B), 
         function(non.mi.interessa) equilibrio(rsalmoni()))
}

#funzione che genra senza assumere ne' normalita' ne' linearita'
n <- NROW(d)
rsalmoni.np <- function() {
  index <- floor(1 + n * runif(n))
  d[index, ]
}

np <- rEconilcappello(rsalmoni.np)
hist(np - equilibrio(d), nclass = 50)

mean(np) - equilibrio(d) #distorsione: stima - vero valore del parametro
#e' piccola. quindi E e' corretto.
sd(np) #valore basso considerando che siamo sull'ordine di grandezza delle 
#centinaia




#stima semi-parametrica

#assumo normalita': approccio semi-parametrico
n <- NROW(d)
res <- residuals(m)
beta <- coef(m)
rsalmoni.sp <- function() {
  index <- floor( 1 + n * runif(n))
  d.new <- d
  d.new[,"giovani"] <- 1/(beta[1]+beta[2]/d[,"adulti"]+res[index] )
  d.new
}

sp <- rEconilcappello(rsalmoni.sp)
hist(sp - equilibrio(d), nclass = 50)

mean(sp) - equilibrio(d) #corretto
sd(sp) #la deviazione standard e' ancora bassa


#stima parametrica

n <- NROW(d)
res <- residuals(m)
beta <- coef(m)
sigma <- sqrt(sum(res^2)/(length(res)-length(beta))) 
rsalmoni.par <- function() {
  d.new <- d
  d.new[,"giovani"] <- 1/(beta[1]+beta[2]/d[,"adulti"]+sigma*rnorm(n))
  d.new
}

par <- rEconilcappello(rsalmoni.par)
hist(par - equilibrio(d), nclass = 50)

mean(par) - equilibrio(d) #corretto
sd(par) #la deviazione standard e' ancora bassa



#ic per i tre procedimenti diversi: alpha = 0.05: 
#errore di stima (theta.cap-theta)
equilibrio(d) - quantile(np - equilibrio(d), c(0.975, 0.025))
equilibrio(d) - quantile(sp - equilibrio(d), c(0.975, 0.025))
equilibrio(d) - quantile(par - equilibrio(d), c(0.975, 0.025))

#ic per i tre procedimenti diversi: alpha = 0.05: 
#stimatore
2 * equilibrio(d) - quantile(np, c(0.975, 0.025))
2 * equilibrio(d) - quantile(sp, c(0.975, 0.025))
2 * equilibrio(d) - quantile(par, c(0.975, 0.025))




#funzione che stima anche la deviazione standard della stima; per farlo si usa 
#il bootstrap (di secondo livello), con B2 = 20 replicazioni di default:
#bootstrap sopra il bootstrap
equilibrio.s <- function(d, rsalmoni, B2=20) {
  beta <- coef(lm(I(1/d[, "giovani"]) ~ I(1/d[,"adulti"]))) 
  E <- (1 - beta[2])/beta[1]
  d.orig <- globalenv()$d
  d <<- d
  sd.E <- sd(rEconilcappello(rsalmoni, B=B2))
  d <<- d.orig
  names(E) <- "E"
  names(sd.E) <- "sd.E"
  return(c(E,sd.E))
}


#In equilibrio.s bisogna fare attenzione al fatto che per definizione 
#rsalmoni.np, rsalmoni.sp e rsalmoni.par usano il valore di d che si trova 
#nell’ambiente globale, e non quello che viene passato come argomento.
rEconilcappello.s <- function(rsalmoni, B=1000) { 
  replicate(B, equilibrio.s(rsalmoni(), rsalmoni))
}



#approccio non parametrico:
sd.E.np <- equilibrio.s(d, rsalmoni.np)[2]
np.s <- rEconilcappello.s(rsalmoni.np)

#approccio semiparametrico:
sd.E.sp <- equilibrio.s(d, rsalmoni.sp)[2]
sp.s <- rEconilcappello.s(rsalmoni.sp)

#approccio parametrico:
sd.E.par <- equilibrio.s(d, rsalmoni.par)[2]
par.s <- rEconilcappello.s(rsalmoni.par)



#applicazione alla popolazione di equilibrio

#approccio non parametrico:
t.np <- (np.s[1,] - equilibrio(d))/np.s[2,]
equilibrio(d) - sd.E.np * quantile(t.np, c(0.975, 0.025))

#approccio semiparametrico:
t.sp <- (sp.s[1,] - equilibrio(d))/sp.s[2,]
equilibrio(d) - sd.E.sp * quantile(t.sp, c(0.975, 0.025))

#approccio parametrico:
t.par <- (par.s[1,] - equilibrio(d))/par.s[2,]
equilibrio(d) - sd.E.par * quantile(t.par, c(0.975, 0.025))





#STIMA DELLA DISTRIBUZIONE DEL COEFFICIENTE DI CORRELAZIONE DI CORRELAZIONE


#genera una coppia di valori X e Y da N(0,1)
rbivnorm <- function(n, rho) { 
  x <- rnorm(n)
  cbind(x, rho * x + sqrt(1 - rho^2) * rnorm(n)) 
}

dati <- rbivnorm(1000, -0.8)
dim(dati)

cor(dati[, 1], dati[, 2])

apply(dati, 2, mean)

apply(dati, 2, sd)

plot(dati[, 1], dati[, 2])


#campione di 50 coppie (X,Y)
n <- 50


#ritorna una stima del coefficiente per B volte per un campione di 50 coppie
rrho <- function(B, n, rho) { 
  un.rho <- function(n) {
    xy <- rbivnorm(n, rho)
    cor(xy[, 1], xy[, 2]) 
  }
  sapply(rep(n, B), un.rho) 
}

#rho <- 0
vero.rho <- 0
B <- 10000
r <- rrho(B, n, vero.rho)
mean(r)
#rho stimato e' quasi zero: bene

sd(r)
#deviazione standard molto bassa dato che l'ordine di grandezza della stima è 
#delle centinaia

hist(r - vero.rho, nclass = 40)
#bene anche dalla distribuzione dell'erorre di stima: centrata in zero, 
#simmetrica e con bassa variabilita'



#rho <- 0.95
vero.rho <- 0.95
r <- rrho(B, n, vero.rho)
mean(r)
#bene: la stima e' molto vicina al vero valore di rho

sd(r)
#molto bene anche la deviazione standard

hist(r - vero.rho, nclass = 40)
#nell'istogramma succede qualosa di inaspetato: non c'e' simmetria.
#questo succede perche' la distribuzione e' tra [-1,1]




#stima bootstrap dell'errore di stima 

#genero gli r*
cor.boot <- function(dati, B = 10000) { 
  n <- nrow(dati)
  r <- rep(0, B) 
  for (i in 1:B) {
    idx <- floor(1 + n * runif(n)). #genero una sola volta l'indice poiche' devo 
    #prendere le coppie di valori in modo tale da non alterare la correlazione 
    #che c'e' nel mondo reale
    r[i] <- cor(dati[idx, 1], dati[idx, 2]) 
  }
  r
}


#MC sopra bootsrap: (MC con B=6 non è un vero e proprio MC pero' concettualmente
#funziona cosi')

#rho <- 0
nvolte <- 6
par(mfrow=c(2,3))
vero.rho <- 0
for ( i in 1:nvolte ) {
  dati <- rbivnorm( n , vero.rho ) 
  r <- cor.boot( dati )
  hist( r - cor(dati[,1],dati[,2]) , nclass=40)
  a <- round(c(cor(dati[,1],dati[,2]),
               quantile(r,c(0.025,0.975))),3) 
  cat("stima: ", a[1],";\t i. c. = [",2*a[1]-a[3]," ,",
      2*a[1]-a[2],"]\n",sep="") 
}


#rho <-  0.95
nvolte <- 6
par(mfrow=c(2,3))
vero.rho <- 0.95
for ( i in 1:nvolte ) {
  dati <- rbivnorm( n , vero.rho ) 
  r <- cor.boot( dati )
  hist( r - cor(dati[,1],dati[,2]) , nclass=40)
  a <- round(c(cor(dati[,1],dati[,2]),
               quantile(r,c(0.025,0.975))),3) 
  cat("stima: ", a[1],";\t i. c. = [",2*a[1]-a[3]," ,",
      2*a[1]-a[2],"]\n",sep="") 
}

















