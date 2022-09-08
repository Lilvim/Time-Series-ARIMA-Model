library(readxl)

install.packages("urca")
library(urca)
library(lmtest)
library(FitAR)

library(astsa)

library(stargazer)

library(lmtest)

library(readr)

library(zoo)

library(forecast)
install.packages('funtimes')
library(funtimes)
library(ellipse)
library(ggplot2)
library(dbplyr)
library(xts)
install.packages('aTSA')
library(aTSA)
install.packages('fBasics')
library(fBasics)
install.packages('FitAR')
library(FitAR)
install.packages('portes')
library(portes)


#Manipulation sur la dataframe
data = read.csv('Bureau/data_hydro.csv', sep = ';')
data$Mois = match(data$Mois, c('Janvier','Février','Mars','Avril','Mai','Juin','Juillet','Août','Septembre','Octobre','Novembre','Décembre'))
data$Jour = '01'
data$date <- as.Date(with(data, paste(Année, Mois, Jour,sep="-")), "%Y-%m-%d")

data$date = as.Date(data$date, tryFormats = c("%Y-%m-%d"))
data$Valeur =as.numeric(paste(data$Valeur))

#Régression pour vérifier la présence d'une trend linéaire
lm = lm(Valeur~date, data = data)
summary(lm)


#On plot les séries
p <- ggplot(data, aes(x=date, y=Valeur)) +
  geom_line() + 
  xlab("Années") + ylab('Indice de production industrielle de savon')

p + theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid",
                                   colour = "black"))

data$logV = log(data$Valeur)

p_log <- ggplot(data, aes(x=date, y=logV)) +
  geom_line() + 
  xlab("Années") + ylab('Logarithme de l\'ndice de production industrielle de savon')
p_log + theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + theme(panel.border = element_blank(),
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),axis.line = element_line(size = 0.5, linetype = "solid",
                                                                               colour = "black"))
#On transforme en série temporelle
Xt.ts<-ts(rev(data$Valeur),start= c(1990,1), frequency=12)
print(Xt.ts)

Gt.ts=log(Xt.ts)


# On trace l'acf de la serie
acf(Xt.ts, main="",60) 

dev.off()

#On decompose Gt 
plot(decompose(Gt.ts)) 
dev.off()

#Différenciciton de la log série et on trace les ACF PACF
Zt.ts=diff(Gt.ts, 1)
plot(Zt.ts)
plot(decompose(Zt.ts))
acf(Zt.ts,20,main="")
pacf(Zt.ts,50)


#On trace la moyenne mobile 

MA<-filter(Zt.ts, filter=array(1/10,dim=10), method = c("convolution"),
           sides = 2, circular = FALSE)
plot(Zt.ts,type='l')
lines(MA,col='red')

#On vérifie que la série n'a plus de trend linéaire
lm <- lm(Zt.ts ~ c(1:length(Zt.ts)))
summary(lm)

#Test de stationnarité
summary(ur.kpss(Zt.ts,type = c('tau'),lags = c('short','long','nil')))
adf.test(Zt.ts, nlag = NULL, output = TRUE)
pp.test(Zt.ts)     



#On calcule le BIC et l'AIC des modèles retenus
pmax = 3
qmax = 2
aic =matrix(0,nrow = pmax+1, ncol = qmax+1) 
bic = matrix(0,nrow = pmax+1, ncol = qmax+1)


for (p in seq(0,pmax)) {
  for (q in seq(0,qmax)) {
    modele=try(arima(Zt.ts,order = c(p,0,q)))
      if (class(modele)=="try-error"){
        next
        }
      else {
        aic[p+1,q+1] =AIC(modele)
        bic[p+1,q+1] =BIC(modele)
        } 
                         }
} 
(aic == min(aic))
(bic == min(bic))


stargazer(as.data.frame((aic == min(aic))),summary = FALSE)
stargazer(as.data.frame((bic == min(bic))),summary = FALSE)
#Juste pour le Latex

#Dans le rapport on choisit le ARMA(0,2)
AR12 <- arima(Zt.ts, order = c(1,0,2),method='ML')
portest(AR12,NREP=100, test = 'LjungBox')


#Le test de Portmanteau
AR02 <- arima(Zt.ts, order = c(0,0,2),method='ML')
portest(AR02,NREP=100, test = 'LjungBox')



#On calcule les pvaleurs des coefficients pour vérifier qu'ils sont bien significatifs
((1-pnorm(abs(AR02$coef)/sqrt(diag(AR02$var.coef))))*2)


#Les résidus sont ils normaux ?
tsdiag(AR02,100)
qqnorm(AR02$residuals)
hist(AR02$residuals)
normality = shapiro.test(AR02$residuals)
normality



#On vérifie les fits des ARIMA


AR112 = arima(Gt.ts, order = c(1,1,2),method='ML')
portest(AR112,NREP=100, test = 'LjungBox')
plot((Gt.ts))
lines(fitted(AR112),col="blue")

AR012 = arima(Gt.ts, order = c(0,1,2),method='ML', include.mean=F)
portest(AR012,NREP=100, test = 'LjungBox')
plot((Gt.ts))
lines(fitted(AR012),col="red")

#prediction region de confiance
sarima.for(tail(Zt.ts,n=20), n.ahead = 3, p = 0,d = 0,q = 2)



