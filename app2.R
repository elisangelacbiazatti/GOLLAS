# case-fatality ratio from 20 countries most affected by COVID-19 worldwide
dados<-c(0.076,0.056,0.04,0.029,0.028,0.028,0.025,0.025,0.024,0.024,0.021,0.021,
         0.016,0.015,0.015,0.013,0.012,0.01,0.009,0.009)

library(AdequacyModel)
library(MASS)
require(GenSA)
require(survival)

n=length(dados)
n

TTT(dados, col="red", lwd=2, grid=TRUE, lty = 2)

boxplot(dados,ylab="")
summary(dados)

#######

fit.sa2<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2])))   
  lower <- c(0,0) 
  upper <- c(10,10)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,maxit=7000,max.time=3))
  return(out[c("value","par","counts")])
}


pdf.GOLLAS<-function(x,alpha,betta){
  g<-1/(pi*sqrt(x-x^2))
  G<-2/pi*asin(sqrt(x))
  fd<-(alpha*betta*g*G^(alpha*betta-1)*(1-G^betta)^(alpha-1))/
    (G^(alpha*betta)+(1-G^betta)^(alpha))^2
  fd
}
set.seed(1729)
fit.sa2(dados,pdf.GOLLAS)

GOLLAS_pdf<-function(par,x){
  alpha<-par[1]
  betta<-par[2]
  g<-1/(pi*sqrt(x-x^2))
  G<-2/pi*asin(sqrt(x))
  (alpha*betta*g*G^(alpha*betta-1)*(1-G^betta)^(alpha-1))/(G^(alpha*betta)+(1-G^betta)^(alpha))^2
  
}
GOLLAS_cdf<-function(par,x){
  alpha<-par[1]
  betta<-par[2]
  
  G<-2/pi*asin(sqrt(x))
  (G^(alpha*betta))/(G^(alpha*betta)+(1-G^betta)^(alpha))
  
}

set.seed(1729)
resultsGOLLAS  = goodness.fit(pdf = GOLLAS_pdf, cdf = GOLLAS_cdf,
                              starts = c(15, 2.045752),
                              data = dados, method = "SANN", domain = c(0, 1),
                              mle = NULL);resultsGOLLAS


###########
beta.pdf=function(x,a,b){
  
  dbeta(x,a, b)
}

set.seed(1729)
fit.sa2(dados,beta.pdf)

cdf_beta=function(par,x){
  a=par[1]
  b=par[2]
  
  pbeta(x,a, b)
}
pdf_beta=function(par,x){
  a=par[1]
  b=par[2]
  
  dbeta(x,a, b)
}
set.seed(1729)
resultsbeta = goodness.fit(pdf = pdf_beta, cdf = cdf_beta,
                           starts = c(0.5973426, 10), data = dados, 
                           method = "SANN", 
                           domain = c(0, 1),
                           mle = NULL);resultsbeta
chutes<-c(0.5973426, 10)
fit.gtnhb<- fitdistr(dados,beta.pdf,start=list(a=chutes[1],b=chutes[2]),control=list(ndeps=c(1e-8,1e-12),maxit=10000))
fit.gtnhb 

#########
#Kw
#########
Kw.pdf<-function(x,a,b){
  fd<-a*b*x^(a-1)*(1-x^a)^(b-1)
}
set.seed(1729)
fit.sa2(dados,Kw.pdf)

pdf_kw=function(par,x){
  a=par[1]
  b=par[2]
  
  a*b*x^(a-1)*(1-x^a)^(b-1)
}

cdf_kw=function(par,x){
  a=par[1]
  b=par[2]
  
  (1-(1-x^a)^b)
}
set.seed(1729)
resultskw = goodness.fit(pdf = pdf_kw, cdf = cdf_kw,
                         starts = c(0.7119176, 10), data = dados, 
                         method = "SANN", domain = c(0, 1),
                         mle = NULL);resultskw
chutes<-c(0.7119176, 10)
fit.gtnhkw<- fitdistr(dados,Kw.pdf,start=list(a=chutes[1],b=chutes[2]),control=list(ndeps=c(1e-6,1e-12),maxit=10000))
fit.gtnhkw

###########
#EAS
#############
EAS.pdf<-function(x,a,b){
  g<-1/(pi*sqrt(x-x^2))
  G<-2/pi*asin(sqrt(x))
  fd<-a*b*(1-G)^(a-1)*(1-(1-G)^a)^(b-1)*g
  fd
}
set.seed(1729)
fit.sa2(dados,EAS.pdf)

EAS_pdf<-function(par,x){
  a<-par[1]
  b<-par[2]
  g<-1/(pi*sqrt(x-x^2))
  G<-2/pi*asin(sqrt(x))
  fd<-a*b*(1-G)^(a-1)*(1-(1-G)^a)^(b-1)*g
  fd
} 
EAS_cdf<-function(par,x){
  a<-par[1]
  b<-par[2]
  #g<-1/(pi*sqrt(x-x^2))
  G<-2/pi*asin(sqrt(x))
  fa<-(1-(1-G)^a)^(b)
  fa
}
set.seed(1729)
resultseas = goodness.fit(pdf = EAS_pdf, cdf = EAS_cdf,
                          starts = c(12, 2.5), data = dados, 
                          method = "SANN", domain = c(0, 1),
                          mle = NULL);resultseas

chutes<-c(12, 2.5)
fit.gtnheas<- fitdistr(dados,EAS.pdf,start=list(a=chutes[1],b=chutes[2]),control=list(ndeps=c(1e-10,1e-10),maxit=10000))
fit.gtnheas
fit.gtnheas$estimate[2]

#graficos
truehist(dados,ylim=c(0,45),nbins = 5,col = "white",ylab="f(x)",xlab = "x")
curve(GOLLAS_pdf(resultsGOLLAS$mle,x),add=TRUE, lwd = 3,lty=1, col="red")
curve(EAS_pdf(fit.gtnheas$estimate,x),add=TRUE, lwd = 3, lty=4, col="yellow")
curve(pdf_beta(fit.gtnhb$estimate,x),add=TRUE, lwd = 3, lty=2, col="blue")
curve(pdf_kw(fit.gtnhkw$estimate,x),add=TRUE, lwd = 3, lty=3, col="green")

legend(0.04,35, legend = c( "GOLLAS","EAS", "Beta","Kw"),
       col = c("red","yellow","blue","green"), lwd = 3,lty=c(1,4,2,3), bty ="n")

#############
#Acumulada
#############
km<- survfit(Surv(dados) ~ 1) #Kaplan-Meier

plot(km$time, 1-km$surv, xlab = "x", ylab="F(x)",lwd=2,type = "l") #curva estimada

curve(GOLLAS_cdf(resultsGOLLAS$mle,x),add=TRUE, lwd = 3,lty=1, col="red")
curve(EAS_cdf(fit.gtnheas$estimate,x),add=TRUE, lwd = 3, lty=4,col="yellow")
curve(cdf_beta(fit.gtnhb$estimate,x),add=TRUE, lwd = 3, lty=2,col="blue")
curve(cdf_kw(fit.gtnhkw$estimate,x),add=TRUE, lwd = 3, lty=3,col="green")

legend(0.04,.6, legend = c( "GOLLAS","EAS", "Beta","Kw"),
       col = c("red","yellow","blue","green"), lwd = 3,lty=c(1,4,2,3), bty ="n")

