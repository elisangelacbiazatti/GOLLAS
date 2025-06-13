# Proportion of severe dengue. Countries and territories of The Americas, 2022
dados<-c(0.0741,0.044,0.0259,0.0245,0.0197,0.0172,
         0.0156,0.0068,0.0068,0.0059, 0.0038,  0.0038, 
         0.0037,  0.0035,  0.0031, 0.0020,0.0008,0.0006,0.0004)

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
                              starts = c(5.1416453, 0.2340499), 
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
                           starts = c(0.3668479, 10), data = dados, 
                           method = "SANN", 
                           domain = c(0, 1),
                           mle = NULL);resultsbeta

chutes<-c(0.3668479, 10)
fit.gtnhb<- fitdistr(dados,beta.pdf,start=list(a=chutes[1],b=chutes[2]),control=list(ndeps=c(1e-6,1e-12),maxit=10000))
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
                         starts = c(0.7737536, 2), data = dados, 
                         method = "SANN", domain = c(0, 1),
                         mle = NULL);resultskw
chutes<-c(0.7737536, 2)
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
  G<-2/pi*asin(sqrt(x))
  fa<-(1-(1-G)^a)^(b)
  fa
}
set.seed(1729)
resultseas = goodness.fit(pdf = EAS_pdf, cdf = EAS_cdf,
                          starts = c(100, 100), data = dados, 
                          method = "SANN", domain = c(0, 1),
                          mle = NULL);resultseas

chutes<-c(100, 100)
fit.gtnheas<- fitdistr(dados,EAS.pdf,start=list(a=chutes[1],b=chutes[2]),control=list(ndeps=c(1e-8,1e-6),maxit=10000))
fit.gtnheas


#graficos
truehist(dados,ylim=c(0,100),nbins = 15,col = "white",ylab="f(x)",xlab = "x")
curve(GOLLAS_pdf(resultsGOLLAS$mle,x),add=TRUE, lwd = 3,lty=1, col="red")
curve(EAS_pdf(fit.gtnheas$estimate,x),add=TRUE, lwd = 3, lty=4, col="yellow")
curve(pdf_beta(fit.gtnhb$estimate,x),add=TRUE, lwd = 3, lty=2, col="blue")
curve(pdf_kw(fit.gtnhkw$estimate,x),add=TRUE, lwd = 3, lty=3, col="green")

legend(0.04,90, legend = c( "GOLLAS","EAS", "Beta","Kw"),
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

legend(0.04,.6, legend = c( "GOLLAS", "EAS","Beta","Kw"),
       col = c("red","yellow","blue","green"), lwd = 3,lty=c(1,4,2,3), bty ="n")

