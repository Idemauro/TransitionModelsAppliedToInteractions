#Program for L4

rm(list=ls(all=TRUE))

#######################################################################################################
# the data-set
setwd("~/Documents/Pesquisa/Livro_Entomo/Programa_discreto")
dados <- read.csv("dados_lagartas_L4.csv", head=TRUE, sep=";", dec=",")
dados=na.omit(dados)
head(dados)
dim(dados)

dados <- dados[,c('instar','planta','resp1','resp2','resp3','resp4','resp5','resp6')]
dados$id <- seq(1,nrow(dados))

# Conversion for factors

dados$instar <- as.factor(dados$instar)
dados$planta <- as.factor(dados$planta)
dados$resp1 <- as.factor(dados$resp1)
dados$resp2 <- as.factor(dados$resp2)
dados$resp3 <- as.factor(dados$resp3)
dados$resp4 <- as.factor(dados$resp4)
dados$resp5 <- as.factor(dados$resp5)
dados$resp6 <- as.factor(dados$resp6)

# Packages
require(MASS)
require(car)
require(VGAM)
require(Matrix)
require(nnet)

head(dados)
#######################################################################################################
# Part 2 - Functions to change structure of data

widetolong<-function(basewide,varswide){
  r1<-reshape(basewide,varying=list(varswide),direction='long',v.names='resp')
  r1$time<-factor(r1$time)
  rownames(r1)<-1:nrow(r1)
  r1<-r1[,-which(names(r1)=='id')]
  r2<-r1[-which(r1$time==levels(r1$time)[1]),]
  r2<-cbind(r2,r1[-which(r1$time==levels(r1$time)[length(levels(r1$time))]),'resp'])
  names(r2)[ncol(r2)]<-'respp';head(r2)
  return(r2)
}

longtowide<-function(baselonga,resplonga,vartime,respprevia,idsuj,resp1){
  indicat1<-names(which(table(baselonga$time)!=0))[1]
  r3<-cbind(reshape(baselonga[,-which(names(baselonga)==respprevia)],idvar = idsuj,v.names = resplonga,
                    timevar = vartime, direction = "wide",sep=""),
            baselonga[,respprevia][which(baselonga[,vartime]==indicat1)])
  names(r3)[ncol(r3)]<-resp1
  return(r3)
}

#######################################################################################################
dadosl <- widetolong(dados,c('resp1','resp2','resp3','resp4','resp5','resp6')) # data in long way.

dadosl$resp <- factor(dadosl$resp)
dadosl$respp <- factor(dadosl$respp)
dadosl$resp=relevel(dadosl$resp,ref="nothing")
dadosl$respp=relevel(dadosl$respp,ref="nothing")
head(dadosl)

# Model suposing stationarity
modE0<- multinom(resp ~ planta, family=multinomial,data=dadosl)
modE <- multinom(resp ~ planta+respp, family=multinomial,data=dadosl)
anova(modE,modE0)


beta0<-as.numeric(c(coef(modE)[1,],coef(modE)[2,])) # estimated vector-stationary model
vcov0<-vcov(modE) # covar matrix-stationary model 
loglik0<-as.numeric(logLik(modE)) # Log-likelihood-stationary model 


dadospred<-with(dadosl,expand.grid(planta=levels(planta),respp=levels(respp)))
# Grid for prediction

preds<-cbind(dadospred,predict(modE,newdata=dadospred, type="probs"))
predsmE<-preds;predsmE # Transtion matrix under stationarity

# Five Models supposing non stationarity

mod12<-multinom(resp2~planta+resp1,family='multinomial',data=dados)


beta12<-as.numeric(c(coef(mod12)[1,],coef(mod12)[2,])) # estimated vector
vcov12<-vcov(mod12) # Covar matrix
loglik12<-as.numeric(logLik(mod12)) # Log-likelihood

mod23<-multinom(resp3~planta+resp2,family='multinomial',data=dados)

beta23<-as.numeric(c(coef(mod23)[1,],coef(mod23)[2,]))
vcov23<-vcov(mod23)
loglik23<-as.numeric(logLik(mod23))


mod34<-multinom(resp4~planta+resp3,family='multinomial',data=dados)


beta34<-as.numeric(c(coef(mod34)[1,],coef(mod34)[2,]))
vcov34<-vcov(mod34)
loglik34<-as.numeric(logLik(mod34))

mod45<-multinom(resp5~planta+resp4,family='multinomial',data=dados)

beta45<-as.numeric(c(coef(mod45)[1,],coef(mod45)[2,]))
vcov45<-vcov(mod45)
loglik45<-as.numeric(logLik(mod45))


mod56<-multinom(resp6~planta+resp5,family='multinomial',data=dados)

beta56<-as.numeric(c(coef(mod56)[1,],coef(mod56)[2,]))
vcov56<-vcov(mod56)
loglik56<-as.numeric(logLik(mod56))

#######################################################################################################
#Likelihood ratio test
################################################################################################
logver <- (-2*(loglik0-sum(loglik12,loglik23,loglik34,loglik45,loglik56))) ### LRT's test statistic 
logver
pchisq(logver, 32, lower.tail = FALSE)

#######################################################################################################
## Wald test (article)
################################################################################################
beta <- c(beta12,beta23,beta34,beta45,beta56)
varcov <- as.matrix(bdiag(vcov12,vcov23,vcov34,vcov45,vcov56))
beta0vec <- rep(beta0,5)

wald <- t(beta-beta0vec)%*%solve(varcov)%*%(beta-beta0vec) ### Wald's test statistic problema na matriz
wald
pchisq(wald, 32, lower.tail = FALSE)


# W1: 08/09 liao(1996)

beta1 <- c(beta12-beta23,beta12-beta34,beta12-beta45,beta12-beta56)
varcov1 <- as.matrix(bdiag(vcov12+vcov23,vcov12+vcov34,vcov12+vcov45,vcov12+vcov56))
dim(varcov1)
dim(vcov12)

wald1 <- t(beta1)%*%solve(varcov1)%*%(beta1) ### Wald's test statistic
wald1
pchisq(wald1, 32, lower.tail = FALSE)


# W3: 12/09 liao(1996)

length(beta1)
dim(varcov1)
dim(vcov12)

#Formando a matrix de covari?ncias mista
C<-matrix(c(0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0),ncol=4,nrow=4,byrow=T)
C
A<-C%x%vcov12 #usando o produto de Kronecker
dim(A)
varcovmod<-varcov1+A 
wald3 <- t(beta1)%*%solve(varcovmod)%*%(beta1) ### Wald's test statistic
wald3
pchisq(wald3, 32, lower.tail = FALSE) 


#######################################################################################################
#Proposed test 3- using the stationary model, with visit term in a predictor
#######################################################################################################
head(dadosl)
levels(dadosl$time)
dadosl$time <- as.factor(as.character(dadosl$time))
summary(dadosl)


# Using multinom function of nnet

modE=multinom(resp~planta+respp, family=multinomial, data=dadosl)
summary(modE)
loglikE<-as.numeric(logLik(modE))


modE1=multinom(resp~planta+respp+time, family=multinomial,data=dadosl)
summary(modE1)
loglikE1<-as.numeric(logLik(modE1))
anova(modE,modE1)


modE2=multinom(resp~planta+respp*time, family=multinomial,data=dadosl)
summary(modE2)
loglikE2<-as.numeric(logLik(modE2))

modE3=multinom(resp~planta*time+respp, family=multinomial,data=dadosl)
summary(modE3)
loglikE3<-as.numeric(logLik(modE3))

modE4=multinom(resp~(planta+respp)*time, family=multinomial,data=dadosl)
summary(modE4)
loglikE4<-as.numeric(logLik(modE4))



anova(modE,modE1,modE2,modE4)
anova(modE,modE1,modE3,modE4)

modE0=multinom(resp~respp, family=multinomial, data=dadosl)
modE=multinom(resp~planta+respp, family=multinomial, data=dadosl)
modE1=multinom(resp~planta+respp+time, family=multinomial,data=dadosl)
anova(modE0,modE)
anova(modE,modE1)



#Selected model: previous response + time

summary(modE1)

#####################################################################################################

###############  PARA FAZER GR?FICO HEAT MAP - CESAR TACO ####################################

dadospred<-with(dadosl,expand.grid(planta=levels(planta),respp=levels(respp),time=levels(time)))
# Grid for prediction

preds<-cbind(dadospred,predict(modE1,newdata=dadospred, type="probs"))
predsmE<-preds;predsmE # Transtion matrix under stationarity

################################################################################################

#########################################################################################################


require(hnp)
hnp(modE1,how.many.out = T, paint.out = T, xlab="theorical quantiles",
    ylab="residuals", main="Fourth instar")
text(locator(1), "Points out of envelope: 1.88%")

#######################################################################################################

#Using VGLM function (vgam) only to check!!! Don't use

modEv=vglm(resp~planta+respp, family=multinomial, data=dadosl)
summary(modEv)
nparam(modEv)

modE1v=vglm(resp~planta+respp+time, family=multinomial, data=dadosl)
summary(modE1v)
nparam(modE1v)

modE2v=vglm(resp~planta+respp*time, family=multinomial,data=dadosl)
summary(modE2v)
nparam(modE2v)

modE3v=vglm(resp~planta*time+respp, family=multinomial,data=dadosl)
summary(modE3v)
nparam(modE3v)

modE4v=vglm(resp~(planta+respp)*time, family=multinomial,data=dadosl)
summary(modE4v)
nparam(modE4v)



# Five Models supposing non stationarity

mod12v<-vglm(resp2~planta+resp1,family='multinomial',data=dados)
mod23v<-vglm(resp3~planta+resp2,family='multinomial',data=dados)
mod34v<-vglm(resp4~planta+resp3,family='multinomial',data=dados)
mod45v<-vglm(resp5~planta+resp4,family='multinomial',data=dados)
mod56v<-vglm(resp6~planta+resp5,family='multinomial',data=dados)

summary(mod12v)
summary(mod23v)
summary(mod34v)
summary(mod45v)
summary(mod56v)

#######################################################################################################

#Exploratory data analysis - frequencies table

dadosa<-dados[dados$planta=='soya',];dadosa
dadosb<-dados[dados$planta=='cotton',];dadosb

head(dadosa)
head(dadosb)


#transition tables for soya
attach(dadosa)
tab1a<-table(resp1,resp2)
tab1a
sum(tab1a)

tab2a<-table(resp2,resp3)
tab2a
sum(tab2a)

tab3a<-table(resp3,resp4)
tab3a
sum(tab3a)

tab4a<-table(resp4,resp5)
tab4a
sum(tab4a)


tab5a<-table(resp5,resp6)
tab5a
sum(tab5a)

# Transition table soya (total)
tabEa<-(tab1a+tab2a+tab3a+tab4a+tab5a)
sum(tab5a)
detach(dadosa)


# Transition tables for cotton
attach(dadosb)
tab1b<-table(resp1,resp2)
tab1b
sum(tab1b)

tab2b<-table(resp2,resp3)
tab2b
sum(tab2b)

tab3b<-table(resp3,resp4)
tab3b
sum(tab3b)

tab4b<-table(resp4,resp5)
tab4b
sum(tab4b)


tab5b<-table(resp5,resp6)
tab5b
sum(tab5b)

#Transition table cotton (total)
tabEb<-(tab1b+tab2b+tab3b+tab4b+tab5b);tabEb
sum(tabEb)

detach(dadosb)

#Exploratory data analysis - frequencies table

dadosa1<-dadosl[dadosl$planta=='soya',];dadosa1
dadosb1<-dadosl[dadosl$planta=='cotton',];dadosb1
attach(dadosa1)
tab_soya<-table(respp,resp);tab_soya
detach(dadosa1)

attach(dadosb1)
tab_cotton<-table(respp,resp);tab_cotton
detach(dadosb1)

