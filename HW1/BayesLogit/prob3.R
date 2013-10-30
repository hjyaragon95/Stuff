library("lattice")
library("coda")
library("MASS")
library("mvtnorm")


#load in data and see the basic of the data
data<-read.table("/Users/aragon95/Google_Drive/STA250/homework/breast_cancer.txt",header=TRUE)
cancer<-data

#format the data
#change M&F value
cancer[,11]=sapply(cancer[,11],as.integer)-1
#standardize
#cancer[,c(1:10)]=sapply(cancer[,c(1:10)], function(x){(x-mean(x))/sd(x)})
cancer$n=rep(1,nrow(cancer))
cancer=cancer[,c(12,1:11)]

#fit logistic regression model
   #fmla<-paste0(names(cancer[,1:11]),sep="+")
glm.cancer<-glm(diagnosis~.,data=cancer[,-1],family=binomial)
summary(glm.cancer)
sigma=cov(cancer[,1:11])



# Fit the Bayesian model:
# set the initial value


# the bayes.logreg function, output  a list{beta_1, beta_2, and v_sq}
cancer.bayes.logreg<-function(m,y,X,beta.0, Sigma.0,psigma,
                       niter, burnin,print.every=1000,
                       retune=100,verbose=TRUE) {
  #intial a dataset to store beta 
  draws<-data.frame()
  #set the datatype to better fit my code
  X=as.matrix(X)
  beta.0=as.matrix(beta.0,nrow=2)
  
  #logit function
  logit<-function(u){
    ratio=exp(u)/(1+exp(u))
  }
  
  #beta.update function: 
  beta.update<-function(beta.cur,v_sq){
    #browser()
    # sample beta.star from the normal distribution with mean=beta.cur and Sigma=v_sq
    beta.cur<-as.matrix(beta.cur)
    beta.star<-rmvnorm(n=1,mean=beta.cur,v_sq)
    beta.star<-t(as.matrix(beta.star))
    #calculate the accept.prob
    pdf.y.beta.cur=sum(dbinom(y,size=m,prob=logit(X%*%beta.cur),log=TRUE))
    pdf.y.beta.star=sum(dbinom(y,size=m,prob=logit(X%*%beta.star),log=TRUE))
    pdf.beta.star=dmvnorm(as.vector(beta.star),mean=beta.0,sigma=Sigma.0,log=TRUE)
    pdf.beta.cur=dmvnorm(as.vector(beta.cur),mean=beta.0,sigma=Sigma.0,log=TRUE)
    accept.prob=min(0,pdf.y.beta.star+pdf.beta.star-(pdf.y.beta.cur+pdf.beta.cur))
    #decide whether to accept beta.star or not
    if(is.na(accept.prob))
      {
         return(beta.cur)
         accept.prob=-1
      }
    else if(log(runif(1))< accept.prob)
         return(beta.star)
    else if(log(runif(1))> accept.prob)
      return(beta.cur)
  }
  #a<-beta.update(beta.0,sigma)
  
  # burnin period and tune
  #set up the initial value for the burnin period
  beta.cur=beta.0
  v_sq=psigma
  v_sq_ls<-list()
  accept=0
  #for loop for the burn in and tuning

  for(i in 1: burnin){
    betaup1<-beta.update(beta.cur,v_sq)
    #count the acceptence time 
    if(!all(betaup1==beta.cur))
      accept=accept+1
    else
      accpet=accept
    #record the updated beta
    draws<- rbind(draws, t(as.matrix(beta.cur<-betaup1))) 
    # update v_sq
    if(i%% retune==0){
      v_sq_ls[[i %/% retune]]=v_sq
      if(accept/retune<0.3 )
        v_sq=v_sq/exp(1)
      else if(accept/retune>0.6)
        v_sq=v_sq*exp(1)
      accept=0
    }
  }
  
 
  ##############iteration period
  accept=0
  for(j in 1:niter){
    betaup1<-beta.update(beta.cur,v_sq)
    #count the accept time
    if(!all(betaup1==beta.cur))
      accept=accept+1
    else
      accpet=accept
    #store the result in the draws
    draws<- rbind(draws, t(beta.cur<-betaup1))
    #print every 1000 iteration is finished
    if(j%%1000==0){cat(j,'finished','\n')}
  }  
  
  return(list(draws,v_sq_ls))
  
}


###############################
# set the original value for bayes.logreg
m=cancer[,1]
y=cancer[,12]
X=as.matrix(cbind(rep(1,nrow(cancer)),cancer[,c(2:11)]))
beta.0=as.matrix(rep(0,11))
Sigma.0=diag(1000,11)
###############################


## run the bayes.logreg function using 
result<-cancer.bayes.logreg(m=cancer[,1],y=cancer[,12],
                     X=as.matrix(cbind(rep(1,nrow(cancer)),cancer[,c(2:11)])),
                     beta.0=as.matrix(rep(0,11)),
                     Sigma.0=diag(1000,11),psigma=10*sigma,
                     niter=10000,burnin=1000)
result1<-cancer.bayes.logreg(m=cancer[,1],y=cancer[,12],
                            X=as.matrix(cbind(rep(1,nrow(cancer)),cancer[,c(2:11)])),
                            beta.0=as.matrix(rep(0,11)),
                            Sigma.0=diag(1000,11),psigma=10*sigma,
                            niter=100000,burnin=1000)

#Trellis plot for mcmc result
par(mfrow=c(3,3))
sapply(result[[1]],function(x){
  plot(mcmc(x,start=5001))
})

sapply(result[[1]][-1],function(x){
  acf(x,plot=F,lag.max=1)
})

# see which covariates is more related tot the cancer diagonosis
cor(cancer[-1])

#acf(result[[1]][10],plot=F,lag.max=1)


# posterie predictive check. 
beta.pred<-sapply(result[[1]],function(x){
            mean(x[7000:11000])
        })
pred.y=as.matrix(cancer[,-12]) %*% as.matrix(beta.pred)
logit<-function(u){
  ratio=exp(u)/(1+exp(u))
  return(ratio)
}
real.y<-logit(cancer[,12])
ss.error<-sum((pred.y-real.y)^2/var(real.y))/10000
head(cbind(pred.y,real.y))








