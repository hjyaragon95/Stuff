#library("gdata")
library("MASS")
library("mvtnorm")
#library("MCMCpack")
#library("debug")
#data1<-read.csv("/Users/aragon95/Stuff/hw1/bayeslogit/data/blr_data_1001.csv")
#data2<-read.csv("/Users/aragon95/Stuff/hw1/bayeslogit/data/blr_pars_1001.csv")
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)


sim_start <- 1000
length.datasets <- 200
if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}
filename=paste0("data/blr_data_",sim_num,".csv")
data1<-read.csv(filename)
#version 2   with 


bayes.logreg<-function(m,y,X,beta.0, Sigma.0.inv,
                       niter=10000,burnin=1000,print.every=1000,
                       retune=100,verbose=TRUE,sim_num) {
  #intial a dataset to store beta 
  draws<-data.frame()
  X=as.matrix(X)
  beta.0=as.matrix(beta.0,nrow=2)
  #logit function
  logit<-function(u){
    ratio=exp(u)/(1+exp(u))
  }
  #browser()
  #beta.update function
  beta.update<-function(beta.cur,v_sq){
    #browser()
    beta.star<-mvrnorm(n=1,mu=beta.cur, Sigma=v_sq)
    pdf.y.beta.cur=sum(dbinom(y,size=m,prob=logit(X%*%beta.cur)/(1+logit(X%*%beta.cur)),log=TRUE))
    pdf.y.beta.star=sum(dbinom(y,size=m,prob=logit(X%*%beta.star)/(1+logit(X%*%beta.star)),log=TRUE))
    pdf.beta.star=dmvnorm(as.vector(beta.star),mean=c(0,0),sigma=diag(1,2),log=TRUE)
    pdf.beta.cur=dmvnorm(as.vector(beta.cur),mean=c(0,0),sigma=diag(1,2),log=TRUE)
    accept.prob=min(0,pdf.y.beta.star+pdf.beta.star-(pdf.y.beta.cur+pdf.beta.cur))
    if(log(runif(1))< accept.prob)
      return(as.matrix(beta.star,nrow=2))
    else 
      return(beta.cur)
  }
  # burnin period and tune
  beta.cur=beta.0
  v_sq=diag(1,2)
  v_sq_ls<-list()
  accept=0
  for(i in 1: burnin){
    betaup1<-beta.update(beta.cur,v_sq)
    if(!all(betaup1==beta.cur))
      accept=accept+1
    else
      accpet=accept
    draws<- rbind(draws, t(beta.cur<-betaup1)) 
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
  
  
  #browser()
  ##############iteration 
  for(j in 1:niter){
    draws<- rbind(draws, t(beta.cur<-beta.update(beta.cur,v_sq)))
    if(j%%1000==0){cat(j,'finished','\n')}
  }  
  
  beta1_q=quantile(draws[-c(1:1000),1],prob=seq(0.01,0.99,0.01))
  beta2_q=quantile(draws[-c(1:1000),2],prob=seq(0.01,0.99,0.01))
  return(cbind(beta1_q,beta2_q))
  write.csv(as.data.frame(cbind(beta1_q,beta2_q)), file=paste0("results/blr_res_",sim_num,".csv"))
  
}


result<-bayes.logreg(m=data1$n,y=data1$y,
                     X=as.matrix(data1[,c(3,4)]),
                     beta.0=as.matrix(c(0,0),nrow=2),
                     Sigma.0.inv=solve(diag(1,2)))
#write.csv(result, file=paste0("results/blr_res_",sim_num,".csv"),col.names=FALSE,row.names=FALSE)
write.table(result, file=paste0("results/blr_res_",sim_num,".csv"),quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")

#write.table(data.frame(data here),file=filename,sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)
