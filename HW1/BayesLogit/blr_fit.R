library("lattice")
library("coda")
library("MASS")
library("mvtnorm")

#set the data reading in Gauss
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



# the bayes.logreg function, output  a list{beta_1, beta_2, and v_sq}
bayes.logreg<-function(m,y,X,beta.0, Sigma.0.inv,
                       niter=10000,burnin=1000,print.every=1000,
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
    # sample beta.star from the normal distribution with mean=beta.cur and Sigma=v_sq
    beta.star<-mvrnorm(n=1,mu=beta.cur, Sigma=v_sq)
    #calculate the accept.prob
    pdf.y.beta.cur=sum(dbinom(y,size=m,prob=logit(X%*%beta.cur),log=TRUE))
    pdf.y.beta.star=sum(dbinom(y,size=m,prob=logit(X%*%beta.star),log=TRUE))
    pdf.beta.star=dmvnorm(as.vector(beta.star),mean=c(0,0),sigma=diag(1,2),log=TRUE)
    pdf.beta.cur=dmvnorm(as.vector(beta.cur),mean=c(0,0),sigma=diag(1,2),log=TRUE)
    accept.prob=min(0,pdf.y.beta.star+pdf.beta.star-(pdf.y.beta.cur+pdf.beta.cur))
    #decide whether to accept beta.star or not
    if(log(runif(1))< accept.prob)
      return(as.matrix(beta.star,nrow=2))
    else 
      return(beta.cur)
  }
  
  
  # burnin period and tune
  #set up the initial value for the burnin period
  beta.cur=beta.0
  v_sq=diag(1,2)
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
  
  return(list(draws[,1],draws[,2],v_sq_ls,accpet/niter))
  
}


## run the bayes.logreg function using 
result<-bayes.logreg(m=data1$n,y=data1$y,
                     X=as.matrix(data1[,c(3,4)]),
                     beta.0=as.matrix(c(0,0),nrow=2),
                     Sigma.0.inv=solve(diag(1,2)))

#Trellis plot for mcmc result
plot(mcmc(result[[1]],start=1001))
plot(mcmc(result[[2]],start=1001))
result[[3]]
result[[4]]
beta1_q=quantile(result[[1]],probs=seq(0.01,0.99,0.01))
beta2_q=quantile(result[[2]],probs=seq(0.01,0.99,0.01))

#write.csv(result, file=paste0("results/blr_res_",sim_num,".csv"),col.names=FALSE,row.names=FALSE)
write.table(cbind(beta1_q,beta2_q), file=paste0("results/blr_res_",sim_num,".csv"),quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")

#write.table(data.frame(data here),file=filename,sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)
