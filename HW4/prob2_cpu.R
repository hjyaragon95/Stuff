library("mvtnorm")
library("RCUDA")
print("End library")

mini<-read.table("mini_data.txt",header=TRUE)
data01<-read.table("data_01.txt",header=TRUE)
data02<-read.table("data_02.txt",header=TRUE)
data03<-read.table("data_03.txt",header=TRUE)
data04<-read.table("data_04.txt",header=TRUE)
data05<-read.table("data_05.txt",header=TRUE)


print("Finished loading Data")


compute_grid <- function(N,sqrt_threads_per_block=16L,grid_nd=1){
    # if...
    # N = 1,000,000
    # => 1954 blocks of 512 threads will suffice
    # => (62 x 32) grid, (512 x 1 x 1) blocks
    # Fix block dims:
    block_dims <- c(as.integer(sqrt_threads_per_block), as.integer(sqrt_threads_per_block), 1L)
    threads_per_block <- prod(block_dims)
    if (grid_nd==1){
      grid_d1 <- as.integer(max(1L,ceiling(N/threads_per_block)))
      grid_d2 <- 1L
    } else {
      grid_d1 <- as.integer(max(1L, floor(sqrt(N/threads_per_block))))
      grid_d2 <- as.integer(ceiling(N/(grid_d1*threads_per_block)))
    }
    grid_dims <- c(grid_d1, grid_d2, 1L)
    return(list("grid_dims"=grid_dims,"block_dims"=block_dims))
}

print("compute_grid finished")

rtruncnorm_cpu<-function(mu,sigma, a,b,n,mode){
  x=c(rep(NA,n))
  #(inf,b)
  for(i in 1:n){
    if(mode[i]==0){
      ai=(-b[i]+mu[i])/sqrt(sigma[i])
      alpha=(ai+sqrt(ai*ai+4))/2
      accept=0
      while(accept==0){
        #sample z
        z=ai+rexp(1,rate=alpha)
        while(z<ai){
          z=ai+rexp(1,rate=alpha)
        }
        #calculate g(x)
        if(ai<alpha){
          gx=exp(-(alpha-z)^2/2)
        }
        else {
          gx=exp(-(alpha-z)^2/2)*exp((alpha-ai)^2/2)
        }
        #decide whether to adopt z 
        u=runif(1,min=0,max=1)
        if(u<gx){
          x[i]=-(sqrt(sigma[i])*z-mu[i])
          accept=1
        }
      }
    }
    
    #(a,inf)
    if(mode[i]==1){
      ai=(a[i]-mu[i])/sqrt(sigma[i])
      alpha=(ai+sqrt(ai*ai+4))/2
      accept=0
      while(accept==0){
        #sample z
        z=ai+rexp(1,rate=alpha)
        while(z<ai){
          z=ai+rexp(1,rate=alpha)
        }
        #calculate g(x)
        if(ai<alpha){
          gx=exp(-(alpha-z)^2/2)
        }
        else {
          gx=exp(-(alpha-z)^2/2)*exp((alpha-ai)^2/2)
        }
        #decide whether to adopt z 
        u=runif(1,min=0,max=1)
        if(u<gx){
          x[i]=sqrt(sigma[i])*z+mu[i]
          accept=1
        }
      }
    }
    #(inf,inf)
    if(mode[i]==2){
      x[i]=rnorm(1,mean=mu[i],sd=sigma[i])
    }
    #(a,b)
    if(mode[i]==3){
      ai=(a[i]-mu[i])/sqrt(sigma[i])
      bi=(b[i]-mu[i])/sqrt(sigma[i])
      accept=0
      while(accept==0){
        z=runif(1,ai,bi);
        if(ai<=0 && bi>=0){
          gx=exp(-z*z/2)
        }
        else if (bi<0) {
          gx=exp((bi*bi-z*z)/2)
        }
        else if(ai>0){
          gx=exp((ai*ai-z*z)/2)
        }
        #decide whether to adopt z 
        u=runif(1,0,1)
        if(u<gx){
          x[i]=sqrt(sigma[i])*z+mu[i]
          accept=1
        }
      }
    }  
  }
  return(x)  
}

print("trunnormal_cpu finished")

probit_mcmc_cpu<-function(
  y,           # vector of length n 
  X,           # (n x p) design matrix
  beta_0,      # (p x 1) prior mean
  Sigma_0_inv, # (p x p) prior precision 
  niter,       # number of post burnin iterations
  burnin)      # number of burnin iterations
{
  beta<-list()
  for(i in 1:(burnin)){
      if(i==1){
        if(i%%100==0)cat(paste0(i,"th job of in burnin has been finished.      "))
        z=rtruncnorm_cpu(mu=X%*%as.matrix(beta_0),sigma=rep(1,length(y)),
                         a=rep(0,length(y)),b=rep(0,length(y)),length(y),mode=y)
      }else{
      z=rtruncnorm_cpu(mu=X%*%as.matrix(unlist(beta[i-1])),sigma=rep(1,length(y)),
                   a=rep(0,length(y)),b=rep(0,length(y)),length(y),mode=y)
      }
      beta[[i]]<-rmvnorm(1,mean=solve(Sigma_0_inv+t(X)%*%X)%*%(Sigma_0_inv%*%beta_0+t(X)%*%as.matrix(z)),
                     sigma=solve(Sigma_0_inv+t(X)%*%X))
  }
  for(i in (burnin+1):(burnin+niter)){
    if(i%%100==0)cat(paste0(i,"th job of in niter has been finished.      "))
    z=rtruncnorm_cpu(mu=X%*%as.matrix(unlist(beta[i-1])),sigma=rep(1,length(y)),
                     a=rep(0,length(y)),b=rep(0,length(y)),length(y),mode=y)
    beta[[i]]<-rmvnorm(1,mean=solve(Sigma_0_inv+t(X)%*%X)%*%(Sigma_0_inv%*%beta_0+t(X)%*%as.matrix(z)),
                                                         sigma=solve(Sigma_0_inv+t(X)%*%X))
  }
  beta.out<-t(matrix(unlist(beta),nrow=8))
  return(beta.out)
}

print("probit_mcmc_cpu finished")
################################
#arguments for probit_mcmc_cpu
############################
#y<-mini[,1]
#X<-as.matrix(mini[,-1])
#beta_0<-as.matrix(rep(0,ncol(X)))
#Sigma_0_inv<-diag(0,ncol(X))
#niter<-500
#burnin<-2000


#time.cpu<-system.time(beta.out<-probit_mcmc_cpu(y,X,beta_0,Sigma_0_inv,niter,burnin))
#print(paste("Time for CPU is ", time.cpu))
#apply(beta.out[-c(1:10000),],2,mean)
#apply(beta.out[-c(1:100),],2,sd)



ptx = "rtruncnorm.ptx"
#We can then load the resulting PTX le with
mod = loadModule(ptx)
#Now that we have the kernel, we can obtain a reference to the kernel with
rtnorm = mod$truncnorm_robert_kernel

probit_mcmc_gpu<-function(
  y,           # vector of length n 
  X,           # (n x p) design matrix
  beta_0,      # (p x 1) prior mean
  Sigma_0_inv, # (p x p) prior precision 
  niter,       # number of post burnin iterations
  burnin,      # number of burnin iterations
  block_dims,  # block dimensions
  grid_dims){    # grid_dimensions
  
  beta<-list()
  N<- as.integer(length(y))
  out<- rep(0.0,N)
  a <- rep(0,N)
  b <- rep(0,N)
  mode<-as.integer(y)
  rng_a<-1L
  rng_b<-2L
  rng_c<-3L
  for(i in 1:(burnin)){
    if(i==1){
      mem <- copyToDevice(out)
      mu<-X%*%as.matrix(beta_0)
      sigma<-rep(1,N)
      .cuda(rtnorm, mem, N, mu, sigma, a, b,mode, rng_a, 
                                          rng_b, rng_c, gridDim=grid_dims, blockDim=block_dims)  
      z = mem[]
    }else{
      if(i%%50==0){cat("ith job of in niter has been finished.      ")}
      mem <- copyToDevice(out)
      mu<-X%*%as.matrix(unlist(beta[i-1]))
      sigma<-rep(1,N)
      .cuda(rtnorm, mem, N, mu, sigma, a, b,mode, rng_a,rng_b, rng_c, gridDim=grid_dims, blockDim=block_dims)  
      z = mem[]
    }
    beta[[i]]<-rmvnorm(1,mean=solve(Sigma_0_inv+t(X)%*%X)%*%(Sigma_0_inv%*%beta_0+t(X)%*%as.matrix(z)),
                       sigma=solve(Sigma_0_inv+t(X)%*%X))
  }
  for(i in (burnin+1):(burnin+niter)){
    if(i%%50==0) {cat("ith job of in niter has been finished.      ")}
    mem <- copyToDevice(out)
    mu<-X%*%as.matrix(unlist(beta[i-1]))
    sigma<-rep(1,N)
    .cuda(rtnorm, mem, N, mu, sigma, a, b,mode, rng_a,rng_b, rng_c, gridDim=grid_dims, blockDim=block_dims)  
    z = mem[]
    beta[[i]]<-rmvnorm(1,mean=solve(Sigma_0_inv+t(X)%*%X)%*%(Sigma_0_inv%*%beta_0+t(X)%*%as.matrix(z)),
                       sigma=solve(Sigma_0_inv+t(X)%*%X))
  }
  beta.out<-t(matrix(unlist(beta),nrow=8))
  return(beta.out)
   
}

print("probit mcmc gpu finished")



#y1<-mini[,1]
#X1<-as.matrix(mini[,-1])
#beta1_0<-as.matrix(rep(0,ncol(X1)))
#Sigma1_0_inv<-diag(0,ncol(X1))
#niter<-500
#burnin<-2000
#grid_dims1<-compute_grid(length(y1))$grid_dims
#block_dims1<-compute_grid(length(y1))$block_dims
#time.cpu1<-system.time(beta1.out.cpu<-probit_mcmc_cpu(y1,X1,beta1_0,Sigma1_0_inv,
#                                                 niter,burnin))[3]
#time.gpu1<-system.time(beta1.out.gpu<-probit_mcmc_gpu(y1,X1,beta1_0,Sigma1_0_inv,
#                                                 niter,burnin,block_dims1,grid_dims1))[3]
#print("data01 finished")
#write.table(beta1.out.cpu,file="data_mini_xcpu.txt",row.names=FALSE)
#write.table(beta1.out.gpu,file="data_mini_xgpu.txt",row.names=FALSE)

y1<-data01[,1]
X1<-as.matrix(data01[,-1])
beta1_0<-as.matrix(rep(0,ncol(X1)))
Sigma1_0_inv<-diag(0,ncol(X1))
niter<-500
burnin<-2000
grid_dims1<-compute_grid(length(y1))$grid_dims
block_dims1<-compute_grid(length(y1))$block_dims
time.cpu1<-system.time(beta1.out.cpu<-probit_mcmc_cpu(y1,X1,beta1_0,Sigma1_0_inv,
                                                 niter,burnin))[3]
time.gpu1<-system.time(beta1.out.gpu<-probit_mcmc_gpu(y1,X1,beta1_0,Sigma1_0_inv,
                                                 niter,burnin,block_dims1,grid_dims1))[3]
print("data01 finished")
write.table(beta1.out.cpu,file="data01_xcpu.txt",row.names=FALSE)
write.table(beta1.out.gpu,file="data01_xgpu.txt",row.names=FALSE)


y2<-data02[,1]
X2<-as.matrix(data02[,-1])
beta2_0<-as.matrix(rep(0,ncol(X2)))
Sigma2_0_inv<-diag(0,ncol(X2))
grid_dims2<-compute_grid(length(y2))$grid_dims
block_dims2<-compute_grid(length(y2))$block_dims
time.cpu2<-system.time(beta2.out.cpu<-probit_mcmc_cpu(y2,X2,beta2_0,Sigma2_0_inv,
                                                      niter,burnin))[3]
time.gpu2<-system.time(beta2.out.gpu<-probit_mcmc_gpu(y2,X2,beta2_0,Sigma2_0_inv,
                                                      niter,burnin,block_dims2,grid_dims2))[3]
print("data02 finished")
write.table(beta2.out.cpu,file="data02_xcpu.txt",row.names=FALSE)
write.table(beta2.out.gpu,file="data02_xgpu.txt",row.names=FALSE)

y3<-data03[,1]
X3<-as.matrix(data03[,-1])
beta3_0<-as.matrix(rep(0,ncol(X3)))
Sigma3_0_inv<-diag(0,ncol(X3))
grid_dims3<-compute_grid(length(y3))$grid_dims
block_dims3<-compute_grid(length(y3))$block_dims
time.cpu3<-system.time(beta3.out.cpu<-probit_mcmc_cpu(y3,X3,beta3_0,Sigma3_0_inv,
                                                      niter,burnin))[3]
time.gpu3<-system.time(beta3.out.gpu<-probit_mcmc_gpu(y3,X3,beta3_0,Sigma3_0_inv,
                                                      niter,burnin,block_dims3,grid_dims3))[3]

print("data03 finished")
write.table(beta3.out.cpu,file="data03_xcpu.txt",row.names=FALSE)
write.table(beta3.out.gpu,file="data03_xgpu.txt",row.names=FALSE)




y4<-data04[,1]
X4<-as.matrix(data04[,-1])
beta4_0<-as.matrix(rep(0,ncol(X4)))
Sigma4_0_inv<-diag(0,ncol(X4))
grid_dims4<-compute_grid(length(y4))$grid_dims
block_dims4<-compute_grid(length(y4))$block_dims
time.cpu4<-system.time(beta4.out.cpu<-probit_mcmc_cpu(y4,X4,beta4_0,Sigma4_0_inv,
                                                      niter,burnin))[3]
time.gpu4<-system.time(beta4.out.gpu<-probit_mcmc_gpu(y4,X4,beta4_0,Sigma4_0_inv,
                                                      niter,burnin,block_dims4,grid_dims4))[3]
print("data04 finished")
write.table(beta4.out.cpu,file="data04_xcpu.txt",row.names=FALSE)
write.table(beta4.out.gpu,file="data04_xgpu.txt",row.names=FALSE)



y5<-data05[,1]
X5<-as.matrix(data05[,-1])
beta5_0<-as.matrix(rep(0,ncol(X5)))
Sigma5_0_inv<-diag(0,ncol(X5))
grid_dims5<-compute_grid(length(y5))$grid_dims
block_dims5<-compute_grid(length(y5))$block_dims
time.cpu5<-system.time(beta5.out.cpu<-probit_mcmc_cpu(y5,X5,beta5_0,Sigma5_0_inv,
                                                      niter,burnin))[3]
time.gpu5<-system.time(beta5.out.gpu<-probit_mcmc_gpu(y5,X5,beta5_0,Sigma5_0_inv,
                                                      niter,burnin,block_dims5,grid_dims5))[3]
print("data05 finished")
write.table(beta5.out.cpu,file="data05_xcpu.txt",row.names=FALSE)
write.table(beta5.out.gpu,file="data05_xgpu.txt",row.names=FALSE)

 write.table(c(time.cpu1,time.gpu1,time.cpu2,time.gpu2,time.cpu3,time.gpu3,time.cpu4,time.gpu4,time.cpu5,time.cpu5), file="time_for_prob2.txt")


