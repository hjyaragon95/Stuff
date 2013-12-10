#!/usr/bin/Rscript
library("RCUDA")

# N = 1,000,000
# => 1954 blocks of 512 threads will suffice
# => (62 x 32) grid, (512 x 1 x 1) blocks


############################
#compile
############################
ptx = "rtruncnorm.ptx"
#We can then load the resulting PTX le with
mod = loadModule("rtruncnorm.ptx")
#Now that we have the kernel, we can obtain a reference to the kernel with
my_kernel = mod$truncnorm_robert_kernel



############################
# Fix block dims:
############################
# Function for computing default grid/block sizes
# Example usage:
# bg <- computed_grid(N=10000)
# grid_dims <- bg$grid_dims
# block_dims <- bg$block_dims
##

compute_grid <- function(N,sqrt_threads_per_block=16L,grid_nd=1)
{
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




############################
#rtruncnorm_cpu()
############################

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
print("all the function is loaded")



#################
#Part C  try TN(2,1;(0,1.5))
#################

#Arguments
N<- 10000L
x<- rep(0.0,N)
mu <- rep(2,N)
sigma <-rep(1,N)
a <- rep(0,N)
b <- rep(1.5,N)
mode<-rep(3L,N)
rng_a<-1L
rng_b<-2L
rng_c<-3L
grid_dims<-compute_grid(N)$grid_dims
block_dims<-compute_grid(N)$block_dims
nthreads <- prod(grid_dims)*prod(block_dims)
grid=FALSE
if (nthreads < N && grid==TRUE ){
  stop("Grid is not large enough...!")
}

print("begin part c gpu calculation")

mem <- copyToDevice(x)
.cuda(my_kernel, mem, N, mu, sigma, a, b,mode, rng_a,rng_b, rng_c, gridDim=grid_dims, blockDim=block_dims)
x_try_gpu<-mem[]
plot(density(x_try_gpu))
hist(x_try_gpu)
print("Part C over")

          
#####################
#Part D  CPU try
#####################
x.try.cpu<-rtruncnorm_cpu(mu,sigma,a,b,N,mode)
plot(density(x.try.cpu))
hist(x.try.cpu)
write.table(cbind(x_try_gpu,x.try.cpu),file="x_cpu_gpu.txt",row.name=FALSE)
print("Part D over")

###############################
#PArt E 
###############################
k=8
Nk<-sapply(c(1:k),function(x){10^x})
gpu.time<-c(rep(NA,k))
cpu.time<-c(rep(NA,k))
for (i in 1:k){
  N<-as.integer(Nk[i])
  x<- rep(0.0,N)
  mu <- rep(2,N)
  sigma <-rep(1,N)
  a <- rep(0,N)
  b <- rep(1.5,N)
  mode<-rep(3L,N) 
  rng_a<-1L
  rng_b<-2L
  rng_c<-3L
  grid_dims<-compute_grid(N)$grid_dims
  block_dims<-compute_grid(N)$block_dims
  nthreads <- prod(grid_dims)*prod(block_dims)
    
  print(paste0("The k for the kernel execution is ", i))
  gpu_copyto_time <- system.time(mem <- copyToDevice(x))[3]
  gpu_kernel_time <- system.time(.cuda(my_kernel, mem, N, mu, sigma, a, b,mode, rng_a,
                                      rng_b, rng_c, gridDim=grid_dims, blockDim=block_dims))[3]
  gpu_copyfrom_time<-system.time(x.gpu<-mem[])[3]
  gpu.time[i]<-gpu_copyto_time+gpu_kernel_time+gpu_copyfrom_time
  cpu.time[i]<-system.time(x.cpu<-rtruncnorm_cpu(mu,sigma,a,b,N,mode))[3]
   
}



png("gpu_cpu_time.png")
plot(cpu.time,x=c(1:k),xlab="n",ylab="time",main="gpu_cpu_time",col="red")
lines(gpu.time,col="green")

dev.off()
write.table(cbind(gpu.time,cpu.time),file="x_cpu_gpu_time.txt",row.names=FALSE)


print("Part e over")
###########################################################
#Part F veryfy (inf, inf works) for cpu and gpu
###########################################################
mem <- copyToDevice(x)
.cuda(my_kernel, mem, N, mu, sigma, a, b,rep(2L,N), rng_a,rng_b, rng_c, gridDim=grid_dims, blockDim=block_dims)
x.try.gpu.inf<-mem[]
plot(density(x.try.gpu.inf))
hist(x.try.gpu.inf)

x.try.cpu.inf<-rtruncnorm_cpu(mu,sigma,a,b,N,rep(2,N))
plot(density(x.try.cpu.inf))
hist(x.try.cpu.inf)
print("Part F over")
write.table(cbind(x.try.gpu.inf,x.try.cpu.inf),file="prob1_F_inf.txt",row.names=FALSE)


###########################################################
#Part G: TN(0,1,-5,-3)
###########################################################

N<- 10000L
x<- rep(0.0,N)
mu <- rep(0,N)
sigma <-rep(1,N)
a <- rep(-5,N)
b <- rep(-3,N)
mode<-rep(3L,N)
rng_a<-1L
rng_b<-2L
rng_c<-3L
grid_dims<-compute_grid(N)$grid_dims
block_dims<-compute_grid(N)$block_dims
nthreads <- prod(grid_dims)*prod(block_dims)
grid=FALSE
if (nthreads < N && grid==TRUE ){
  stop("Grid is not large enough...!")
}


mem <- copyToDevice(x)
.cuda(my_kernel, mem, N, mu, sigma, a, b,mode, rng_a,rng_b, rng_c, gridDim=grid_dims, blockDim=block_dims)
x.try2.gpu<-mem[]
plot(density(x.try2.gpu))
hist(x.try2.gpu)

x.try2.cpu<-rtruncnorm_cpu(mu,sigma,a,b,N,mode)
plot(density(x.try2.cpu))
hist(x.try2.cpu)
write.table(cbind(x.try2.gpu,x.try2.cpu),file="prob1_g.txt",row.names=FALSE)
print("Part G over")





