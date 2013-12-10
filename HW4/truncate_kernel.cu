#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{  
  __global__ void
  truncnorm_robert_kernel(
    
    float *x,      // Vector to contain returned samples 
    int n,         // Number of samples to retur
    float *mu,     // Vector of mu_s
    float *sigma,  // Vector of sigma_s
    float *a,      // Vector of lower-truncation values
    float *b,      // Vector of upper-truncation values
    int *mode,      // 1:(inf,inf); 2:(inf,b); 3:(a,inf);4:(a,b)
    int rng_a,     //input parameters as RNG seed constant
    int rng_b,     //input parameters as RNG seed constant
    int rng_c)      //input parameters as RNG seed constant
{
    int myblock = blockIdx.x + blockIdx.y * gridDim.x;
    int blocksize = blockDim.x * blockDim.y * blockDim.z;
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) +threadIdx.y*blockDim.x + threadIdx.x;
    int idx = myblock * blocksize + subthread;
    int nthreads <- prod(grid_dims)*prod(block_dims)
    double e = 4.71828183; double ai,bi,y,alpha,u,gx,z,t2; int accept;
    
    // if(n>=nthreads){idx=idx-nthreads}
    if(idx<n){
      if(mode[idx]==2){
        // set up the random number generating
        curandState rng;
	curand_init(rng_a+idx*rng_b,rng_c,0,&rng);
        x[idx]=mu[idx]+sqrt(sigma[idx])*curand_normal(&rng);
      }
      
      else if (mode[idx]==0){
        // set up the random number generating
        curandState rng;curand_init(rng_a+idx*rng_b,rng_c,0,&rng);
        //mirrow
         ai=(-b[idx]+mu[idx])/sqrt(sigma[idx]);
	 accept=0;
        //Then sample the truncated normal
        while(accept==0){
          //generate alpha
          alpha=(ai+sqrt(ai*ai+4.0))/2.0;
          //generate the z using inverse CDF
           y=curand_uniform(&rng); z=ai-1/alpha*log(1-y); t2=0;
          while(z<ai){
            t2=t2+1;curand_init(rng_a+idx*rng_b,rng_c+t2,0,&rng);
            y=curand_uniform(&rng);z=ai-1/alpha*log(1-y);  
          }
          //calculate g(x)
          if(ai<alpha){
	    gx=pow(e,-pow(alpha-z,2)/2);
          }
          else {
            gx=pow(e,-pow(ai-alpha,2)/2-pow(alpha-z,2)/2);
          }
          //decide whether to adopt z 
          u=curand_uniform(&rng);
          if(u<gx){
            accept=1;x[idx]=sqrt(sigma[idx])*z-mu[idx];x[idx]=-x[idx];
          }
        }
      }
      // mode==3:(a,inf)
      else if (mode[idx]==1){
        // set up the random number generating
        curandState rng;curand_init(rng_a+idx*rng_b,rng_c,0,&rng);
        //Then sample the truncated normal
        accept=0;ai=(a[idx]-mu[idx])/sqrt(sigma[idx]);
        while(accept==0){
          //generate alpha
         alpha=(ai+sqrt(ai*ai+4.0))/2.0;
          //generate the z using inverse CDF
          y=curand_uniform(&rng);z=ai-1/alpha*log(1-y);t2=0;
          while(z<ai){
            t2=t2+1;curand_init(rng_a+idx*rng_b,rng_c+t2,0,&rng);
            y=curand_uniform(&rng);z=ai-1/alpha*log(1-y);  
          }
          //calculate g(x)
          if(ai<alpha){
            gx=pow(e,(-pow((alpha-z),2)/2));
          }
          else {
            gx=pow(e,(pow(ai-alpha,2)/2))*pow(e,(-pow(alpha-z,2)/2));
          }
          //decide whether to adopt z 
          u=curand_uniform(&rng);
          if(u<gx){
            x[idx]=z;accept=1;x[idx]=sqrt(sigma[idx])*x[idx]+mu[idx];
          }
        } 
      }

      else if (mode[idx]==3){
        // set up the random number generating
        curandState rng;curand_init(rng_a+idx*rng_b,rng_c,0,&rng);
        ai=(a[idx]-mu[idx])/sqrt(sigma[idx]); bi=(b[idx]-mu[idx])/sqrt(sigma[idx]);
        //Then sample the truncated normal
        accept=0;
        while(accept==0){
          //generate the z using uniform
          z=ai+(bi-ai)*curand_uniform(&rng);
          //calculate g(x)
          if(ai<=0 && bi>=0){
            gx=pow(e,-z*z/2);
          }
          else if (bi<0) {
            gx=pow(e,((bi*bi-z*z)/2));
          }
          else if(ai>0){
            gx=pow(e,((ai*ai-z*z)/2));
          }
          //decide whether to adopt z 
          u=curand_uniform(&rng);
          if(u<gx){
            x[idx]=z;accept=1;x[idx]=sqrt(sigma[idx])*x[idx]+mu[idx];
            }
        } 
      }
    }  
    return;
  }  
  
} // END extern "C"




