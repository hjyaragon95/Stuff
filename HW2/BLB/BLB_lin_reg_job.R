
#============================== Run the simulation study ==============================#

# Load packages:
mini <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:
  if((sim_num-1000)%%50!=0){
    s_index=(sim_num-1000)%/%50+1
  } else {
    s_index=(sim_num-1000)%/%50
  }
r_index=(sim_num-1000)-(s_index-1)*50

#============================== Run the simulation study ==============================#

# Load packages:
library("BH")
library("bigmemory.sri")
library("bigmemory")
library("biganalytics")
library("bigtabulate")
 

# I/O specifications:
datapath <- "/home/pdbaines/data"
outpath <- "output/"

# mini or full?
if (mini){
        rootfilename <- "blb_lin_reg_mini"
} else {
         rootfilename <- "blb_lin_reg_data"
}

# Filenames:
# Set up I/O stuff:
# Attach big.matrix :
# Remaining BLB specs:
# Extract the subset:
# Reset simulation seed:
# Bootstrap dataset:
# Fit lm:
# Output file:
# Save estimates to file:


#########################
#load in data
#########################
# in python open("/home/pdbaines/data/blb_lin_reg_data.txt","rb")
#install.packages("bigmemory")
#data.desc<-read.csv("/Users/aragon95/Google_Drive/STA250/homework/homework2/blb_lin_reg_mini.txt",
#               header=FALSE)
#data.desc<-attach.big.matrix("/Users/aragon95/Google_Drive/STA250/homework/homework2/blb_lin_reg_mini.desc","rb")
#data.desc<-attach.big.matrix("/home/pdbaines/data/blb_lin_reg_data_mini.desc","rb")

filepath<-paste0(datapath,"/",rootfilename,".desc")
data.desc<-attach.big.matrix(filepath)

################
# Bagging and bootstrap
################

#bagging function, return the index being sampled with nrow=n
bag.fun<-function(dat,b){
    ind<-sample(1:nrow(dat),size=b,replace=FALSE)
    return(ind)
  }


#bagging and bootstrap function
blb.fun<-function(dat,b,m){
  bag.ind<-bag.fun(dat,b)
  boot.ind<-as.vector(rmultinom(1,size=m,prob=1/b*as.vector(rep(1,b))))
  
  #lm regression
  fmla<-paste0("V",ncol(dat),"~.")
  #tmp.df$boot.ind <- boot.ind
  reg1<-lm(formula(fmla),data=as.data.frame(dat[bag.ind,]),weight=boot.ind)
  return(reg1)
  }


#pre specification
#########################
dat<-data.desc
b<-(nrow(data.desc))^0.7
m<-nrow(data.desc)
############################

result<-blb.fun(dat,b,m)
outfile = paste0("output/","coef_",sprintf("%02d",s_index),
                 "_",sprintf("%02d",r_index),".txt")
write.table(result$coeff,outfile,row.names=FALSE,col.names=FALSE)

  

