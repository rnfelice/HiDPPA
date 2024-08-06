library(phylopath)
library(geomorph)
library(geiger)
library(MASS)
library(tidyverse)
library(mvMORPH)
source("functions.R")

#simulate tree
ntaxa <-200
tree<-sim.bdtree(b=1,d=0,n=200)



#Now define the stength of the links between nodes for the simuation
b13 <- 0.65
b23 <- 0.48
b3y <- 0.55


#Define the number of dimensions for each node:
dim_x1 <-1
dim_x2 <-20
dim_x3 <-20
dim_y <-1

#define DAG models:
candidates <- define_model_set(
  A = c(y ~ x3, x3 ~ x1, x3 ~ x2),
  B = c(y ~ x3, x3 ~ x2, x2 ~ x1),
  C = c(y ~ x2, x2 ~ x1, x2 ~ x3),
  D = c(y ~ x1, x1 ~ x2, x1 ~ x3)
)

#define the number of simualtions to run
nsim<-10
#define an empty vector for recording model specification rate
hits<-c()

for (i in 1:nsim){
#simulate the data. this will have to be repeated a 1000x later
mean_x1 <- 10
x1<-mvrnorm(n=dim_x1,Sigma=1*vcv(tree),  mu=rep(mean_x1,length(tree$tip)))

mean_x2 <- 0
x2<-mvrnorm(n=dim_x2,Sigma=1*vcv(tree),  mu=rep(mean_x2,length(tree$tip)))

mean_x3 <- 0
x3<-(as.vector(x1)*b13)+(as.matrix(x2)*b23)+(as.matrix(mvrnorm(n=dim_x3,Sigma=1*vcv(tree),  mu=rep(mean_x3,length(tree$tip)))))

mean_y <- 10
y<-colSums(as.matrix(x3)*b3y)+(as.vector(mvrnorm(n=dim_y,Sigma=1*vcv(tree),  mu=rep(mean_y,length(tree$tip)))))

list1 <- list(y = as.matrix(y), x1 = as.matrix(x1), x2 = t(x2), x3 = t(x3))
rm(y, x1, x2, x3)


p <- HiDPPA(candidates, list1, tree)
summ1<-hiDsummary(p)
bestmodel<-rownames(summ1)[which.min(summ1$delta_CICc)]
hits<-c(hits, bestmodel)
}
write.csv(as.matrix(hits), file = paste0("testresults_",ntaxa,"sp_",dim_x2,"p.csv"))