source("logistic_GEM.R")
library(tidyverse)
library(parallel)

## Start at the ecological equilibrium given the trait value
bmax_seq <- seq(0.5,5.5,0.5)
## choose the CV so that the trait SD remains constant at 0.2, even as bmax increases
## note: traitsd = traitmean*traitcv
cv_seq <- 0.2/bmax_seq
s <- 1/6
bs <- ds <- 0.0125
N0_seq <- sapply(bmax_seq, function(b) (b-s*b^2)/(bs+ds))
output1 <- vector(mode='list', length=length(bmax_seq))
for (i in 1:length(bmax_seq)) {
    print(i)
    tin = Sys.time()
  set.seed(1343421)
  seeds <- floor(runif(100, 1, 1e5))
  mclapply(seeds,
           function(seed) logistic_GEM(seed=seed, dt=1, tmax=1000, N0=N0_seq[i], traitmean=bmax_seq[i], traitcv=cv_seq[i], h2=0.75, bs=bs, ds=ds, slope=s),
           mc.cores=8
  ) -> out
tout = Sys.time()
print(tout-tin)
  ## create a dataframe storing the results (the median population size across the replicates and the median mean trait)
  output1[[i]] <- data.frame(time=0:1000,
                             N=lapply(out, function(l) lapply(l$traits, function(j) ifelse(length(j)>0,length(j),NA)) %>% unlist) %>% do.call("cbind.data.frame",.) %>% apply(.,1,function(r) median(r,na.rm=TRUE)),
                             b=lapply(out, function(l) lapply(l$traits, function(t) mean(t,na.rm=TRUE)) %>% unlist) %>% do.call("cbind.data.frame",.) %>% apply(.,1,function(r) median(r,na.rm=TRUE)),
                             init=paste("b0=",signif(bmax_seq[i],3),"\n","N0=",signif(N0_seq[i],3),sep="")) %>%
    mutate(.,
           QGselGrad=1-2*s*b,
           LRSselGrad=(ds*N-b*(b-2*bs*N)*s)/(ds*N+b^2*s)^2)

}
output1 %>% do.call("rbind.data.frame",.) -> out1
saveRDS(out1, file="Logistic_GEM_bs=ds=0125_varying_b0_and_N0.RDS")


bmax_seq <- seq(0.5,5.5,0.5)
## choose the CV so that the trait SD remains constant at 0.2, even as bmax increases
## note: traitsd = traitmean*traitcv
cv_seq <- 0.2/bmax_seq
s <- 1/6
bs <- ds <- 0.025
N0_seq <- sapply(bmax_seq, function(b) (b-s*b^2)/(bs+ds))
output1 <- vector(mode='list', length=length(bmax_seq))
for (i in 1:length(bmax_seq)) {
    print(i)
    tin = Sys.time()
  set.seed(1343421)
  seeds <- floor(runif(100, 1, 1e5))
  mclapply(seeds,
           function(seed) logistic_GEM(seed=seed, dt=1, tmax=1000, N0=N0_seq[i], traitmean=bmax_seq[i], traitcv=cv_seq[i], h2=0.75, bs=bs, ds=ds, slope=s),
           mc.cores=8
  ) -> out
tout = Sys.time()
print(tout-tin)
  ## create a dataframe storing the results (the median population size across the replicates and the median mean trait)
  output1[[i]] <- data.frame(time=0:1000,
                             N=lapply(out, function(l) lapply(l$traits, function(j) ifelse(length(j)>0,length(j),NA)) %>% unlist) %>% do.call("cbind.data.frame",.) %>% apply(.,1,function(r) median(r,na.rm=TRUE)),
                             b=lapply(out, function(l) lapply(l$traits, function(t) mean(t,na.rm=TRUE)) %>% unlist) %>% do.call("cbind.data.frame",.) %>% apply(.,1,function(r) median(r,na.rm=TRUE)),
                             init=paste("b0=",signif(bmax_seq[i],3),"\n","N0=",signif(N0_seq[i],3),sep="")) %>%
    mutate(.,
           QGselGrad=1-2*s*b,
           LRSselGrad=(ds*N-b*(b-2*bs*N)*s)/(ds*N+b^2*s)^2)

}
output1 %>% do.call("rbind.data.frame",.) -> out1
saveRDS(out1, file="Logistic_GEM_bs=ds=025_varying_b0_and_N0.RDS")

