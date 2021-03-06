## The purpose of this script is to compare logistic_GEM, which just outputs the trait distribution every dt time steps, to logistic_GEM_storage, which stores every individual ever born, recording its time of birth, time of death, and parent's traits. This storage allows us to do things like look at parent-offspring regressions and, more importantly, to compute the fitness of every individual ever born to look at the relationship between traits and fitness.
seed <- 10110
tmax <- 200
N0 <- 5
traitmean <- 1.8
traitcv <- 0.3
h2 <- 0.75
bs <- ds <- 0.02
slope <- 0.3/1.8^2

source("logistic_GEM.R")
out1 <- logistic_GEM(seed, dt=0.1, tmax, N0, traitmean, traitcv, h2, bs, ds, slope)

source("logistic_GEM_with_storage.R")
out2 <- logistic_GEM_storage(seed, tmax, N0, traitmean, traitcv, h2, bs, ds, slope)

## Compare the dynamics
time.seq <- seq(0, tmax, 0.1)
sapply(1:length(time.seq),
       function(t) {
           ## who is alive?
           if (t < length(time.seq))
               alive <- subset(out2, tBirth <= time.seq[t] & tDeath > time.seq[t])
           else alive <- subset(out2, is.na(tDeath))
           ## at the current moment, how big is the population size and what are the mean and SD of the current trait distribution?
           c(nrow(alive), mean(alive$trait), sd(alive$trait))
       }) %>% t %>% as.data.frame %>% mutate(., time=time.seq) -> traj2
colnames(traj2)[1:3] <- c("N","mean","sd")

## Plot the two
plot(time.seq, traj2$N, type='l')
lines(time.seq, lapply(out1[[2]], length) %>% unlist, col=2)
