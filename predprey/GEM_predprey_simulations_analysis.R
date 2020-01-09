## mostly for plotting results from GEM_predprey_simulations.R and GEM_predprey_w_cull_simulations.R
library(magrittr)
library(plyr)
library(tidyr)
library(ggplot2)

paramsets <- vector(mode='list', length=5)
paramsets[[1]] <- c(K=250, a0=0.001, e=3.5, m=2, d=1)
paramsets[[2]] <- c(K=200, a0=0.001, e=2.5, m=1, d=1)
paramsets[[3]] <- c(K=150, a0=0.001, e=1.5, m=0.5, d=1)
paramsets[[4]] <- c(K=100, a0=0.001, e=0.5, m=0.1, d=1)
paramsets[[5]] <- c(K=50, a0=0.001, e=0.2, m=0.02, d=1)

for (i in 1:5) {
    params <- paramsets[[i]]
    assign(paste0("out",i), readRDS(file=paste0("predprey_GEM_K=",unname(params["K"]),"_a0=0.001_e=",unname(params["e"]),"_m=",unname(params["m"]),"_d=1.RDS")))
}

for (j in 1:5) {
     out <- get(paste0("out",j))
     lapply(1:length(out),
       function(i)
           data.frame(time=seq(1:length(out[[i]]$traits))/10-0.1,
                      C = lapply(out[[i]]$states, function(x) x[1]) %>% unlist,
                      P = lapply(out[[i]]$states, function(x) x[2]) %>% unlist,
                      med_b = lapply(out[[i]]$traits, median) %>% unlist,
                      sd_b = lapply(out[[i]]$traits, sd) %>% unlist,
                      rep=i)
            ) %>% do.call(rbind, .) -> out
     params <- paramsets[[j]]
     ## compute the QG expectation for C, P, and b
     Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
     Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
     beq <- unname(2*params["d"])


     png(file=paste0("GEM_predprey_simulations_figure_",j,".png"), height=5, width=5, units='in', res=300)
     par(mfrow=c(2,2), mar=c(2,4,0.5,0.5), oma=rep(0.5,4))
     plot.new()
     plot.window(xlim=range(out$time), ylim=range(out$C))
     axis(1); axis(2); box('plot')
     for (i in 1:max(out$rep))
         with(subset(out, rep==i), lines(time,C,col=gray(0.7),lwd=0.75))
     points(unique(out$time), sapply(unique(out$time), function(t) subset(out, time==t)$C %>% mean), col=4, lwd=2, pch=21, cex=0.25)
     abline(h=Ceq, lwd=1.6)
     mtext(side=2, line=2.5, "Prey abundance")

     plot.new()
     plot.window(xlim=range(out$time), ylim=range(out$P))
     axis(1); axis(2); box('plot')
     for (i in 1:max(out$rep))
         with(subset(out, rep==i), lines(time,P,col=gray(0.7),lwd=0.75))
     points(out$time %>% unique, sapply(out$time %>% unique, function(t) subset(out, time==t)$P %>% mean), col=4, lwd=2, pch=21, cex=0.25)
     abline(h=Peq, lwd=1.6)
     mtext(side=2, line=2.5, "Predator abundance")

     plot.new()
     plot.window(xlim=range(out$time), ylim=range(out$med_b, na.rm=T))
     axis(1); axis(2); box('plot')
     for (i in 1:max(out$rep))
         with(subset(out, rep==i), lines(time,med_b,col=gray(0.7),lwd=0.75))
     points(out$time %>% unique, sapply(out$time %>% unique, function(t) subset(out, time==t)$med_b%>% mean), col=4, lwd=2, pch=21, cex=0.25)
     mtext(side=2, line=2.5, "Median birth rate")
     abline(h=2, lwd=1.6)
     ## compute the tea using the average C and P over the last 50 time steps (across all reps)
     meanP <- apply(subset(out, time > 150), 2, mean)["P"]
     meanC <- apply(subset(out, time > 150), 2, mean)["C"]
     abline(h=(params["K"]-meanC)/(2*params["a0"]*params["K"]*meanP), col=2, lwd=1.6)
     legend(x='topright', c("ESS","TEA"), lwd=1.2, col=c(1,2), bty='n')

     plot.new()
     plot.window(xlim=range(out$time), ylim=range(out$sd_b, na.rm=T))
     axis(1); axis(2); box('plot')
     for (i in 1:max(out$rep))
         with(subset(out, rep==i), lines(time,sd_b,col=gray(0.7),lwd=0.75))
     points(out$time %>% unique, sapply(out$time %>% unique, function(t) subset(out, time==t)$sd_b%>% mean), col=4, lwd=2, pch=21, cex=0.25)
     mtext(side=2, line=2.5, "Trait SD")
     dev.off()
 }

## MODEL WITH CULLING CULLING CULLING
params <- c(K=250, a0=0.001, e=3.5, m=2, d=1)
out <- readRDS(file=paste0("predprey_GEM_cull_K=",unname(params["K"]),"_a0=0.001_e=",unname(params["e"]),"_m=",unname(params["m"]),"_d=1.RDS"))
lapply(1:length(out),
       function(i)
           data.frame(time=seq(1:length(out[[i]]$traits))/10-0.1,
                      C = lapply(out[[i]]$states, function(x) x[1]) %>% unlist,
                      P = lapply(out[[i]]$states, function(x) x[2]) %>% unlist,
                      med_b = lapply(out[[i]]$traits, median) %>% unlist,
                      sd_b = lapply(out[[i]]$traits, sd) %>% unlist,
                      rep=i)
       ) %>% do.call(rbind, .) -> out
mutate(out,
       tea=(params["K"]-C)/(2*params["a0"]*params["K"]*P)) -> out

## subset to just plot those that survived until the end
goodreps = subset(out, time > 499.9)$rep
out2 <- out
out <- subset(out2, rep%in%goodreps)

meanC <- sapply(unique(out$time), function(t) subset(out, time==t)$C %>% mean)
meanP <- sapply(out$time %>% unique, function(t) subset(out, time==t)$P %>% mean)
meanB <- sapply(out$time %>% unique, function(t) subset(out, time==t)$med_b%>% mean)
meanTEA <- sapply(out$time %>% unique, function(t) subset(out, time==t)$tea%>% mean)
meanV <- sapply(out$time %>% unique, function(t) subset(out, time==t)$sd_b%>% mean)

png(file=paste0("GEM_predprey_cull_simulations_figure_1.png"), height=5, width=5, units='in', res=300)
par(mfrow=c(2,2), mar=c(2,4,0.5,0.5), oma=rep(0.5,4))
plot.new()
plot.window(xlim=range(out$time), ylim=range(out$C))
axis(1); axis(2); box('plot')
for (i in 1:max(out$rep))
    with(subset(out, rep==i), lines(time,C,col=gray(0.7),lwd=0.75))
points(unique(out$time), meanC, col=4, pch=21, cex=0.25)
mtext(side=2, line=2.5, "Prey abundance")

plot.new()
plot.window(xlim=range(out$time), ylim=range(out$P))
axis(1); axis(2); box('plot')
for (i in 1:max(out$rep))
    with(subset(out, rep==i), lines(time,P,col=gray(0.7),lwd=0.75))
points(out$time %>% unique, meanP, col=4, pch=21, cex=0.25)
mtext(side=2, line=2.5, "Predator abundance")

plot.new()
plot.window(xlim=range(out$time), ylim=range(out$med_b, na.rm=T))
axis(1); axis(2); box('plot')
for (i in 1:max(out$rep))
    with(subset(out, rep==i), lines(time,med_b,col=gray(0.7),lwd=0.75))
points(out$time %>% unique, meanB, col=4, pch=21, cex=0.25)
mtext(side=2, line=2.5, "Median birth rate")
## add the TEA by averaging the final
abline(h=median(meanTEA), col=3, lwd=2)
abline(h=2, col=1, lwd=2)
legend(x='topleft', c('TEA','ESS'), lwd=2, col=c(3,1), bty='n')

plot.new()
plot.window(xlim=range(out$time), ylim=range(out$sd_b, na.rm=T))
axis(1); axis(2); box('plot')
for (i in 1:max(out$rep))
    with(subset(out, rep==i), lines(time,sd_b,col=gray(0.7),lwd=0.75))
points(out$time %>% unique, meanV, col=4, pch=21, cex=0.25)
mtext(side=2, line=2.5, "Trait SD")
dev.off()

## OLDER STUFF
for (i in 1:8) {
    this.a0 <- 0.0005*i
    assign(paste0("out",i), readRDS(file=paste0("predprey_GEM_K=500_a0=",this.a0,"_e=5_m=4_d=1.RDS")))
}


## process this "out" file
## it is a list, where each element of the list is a single replicate run of the parameter set
for (i in 1:8) {
    out <- get(paste0("out",i))
    lapply(1:length(out),
       function(i)
           data.frame(time=seq(1:length(out[[i]]$traits))/10-0.1,
                      C = lapply(out[[i]]$states, function(x) x[1]) %>% unlist,
                      P = lapply(out[[i]]$states, function(x) x[2]) %>% unlist,
                      med_b = lapply(out[[i]]$traits, median) %>% unlist,
                      sd_b = lapply(out[[i]]$traits, sd) %>% unlist,
                      rep=i)
       ) %>% do.call(rbind, .) -> out

    png(file=paste0("GEM_predprey_simulations_figure_",i,".png"), height=5, width=5, units='in', res=300)
    par(mfrow=c(2,2), mar=c(2,4,0.5,0.5), oma=rep(0.5,4))
    plot.new()
    plot.window(xlim=c(0, 200), ylim=c(0,500))
    axis(1); axis(2); box('plot')
    for (i in 1:max(out$rep))
        with(subset(out, rep==i), lines(time,C,col=gray(0.7),lwd=0.75))
    points(seq(0,200,0.1), sapply(seq(0,200,0.1), function(t) subset(out, time==t)$C %>% mean), col=4, lwd=2, pch=21, cex=0.25)
    mtext(side=2, line=2.5, "Prey abundance")

    plot.new()
    plot.window(xlim=c(0, 200), ylim=c(0,500))
    axis(1); axis(2); box('plot')
    for (i in 1:max(out$rep))
        with(subset(out, rep==i), lines(time,P,col=gray(0.7),lwd=0.75))
    points(seq(0,200,0.1), sapply(seq(0,200,0.1), function(t) subset(out, time==t)$P %>% mean), col=4, lwd=2, pch=21, cex=0.25)
    mtext(side=2, line=2.5, "Predator abundance")

    plot.new()
    plot.window(xlim=c(0, 200), ylim=c(1.8,4))
    axis(1); axis(2); box('plot')
    for (i in 1:max(out$rep))
        with(subset(out, rep==i), lines(time,med_b,col=gray(0.7),lwd=0.75))
    points(seq(0,200,0.1), sapply(seq(0,200,0.1), function(t) subset(out, time==t)$med_b%>% mean), col=4, lwd=2, pch=21, cex=0.25)
    mtext(side=2, line=2.5, "Median birth rate")

    plot.new()
    plot.window(xlim=c(0, 200), ylim=c(0,1))
    axis(1); axis(2); box('plot')
    for (i in 1:max(out$rep))
        with(subset(out, rep==i), lines(time,sd_b,col=gray(0.7),lwd=0.75))
    points(seq(0,200,0.1), sapply(seq(0,200,0.1), function(t) subset(out, time==t)$sd_b%>% mean), col=4, lwd=2, pch=21, cex=0.25)
    mtext(side=2, line=2.5, "Trait SD")
    dev.off()

    ## average over the last 20 time steps
    subset(out, time > 180) %>% apply(., 2, mean)

}




    ## one more processing step to make this very easy to plot using ggplot
    out %>% gather(., key="variable", value="value", 2:5) -> outgg

    ggplot(outgg, aes(x=time, y=value, group=rep)) +
        facet_wrap(~variable, scales="free_y") +
        geom_line() +
            theme_bw()
