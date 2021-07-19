load("SiteInfo.RData")
source("Function_get_cluster_lengths.R")
source("Function_perfect_match_trials_july5th.R")
sites.selected <- raw.sites[1:1e4, ]
cluster.lengths <- get_cluster_lengths(sites.selected,min.size =5, max.size = 15)

set.seed(123)

mu1.a <- 1.5

mu1.b <- 3


lst.trials.a <- list()
for (i in 1:10){
  lst.trials.a[[i]] <- perfect_match_trials_july5th(mu0 = 0, mu1 = 1.5, sample.N = 10, B =  250, delta = 1.5, cluster.lengths = cluster.lengths)
}

lst.eval.a <- lapply(lst.trials.a, function(lst.of.multiple.trials.of.same.dmr) {
 lst.tabs <- lst.of.multiple.trials.of.same.dmr[[1]]
 real.index <- lst.of.multiple.trials.of.same.dmr[[2]]
 num.real.dmr <- length(real.index)
 num.clusters <- length(cluster.lengths)
 
 lst.eval <- lapply(lst.tabs, function(df.single.run) {
   
   df.identified <- df.single.run[df.single.run$fwerArea <0.05 & (df.single.run$L > 0.5 *df.single.run$clusterL), ]
   num.true.positive <- sum(df.identified$cluster %in% real.index)
   num.false.positive <- nrow(df.identified) - num.true.positive
   tpr <- num.true.positive / num.real.dmr
   return( data.frame(num.true.positive, num.false.positive, tpr))
 })
 df.eval.stats <- do.call(rbind, lst.eval)
 return(df.eval.stats) 
})

df.eval.a <- do.call(rbind, lst.eval.a)



save(lst.trials.a, df.eval.a, file = "lst_tabs_July5th.RData")

# test for FDR with 100 round of paired sample without DMR
set.seed(123)
lst.null.test <-  perfect_match_trials_july5th(mu0 = 0, mu1 = 0, sample.N = 100, B =  250, num.real.dmr = 0, cluster.lengths = cluster.lengths)
df.null.tab <- do.call(rbind,lst.null.test[[1]])
sum(df.null.tab$p.valueArea <0.05/250)
save(lst.null.test, df.null.tab,file = "null_tests_july7th.RData")

# test for the case where only difference in the variance is present
set.seed(123)
lst.trials.b.var <- list()
for (i in 1:10){
   lst.trials.b.var[[i]] <- perfect_match_trials_july5th(mu0 = 0, mu1 = 0, sample.N = 10, B =  250, delta = 1.5, cluster.lengths = cluster.lengths)
}
save(lst.trials.b.var, file = "variance_only_tests_july8th.RData")

#test the case with different cluster size distribution

set.seed(123)
cluster.lengths.3 <- get_cluster_lengths(sites.selected,min.size =3, max.size = 15)
lst.trials.a.3 <- list()
for (i in 1:10){
   lst.trials.a[[i]] <- perfect_match_trials_july5th(mu0 = 0, mu1 = 1.5, sample.N = 10, B =  250, delta = 1.5, cluster.lengths = cluster.lengths.3)
}


save(lst.trials.a.3, file = "simulation_a_3_july9th.RData")

