#anaysis trial.a-----
#load the file lst_tabs_July5th.RData

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

# evaluation of the test for null cases

#lst.null.test <-  perfect_match_trials_july5th(mu0 = 0, mu1 = 0, sample.N = 100, B =  250, num.real.dmr = 0, cluster.lengths = cluster.lengths)
#load the null_tests_july7th.RData

df.null.tab <- do.call(rbind,lst.null.test[[1]])
sum(df.null.tab$p.valueArea <0.05/250)

#evaluation of the tests where only variance difference present

lst.eval.variance.only <- lapply(lst.trials.b.var, function(lst.of.multiple.trials.of.same.dmr) {
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

df.variance.only <- do.call(rbind, lst.eval.variance.only)

colMeans(df.variance.only)

sum(df.variance.only$num.false.positive >0)

sample.var.only.trials <- lst.trials.b.var[[1]]
sample.var.only.data.generated <- sample.var.only.trials[[3]]
sample.var.only.normal <- sample.var.only.data.generated[,1]
summary(sample.var.only.normal)
sample.var.only.cancer <- sample.var.only.data.generated[,31]
summary(sample.var.only.cancer)
plot(sample.var.only.normal, col = "black",type = "l")
par(new = TRUE)
plot(sample.var.only.cancer,col = "red",type= "l")

real.index.var.only <- lst.trials.b.var[[3]]

inds <- clusters %in% real.index.var.only


# analysis the performance of the algorithm with different cluster length vary from 3 to 15

load("SiteInfo.RData")
source("Function_get_cluster_lengths.R")
source("Function_perfect_match_trials_july5th.R")
sites.selected <- raw.sites[1:1e4, ]
cluster.lengths <- get_cluster_lengths(sites.selected,min.size =3, max.size = 15)


lst.eval.a.3 <- lapply(lst.trials.a.3, function(lst.of.multiple.trials.of.same.dmr) {
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

df.eval.a.3 <- do.call(rbind, lst.eval.a.3)




