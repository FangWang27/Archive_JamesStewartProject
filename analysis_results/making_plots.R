
library(ggplot2)



# load Wang's sample and run their code to produce figure-1.
map.info <- sample.data[1:300,3]

#figure-1
#ellustrative example for DMRs in the sample data set provided by Wang's paper

#run Wang's code to obtain necessary data
# code to produce the illustrative example of DMRs in the real data
start <- tab$start[tab$fwerArea < 0.05& tab$end < 1.5e6 &tab$start > 0 ]
end <- tab$end[tab$fwerArea < 0.05& tab$end < 1.5e6  &tab$start > 0 ]


df.dmr.region <- data.frame(start,end,group = seq_along(start))

mean.cancer <- rowMeans(sample.data[,c(34:63)])
mean.normal <- rowMeans(sample.data[,c(4:33)])
df.illustration <- data.frame(map.info,mean.normal[1:300],mean.cancer[1:300])


names(df.illustration) <- c("MapInfo", "normal", "cancer")

ggplot(df.illustration, aes(x = map.info, color = sample.type)) +
  geom_rect(df.dmr.region, inherit.aes = FALSE, 
            mapping = aes(xmin = start, xmax = end, ymin = -Inf, 
                          ymax = Inf,
                          group = group, fill= "DMRs identifeid"),
            color="transparent",alpha = 0.3) +
  xlab("genetic distance") + ylab("DNA methylation level") +
  geom_point( aes(y = cancer, col = "cancer")) +
  geom_point( aes(y = normal, col = "normal")) +
  scale_fill_manual("",values = "green")

#clear the workspace
load("SiteInfo.RData")
source("Function_get_cluster_lengths.R")
source("Function_perfect_match_trials_july5th.R")
sites.selected <- raw.sites[1:1e4, ]
cluster.lengths <- get_cluster_lengths(sites.selected,min.size =5, max.size = 15)
clusters <- rep(1:length(cluster.lengths), times = cluster.lengths)


#figure-2--------
#cluster_lengths
 # a) the distribution of genetic distanced
 df.genetic.distance <- data.frame(gene.distance <- sites.selected$MAPINFO)
 
plot2.a <- 
  ggplot(df.genetic.distance, aes(gene.distance)) +
   geom_histogram(bins=100,color = "white") +
   xlab("CpG sites position") +
   ylab("Frequency")

 
# b1) the disribution of the length of clusters of sites selected
cluster.lengths.raw <- 
  data.frame(cluster.lengths = get_cluster_lengths(sites.selected, maxGap = 300,min.size = -Inf,max.size = Inf))

plot2.b1<- 
  ggplot(cluster.lengths.raw, aes(cluster.lengths)) +
  geom_bar(fill = "darkblue") +
 xlim(c(0,26)) +
 xlab("number of CpG sites in each cluster") +
 ylab("Frequency")

# b2) the distibution of the length of the clusters after removing clusters with size less than 5 and greater than 15
  
cluster.lengths <- 
  data.frame(cluster.lengths = get_cluster_lengths(sites.selected, maxGap = 300,min.size = 5,max.size = 15))

plot2.b2 <-
ggplot(cluster.lengths, aes(cluster.lengths)) +
  geom_bar(fill = "darkblue") +
  xlim(c(4,16)) +
  ylim(c(0,150)) +
  geom_text(stat = 'count',aes(label = ..count..), vjust = -1) +
  xlab("number of CpG sites in each cluster") +
  ylab("Frequency")

# b3) the distibution of the length of the clusters after removing clusters with size less than 5 and greater than 15
  
cluster.lengths <- 
  data.frame(cluster.lengths = get_cluster_lengths(sites.selected, maxGap = 300,min.size = 3,max.size = 15))
plo2.b3 <-
ggplot(cluster.lengths, aes(cluster.lengths)) +
  geom_bar(fill = "darkblue") +
  xlim(c(2,16)) +
  ylim(c(0,400)) +
  geom_text(stat = 'count',aes(label = ..count..), vjust = -1) +
  xlab("number of CpG sites in each cluster") +
  ylab("Frequency")


figure2 <- cowplot::plot_grid(plot2.a,plot2.b1,plot2.b2,plo2.b3,
                   labels = c("A","B1","B2","B3"),
                   nrow = 2,ncol = 2
                   ,scale = 1,vjust = 2.5,hjust = 0)

figure2

#figure-3 --------------
# distribution of sample generated under null hypothesis,
# load the data from null_tests_july7th.RData
seq.bin <- seq(-0.015,0.015,along.with = null.samples[,1])
upper <- seq.bin[-1]
lower <- seq.bin[-length(seq.bin)]
theoretical.norm <- (pnorm(upper, 0, 1/sqrt(6000)) - pnorm(lower, 0, 1/sqrt(6000)))

null.samples <- lst.null.test[[3]]
df.null.dis <- data.frame(
  combined.mean.of.sites = rowMeans(null.samples),
  first.pair.1 = null.samples[,1],
  first.pair.2 = null.samples[,2],
  theoretical.norm = c(0,theoretical.norm)
)

#a)characteristic of samples generated

#a) # histgram of mean of Mval of each sites
plot3.a <- ggplot(df.null.dis, aes(x = combined.mean.of.sites)) +
  geom_histogram(aes(y = ..count../sum(..count..)), bins  = 100,color = "white") +
  geom_line(aes(x = seq.bin,y = theoretical.norm)) +
  xlab("mean methylation level (M value) across each CpG sites") +
  ylab("Density")

# b) plot of mean Mval of each sitess

plot3.b <- ggplot(df.null.dis, aes(x=seq_along(combined.mean.of.sites))) +
  geom_point( aes(y = combined.mean.of.sites,color = combined.mean.of.sites ), show.legend = FALSE) +
  xlab("index of CpG sites") +
  ylab("methylation level (M value)")

#c) plot of two paired null samples

plot3.c <- ggplot(df.null.dis, aes(x=seq_along(combined.mean.of.sites))) +
  geom_point( aes(y = first.pair.1,color = "sample1"),alpha = 0.3) +
  geom_point( aes(y = first.pair.2,color = "paired sample of sample1" ),alpha = 0.3) +
  xlab("index of CpG sites") +
  ylab("methylation level (M value)") +
  theme(legend.title=element_blank())

#d) correlation between paried sample
plot3.d <- ggplot(plot_null_samples, aes(first.pair.1,first.pair.2)) +
  geom_jitter( aes(x = first.pair.1, y = first.pair.2),alpha = 0.3) +
  xlab("methylation level of sample 1") +
  ylab("methylation level of paired sample 1")


figure3 <- ggpubr::ggarrange( ggpubr::ggarrange(plot3.a, plot3.b, ncol = 2,labels = c("A","B")),
                              ggpubr::ggarrange(plot3.c, plot3.d, labels = c("C1","C2"))
                              ,nrow = 2)
figure3


#4b) characteristic of p-values and fwer

#4b1) histogram of fwer

plot4.a <- ggplot(df.null.tab, aes(fwerArea)) +
  geom_histogram(aes(y = ..count../sum(..count..)),bins = 200,color = "white") +
  xlab("family wise error rate") +
  ylab("Proportion")+
  geom_vline(aes(xintercept= c(0.05)),
             colour="red")
  
#b2) histogram of p-value

plot4.b <- ggplot(df.null.tab, aes(p.valueArea)) +
 geom_histogram(aes(y = ..count../sum(..count..)),bins = 200,color = "white") +
    geom_vline(aes(xintercept= c(0.05)),
               colour="red")+
  xlab("p value") +
  ylab("Proportion")
  
figure4 <- ggpubr::ggarrange(plot4.a,plot4.b, nrow = 2,labels = c("A","B"))



# figure 5) the disribution under alternative hypothesis----
# mu1 = 1.5, delta = 1.5
# need to load lst_tabs_July5th.RDataa
# use the first trial for illustration
lst.alt.example <- lst.trials.a[[1]]
real.index <- lst.alt.example[[2]]
dmr.sites <- clusters %in% real.index
generated.samples <- lst.alt.example[[3]]
normal.example <- generated.samples[,31]
cancer.example <- generated.samples[,1]



#5a1)historam of site-level mean of normal sample and cancer sample ------
avg.normal.samples <- rowMeans(generated.samples[,1:30])
avg.cancer.samples <- rowMeans(generated.samples[,31:60])

df.mean.normal.samples <- 
  data.frame( methylation.level= avg.normal.samples, sample.type = "normal")
df.mean.cancer.samples <- 
  data.frame( methylation.level= avg.cancer.samples, sample.type = "cancer")

df.mean.paired.samples <- rbind(df.mean.normal.samples,df.mean.cancer.samples)

plot5.a1 <- 
  ggplot(df.mean.paired.samples, aes(methylation.level, fill = sample.type)) +
  geom_histogram(alpha = 0.5, binwidth = 0.05, position = 'identity') +
  xlab("methylation level (M value)") +
  ylab("Frequency")


# 5a2) dot plot of mean with lable of DMRs simulated--------

inds <- diff( c(0, dmr.sites))
dmr.start <- seq_along(inds)[inds == 1]
dmr.end <- seq_along(inds)[inds == - 1]
dmr.rect <- data.frame(start = dmr.start, end =dmr.end, group =  seq_along(dmr.start) )

df.mean.paired.samples <- data.frame(avg.normal.samples, avg.cancer.samples)
plot5.a2 <-
  ggplot(df.mean.paired.samples, aes(x = seq_along(avg.normal.samples))) +
  geom_point( aes(y = avg.normal.samples, col = "cancer"),alpha = 0.5) +
  geom_point( aes(y = avg.cancer.samples, col = "normal"), alpha = 0.5) +
  geom_rect(dmr.rect, inherit.aes=FALSE,
            mapping = aes(xmin = dmr.start, xmax = dmr.end, ymin = -Inf, ymax = Inf, 
                          group = group, fill= "DMRs simulated"),
            color="transparent", alpha=0.3) +
  xlab("index of CpG sites") +
  ylab("mean methylation level") +
scale_fill_manual("",values = "green")+
theme(legend.title=element_blank())


# 5b) dot plot of the first paired sample--------
  df.single.paired.sample <- data.frame(normal.example, cancer.example)
plot5.b <- 
ggplot(df.single.paired.sample, aes(x = seq_along(normal.example))) +
  geom_point( aes(y = cancer.example, col = "cancer"),alpha = 0.5) +
  geom_point( aes(y = normal.example, col = "normal"), alpha = 0.5) +
  geom_rect(dmr.rect, inherit.aes=FALSE,
            mapping = aes(xmin = dmr.start, xmax = dmr.end, ymin = -Inf, ymax = Inf, 
                          group = group, fill= "DMRs simulated"),
            color="transparent", alpha=0.3) +
  xlab("index of CpG sites") +
  ylab("methylation level") +
  scale_fill_manual("",values = "green") +
theme(legend.title=element_blank())



#5c1) distribution of tpr,ppv and fdr ----


plot5.c1 <-
  ggplot(df.eval.a) +
  geom_bar(aes(x = ppv) ,fill = "darkblue") +
  geom_text(stat = 'count',aes(ppv,label = ..count..), vjust = -.5) +
  ylim(0,100) +
  xlab("Positive predictive value") +
  ylab("Frequency ")

plot5.c2 <-
  ggplot(df.eval.a) +
  geom_bar(aes(x = tpr) ,fill = "darkblue") +
  geom_text(stat = 'count',aes(tpr,label = ..count..), vjust = -.5) +
  ylim(0,110) +
  xlab("True positive rate") +
  ylab("Frequency")


#simulation with difference in variance only

figure5 <-  ggpubr::ggarrange(
  ggpubr::ggarrange(plot5.a1,plot5.a2,labels = c("A1","A2"),ncol = 2),
  ggpubr::ggarrange(plot5.b,labels = c("B")),
  ggpubr::ggarrange(plot5.c1,plot5.c2,labels = c("C1","C2")),
  nrow = 3
)


#anaysis trial.a-----




lst.eval.a <- lapply(lst.trials.a, function(lst.of.multiple.trials.of.same.dmr) {
  lst.tabs <- lst.of.multiple.trials.of.same.dmr[[1]]
  real.index <- lst.of.multiple.trials.of.same.dmr[[2]]
  num.real.dmr <- length(real.index)
  num.clusters <- length(cluster.lengths)
  
  lst.eval <- lapply(lst.tabs, function(df.single.run) {
    
    df.identified <- df.single.run[df.single.run$p.valueArea <0.05 & (df.single.run$L > 0.5 *df.single.run$clusterL), ]
    num.true.positive <- sum(df.identified$cluster %in% real.index)
    num.false.positive <- nrow(df.identified) - num.true.positive
    tpr <- num.true.positive / num.real.dmr
    ppv <- num.true.positive / nrow(df.identified)
    return( data.frame(num.true.positive, num.false.positive, tpr,ppv))
  })
  df.eval.stats <- do.call(rbind, lst.eval)
  return(df.eval.stats) 
})

df.eval.a <- do.call(rbind, lst.eval.a)
