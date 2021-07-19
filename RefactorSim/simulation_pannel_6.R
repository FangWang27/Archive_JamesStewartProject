#data created : July 4th
# simulating samples with different in mean

load("SiteInfo.RData")
source("Function_get_cluster_lengths_July4th.R")
source("Function_dmrIdentifier_opt.R")
sites.selected <- raw.sites[1:5e3, ]
cluster.lengths <- get_cluster_lengths(sites.selected)

set.seed(123)

mu1.a <- 1.5

mu1.b <- 3

lst.tab.a <- perfect_match_trials_july4th(mu0 = 0, mu1 = mu1.a, sample.N = 100, B =  1, cluster.lengths = cluster.lengths)
save(lst.tab.a, file = "lst_tabs_July4th")
lst.tab.b <- perfect_match_trials_july4th(mu0 = 0, mu1 = mu1.b, sample.N = 100, B = 1, cluster.lengths = cluster.lengths)
save(lst.tab.a,lst.tab.b, file = "lst_tabs_July4th")


