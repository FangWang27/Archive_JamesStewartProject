# This file contains R code that generate and evaluate samples for a sequence of samples
load("SiteInfo.RData")
source("Function_get_cluster_lengths.R")
source("Function_perfect_match_trials.R")
source("Function_perfect_match_test_evaluator.R")
sites.selected <- raw.sites[1:5e3, ]
cluster.lengths <- get_cluster_lengths(sites.selected)

set.set.seed(123)

seq.mu <- seq(0, 4, 0.25)

lst.tabs.mu <- list()
lst.tabs.delta <- list()

for (i in seq_along(seq.mu)) {
    mu.1 <- seq.mu[i]
    lst.tabs.mu[[i]] <- perfect_match_trials(
        mu0 = 0,
        mu1 = mu.1,
        num.real.dmr = 10,
        delta = 1,
        num.pairs = 100,
        sample.N = 20,
        B = 250,
        cluster.lengths = cluster.lengths
    )
    save(lst.tabs.mu, file = "Trial_4_mu_tabs.RData")
}

seq.delta.below <- seq(0.1, 0.9, 0.1)
seq.delta.above <- seq(1, 3.5, 0.25)

seq.delta <- c(seq.delta.below, seq.delta.above)

for (i in seq_along(seq.delta)) {
    delta <- seq.delta[i]
    lst.tabs.delta[[i]] <- perfect_match_trials(
        mu0 = 0,
        mu1 = 3,
        num.real.dmr = 10,
        delta = delta,
        num.pairs = 100,
        sample.N = 20,
        B = 250,
        cluster.lengths = cluster.lengths
    )
    save(lst.tabs.delta, file = "Trial_4_delta_tabs.RData")
}


#analysis the tabs generated
lst.of.tab.combined.mu <- lapply(lst.tabs.mu, function(lst.tab) perfect_match_test_evaluator(lst.tab,0.05))
df.eval.mu.combined <- do.call(rbind,lst.of.tab.combined.mu)

lst.of.tab.combined.delta <- lapply(lst.tabs.delta, function(lst.tab) perfect_match_test_evaluator(lst.tab,0.05))
df.eval.delta.combined <- do.call(rbind,lst.of.tab.combined.delta)

