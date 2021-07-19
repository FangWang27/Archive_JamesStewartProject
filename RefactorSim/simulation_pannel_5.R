load("SiteInfo.RData")
source("Function_get_cluster_lengths.R")
source("Function_dmrIdentifier_opt.R")
sites.selected <- raw.sites[1:5e3, ]
cluster.lengths <- get_cluster_lengths(sites.selected)

set.seed(123)

mu0 <- 0
sigma <- 0.25; rho <- 0.5
a <- 1; b <- 1;
num.pairs <- 30; sample.N <- 100; B <- 250

clusters <- rep(1:length(cluster.lengths), times = cluster.lengths) 

num.sites <- sum(cluster.lengths)

sigma_generator <- function(l, rho, sigma) {
  # function that calculate the squared variance-covariance matrix
  # l : dimension of the matrix
  # sigma,rho : parameter specified
  mat_rho <- rho * diag(l)
  mat_row <- row(mat_rho)
  mat_col <- col(mat_rho)
  mat_Sigma <- sigma * rho^(abs(mat_row - mat_col))
  return(mat_Sigma)
}

# create a dictionary to store Sigmas of all possible size. dict.Sigma[[l]] wil return a Sigma matrix of size l
dict.Sigma <- list()
for (length in 1:max(cluster.lengths)) {
  dict.Sigma[[length]] <- sigma_generator(length, rho, sigma)
}

# using the information about sizes of various clusters to gnenerate a block diagonal matrix

# the list all sub-matrixes
lst.Sigma <- lapply(cluster.lengths, function(length) {
  sub.matrix.of.Sigma <-  dict.Sigma[[length]]
  return(sub.matrix.of.Sigma)})
mat.Sigma <- Matrix::bdiag(lst.Sigma) # combine the list of sigma matrixes into a block diagonal matrix

# generating mvn samples of mean 0 and covariance matrix given by the combined Sigma matrix
mvsample <- MASS::mvrnorm(
  n = num.pairs * sample.N * 2,
  mu = rep(0, num.sites),
  Sigma = mat.Sigma
)

# generating variable z
mat.z <- (rbeta(n = num.pairs * num.sites, a, b))

null.samples <- (mu0 + t(mvsample)) * mat.z

null.tabs <- list()
for (index in 1:sample.N) {
  begin <- 1 + (index - 1) * num.pairs * 2 # samples for this trial begin
  end <- index * num.pairs * 2 # samples for this trial end
  null.sample <- null.samples[, begin:end]
  null.tabs[[index]] <- dmrIdentifier_opt(null.sample, clusters, num.pairs, B, k = 5, cluster.lengths)
}


#extract one sample for analysis

#type 1 error 

lst.type1.error <- lapply(null.tabs_july4, function(tab) {
  df.indentified <- tab[tab$fwerArea < 0.05 && tab$L == tab$clusterL,]
  is.type1.error <- ifelse( nrow(df.indentified)>0,TRUE,FALSE)
})

type1.error.val <- mean(unlist(lst.type1.error))


# test with t-test
pval.t.test <- numeric()
for (index in 1:sample.N) {
  begin <- 1 + (index - 1) * num.pairs * 2 # samples for this trial begin
  end <- index * num.pairs * 2 # samples for this trial end
  pval.t.test[index] <- t.test(x = null.samples[,begin: (begin + num.pairs-1)], y = null.samples[,(begin + num.pairs) : end],paired = TRUE)$p.value
}

sum(pval.t.test <0.05)

# test with p.value 


