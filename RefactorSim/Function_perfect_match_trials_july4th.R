# this file contains the function perfect_match_trials that takes necessary parameters then
# simulating samples and returns two list of data.frames and a numeric vector indicating the
#cluster number of the real DMR simulated.

source("Function_dmrIdentifier_opt.R")
perfect_match_trials_july4th <- function(
mu0 = 0, 
mu1 = 1,
num.real.dmr = 10,
delta = 1,
sigma = 0.25, rho = 0.5,
a = 1, b = 1,
num.pairs = 30, sample.N = 10, B = 250,
cluster.lengths = NULL
){


real.index <-sample(1: length(cluster.lengths), num.real.dmr)
library(tictoc)
# generating clusters from cluster_lengths an
clusters <- rep(1:length(cluster.lengths), times = cluster.lengths) 

tic("the total running time for simulating samples and finding DMRs")
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

# using the information about sizes of various clusters to gnenerate a block diggonal matrix

# the list all sub-matrixes
lst.Sigma <- lapply(cluster.lengths, function(length) {
  sub.matrix.of.Sigma <-  dict.Sigma[[length]]
  return(sub.matrix.of.Sigma)})
mat.Sigma <- Matrix::bdiag(lst.Sigma) # combine the list of sigma matrixes into a block diagonal matrix

# generating mvn samples of mean 0 and covariance matrix given by the combined Sigma matrix
mvsample.Ns <- MASS::mvrnorm(
  n = num.pairs * sample.N * 2,
  mu = rep(0, num.sites),
  Sigma = mat.Sigma
)

# generating variable z
mat.z <- (rbeta(n = num.pairs * num.sites, a, b))


# Since the matrix Delta is just a multiple a diagonal matrix, hence we may apply the corresponding linera transformation
# to the reading of normal sites to generate the sample for the DMR sites

# find the dmr regions across the sites

region.dmr <- clusters %in% real.index # a logical vector indicating whether a site is in a dmr cluster

# construct the vector used to perform the linear transformation on the DMR sites.

# construct the vector that control the spread of the samples 

# Here vec.delta1 is a vector that control the spread of a single sample consisting of all sites.
#By replicating it (num.pairs) times, we obtain a long vector such that controls the spread of a 
# matrix with (num.pairs) columns, (num.sites) rows. For the matrix being multiplied with, each column
# follows MVN(0, \Sigma) disribution. 
vec.delta0 <- rep(1, num.sites * num.pairs) # no changes for the normal tissues
vec.delta1 <- rep(1, num.sites)
vec.delta1[region.dmr] <- delta
vec.delta1 <- rep(vec.delta1, num.pairs)
vec.delta <- c(vec.delta0, vec.delta1)

# construct the vector that control the location of mean of the sample in the similar way

vec.mu0 <- rep(mu0, num.sites * num.pairs)
vec.mu1 <- rep(mu0, num.sites)
vec.mu1[region.dmr] <- mu1
vec.mu1 <- rep(vec.mu1, num.pairs)
vec.mu <- c(vec.mu0, vec.mu1)

# generating actual paired samples by performing the linear transformation on the multivaraitate normal samples
alt.samples <- (vec.mu + t(mvsample.Ns) * vec.delta) * mat.z

# applying Wang's algrithm
alt.tabs <- list()
sites <- rep(1, num.sites)
tic("time for identifying DMRs")
for (index in 1:sample.N) {
  begin <- 1 + (index - 1) * num.pairs * 2 # samples for this trial begin
  end <- index * num.pairs * 2 # samples for this trial end
  alt.sample <- alt.samples[, begin:end]
  alt.tabs[[index]] <- dmrIdentifier_opt(alt.sample, clusters, num.pairs, B, k = 5, cluster.lengths)
}
toc()
toc()
return( list( alt.tabs = alt.tabs, real.index = real.index))
}
