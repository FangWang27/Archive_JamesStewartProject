# code for generating case-control sample, where the sub-population structure is present

case.p1 <- 0.5 # proportion of sub-population 1 in the case sample
control.p1 <- 0.5 # proportion of sub-population1 in the control sample

mu.control.1 <- 0 # parameter mu of sub-population 1 for control samples
mu.case.1 <- 1 # parameter mu of sub-population 1 for the case samples

mu.control.2 <- 0 # parameter mu of sub-population 2 for control samples
mu.case.2 <- 1  # parameter mu of sub-population 2 for the case samples

delta.1 <- 1 # parameter delta for sub-population 1

delta.2 <- 1 # parameter delta for sub-population 2

rho <- 0.5; sigma <-0.25

n.case <- 30 # sample size for case

n.control <- 30 # sample size for control

sample.N <- 20 # how many group of sample should be generated

num.real.dmr <- 10


real.index <-sample(1: length(cluster.lengths), num.real.dmr)
# generating clusters from cluster_lengths an
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

# using the information about sizes of various clusters to gnenerate a block diggonal matrix

# the list all sub-matrixes
lst.Sigma <- lapply(cluster.lengths, function(length) {
  sub.matrix.of.Sigma <-  dict.Sigma[[length]]
  return(sub.matrix.of.Sigma)})
mat.Sigma <- Matrix::bdiag(lst.Sigma) # combine the list of sigma matrixes into a block diagonal matrix

# generating mvn samples of mean 0 and covariance matrix given by the combined Sigma matrix
mvsample <- MASS::mvrnorm(
  n = sample.N * (n.control + n.case),
  mu = rep(0, num.sites),
  Sigma = mat.Sigma
)


# generating mu and delta that used for linear transformation for case

n.case.1 <- round(n.case * case.p1) # how many samples of case are from sub-population 1
n.case.2 <- n.case - n.case.1 # how many case samples are from sub-population 2
region.dmr <- clusters %in% real.index # a logical vector indicating whether a site is in a dmr cluster

mu.single.case.1 <- rep(mu.control.1, num.sites)
mu.single.case.1[region.dmr] <- mu.case.1
mu.single.case.2 <- rep(mu.control.2, num.sites)
mu.single.case.2[region.dmr] <- mu.case.2

delta.single.case.1 <- rep(1, num.sites)
delta.single.case.1[region.dmr] <- delta.1
delta.single.case.2 <- rep(1, num.sites)
delta.single.case.2[region.dmr] <- delta.2

delta.case <- c( rep(delta.single.case.1, n.case.1), rep(delta.single.case.2, n.case.2))
mu.case <- c( rep(mu.single.case.1, n.case.1), rep(mu.single.case.2, n.case.2))


sample.case <- t(mvsample[1:(sample.N * n.case),])

sample.case <- mu.case + sample.case * delta.case

# generating sample for control
n.control.1 <- round(n.control * control.p1) # how many samples of case are from sub-population 1
n.control.2 <- n.control - n.control.1 # how many case samples are from sub-population 2

mu.single.control.1 <- rep( mu.control.1, num.sites)
mu.single.control.2 <- rep( mu.control.2, num.sites)
mu.control <- c( rep(mu.single.case.1, n.control.1), rep(mu.single.case.2, n.control.2) )

sample.control <- t(mvsample[ -(1:(sample.N * n.case)),])
sample.control <- mu.control + sample.control

