# this file contains the file that construct SimulationTrial


setClass(
  Class = "SimulationTrial",
  representation(
    site.info = "data.frame", # number of clusters and index of clusters
    num.observed = "numeric", # number of pairs(of patients) generate in the each run
    null.parameter = "numeric", # parameters used for null distribution
    alt.parameter = "numeric", # parameters used for alternative parameter
    global.parameter = "numeric",
    SEED = "numeric", # seed used for sample for the real index
    real.index = "numeric", # the parameters that remain the same for all groups
    dmr.num = "numeric" # real dmr index
  ),
  prototype(
    site.info = data.frame(length = numeric(), index = numeric()),
    num.observed = numeric(),
    null.parameter = numeric(),
    alt.parameter = numeric(),
    global.parameter = c(sigma = 0.3, rho = 0.5, a = 1, b = 1),
    SEED = numeric(),
    real.index = numeric(),
    dmr.num = 10
  )
)

setClass(
  Class = "PerfectMatchSimulationTrial",
  prototype = prototype(
    null.parameter = c(mu0 = 0, delta0 = 1),
    alt.parameter = c(mu1 = 1, delta1 = 1),
    global.parameter = c(sigma = 0.3, rho = 0.5, a = 1, b = 1),
    SEED = numeric(),
    real.index = numeric()
  ), contains = "SimulationTrial"
)



setGeneric(
  name = "ST.generate.sample", # method that generate single paired sampled using parameters specified. The sample are returned as a list of
  def = function(object) {
    standardGeneric("ST.generate.sample")
  }
)
setMethod(
  f = "ST.generate.sample", signature = "PerfectMatchSimulationTrial",
  definition = function(object) {
    site.info <- object@site.info
    null.parameter <- object@null.parameter
    alt.parameter <- object@alt.parameter
    global.parameter <- object@global.parameter
    rho <- global.parameter["rho"]
    sigma <- global.parameter["sigma"]
    cl.length <- site.info$length
    cl.index <- site.info$index
    mu0 <- null.parameter["mu0"]
    delta0 <- null.parameter["delta0"]
    mu1 <- alt.parameter["mu1"]
    delta1 <- alt.parameter["delta1"]

    sigma_generator <- function(l, rho, sigma) {
      # calculate the squared variance-covariance matrix
      # l : dimension of the matrix
      # sigma,rho : parameter specified
      mat_rho <- rho * diag(l)
      mat_row <- row(mat_rho)
      mat_col <- col(mat_rho)
      mat_Sigma <- sigma * rho^(abs(mat_row - mat_col))
      return(mat_Sigma)
    }

    condi_MVN_generator <- function(lst.para, Y = 0) {
      # takes a list of parameters for the model and returns the list of simulated sample for X|Y=0 or X|Y =1
      # parameter:
      # lst.para : a list of named list of parameters, each sublist is of the form list(mu, Delta, Sigma)
      # Y : value for Y, can only be 0 or 1,default set Y = 0
      lst.mvn <- lapply(lst.para, function(lst) {
        if (!(Y == 1 | Y == 0)) {
          stop(" Y can only be 0 or 1")
        }
        mu <- lst$mu
        Delta <- (lst$Delta)
        Sigma <- lst$Sigma
        if (Y == 1) {
          Sigma <- t(Delta) %*% Sigma %*% Delta
        }
        return(MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma))
      })
      return(unlist(lst.mvn))
    }

    sample.index <- object@real.index

    # list contains Parameter for H0
    lst.null.par <- lst.simulated.par <- lapply(cl.length, function(l) {
      # function that takes a vector of lengths of each sites and generate list contains
      # mu, Delta and Sigma for each
      list(mu = rep(mu0, l), Delta = diag(rep(delta0, l)), Sigma = sigma_generator(l = l, rho = rho, sigma = sigma))
    })
    # list contains parameters for H1
    lst.simulated.par[sample.index] <- lapply(sample.index, function(i) {
      unlist(lapply(cl.length[i], function(l) {
        Delta <- diag(rep(delta1, l))
        list(mu = rep(mu1, l), Delta = Delta, Sigma = sigma_generator(l = l, rho = rho, sigma = sigma))
      }), recursive = FALSE)
    })

    z <- rbeta(n = sum(cl.length), shape1 = global.parameter["a"], shape2 = global.parameter["b"])
    list.generated <- list(
      H0 = list(
        X0 = sqrt(z) * condi_MVN_generator(lst.null.par, Y = 0),
        X1 = sqrt(z) * condi_MVN_generator(lst.null.par, Y = 1)
      ),

      H1 = list(
        X0 = sqrt(z) * condi_MVN_generator(lst.null.par, Y = 0),
        X1 = sqrt(z) * condi_MVN_generator(lst.simulated.par, Y = 1)
      )
    )

    list.generated$DMRIndex <- sample.index

    return(list.generated)
  }
)


setGeneric(
  name = "ST.generate.tab",
  def = function(object, ...) {
    standardGeneric("ST.generate.tab")
  }
)
setMethod(
  f = "ST.generate.tab", signature = "PerfectMatchSimulationTrial",

  definition = function(object, B = 250) {
    site.info <- object@site.info
    cluster.rle <- list(lengths = site.info$length, values = site.info$index)
    cluster <- inverse.rle(cluster.rle)
    sites <- rep(1, length(cluster))
    patient.num <- object@num.observed
    df.sample.H0 <- data.frame(matrix(rep(0, patient.num * 2 * length(sites)), ncol = patient.num * 2))
    colnames(df.sample.H0) <- paste("sample", c(1:patient.num, 1:patient.num), c(rep("normal", patient.num), rep("cancer", patient.num)), sep = "")
    df.sample.H1 <- df.sample.H0
    for (i in 1:patient.num) {
      sample.list <- ST.generate.sample(object)
      H0 <- sample.list$H0
      H1 <- sample.list$H1
      df.sample.H0[, c(i, i + patient.num)] <- as.data.frame(H0[c("X0", "X1")])
      df.sample.H1[, c(i, i + patient.num)] <- as.data.frame(H1[c("X0", "X1")])
    }


    dmrIdentifier <- function(methyl.MVal, cluster, sites, patient.num = 30, B = 250) {
      # The function obtained by modifying author's code
      require(bumphunter)
      methyl.MVal <- data.matrix(methyl.MVal)
      num.observeds <- patient.num
      case.ind <- 1:patient.num
      ctrl.ind <- (patient.num + 1):(2 * patient.num)

      ##### Define function to calculate p-value from one-sided Pitman-Morgan test
      pitman.var.grt <- function(x, y, ratio = 1) {
        xok <- yok <- complete.cases(x, y)
        x <- x[xok]
        y <- y[yok]
        n <- length(x)
        df <- n - 2
        r <- cor(x, y)
        Var1 <- var(x)
        Var2 <- var(y)
        w <- var(x) / var(y)
        stat.t <- ((w - ratio) * sqrt(n - 2)) / sqrt(4 * (1 - r^2) * w * ratio)
        pval <- pt(-stat.t, df = df)
        return(list(statistic = stat.t, p.value = pval))
      }


      ##### Calculate combined signal scores for observed data

      rawStat <- sapply(1:nrow(methyl.MVal), function(i) {
        res.ttest <- t.test(methyl.MVal[i, case.ind], methyl.MVal[i, ctrl.ind], paired = TRUE)
        t.stat <- res.ttest$statistic

        pval.mean <- res.ttest$p.value
        pval.var <- pitman.var.grt(methyl.MVal[i, case.ind], methyl.MVal[i, ctrl.ind])$p.value
        if (pval.mean < 0.5 || pval.var < 0.5) { # so at least one of pval.mean, pval.var is less that 0.5
          return(c(pval.mean, pval.var, NA, NA, NA, NA))
        } else {
          Qm <- qnorm(pval.mean, lower.tail = FALSE)
          Qv <- qnorm(pval.var, lower.tail = FALSE)
          Qm <- Qm[Qm >= 0]
          Qv <- 0
          ratio <- Qv / (Qm + Qv)
          return(c(pval.mean, pval.var, Qm, Qv, ratio, t.stat))
        }
      })
      lambda <- mean(rawStat[5, ], na.rm = T)
      rawStat.wt <- lambda * rawStat[3, ] + (1 - lambda) * rawStat[4, ]
      rawStat.wt.sign <- rawStat.wt * sign(rawStat[6, ])

      stat <- bumphunter::runmedByCluster(y = rawStat.wt.sign, x = NULL, cluster = cluster, weights = NULL, k = 5, endrule = "constant", verbose = FALSE)
      Index <- which(stat$smoothed)
      stat <- stat$fitted


      ##### Calculate combined signal scores for 250 permutation datasets
      set.seed(123456)
      perm.sign <- matrix(sample(c(-1, 1), B * num.observeds, replace = TRUE), ncol = B)

      permStat <- sapply(1:B, function(i) {
        shuff.ind <- which(perm.sign[, i] == -1)

        perm.case.ind <- case.ind
        perm.case.ind[shuff.ind] <- case.ind[shuff.ind] + num.observeds

        perm.ctrl.ind <- ctrl.ind
        perm.ctrl.ind[shuff.ind] <- ctrl.ind[shuff.ind] - num.observeds
        permRawStat <- sapply(1:nrow(methyl.MVal), function(i) {
          res.ttest <- t.test(methyl.MVal[i, perm.case.ind], methyl.MVal[i, perm.ctrl.ind], paired = TRUE)
          t.stat <- res.ttest$statistic

          pval.mean <- res.ttest$p.value
          pval.var <- pitman.var.grt(methyl.MVal[i, perm.case.ind], methyl.MVal[i, perm.ctrl.ind])$p.value
          if (pval.mean < 0.5 || pval.var < 0.5) { # so at least one of pval.mean, pval.var is less that 0.5
            return(c(pval.mean, pval.var, NA, NA, NA, NA))
          } else {
            Qm <- qnorm(pval.mean, lower.tail = FALSE)
            Qv <- qnorm(pval.var, lower.tail = FALSE)

            Qm[Qm < 0] <- 0
            Qv[Qv < 0] <- 0
            ratio <- Qv / (Qm + Qv)
            return(c(pval.mean, pval.var, Qm, Qv, ratio, t.stat))
          }
        })

        perm.lambda <- mean(permRawStat[5, ], na.rm = T)
        permRawStat.wt <- perm.lambda * permRawStat[3, ] + (1 - perm.lambda) * permRawStat[4, ]
        permRawStat.wt.sign <- permRawStat.wt * sign(permRawStat[6, ])

        permStat <- bumphunter::runmedByCluster(y = permRawStat.wt.sign, x = NULL, cluster = cluster, weights = NULL, k = 5, endrule = "constant", verbose = FALSE)
        permStat <- permStat$fitted
        return(permStat)
      })


      ##### Get threshold as the 99th percentile of genome-wide smoothed signal scores
      cutoff <- quantile(abs(permStat), 0.99, na.rm = T)

      # browser()
      ##### Identify candidate DMRs
      pos.gen <- 1:length(stat) # fake positions
      tab <- bumphunter::regionFinder(x = stat, chr = sites, pos = pos.gen, cluster = cluster, cutoff = cutoff, ind = Index, assumeSorted = TRUE, verbose = FALSE)


      ##### Calculate p-values for candidate DMRs
      if (nrow(tab) > 0) {
        nulltabs <- lapply(1:B, function(i) {
          bumphunter::regionFinder(x = permStat[, i], chr = sites, pos = pos.gen, cluster = cluster, cutoff = cutoff, ind = Index, assumeSorted = TRUE, verbose = FALSE)
        })

        A <- as.list(rep(0, B))
        for (i in 1:B) {
          nulltab <- nulltabs[[i]]
          if (nrow(nulltab) > 0) {
            A[[i]] <- nulltab$area
          }
        }

        Avalue <- matrix(tab$area, ncol = 1)
        tots <- sapply(1:nrow(Avalue), function(i) {
          return(sapply(seq(along = A), function(j) {
            sum((A[[j]] >= Avalue[i]) | (abs(A[[j]] - Avalue[i]) <= sqrt(.Machine$double.eps)))
          }))
        })
        if (is.vector(tots)) {
          tots <- matrix(tots, nrow = 1)
        }
        tab$p.valueArea <- colSums(tots) / sum(sapply(nulltabs, nrow)) # nrow(nulltabs) corresponding to how many regions generate at each iteration, and summing up corresponding to the denominator.
        tab$fwerArea <- colMeans(tots > 0) # the family wise error rate is given by
        tab <- tab[order(tab$fwer, -tab$area), ]
      }
      return(tab)
    }
    tab.H0 <- dmrIdentifier(df.sample.H0, cluster, sites, patient.num, B)
    tab.H1 <- dmrIdentifier(df.sample.H1, cluster, sites, patient.num, B)
    tab <- lapply(list(df.sample.H0, df.sample.H1), function(df.sample) {
      dmrIdentifier(df.sample, cluster, sites, patient.num, B)
    })
    tab[[3]] <- B
    names(tab) <- c("H0", "H1", "permutation")
    return(tab)
  }
)


setGeneric(
  name = "ST.set.real.index",
  def = function(object) {
    standardGeneric("ST.set.real.index")
  }
)
setMethod(
  f = "ST.set.real.index", signature = "SimulationTrial",
  definition = function(object) {
    if (length(object@SEED) != 0) {
      set.seed(object@SEED)
    }
    site.info <- object@site.info
    # browser()
    object@real.index <- sample(site.info$index, object@dmr.num)
    return(object)
  }
)



ST.get.site.info <- function(sites, maxGap = 1000, min.size = 3) {
  # takes a sites and maxGap gives characteristic of clusters
  # define with given maxGap and clusters smaller than min.size is discarded
  require("bumphunter")
  cl <- bumphunter::clusterMaker(
    chr = sites$Chromosome_36,
    pos = sites$MAPINFO,
    maxGap = maxGap
  )
  cl <- sort(cl)
  cl.lst <- rle(cl)
  lengths <- cl.lst$lengths[cl.lst$lengths >= min.size]
  index <- (1:length(lengths)) + 1
  # browser()
  return(data.frame(length = lengths, index = index))
}

setGeneric(
  name = "ST.trial.evaluation",
  def = function(object, tab, p.val = 0.05) { # function that calculate the false postie rate and true negative rate under a given p-value
    standardGeneric("ST.trial.evaluation")
  }
)
setMethod(
  f = "ST.trial.evaluation", signature = "PerfectMatchSimulationTrial",
  definition = function(object, tab, p.val = 0.05) {
    real.index <- object@real.index
    dmr.num <- object@dmr.num
    site.info <- object@site.info
    tab1 <- tab$H1
    tab1 <- tab1[tab1$p.valueArea < p.val, ]
    tab0 <- tab$H0
    tab0 <- tab0[tab0$p.valueArea < p.val, ]
    cl.identified <- tab1$cluster
    cl.total <- length(site.info$index)
    # get the true positive rate (TPR)
    if (nrow(tab1) == 0) {
      TPR <- 0
    } else {
      TPR <- sum(cl.identified %in% real.index) / dmr.num
    }
    # get the false positive rate (FPR)
    FPR <- (length(cl.identified) - sum(cl.identified %in% real.index)) / (cl.total - dmr.num)

    evaluation <- data.frame(
      "p-value" = p.val,
      "permutation" = tab$permutation,
      "TPR" = TPR,
      "FPR" = FPR,
      "FP.H0" = nrow(tab0),
      t(object@null.parameter), t(object@alt.parameter), t(object@global.parameter)
    )
    return(evaluation)
  }
)