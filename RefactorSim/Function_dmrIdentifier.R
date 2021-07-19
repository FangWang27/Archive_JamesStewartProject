dmrIdentifier_opt <- function(methyl.MVal, cluster, sites, patient.num = 30, B = 250) {
  # The function obtained by modifying author's code
  # note that the modified function took the transpose of what Wang's original dataform
  require(bumphunter)
  methyl.MVal <- t(methyl.MVal)
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

    Qm <- qnorm(pval.mean, lower.tail = FALSE)
    Qv <- qnorm(pval.var, lower.tail = FALSE)
    Qm[Qm < 0] <- 0
    Qv[Qv < 0] <- 0
    ratio <- Qv / (Qm + Qv)
    return(c(pval.mean, pval.var, Qm, Qv, ratio, t.stat))
  })
  # browser()
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

      Qm <- qnorm(pval.mean, lower.tail = FALSE)
      Qv <- qnorm(pval.var, lower.tail = FALSE)

      Qm[Qm < 0] <- 0
      Qv[Qv < 0] <- 0
      ratio <- Qv / (Qm + Qv)
      return(c(pval.mean, pval.var, Qm, Qv, ratio, t.stat))
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
    nulltabs = lapply(1:B, function(i) {
      bumphunter::regionFinder(x = permStat[, i], chr = sites, pos = pos.gen, cluster = cluster, cutoff = cutoff, ind = Index, assumeSorted = TRUE, verbose = FALSE)
    })

    A = as.list(rep(0, B))
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