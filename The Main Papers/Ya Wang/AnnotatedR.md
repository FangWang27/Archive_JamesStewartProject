# The annotated R code

## Introduction

This MD file serves as a translation for the R code used in paper by Wang. The code is preserved in the exact the same
way from R file given by the author.

## Actual code and explanation

### Environment setup

The working directory is set up and the necessary package is installed here. The package `bumphunter` is necessary for setting up the DMRs.

```R
require(bumphunter)

#####Set working directory (PLEASE REPLACE IT WITH YOUR OWN PATH)
setwd("D:/R_tutorial")


#####Load sample data
load("Sample_data.RData")

```

### Data preparation

The author subset the data obtained. The balance case-control experiment is assumed, here we have each group has $30$ individual. The function `data.matrix` is a function that transform the data from a `dataframe` to `Numeric Matrix`.

The data is organized in the way that the first 30 columns are the reading from cancer tissue and the 31st to 60th columns are from normal tissue.

```R
#####Prepare data for analysis
methyl.MVal=sample.data[,4:63]
methyl.MVal=data.matrix(methyl.MVal)

num.pairs=30
case.ind=1:30
ctrl.ind=31:60

```

The author using `::` operator to extract the `clusterMaker` from package `bumphunter`. By the documentation of `bumphunter`, the function `clusterMaker` group sites within the genomic distance specified by `maxGap`. The distance is given my the `MAPINFO` comes with the data and `CHR` is the chromosome number that coded as `factor`.

```R
#####Assign CpG sites into clusters
cluster=bumphunter::clusterMaker(sample.data$CHR, sample.data$MAPINFO, assumeSorted=TRUE, maxGap=1000)

```

### Function perform the test

The function below calculate p-value from one-sided Pitman-Morgan test. The function `complete.cases` takes a sequence of vectors/matrices
a logical vector A logical vector specifying which observations/rows have no missing values across the entire sequence.

This theoretical basis for this function comes from paper *A NOTE ON NORMAL CORRELATION* written by E.J.G Pitman in 1939 in *Biomerika*, We present the result used below.

#### one-sided Pitman-Morgan test

Let $X$ and $Y$ be normally correlated random variables with variance $\sigma_1^2$ and $\sigma_2^2$; $\mathbf{x} = (x_1,\ldots,x_n)$ and $\mathbf{y} = (y_1,\ldots,y_n)$ are paired sample obtained. Let

$$
\omega = \frac{\sigma_1^2}{\sigma_2^2} \qquad \text{and}\qquad  w = \frac{ \sum (x_i - \overline{\mathbf{x}})^2}{ \sum (y_i - \overline{\mathbf{y}})^2}.
$$

Then the random variable

$$
T = \frac{(w - \omega)\sqrt{n-2}}{\sqrt{4(1-r^2)w\omega}}
$$
follows $T_{n-2}$ distribution. One can using this to test the hypothesis $\sigma_1^2 \neq \sigma_2^2$ using null hypothesis $H_0 : w = 1$ and $H_1 : w > 1$.

```R
#####Define function to calculate p-value from one-sided Pitman-Morgan test
pitman.var.grt=function(x,y,ratio=1){
  xok <- yok <- complete.cases(x, y)
  x <- x[xok]
  y <- y[yok]
  n <- length(x)
  df = n - 2
  r <- cor(x, y)
  Var1 <- var(x)
  Var2 <- var(y)
  w <- var(x)/var(y)
  stat.t <- ((w - ratio) * sqrt(n - 2))/sqrt(4 * (1 - r^2) * w * ratio)  
  pval <- pt(-stat.t, df = df)
  return(list(statistic=stat.t, p.value=pval))
}
```

### Step 1 : Define the site-level mean and combined signal score

The score $S_i$ is given by

$$
S_i = \frac{|{T_{m_i}}|}{T_{m_i}} (\lambda m_i + (1 - \lambda)v_i)
$$

where $m_i = \mathbf{\Phi}^{-1}(1-p_{m_i})$ and $v_i = \mathbf{\Phi}^{-1}(1-p_{v_i})$ and $p_{m_i}$ and $p_{v_i}$ are $p$-valueP-values from the
two-sided paired $t$-test testing if the mean methylation measures
are the same between tumor and normal-adjacent tissues
and from the one-sided Pitman–Morgan test testing if the variance
of the methylation measures in tumor tissues is greater
than that in normal-adjacent tissues at CpG site $i$, respectively.

#### Calculate the necessary stats for each site

Since the first and last $30$ columns are reading from cancer and normal tissue respectively, the code below calculated the necessary information at each site. The variable `t.stat` from  `res.ttest` is the the two sample paired $t$ test for the methylation for entire cancer and normal tissue; `pval.mean` and `pval.var` are the $p$-value for the mean and variance for testing the mean and variance, corresponding to $p_{m_i}$ and $p_{v_i}$. Then `Qm` and `Qv` is corresponding to $m_i$ and $v_i$ respectively.

```R
#####Calculate combined signal scores for observed data
rawStat = sapply(1:nrow(methyl.MVal), function(i){
  res.ttest=t.test(methyl.MVal[i,case.ind], methyl.MVal[i,ctrl.ind], paired=TRUE)
  t.stat=res.ttest$statistic
  
  pval.mean=res.ttest$p.value
  pval.var=pitman.var.grt(methyl.MVal[i,case.ind], methyl.MVal[i,ctrl.ind])$p.value
  
  Qm=qnorm(pval.mean, lower.tail = FALSE)
  Qv=qnorm(pval.var, lower.tail = FALSE)
```

In the rest of the code the author set the site with negative $m_i$ and $v_i$ to $0$ and remove the site with $m_i = v_i = 0$ latter, since `ratio` will be `NA` in such case. The `ratio` is corresponding to the site-level scaling parameter $\lambda_i = v_i/(m_i + v_i)$.

```R  
  Qm[Qm<0]=0
  Qv[Qv<0]=0
  ratio=Qv/(Qm+Qv)
  return(c(pval.mean, pval.var, Qm, Qv, ratio, t.stat))
})
```

In the following code the site-level score $S_i$ is given by `rawStat.wt.sign`.

```R

lambda=mean(rawStat[5,],na.rm=T)
rawStat.wt=lambda*rawStat[3,]+(1-lambda)*rawStat[4,]
rawStat.wt.sign=rawStat.wt*sign(rawStat[6,])
```

The function `runmedByCluster` is the function using running median to smooth each cluster of genomic locations, where argument `k` is the window size of the running median. Returned variable `fitted` contains smoothed data and `smoothed` is a boolean vector indicating whether a given position was smoothed. The variable `stat` is the vector contains $\tilde{S_i}$, the smoothed site-level score.

The function `which` takes a boolean vector and returns the indexes of `TRUE` by default.

```R
stat=bumphunter::runmedByCluster(y=rawStat.wt.sign, x=NULL, cluster=cluster, weights=NULL, k=5, endrule="constant", verbose=FALSE)
Index=which(stat$smoothed)
stat=stat$fitted
```

### Step 3,4: Using permutation test to exam the significant of DMRs

Here the author generate a $30$ by 250 matrix called `perm.sign`, where th `i`th column is the result of shuffling `i` th times of all samples.

```R
#####Calculate combined signal scores for 250 permutation datasets
B=250
set.seed(123456)
perm.sign=sapply(1:B, function(i){
  sample(c(-1,1), num.pairs, replace=TRUE)
})
```

The following function perform the actual permutation job, and we will break it to parts to explore each line's specific role.

#### Permutation  on  the indexes

The `shuff.ind` is a vector indicating which position among the paired group we are going to swap.

```R
permStat=sapply(1:B, function(i){
  shuff.ind=which(perm.sign[,i]==-1)
```

Since the data is given in the way that first `i` th column is the paired sample with `i + num.pairs` th column, and therefore the following code shuffled the samples for the groups specified by `shuff.ind`.

```R
  perm.case.ind=case.ind
  perm.case.ind[shuff.ind]=case.ind[shuff.ind]+num.pairs

  perm.ctrl.ind=ctrl.ind
  perm.ctrl.ind[shuff.ind]=ctrl.ind[shuff.ind]-num.pairs
```

#### Calculating the statistics after the permutation

The function `permRawStat` calculated all the necessary statistics for the sample being shuffled. It's essentially the same function as `rawStat`, the only difference is `rawStat` using the fact that the data is grouped in the way that first `num.pairs` and the second `num.pairs` are paired sample wheres `permRawStat` uses `perm.case.ind` generated ealier to specify the sample

```R
  permRawStat = sapply(1:nrow(methyl.MVal), function(i){
    res.ttest=t.test(methyl.MVal[i,perm.case.ind], methyl.MVal[i,perm.ctrl.ind], paired=TRUE)
    t.stat=res.ttest$statistic

    pval.mean=res.ttest$p.value
    pval.var=pitman.var.grt(methyl.MVal[i,perm.case.ind], methyl.MVal[i,perm.ctrl.ind])$p.value

    Qm=qnorm(pval.mean, lower.tail = FALSE)
    Qv=qnorm(pval.var, lower.tail = FALSE)

    Qm[Qm<0]=0
    Qv[Qv<0]=0
    ratio=Qv/(Qm+Qv)
    return(c(pval.mean, pval.var, Qm, Qv, ratio, t.stat))
  })
  ```

The following code performed the same procedure to as before to obtain $\tilde{S_i}$ for the shuffled samples.

  ```R
  perm.lambda=mean(permRawStat[5,],na.rm=T)
  permRawStat.wt=perm.lambda*permRawStat[3,]+(1-perm.lambda)*permRawStat[4,]
  permRawStat.wt.sign=permRawStat.wt*sign(permRawStat[6,])
  
  permStat=bumphunter::runmedByCluster(y=permRawStat.wt.sign, x=NULL, cluster=cluster, weights=NULL, k=5, endrule="constant", verbose=FALSE)
  permStat=permStat$fitted
  return(permStat)
})

```

### Identify candidate DMRs

The code below identifies the candidates DMRs . According to the documentation the function `regionFinder` finds regions for which a numeric vector is above (or below) predefined thresholds. Note that the threshold `cutoff` is estimated from the permutation test while the candidate DMRs `tab` is found from the observed data. In this case, `cutoff` is corresponding to $k$.

```R
#####Get threshold as the 99th percentile of genome-wide smoothed signal scores  
cutoff=quantile(abs(permStat), 0.99, na.rm=T)


#####Identify candidate DMRs 
tab=bumphunter::regionFinder(x=stat, chr=sample.data$CHR, pos=sample.data$MAPINFO, cluster=cluster, cutoff=cutoff, ind=Index, assumeSorted=TRUE, verbose=FALSE)
```

### Assess significance of candidate DMRs

Here I tried to give an explanation for the motivation for the $p$-value they defined. Under $H_0$, there is no difference between the the methylation in two group of samples. Hence, under $H_0$, for the $j$th candidates DMR we found, its score $A_j$ should not be different from another DMR obtained from samples being shuffled randomly. Then we would expect the event $\{ A_{\text{perm}_g,t} > A_j\}$ happen for fairly large number of times, where $g$ is the index for permuatation and $t$ is index for the $t$ th DMR generated from that permutation. Then $P_j$ defined by the author makes sense, the denominator is the normalizing constant to make $0 \leq P_j \leq 1$.

The function `regionFinder` summarized the characteristics of candidates DMRs observed and generated in `nulltabs` and `tab`. In particular, it follows from the author's definition the combined score of a certain region can be convenient be calculated as the area of the corresponding "bump". Then `A` and `Avalue` store the score of observed and generated DMR respectively.

```R
#####Calculate p-values for candidate DMRs
if(nrow(tab)>0){
  nulltabs=lapply(1:B, function(i){
    bumphunter::regionFinder(x=permStat[ ,i], chr=sample.data$CHR, pos=sample.data$MAPINFO, cluster=cluster, cutoff=cutoff, ind=Index, assumeSorted=TRUE, verbose=FALSE)
  })
  
  A=as.list(rep(0, B))
  for (i in 1:B) {
    nulltab=nulltabs[[i]]
    if (nrow(nulltab) > 0) {
      A[[i]]=nulltab$area
    }
  }

  Avalue <- matrix(tab$area, ncol = 1)
```

Here `tots` serves as the summation of the indicator functions in the numerator.

```R
  tots=sapply(1:nrow(Avalue), function(i){ # For ith observed DMR candidates
    return(sapply(seq(along = A), function(j) { # For the jth permutation
      sum((A[[j]]>=Avalue[i]) | (abs(A[[j]]-Avalue[i])<=sqrt(.Machine$double.eps))) # sum of Indicator functions defined
    }))
  })
```

The following code is used to calculate the $p$-value for each site

```R
  if (is.vector(tots)) {
    tots <- matrix(tots, nrow = 1)
  }
  tab$p.valueArea=colSums(tots)/sum(sapply(nulltabs, nrow)) # nrow(nulltabs) corresponding to how many regions generate at each iteration, and summing up corresponding to the denominator.
  tab$fwerArea=colMeans(tots > 0) # the family wise error rate is given by 
  tab=tab[order(tab$fwer, -tab$area), ]
}
```