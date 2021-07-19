	\name{c_score}
\alias{c_score}
\title{
Computing the mean and variance combined signal scores
}
\description{
A function computing gene-level scores combining both signals in mean and variance differences. 
}
\usage{
c_score(mean_sig, var_sig, lambda = "default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mean_sig}{
A dataframe describing mean signals for every CpG sites. Three columns are needed in this dataset, which are 'cpg', 'gene' and 'mean_p'. 'cpg' column has names of CpG sites. 'gene' column has has gene annotations for CpG sites. 'mean_p' column has p values of mean differences for CpG sites.
}
  \item{var_sig}{
A dataframe describing variance signals for every CpG sites. Three columns are needed in this dataset, which are 'cpg', 'gene' and 'var_p'. 'cpg' column has names of CpG sites. 'gene' column has has gene annotations for CpG sites. 'var_p' column has p values of variance differences for CpG sites.
}
  \item{lambda}{
A parameter to balance mean and variance scores. NEpiC computes it as default if it is not otherwise specified by user. If lambda = 1, NEpiC will search dense modules using only mean signals. If lambda = 0, NEpiC will search dense modules using only variance signals. Users can also specify any float between 0 and 1 for lambda if they have any additional prior biology information.
}
}
\value{
c_score returns a dataframe describing the combined gene-level scores, which can be used directly in c_dms.
}
\references{
NEpiC: a Network-assisted algorithm for Epigenetic studies using mean and variance Combined signals. Peifeng Ruan, Jing Shen, Regina M. Santella, Shuigeng Zhou, Shuang Wang.(2015)
}

\examples{
#c_score(mean_sig, var_sig, lambda = "default")
}