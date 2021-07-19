\name{candidate_genes}
\alias{candidate_genes}

\title{
Chosing candidate genes.
}
\description{
Chosing genes in top ranked modules as candidate genes for outcomes. Further prioritized candidate genes can be chosen based on percentages selected by top ranked modules selected if prioritizing = T.
}
\usage{
candidate_genes(res, pct = 0.01, prioritizing)
}

\arguments{
  \item{res}{
The result for c_dms.
}
  \item{pct}{
A float which determines the number of candidate genes. For example, if pct is 0.01 as default setting, the genes in top 1 percent modules will be chosen as the candidate genes.
}
  \item{prioritizing}{
To further prioritize the candidate genes or not.
}
}

\value{
A list for candidate genes or prioritized candidate genes.
}
\references{
NEpiC: a Network-assisted algorithm for Epigenetic studies using mean and variance Combined signals Peifeng Ruan, Jing Shen, Regina M. Santella, Shuigeng Zhou, Shuang Wang.(2015)
}

\examples{
# candidate_genes(res, pct = 0.01, prioritizing = T)
}