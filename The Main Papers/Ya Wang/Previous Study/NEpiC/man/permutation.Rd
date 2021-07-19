\name{permutation}
\alias{permutation}

\title{
Validating identified modules using permutations  
}
\description{
 A function shuffling the tumor/adjacent normal labels within a tumor/adjacent normal pair and repeating the permutation procedures to calculate the permutated module scores for the identified top 1 percentage modules.
}
\usage{
permutation(normal_data, tumor_data, annotation, res, n = 100)
}
\arguments{
  \item{normal_data}{
The Illumina 450k methylation array results for normal samples. The rows represent the methylation sites and the columns represent different samples.
}
  \item{tumor_data}{
The Illumina 450k methylation array results for adjacent tumor paire samples. The rows represent the methylation sites and the columns represent different samples.
}
  \item{annotation}{
The annotation file for Illumina 450k methylation array data. The first column is the methylation sites' names, whose column name is "cpg" and the second column is the methylation sites' gene annotations, whose column name is "gene"
}
  \item{res}{
The result of the c_dms, which contains the information of identified modules.
}
  \item{n}{
The numbers of permutations.
}
}
\value{
A dataset contains identified modules' names, identified modules' score, identified modules' permutated scores and the last column contains the permutation p values for identified modules.
}
\references{
Peifeng Ruan, Jing Shen, Regina M. Santella, Shuigeng Zhou, Shuang Wang.(2015). NEpiC: a Network-assisted algorithm for Epigenetic studies using mean and variance Combined signals.
}
\examples{
#perm_result = function (normal_data, tumor_data, annotation, res, n = 100) 
}