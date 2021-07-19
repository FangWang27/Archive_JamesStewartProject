\name{c_dms}
\alias{c_dms}

\title{
Combined signals dense modules searching algorithm.
}
\description{
c_dms uses the combined signals in both mean and variance differences in methylation measures between case and control groups to search for differentially methylated sub-networks (modules) using a protein-protein interaction (PPI) network. It returns ordered dense modules. }
\usage{
c_dms(network, gene2weight, d = 1, r = 0.1)
}
\arguments{
  \item{network}{
A dataframe describing the edges of s PPI network. Two columns are needed in this dataset: interactor A and interactor B.
}
  \item{gene2weight}{
A dataframe describing combined gene-level scores. The first column of this dataframe, named 'gene', contains the gene names, which should use the same annotation as those used in the network. The second column 'weight' contains gene-level combined scores.
}
  \item{d}{
An integer determining the order of neighbour genes to be searched in NEpiC algorithm. We suggest to set d=1 for NEpiC.
}
  \item{r}{
A number determining the threshold for including a new gene into the module.
}
}
\value{
A list containing several statistics for modules, including genes in modules, module scores before and after normalization. 
}
\references{
Peifeng Ruan, Jing Shen, Regina M. Santella, Shuigeng Zhou, Shuang Wang.(2015). NEpiC: a Network-assisted algorithm for Epigenetic studies using mean and variance Combined signals.\cr
Peilin Jia, Siyuan Zheng, Jirong Rong, Wei Zheng, Zhongming Zhao. (2011) Bioformatics. dmGWAS: dense module searching for genome-wide association studies in protein-protein interaction networks.
}
\examples{
#res=c_dms(network, gene2weight, d = 1, r = 0.1) 
}

 