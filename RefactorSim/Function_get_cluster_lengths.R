get_cluster_lengths <- function(sites, maxGap = 1000, min.size = 3, max.size = 15) {
  # the function takes a sites file and some necessary parameters for generating
  # clusters by clusterMaker. It returns a numeric vector that record the size of each cluster. More
  # specifically, the ith position of the vector records how many CpG sites does Cluster i has.

  # sites : an file with MAPINFO of sites on a chromosome
  # maxGap: the parameter pass to clusterMaker to generate clusters
  # min.size: clusters with sites less than this will be discard
  # max.size: clusters with more sites than this will be discard
  require("bumphunter")

  cl <- bumphunter::clusterMaker(
    chr = sites$Chromosome_36,
    pos = sites$MAPINFO,
    maxGap = maxGap
  )
  cl <- sort(cl)
  cl.lst <- rle(cl)
  lengths <- cl.lst$lengths[cl.lst$lengths >= min.size & cl.lst$lengths <= max.size ]
  lengths <- sort(lengths)
  return(lengths)
}