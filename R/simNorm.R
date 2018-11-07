# Generate data from Gaussian model
# n: sample size
# mu1: mean of group 1
# mu2: mean of group 2
# sd1: standard deviation of group 1
# sd2: standard deviation of group 2
# OUTPUT:
# a data frame with abundance and group
# Xu Qiao
# Created: 22nd, Sep, 2017
# Last Modification: 3rd, Jan, 2018
#' @importFrom  stats rnorm
#'
simNorm <- function(n1, n2, mu1, mu2, sd1, sd2) {
  if(mu1 < 1) mu1 <- 1 # mean count lower than 1 is not allowed
  if(mu2 < 1) mu2 <- 1
  while(1) { #only non-negative values are allowed
    data.group1 <- rnorm(n=n1, mean=mu1, sd=sd1)
    data.group2 <- rnorm(n=n2, mean=mu2, sd=sd2)
    neg.g1 <- sum(data.group1 < 0)
    neg.g2 <- sum(data.group2 < 0)
    if(sum(neg.g1, neg.g2) == 0) break()
    }
  data <- c(data.group1, data.group2)
  group.vec <- rep(c("A", "B"), c(n1,n2))
  abds.btw.2groups <- data.frame(data=data, group=group.vec)
  return(abds.btw.2groups)
}
