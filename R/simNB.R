# Generate data from NB model
# between two selected groups
# n1, n2: number of replicates of each group
# mu1, mu2: mean of each group
# theta1, theta2: dispersion of each group
#
# simNB(5, 5, 300,400,0.3,0.2) to see an example simulated reads
# Author: Xu Qiao
# Created: 19th, Sep, 2017
# Last Modifed: 20th, Dec, 2017
#'@importFrom stats rnbinom var
simNB <- function(n1, n2,
                  mu1, mu2,
                  theta1, theta2) {
  if(mu1 < 1) mu1 <- 1 # mean count lower than 1 is not allowed
  if(mu2 < 1) mu2 <- 1
  # cat(sprintf("\rn1: %s, n2: %s, mu1: %s, mu2:%s, theta1: %s, theta: %s",
  #             n1, n2, mu1, mu2, theta1, theta2))
  while(1) {
    data.group1 <- rnbinom(n=n1, mu=mu1, size=1/theta1)
    data.group2 <- rnbinom(n=n2, mu=mu2, size=1/theta2)

    # avoid zero/negative data from both groups
    if(var(data.group1) > 0 &
       var(data.group2) > 0 &
       (sum(data.group1 < 0) + sum(data.group2 < 0)) == 0) break()
  }
  cts <- c(data.group1, data.group2)
  group.vec <- rep(c("A", "B"), c(n1,n2))
  data.btw.2groups <- data.frame(data=cts, group=group.vec)
  return(data.btw.2groups)
}
