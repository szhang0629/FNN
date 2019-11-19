# .polyno
fct <- function(G, pos, index = NULL, n = 20,
                    np = 10, nq = 10, iscale = T, align = T){
  set.seed(index)
  cf <- runif(n, -2, 2)
  ef <- runif(n, 1/3, 3)
  pos0 <- c(0, sort(runif(np)), 1)
  loc <- c(0, sort(runif(nq)), 1)
  bb1 <- create.bspline.basis(pos0, norder = 4)
  bb2 <- create.bspline.basis(loc, norder = 4)
  n1 <- bb1$nbasis
  n2 <- bb2$nbasis
  sample1 <- sample(1:n1, n, replace = T)
  sample2.1 <- sample(1:n2, n2, replace = F)
  sample2.2 <- sample(1:n2, n - n2, replace = T)
  sample2 <- c(sample2.1, sample2.2)
  B1 <- eval.basis(pos, bb1)
  if (align) {
    func <- function(x = 0.5, i = 1:nrow(G)) {
      G. <- G[i, ,drop = FALSE]
      y <- 0
      for (j in 1:n) {
        y <- y + cf[j] * (G.^ef[j] %*% B1[, sample1[j], drop = F]) %*%
          t(eval.basis(x, bb2)[, sample2[j], drop = F])
      }
      return(y)
    }
  } else {
    func <- function(x = rep(0.5, nrow(G)), i = 1:nrow(G)) {
      G. <- G[i, ,drop = FALSE]
      y <- 0
      for (j in 1:n) {
        y <- y + cf[j] * (G.^ef[j] %*% B1[, sample1[j], drop = F]) *
          eval.basis(x, bb2)[, sample2[j], drop = F]
      }
      return(y)
    }
  }
  force(func)
  func <- scale.func(func, iscale, align)
  return(func)
}