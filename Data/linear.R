# .linear 
fct <- function(G, pos, index = NULL, n = 20, iscale = T, align = T){
  set.seed(index)
  n <- ceiling(n/2)
  cf1 <- runif(n, -12, 12)
  cf2 <- runif(n, -2, 2)
  af <- runif(n, -12, 12)
  caf <- runif(n, -pi, pi)
  bf <- runif(n, -12, 12)
  cbf <- runif(n, -pi, pi)
  loc1 <- c(0, sort(runif(10)), 1)
  loc2 <- c(0, sort(runif(10)), 1)
  bb1 <- create.bspline.basis(loc1, norder = 4)
  bb2 <- create.bspline.basis(loc2, norder = 4)
  n1 <- bb1$nbasis
  n2 <- bb2$nbasis
  sample1 <- sample(1:n1, n, replace = T)
  sample2 <- sample(1:n2, n, replace = T)
  B1 <- eval.basis(pos, bb1)
  if (align) {
    func <- function(x = 0.5, i = 1:nrow(G)) {
      G. <- G[i, ,drop = FALSE]
      y <- 0
      for (i in 1:n) {
        y <- y + cf1[i] * (G. %*% t(t(cos(af[i] * pos + caf[i]))) %*%
                             t(cos(bf[i] *  x + cbf[i]))) + 
          cf2[i] * (G. %*% B1[, sample1[i], drop = F]) %*%
          t(eval.basis(x, bb2)[, sample2[i], drop = F])
      }
      return(y)
    }
  } else {
    func <- function(x = rep(0.5, nrow(G)), i = 1:nrow(G)) {
      G. <- G[i, ,drop = FALSE]
      y <- 0
      for (i in 1:n) {
        y <- y + cf1[i] * (G. %*% t(t(cos(af[i] * pos + caf[i]))) * 
                             t(t(cos(bf[i] *  x + cbf[i])))) + 
          cf2[i] * (G. %*% B1[, sample1[i], drop = F]) *
          eval.basis(x, bb2)[, sample2[i], drop = F]
      }
      return(y)
    }
  }
  force(func)
  func <- scale.func(func, iscale, align)
  return(func)
}