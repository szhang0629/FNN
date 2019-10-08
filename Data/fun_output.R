Fun.Output <- function(fct.name = "polyno") {
  fct <- switch(fct.name, polyno = .polyno, logist = .logist, linear = .linear)
  force(fct)
  return(fct)
}
.polyno <- function(G, pos, index = NULL, n = 20,
                       np = 10, nq = 10, iscale = T){
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
  func <- function(x = rep(0.5, nrow(G)), i = 1:nrow(G)) {
    G. <- G[i, ,drop = FALSE]
    y <- 0
    for (j in 1:n) {
      y <- y + cf[j] * (G.^ef[j] %*% B1[, sample1[j], drop = F]) *
        eval.basis(x, bb2)[, sample2[j], drop = F]
    }
    return(y)
  }
  force(func)
  func <- scale.func(func, iscale)
  return(func)
}
.logist <- function(G, pos, index = NULL, n = 20, iscale = T){
  set.seed(index)
  cf <- runif(n, -2, 2)
  ef <- runif(n, -10, 10)
  af <- runif(n, 0, 1)
  caf <- runif(n, -pi, pi)
  bf <- runif(n, 0, 1)
  cbf <- runif(n, -pi, pi)
  func <- function(x, i = 1:nrow(G)) {
    G. <- G[i, ,drop = FALSE]
    y <- 0
    for (j in 1:n) {
      y <- y + cf[j] * ((1/(1 + exp(ef[j] * G.))) %*%
                          t(t(cos(af[j] * pos + caf[j]))) *
                          t(t(cos(bf[j] * x + cbf[j]))))
    }
    return(y)
  }
  force(func)
  func <- scale.func(func, iscale)
  return(func)
}
.linear <- function(G, pos, index = NULL, n = 20, iscale = T){
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
  func <- function(x, i = 1:nrow(G)) {
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
  force(func)
  func <- scale.func(func, iscale)
  return(func)
}