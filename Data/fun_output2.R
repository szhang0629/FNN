Fun.Output.2 <- function(G, pos, index = NULL, n = 100,
                       np = 10, nq = 10, iscale = T){
  set.seed(index)
  cf <- runif(n, -2, 2)
  ef <- runif(n, 1/3, 3)
  pos0 <- c(0, sort(runif(np)), 1)
  loc1 <- c(0, sort(runif(nq)), 1)
  loc2 <- c(0, sort(runif(nq)), 1)
  bb0 <- create.bspline.basis(pos0, norder = 4)
  bb1 <- create.bspline.basis(loc1, norder = 4)
  bb2 <- create.bspline.basis(loc2, norder = 4)
  sample0 <- sample(1:bb0$nbasis, n, replace = T)
  sample1 <- sample(1:bb1$nbasis, n, replace = T)
  sample2 <- sample(1:bb2$nbasis, n, replace = T)
  B0 <- eval.basis(pos, bb0)
  # force(G)
  func <- function(x1, x2, i = 1:nrow(G)) {
    G. <- G[i, ,drop = FALSE]
    y <- 0
    for (j in 1:n) {
      y <- y + cf[j] * (G.^ef[j] %*% B0[, sample0[j], drop = F]) %*%
        (t(eval.basis(x1, bb1)[, sample1[j], drop = F]) *
           t(eval.basis(x2, bb2)[, sample2[j], drop = F]))
    }
    return(y)
  }
  force(func)
  func <- scale.func.2(func, iscale)
  return(func)
}
scale.func.2 <- function(f, iscale = T){
  if (iscale) {
    points <- (1:10 - 0.5) / 10
    Y.fun <- f(rep(points, 10), rep(points, each  = 10))
    Y.mean <- mean(Y.fun)
    Y.sd <- sqrt(cost.(Y.fun, rep.row(colMeans(Y.fun), nrow(Y.fun))))
    force(Y.mean)
    force(Y.sd)
    g <- function(x, i) {
      y <- (f(x, i) - Y.mean) / Y.sd
      return(y)
    }
    return(g)
  } else return(f)
}
gData.2 <- function(index, noise = 0.1, n = 200, p = 500, loc1, loc2){
  set.seed(index)
  Data. <- Data.Input(n = n, p = p)
  f <- Fun.Output.2(Data.$G, Data.$pos, index)
  Y <- f(loc1, loc2)
  Y <- Y + array(rnorm(prod(dim(Y))), dim = dim(Y))*noise
  Data <- list(G = Data.$G, f = f, Y = Y, pos = Data.$pos)
  return(Data)
}