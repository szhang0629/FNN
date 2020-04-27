Fun.Output.2 <- function(G, pos, index = NULL, n = 100, iscale = T){
  set.seed(index)
  cf <- runif(n, -2, 2)
  ef <- runif(n, 1/3, 3)
  pos0 <- (0:7)/8 + runif(8, 0.025, 0.1)
  loc1 <- (0:7)/8 + runif(8, 0.025, 0.1)
  loc2 <- (0:7)/8 + runif(8, 0.025, 0.1)
  bb0 <- create.bspline.basis(rangeval = 0:1, breaks = pos0, norder = 4)
  bb1 <- create.bspline.basis(rangeval = 0:1, breaks = loc1, norder = 4)
  bb2 <- create.bspline.basis(rangeval = 0:1, breaks = loc2, norder = 4)
  sample0 <- sample(1:bb0$nbasis, n, replace = T)
  sample1 <- sample(1:bb1$nbasis, n, replace = T)
  sample2 <- sample(1:bb2$nbasis, n, replace = T)
  B0 <- eval.basis(pos, bb0)
  # func <- function(x1 = rep(0.5, nrow(G)), x2 = rep(0.5, nrow(G)), 
  #                  i = 1:nrow(G)) {
  #   G. <- G[i, ,drop = FALSE]
  #   y <- 0
  #   for (j in 1:n) {
  #     y <- y + cf[j] * (G.^ef[j] %*% B0[, sample0[j], drop = F]) *
  #       (eval.basis(x1, bb1)[, sample1[j], drop = F] *
  #          eval.basis(x2, bb2)[, sample2[j], drop = F])
  #   }
  #   return(y)
  # }
  func <- function(x1 = 0.5, x2 = 0.5, i = 1:nrow(G)) {
    G. <- G[i, ,drop = FALSE]
    y <- 0
    for (j in 1:n) {
      y <- y + cf[j] * (G.^ef[j] %*% B0[, sample0[j], drop = F]) %*%
        t(eval.basis(x1, bb1)[, sample1[j], drop = F] *
           eval.basis(x2, bb2)[, sample2[j], drop = F])
    }
    return(y)
  }
  force(func)
  func <- scale.func.2(func, iscale)
  return(func)
}
scale.func.2 <- function(f, iscale = T){
  if (iscale) {
    n <- length(f())
    points <- (1:10 - 0.5) / 10
    Y.fun <- f(points, points)
    Y.mean <- mean(Y.fun)
    Y.sd <- sqrt(cost.(Y.fun, rep.row(colMeans(Y.fun), nrow(Y.fun))))
    force(Y.mean)
    force(Y.sd)
    g <- function(x, y, i) {
      y <- (f(x, y, i) - Y.mean) / Y.sd
      return(y)
    }
    return(g)
  } else return(f)
}
gData.2 <- function(index, noise = 0.1, n = 200, p = 500, loc){
  set.seed(index)
  Data. <- Data.Input(n = n, p = p)
  f <- Fun.Output.2(Data.$G, Data.$pos, index)
  Y <- f(loc$loc1, loc$loc2)
  Y <- Y + array(rnorm(prod(dim(Y))), dim = dim(Y))*noise
  Data <- list(G = Data.$G, Y = Y, pos = Data.$pos)
  return(Data)
}