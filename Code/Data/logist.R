# .logist
fct <- function(G, pos, index = NULL, n = 20, iscale = T, align = T){
  set.seed(index)
  cf <- runif(n, -2, 2)
  ef <- runif(n, -10, 10)
  af <- runif(n, 0, 1)
  caf <- runif(n, -pi, pi)
  bf <- runif(n, 0, 1)
  cbf <- runif(n, -pi, pi)
  if (align) {
    func <- function(x = 0.5, i = 1:nrow(G)) {
      G. <- G[i, ,drop = FALSE]
      y <- 0
      for (j in 1:n) {
        y <- y + cf[j] * ((1/(1 + exp(ef[j] * G.))) %*%
                            t(t(cos(af[j] * pos + caf[j]))) %*%
                            t(cos(bf[j] * x + cbf[j])))
      }
      return(y)
    }
  } else {
    func <- function(x = rep(0.5, nrow(G)), i = 1:nrow(G)) {
      G. <- G[i, ,drop = FALSE]
      y <- 0
      for (j in 1:n) {
        y <- y + cf[j] * ((1/(1 + exp(ef[j] * G.))) %*%
                            t(t(cos(af[j] * pos + caf[j]))) *
                            t(t(cos(bf[j] * x + cbf[j]))))
      }
      return(y)
    }
  }
  force(func)
  func <- scale.func(func, iscale, align)
  return(func)
}