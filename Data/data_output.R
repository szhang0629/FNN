Data.Output <- function(f, loc, noise = 0.1){
  force(f)
  if (is.list(loc)) {
    Y <- list()
    for (i in 1:length(loc))
      Y[[i]] <- f((loc[[i]]), i) + rnorm(length(loc[[i]])) * noise
  } else {
    if (is.numeric(loc))
      Y <- f(loc)
    Y <- Y + array(rnorm(prod(dim(Y))), dim = dim(Y))*noise
  }
  return(Y)
}
scale.func <- function(f, iscale = T){
  if (iscale) {
    Y.fun <- f((1:100 - 0.5) / 100)
    Y.mean <- mean(Y.fun)
    Y.sd <- sqrt(cost.(Y.fun, rep.row(colMeans(Y.fun), nrow(Y.fun))))
    force(Y.mean)
    force(Y.sd)
    g <- function(x, i = 1:nrow(Y.fun)) {
      y <- (f(x, i) - Y.mean) / Y.sd
      return(y)
    }
    return(g)
  } else return(f)
}