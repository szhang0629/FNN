scale.func <- function(f, iscale = T, align){
  if (iscale) {
    n <- length(f())
    if (!align) {
      Y.fun <- f(rep((1:100 - 0.5) / 100, n), rep(1:n, each = 100))
      Y.mean <- mean(Y.fun)
      Y.sd <- sqrt(cost.(Y.fun, rep.row(colMeans(Y.fun), nrow(Y.fun))))
      force(Y.mean)
      force(Y.sd)
      g <- function(x, i = 1:nrow(Y.fun)) {
        y <- (f(x, i) - Y.mean) / Y.sd
        return(y)
      }
    } else {
      Y.fun <- f((1:100 - 0.5) / 100)
      Y.mean <- mean(Y.fun)
      Y.sd <- sqrt(cost.(Y.fun, rep.row(colMeans(Y.fun), nrow(Y.fun))))
      force(Y.mean)
      force(Y.sd)
      g <- function(x, i = 1:nrow(Y.fun)) {
        y <- (f(x, i) - Y.mean) / Y.sd
        return(y)
      }
    }
    return(g)
  } else return(f)
}
Data.Output <- function(f, loc, noise = 0.1){
  force(f)
  if (is.data.frame(loc)) {
    Y <- c(f(loc$loc, loc$PTID) + rnorm(length(loc$loc)) * noise)
  }  else {
    if (is.numeric(loc)) {
      Y <- f(loc)
      Y <- Y + array(rnorm(prod(dim(Y))), dim = dim(Y))*noise
    }
  }
  return(Y)
}