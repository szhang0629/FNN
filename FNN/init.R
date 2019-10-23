init <- function(fBases, Ac = NULL, para.ub = 2, seed = NULL, X.col = 0){
  if (is.null(X.col)) X.col <- 0
  if (is.null(Ac)) {
    K <- init.(fBases[[1]]) + X.col
    para.lb <- -para.ub
    if (!is.null(seed))
      set.seed(seed)
    for (i in 2:length(fBases)) {
      K0 <- K
      K <- init.(fBases[[i]])
      A. <- matrix(runif(K*K0, para.lb, para.ub), ncol = K)
      c. <- rep(0, K)
      Ac. <- rbind(c., A.)
      Ac[[length(Ac) + 1]] <- Ac.
    }
  }
  return(Ac)
}
init. <- function(x) {
  if (is.null(x))
    return(0)
  if (is.basis(x))
    return(x$nbasis)
  if (is.numeric(x)) {
    if (is.matrix(x))
      return(ncol(x))
    else if (length(x) == 1)
      return(x)
    else
      return(length(x))
  }
  if (is.list(x))
    return(ncol(x[[1]]))
}
int.fun <- function(fun, int = NULL) {
  if (is.null(int)) {
    if (is.list(fun) && !is.basis(fun)) {
      result <- list()
      for (i in 1:length(fun))
        result[[i]] <- int.fun(fun[[i]])
      return(result)
    } else {
      if (is.basis(fun)) return(TRUE)
      else return(FALSE)
    }
  } else
    return(int)
}
cbb <- function(pos, nbasis = NA, norder= NA, ratio = NA) {
  pos <- sort(pos)
  if (is.na(nbasis))
    nbasis <- length(pos) * ratio
  if (is.na(norder))
    norder = round(sqrt(nbasis))
  breaks. <- ceiling(seq(1, length(pos), length.out = nbasis - norder + 1))
  breaks. <- pos[breaks.]
  breaks <- (c(0, breaks.) + c(breaks., 1))/2
  return(create.bspline.basis(breaks = breaks, norder = norder))
}
sigmoid <- function(x){
  (1 + exp(-x))^(-1)
}
linear <- function(x){
  x
}
relu <- function(x){
  (abs(x) + x) / 2
}