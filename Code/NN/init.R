init <- function(fBases, Ac = NULL, para.ub = 1, seed = NULL){
  if (is.null(Ac)) {
    K <- init.(fBases[[1]])
    if (!is.null(seed))
      set.seed(seed)
    for (i in 2:length(fBases)) {
      K0 <- K
      K <- init.(fBases[[i]])
      A. <- matrix(rnorm(K*K0, sd = para.ub), ncol = K)
      c. <- rep(0, K)
      Ac. <- rbind(c., A.)
      Ac[[length(Ac) + 1]] <- Ac.
    }
  }
  return(Ac)
}
init. <- function(x) {
  if (is.null(x)) return(0)
  if (is.numeric(x)) {
    if (is.matrix(x)) return(ncol(x))
    else if (length(x) == 1) return(x)
    else return(length(x))
  }
}
