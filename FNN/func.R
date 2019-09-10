integ <- function(D, B, int = F){
  if (is.null(D))
    return(NULL)
  if (int) return(as.matrix(D) %*% B / nrow(B))
  else return(as.matrix(D) %*% B)
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
cbb <- function(norder, pos, nbasis = NA, ratio = NA) {
  pos <- sort(pos)
  if (is.na(nbasis))
    nbasis <- length(pos) * ratio
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