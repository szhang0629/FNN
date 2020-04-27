f2m <- function(bb, loc = NA) {
  if (is.matrix(bb)) # B1
    return(bb)
  if (is.basis(bb)) { # B2
    if (is.atomic(loc) && length(loc) == 1 && is.na(loc)) 
      B <- f2m.(bb, ((1:(5*bb$nbasis)) - 0.5)/(5*bb$nbasis))
    else {
      B <- f2m.(bb, loc)
    }
    return(B)
  }
  if (is.numeric(bb) && length(bb) == 1) # B3
    return(diag(nrow = bb, ncol = bb))
  if (is.list(bb)) {
    if ("bb0" %in% names(bb))
      return(as.matrix(bdiag(f2m(bb[[1]], loc$pos1), f2m(bb[[2]], loc$pos2))))
    return(f2m(bb[[1]], loc$loc1)[, rep(1:bb[[1]]$nbasis, bb[[2]]$nbasis)] * 
      f2m(bb[[2]], loc$loc2)[, rep(1:bb[[2]]$nbasis, each = bb[[1]]$nbasis)])
  }
  if (is.null(bb))
    return(NULL)
}
c.fun <- function(bb) {
  return(bb$nbasis/sum(diag(bsplinepen(bb, 0))))
  # return(bb$nbasis)
  # return(1)
}
f2m. <- function(bb, loc) {
  if (!is.matrix(loc)) {
    return(eval.basis(loc, bb) * sqrt(c.fun(bb)))
  } else if (ncol(loc) == 1)
    B <- f2m.(bb, c(loc))
  else B <- (f2m.(bb, loc[, "loc"]) - f2m.(bb, loc[, "loc0"]))
  return(B)
}
Bf2m <- function(fBases, pos = NA, loc = NA){
  Bases <- list()
  n <- length(fBases)
  Bases[[1]] <- f2m(fBases[[1]], pos)
  for (i in 2:(n - 1)) 
    Bases[[i]] <- f2m(fBases[[i]])
  Bases[[n]] <- f2m(fBases[[n]], loc)
  return(Bases = Bases)
}
integ <- function(D, B, int = T){
  if (is.null(D)) return(NULL)
  if (int) return(as.matrix(D) %*% B / nrow(B))
  else return(as.matrix(D) %*% B)
}