f2m <- function(bb, loc = NA) {
  if (is.matrix(bb)) # B1
    return(bb)
  if (is.basis(bb)) { # B2
    if (is.atomic(loc) && length(loc) == 1 &&is.na(loc)) 
      B <- f2m(bb, ((1:(3*bb$nbasis)) - 0.5)/(3*bb$nbasis))
    else if (!is.list(loc)) {
      # P0 <- bsplinepen(bb, 0)
      # coef <- sqrt(sum(diag(P0 %*% P0)))
      # B <- eval.basis(loc, bb) / coef
      B <- eval.basis(loc, bb) * sqrt(bb$nbasis)
    } else {
      B <- list()
      for (i in 1:length(loc))
        B[[i]] <- f2m(bb, loc[[i]])
    }
    return(B)
  }
  if (is.numeric(bb) && length(bb) == 1) # B3
    return(diag(nrow = bb, ncol = bb))
  if (is.null(bb))
    return(NULL)
  if (is.list(bb))
    return(bb)
}
Bf2m <- function(fBases, pos = NA, loc = NA){
  Bases <- list()
  n <- length(fBases)
  Bases[[1]] <- f2m(fBases[[1]], pos)
  for (i in 2:(n - 1)) {
    if ((is.list(fBases[[i]])) && (!is.basis(fBases[[i]]))) {
      BB <- fBases[[i]]
      m <- length(BB)
      B <- list()
      for (j in 1:m) 
        B[[j]] <- f2m(BB[[j]])
      Bases[[i]] <- as.matrix(bdiag(B))
    } else
      Bases[[i]] <- f2m(fBases[[i]])
  }
  Bases[[n]] <- f2m(fBases[[n]], loc)
  return(Bases = Bases)
}