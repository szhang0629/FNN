f2m <- function(bb = .bb, loc = NA) {
  if (is.matrix(bb)) # B1
    return(bb)
  if (is.basis(bb)) { # B2
    if (is.atomic(loc) && length(loc) == 1 &&is.na(loc)) 
      B <- f2m.(bb, ((1:(3*bb$nbasis)) - 0.5)/(3*bb$nbasis))
    else {
      B <- f2m.(bb, loc)
    }
    return(B)
  }
  if (is.numeric(bb) && length(bb) == 1) # B3
    return(diag(nrow = bb, ncol = bb))
  if (is.null(bb))
    return(NULL)
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
