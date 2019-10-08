f2m <- function(bb, loc = NA) {
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
# pen.Bases <- function(Bases.) {
#   result <- list()
#   for (i in 2:length(Bases.))
#     result[[i - 1]] <- pen.fun(Bases.[[i]])
#   return(result)
# }
pen.Bases <- function(Bases.) {
  result <- list()
  for (i in 1:(length(Bases.)-1)) {
    P2 <- pen.fun(Bases.[[i+1]])
    P0 <- pen.fun(Bases.[[i+1]], 0)
    P2. <- pen.fun(Bases.[[i]])
    P0. <- pen.fun(Bases.[[i]], 0)
    PA <- P0 %x% P2. + P2 %x% P0.
    result[[i]] <- list(Pc = P2, PA = PA)
  }
  return(result)
}
# pen.Bases <- function(Bases.) {
#   result <- list()
#   for (i in 1:(length(Bases.)-1)) {
#     P2 <- pen.fun(Bases.[[i+1]])
#     P0 <- pen.fun(Bases.[[i+1]], 0)
#     P2. <- pen.fun(Bases.[[i]])
#     P0. <- pen.fun(Bases.[[i]], 0)
#     J <- matrix(1, ncol = ncol(P2.), nrow = nrow(P2.))
#     Q0 <- J %*% P2. + P2. %*% J
#     Q2 <- J %*% P0. + P0. %*% J
#     result[[i]] <- list(P2 = P2, P0 = P0, Q2 = Q2, Q0 = Q0)
#   }
#   return(result)
# }