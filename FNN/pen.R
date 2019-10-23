RMS <- function(g, epsilon = 1e-6){
  return(sqrt(g + epsilon))
}
adadelta <- function(output, grads, lambda = 0, P = NULL, P0 = NULL, 
                     rho = 0.95) {
  eg <- list()
  ex <- list()
  Ac <- output$Ac
  eg. <- output$eg.
  ex. <- output$ex.
  for (i in 1:length(grads)) {
    grads. <- pen.(Ac[[i]], grads[[i]], lambda, P[[i]])
    # grads. <- pen.(Ac[[i]], grads[[i]], lambda, P[i:(i + 1)])
    eg[[i]] <- rho * eg.[[i]] + (1 - rho) * (grads.)^2
    dAc <- -RMS(ex.[[i]]) / RMS(eg[[i]]) * grads.
    # dAc <- - 10 * grads.
    ex[[i]] <- rho * ex.[[i]] + (1 - rho) * dAc^2
    Ac[[i]] <- Ac[[i]] + dAc
  }
  return(list(Ac = Ac, eg. = eg, ex. = ex))
}
## ------matrix smooth------
pen.Bases <- function(Bases.) {
  result <- list()
  for (i in 1:(length(Bases.) - 1)) {
    P2 <- pen.fun(Bases.[[i + 1]])
    P0 <- pen.fun(Bases.[[i + 1]], 0)
    P2. <- pen.fun(Bases.[[i]])
    P0. <- pen.fun(Bases.[[i]], 0)
    J <- matrix(1, ncol = ncol(P2), nrow = nrow(P2))
    Q0 <- (J %*% P2 + P2 %*% J) / 2
    Q2 <- (J %*% P0 + P0 %*% J) / 2
    result[[i]] <- list(P2 = P2., P0 = P0., Q2 = Q2, Q0 = Q0, Pc = P2)
  }
  return(result)
}
pen. <- function(Ac, dAc, lambda = 0, P = 1) {
  if (lambda == 0) return(dAc)
  else if (length(P) < 2) return(dAc + lambda * Ac)
  else {
    n.w <- nrow(P$P2)
    n <- nrow(dAc)
    row.A <- (n - n.w + 1):n
    row.c <- 1:(n - n.w)
    dc. <- subs(Ac, row.c) %*% P$Pc
    dA. <- P$P0 %*% subs(Ac, row.A) %*% P$Q0 + P$P2 %*% subs(Ac, row.A) %*% P$Q2
    return(dAc + lambda * rbind(dc., dA.))
  }
}
## ------smooth vectorize------
# pen.Bases <- function(Bases.) {
#   result <- list()
#   for (i in 1:(length(Bases.) - 1)) {
#     P2 <- pen.fun(Bases.[[i + 1]])
#     P0 <- pen.fun(Bases.[[i + 1]], 0)
#     P2. <- pen.fun(Bases.[[i]])
#     P0. <- pen.fun(Bases.[[i]], 0)
#     PA <- (P0 %x% P2. + P2 %x% P0.)
#     result[[i]] <- list(Pc = P2, PA = PA)
#   }
#   return(result)
# }
# pen. <- function(Ac, dAc, lambda = 0, P = 1) {
#   if (lambda == 0) return(dAc)
#   else if (length(P) < 2) return(dAc + lambda * Ac)
#   else {
#     n.w <- nrow(P$PA)/ncol(dAc)
#     n <- nrow(dAc)
#     row.A <- (n - n.w + 1):n
#     row.c <- 1:(n - n.w)
#     dc. <- subs(Ac, row.c) %*% P$Pc
#     dA. <- matrix(as.vector(subs(Ac, row.A)) %*% P$PA, nrow = n.w)
#     return(dAc + lambda * rbind(dc., dA.))
#   }
# }
## ------Smooth Output------
# pen.Bases <- function(Bases.) {
#   result <- list()
#   for (i in 2:length(Bases.))
#     result[[i - 1]] <- pen.fun(Bases.[[i]])
#   return(result)
# }
# pen. <- function(Ac, dAc, lambda = 0, P = 1) {
#   if (lambda == 0) return(dAc)
#   else if (length(P) < 2) return(dAc + lambda * Ac)
#   else return(dAc + lambda * Ac %*% P)
# }