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
    # grads. <- pen.(Ac[[i]], grads[[i]], lambda, P[i:(i+1)])
    eg[[i]] <- rho * eg.[[i]] + (1 - rho) * (grads.)^2
    dAc <- -RMS(ex.[[i]]) / RMS(eg[[i]]) * grads.
    # dAc <- - 10 * grads.
    ex[[i]] <- rho * ex.[[i]] + (1 - rho) * dAc^2
    Ac[[i]] <- Ac[[i]] + dAc
  }
  return(list(Ac = Ac, eg. = eg, ex. = ex))
}
# pen. <- function(Ac, dAc, lambda = 0, P = 1) {
#   if (lambda == 0) return(dAc)
#   else if (length(P) < 2) return(dAc + lambda * Ac)
#   else return (dAc + lambda * Ac %*% P)
# }
pen. <- function(Ac, dAc, lambda = 0, P = 1) {
  if (lambda == 0) return(dAc)
  else if (length(P) < 2) return(dAc + lambda * Ac)
  else {
    n.w <- nrow(P$PA)/ncol(dAc)
    n <- nrow(dAc)
    row.A <- (n - n.w + 1):n
    row.c <- 1:(n - n.w)
    dc. <- subs(Ac, row.c) %*% P$Pc
    dA. <- matrix(as.vector(subs(Ac, row.A)) %*% P$PA, nrow = n.w)
    return(dAc + lambda * rbind(dc., dA.))
  }
}
# pen. <- function(Ac, dAc, lambda = 0, P = 1) {
#   if (lambda == 0) return(dAc)
#   else if (length(P) < 2) return(dAc + lambda * Ac)
#   else {
#     n.w <- nrow(P$Q2)
#     n <- nrow(dAc)
#     row.A <- (n - n.w + 1):n
#     row.c <- 1:(n - n.w)
#     dc. <- subs(Ac, row.c) %*% P$P2
#     dA. <- P$Q0 %*% subs(Ac, row.A) %*% P$P0 + P$Q2 %*% subs(Ac, row.A) %*% P$P2
#     return (dAc + lambda * rbind(dc., dA.))
#   }
# }