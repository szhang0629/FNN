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
#    grads. <- pen.(Ac[[i]], grads[[i]], lambda, P[[i]])
    grads. <- pen.(Ac[[i]], grads[[i]], lambda, P[i:(i+1)])
    eg[[i]] <- rho * eg.[[i]] + (1 - rho) * (grads.)^2
    dAc <- - RMS(ex.[[i]]) / RMS(eg[[i]]) * grads.
    ex[[i]] <- rho * ex.[[i]] + (1 - rho) * dAc^2
    Ac[[i]] <- Ac[[i]] + dAc
  }
  return(list(Ac = Ac, eg. = eg, ex. = ex))
}
#pen. <- function(Ac, dAc, lambda = 0, P = 1) {
#  if (lambda == 0) return(dAc)
#  else if (length(P) == 1) return(dAc + lambda * P * Ac)
#  else return (dAc + lambda * Ac %*% P)
#}
#pen.Bases <- function(Bases.) {
#  result <- list()
#  for (i in 2:length(Bases.))
#    result[[i - 1]] <- pen.fun(Bases.[[i]])
#  return(result)
#}
pen.fun <- function(fun, order = 2) {
  if (is.basis(fun)) {
    if (order == 2)
      return(bsplinepen(fun, order) / 1e9) 
      #return(fun$nbasis * bsplinepen(fun, order) / 1e9)
    else return(fun$nbasis * bsplinepen(fun, order))
    # return(bsplinepen(fun, order))
  } else return(1)
}
pen. <- function(Ac, dAc, lambda = 0, P) {
  if (lambda == 0) return(dAc)
  else if (length(P[[2]]$P2) == 1) return(dAc + lambda * Ac)
  else {
    JJ <- matrix(1, ncol = ncol(P[[2]]$P2), nrow = nrow(P[[2]]$P2))
    dAc[-1, ] <- dAc[-1, ] + lambda * (P[[1]]$P0 %*% Ac[-1, ]) %*% 
      (P[[2]]$P2 %*% JJ + JJ %*% P[[2]]$P2)/2 + 
      lambda * (P[[1]]$P2 %*% Ac[-1, ]) %*% 
      (P[[2]]$P0 %*% JJ + JJ %*% P[[2]]$P0)/2
    dAc[1, ] <- dAc[1, ] + lambda * Ac[1, ] %*% P[[2]]$P2
    return(dAc)
  }
}
pen.Bases <- function(Bases.) {
  P <- list()
  for (i in 1:length(Bases.)) 
    P[[i]] <- list(P2 = pen.fun(Bases.[[i]]), P0 = pen.fun(Bases.[[i]], 0))
  return(P)
}
