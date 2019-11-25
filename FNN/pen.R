RMS <- function(g, epsilon = 1e-6){
  return(sqrt(g + epsilon))
}
adadelta <- function(output, grads.g, grads.p, rho = 0.95) {
  for (i in 1:length(grads.g)) {
    grads. <- grads.g[[i]] + grads.p[[i]]
    output$eg[[i]] <- rho * output$eg[[i]] + (1 - rho) * (grads.)^2
    dAc <- -RMS(output$ex[[i]]) / RMS(output$eg[[i]]) * grads.
    # dAc <- -0.03 * grads.
    output$ex[[i]] <- rho * output$ex[[i]] + (1 - rho) * dAc^2
    output$Ac[[i]] <- output$Ac[[i]] + dAc
  }
  return(output)
}
## ------matrix smooth------
pen.Bases <- function(Bases.) {
  result <- list()
  for (i in 1:(length(Bases.))) {
    P2 <- pen.fun(Bases.[[i]])
    P0 <- pen.fun(Bases.[[i]], 0)
    result[[i]] <- list(P2 = P2, P0 = P0)
  }
  return(result)
}
pen <- function(Ac, P, lambda) {
  grads <- list()
  if (length(lambda) == 1) 
    lambda <- rep(lambda, length(P))
  for (i in 1:length(Ac)) 
    grads[[i]] <- pen.(Ac[[i]], P[i:(i + 1)], lambda[i:(i + 1)])
  return(grads)
}
pen. <- function(Ac, P, lambda) {
  n.w <- nrow(P[[1]]$P0)
  n <- nrow(Ac)
  row.A <- (n - n.w + 1):n
  row.c <- 1:(n - n.w)
  dc. <- 2 * lambda[2] * subs(Ac, row.c) %*% P[[2]]$P2
  dA. <- 2 * lambda[2] * P[[1]]$P0 %*% subs(Ac, row.A) %*% P[[2]]$P2 +
    2 * lambda[1] * P[[1]]$P2 %*% subs(Ac, row.A) %*% P[[2]]$P0
  return(rbind(dc., dA.))
  # return(2 * lambda[2] * Ac %*% P[[2]]$P2)
}
Error.p <- function(Ac, P, lambda) {
  result <- rep(0, length(Ac))
  if (length(lambda) == 1) 
    lambda <- rep(lambda, length(P))
  for (i in 1:length(Ac)) 
    result[i] <- sum(Error.p.(Ac[[i]], P[i:(i + 1)], lambda[i:(i + 1)]))
  return(result)
}
Error.p. <- function(Ac, P, lambda) {
  n.w <- nrow(P[[1]]$P0)
  n <- nrow(Ac)
  row.A <- (n - n.w + 1):n
  row.c <- 1:(n - n.w)
  PA1 <- lambda[1] * sum(diag(P[[1]]$P0 %*% subs(Ac, row.A) %*% P[[2]]$P2 %*%
                                t(subs(Ac, row.A))))
  PA0 <- lambda[2] * sum(diag(P[[2]]$P0 %*% t(subs(Ac, row.A)) %*%
                                P[[1]]$P2 %*% subs(Ac, row.A)))
  Pc <- lambda[2] * sum(diag(subs(Ac, row.c) %*% P[[2]]$P2
                             %*% t(subs(Ac, row.c))))
  return(c(PA1, PA0, Pc))
  # return(lambda[2] * Ac %*% P[[2]]$P2 %*% t(Ac))
}
pen.fun <- function(bb, order = 2) {
  if (is.basis(bb)) {
    if (order == 2)
      return(c.fun(bb)*bsplinepen(bb, order)/1e6)
    else 
      return(c.fun(bb)*bsplinepen(bb, order))
  } else if (is.list(bb)) {
    if ("bb0" %in% names(bb)) {
      if (order == 2) 
        return(as.matrix(bdiag(pen.fun(bb[[1]], 2), pen.fun(bb[[2]], 2))))
      else 
        return(as.matrix(bdiag(pen.fun(bb[[1]], 0), pen.fun(bb[[2]], 0))))
    }
    if (order == 2) 
      return(pen.fun(bb[[2]], 2) %x% pen.fun(bb[[1]], 0) +
               pen.fun(bb[[2]], 0) %x% pen.fun(bb[[1]], 2))
    else 
      return(pen.fun(bb[[1]], order) %x% pen.fun(bb[[2]], order))
  } else return(1)
}