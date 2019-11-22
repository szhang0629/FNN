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
  for (i in 1:(length(Bases.) - 1)) {
    P2 <- pen.fun(Bases.[[i + 1]])
    P0 <- pen.fun(Bases.[[i + 1]], 0)
    P2. <- pen.fun(Bases.[[i]])
    P0. <- pen.fun(Bases.[[i]], 0)
    result[[i]] <- list(P2 = P2., P0 = P0., Q2 = P0, Q0 = P2, Pc = P2)
  }
  return(result)
}
pen <- function(Ac, P, lambda) {
  grads <- list()
  if (length(lambda) == 1) {
    for (i in 1:length(Ac)) 
      grads[[i]] <- pen.(Ac[[i]], P[[i]], c(lambda, lambda))
  } else {
    for (i in 1:length(Ac)) 
      grads[[i]] <- pen.(Ac[[i]], P[[i]], lambda[i:(i + 1)])
  }
  
  return(grads)
}
# pen. <- function(Ac, lambda) {
#   n.w <- nrow(.P2)
#   n <- nrow(Ac)
#   row.A <- (n - n.w + 1):n
#   row.c <- 1:(n - n.w)
#   dc. <- subs(Ac, row.c) %*% .P2*lambda[2]*2
#   dA. <- .P0 %*% subs(Ac, row.A) %*% .P2 * lambda[2]*2 +
#     .P2 %*% subs(Ac, row.A) %*% .P0 * lambda[1]*2
#   dAc. <- rbind(dc., dA.)
#   return(dAc.)
# }
pen. <- function(Ac, P, lambda) {
  if (length(P) < 2) return(lambda * Ac)
  else {
    n.w <- nrow(P$P2)
    n <- nrow(Ac)
    row.A <- (n - n.w + 1):n
    row.c <- 1:(n - n.w)
    dc. <- subs(Ac, row.c) %*% P$Pc * lambda[2] * 2
    dA. <- 2 * lambda[2] * P$P0 %*% subs(Ac, row.A) %*% P$Q0 + 
      2 * lambda[1] * P$P2 %*% subs(Ac, row.A) %*% P$Q2
    return(rbind(dc., dA.))
  }
}
Error.p <- function(Ac, P, lambda) {
  result <- rep(0, length(Ac))
  for (i in 1:length(Ac)) {
    if (length(lambda) == 1)
      result[i] <- sum(Error.p.(Ac[[i]], P[[i]], c(lambda, lambda)))
    else
      result[i] <- sum(Error.p.(Ac[[i]], P[[i]], lambda[i:(i +1)]))
  }
  return(result)
}
Error.p. <- function(Ac, P, lambda) {
  n.w <- nrow(P$P2)
  n <- nrow(Ac)
  row.A <- (n - n.w + 1):n
  row.c <- 1:(n - n.w)
  PA1 <- lambda[1] * sum(diag(subs(Ac, row.A) %*% P$Q2 %*% t(subs(Ac, row.A)) %*% P$P2))
  PA0 <- lambda[2] * sum(diag(P$Q0 %*% t(subs(Ac, row.A)) %*% 
                                P$P0 %*% subs(Ac, row.A)))
  Pc <- lambda[2] * sum(diag(subs(Ac, row.c) %*% P$Pc %*% t(subs(Ac, row.c))))
  return(c(PA1, PA0, Pc))
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
        return(as.matrix(bdiag(bsplinepen(bb[[1]], 2), bsplinepen(bb[[2]], 2))))
      else 
        return(as.matrix(bdiag(bsplinepen(bb[[1]], 0), bsplinepen(bb[[2]], 0))))
    }
    if (order == 2) 
      return(bsplinepen(bb[[2]], 2) %x% bsplinepen(bb[[1]], 0)/1e6 +
               bsplinepen(bb[[2]], 0) %x% bsplinepen(bb[[1]], 2)/1e6)
    else 
      return(bsplinepen(bb[[1]], order) %x% bsplinepen(bb[[2]], order))
  } else return(1)
}