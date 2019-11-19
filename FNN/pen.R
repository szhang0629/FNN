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
  for (i in 1:length(Ac)) 
    grads[[i]] <- pen.(Ac[[i]], P[[i]], lambda)
  return(grads)
}
pen. <- function(Ac, P, lambda) {
  if (lambda == 0) return(0)
  else if (length(P) < 2) return(lambda * Ac)
  else {
    n.w <- nrow(P$P2)
    n <- nrow(Ac)
    row.A <- (n - n.w + 1):n
    row.c <- 1:(n - n.w)
    dc. <- subs(Ac, row.c) %*% P$Pc * 2
    dA. <- 2 * P$P0 %*% subs(Ac, row.A) %*% P$Q0 + 
      2 * P$P2 %*% subs(Ac, row.A) %*% P$Q2
    return(lambda * rbind(dc., dA.))
  }
}
Error.p <- function(Ac, P, lambda) {
  result <- rep(0, length(Ac))
  for (i in 1:length(Ac)) {
    if (lambda == 0) result[i] <- 0
    else result[i] <- sum(Error.p.(Ac[[i]], P[[i]])*lambda)
  }
  return(result)
}
Error.p. <- function(Ac, P) {
  n.w <- nrow(P$P2)
  n <- nrow(Ac)
  row.A <- (n - n.w + 1):n
  row.c <- 1:(n - n.w)
  PA1 <- sum(diag(subs(Ac, row.A) %*% P$Q2 %*% t(subs(Ac, row.A)) %*% P$P2))
  PA0 <- sum(diag(P$Q0 %*% t(subs(Ac, row.A)) %*% P$P0 %*% subs(Ac, row.A)))
  Pc <- sum(diag(subs(Ac, row.c) %*% P$Pc %*% t(subs(Ac, row.c))))
  return(c(PA1, PA0, Pc))
}
pen.fun <- function(bb, order = 2) {
  if (is.basis(bb)) {
    if (order == 2)
      return(bb$nbasis*bsplinepen(bb, order)/1e7)
      # return(bsplinepen(bb, order)/1e7)
    else 
      # return(bsplinepen(bb, order))
      return(bb$nbasis*bsplinepen(bb, order))
  } else if (is.list(bb)) {
    if ("bb0" %in% names(bb)) {
      if (order == 2) 
        return(as.matrix(bdiag(bsplinepen(bb[[1]], 2), bsplinepen(bb[[2]], 2))))
      else 
        return(as.matrix(bdiag(bsplinepen(bb[[1]], 0), bsplinepen(bb[[2]], 0))))
    }
    if (order == 2) 
      return(bsplinepen(bb[[2]], 2) %x% bsplinepen(bb[[1]], 0)/1e7 +
               bsplinepen(bb[[2]], 0) %x% bsplinepen(bb[[1]], 2)/1e7)
    else 
      return(bsplinepen(bb[[1]], order) %x% bsplinepen(bb[[2]], order))
  } else return(1)
}