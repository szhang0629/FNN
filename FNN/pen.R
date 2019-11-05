RMS <- function(g, epsilon = 1e-6){
  return(sqrt(g + epsilon))
}
adadelta <- function(output, grads, lambda = 0, P = NULL, P0 = NULL, 
                     rho = 0.95) {
  for (i in 1:length(grads)) {
    grads. <- pen.(output$Ac[[i]], grads[[i]], lambda, P[[i]])
    # grads. <- pen.(Ac[[i]], grads[[i]], lambda, P[i:(i + 1)])
    output$eg.[[i]] <- rho * output$eg.[[i]] + (1 - rho) * (grads.)^2
    dAc <- -RMS(output$ex.[[i]]) / RMS(output$eg.[[i]]) * grads.
    # dAc <- - 10 * grads.
    output$ex.[[i]] <- rho * output$ex.[[i]] + (1 - rho) * dAc^2
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
pen. <- function(Ac, dAc, lambda = 0, P = 1) {
  if (lambda == 0) return(dAc)
  else if (length(P) < 2) return(dAc + lambda * Ac)
  else {
    n.w <- nrow(P$P2)
    n <- nrow(dAc)
    row.A <- (n - n.w + 1):n
    row.c <- 1:(n - n.w)
    dc. <- subs(Ac, row.c) %*% P$Pc * 2
    dA. <- 2 * P$P0 %*% subs(Ac, row.A) %*% P$Q0 + 
      2 * P$P2 %*% subs(Ac, row.A) %*% P$Q2
    return(dAc + lambda * rbind(dc., dA.))
  }
}
Error.p <- function(Ac, P, lambda) {
  result <- rep(0, length(Ac))
  for (i in 1:length(Ac))
    result[i] <- sum(Error.p.(Ac[[i]], P[[i]])*lambda)
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
      # return(bb$nbasis*bsplinepen(bb, order)/1e7)
      return(bsplinepen(bb, order)/1e7)
    else 
      return(bsplinepen(bb, order))
    # return(bb$nbasis*bsplinepen(bb, order))
  } else return(1)
}