RMS <- function(g, epsilon = 1e-6){
  return(sqrt(g + epsilon))
}
adadelta <- function(output, grads, Pc, PA, rho = 0.95) {
  eg <- list()
  ex <- list()
  Ac <- output$Ac
  eg. <- output$eg.
  ex. <- output$ex.
  for (i in 1:length(grads)) {
    # grads. <- pen.(Ac[[i]], grads[[i]], lambda, P[[i]])
    grads. <- pen.(Ac[[i]], grads[[i]], Pc[[i]], PA[[i]])
    eg[[i]] <- rho * eg.[[i]] + (1 - rho) * (grads.)^2
    dAc <- -RMS(ex.[[i]]) / RMS(eg[[i]]) * grads.
    ex[[i]] <- rho * ex.[[i]] + (1 - rho) * dAc^2
    Ac[[i]] <- Ac[[i]] + dAc
  }
  return(list(Ac = Ac, eg. = eg, ex. = ex))
}
pen. <- function(Ac, dAc, Pc, PA) {
  # if (lambda[[2]] == 0) return(dAc)
  # # else if (length(P) < 2) return(dAc + lambda * Ac)
  # else {
    n <- nrow(dAc)
    row.A <- (n - nrow(Pc) + 1):n
    row.c <- 1:(n - nrow(Pc))
    dc. <- subs(Ac, row.c) %*% Pc
    dA. <- matrix(as.vector(subs(Ac, row.A)) %*% PA, nrow = nrow(Pc))
    return(dAc + rbind(dc., dA.))
  # }
}