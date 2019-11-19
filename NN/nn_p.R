NN.p. <- function(Y, E, Bases, A, lambda) {
  if (is.data.frame(Y)) Y <- Y$Y
  
  Ac <- init(Bases, seed = 1)
  e0 <- init(Bases, para.ub = 0)
  output <- list(Ac = Ac, eg = e0, ex = e0)
  A.prime <- lapply(A, Deriv)
  j <- 0
  cache <- list(k = 0, error. = Inf)
  while (cache$k < 30) { 
    if (j %% 10 == 0) {
      error.g <- sum((Y - pred(output$Ac, E, A))^2)
      error.p <- Error.p(output$Ac, lambda)
      error <- sum(error.g, error.p)
      if (j %% 100 == 0) {
        cat(error.g, error.p, j)
        cat("\n")
      }
      # if (error < (cache$error. * 0.999))
      if (error < (cache$error.))
        cache <- list(Ac. = output$Ac, j. = j, k = 0, error. = error)
      else cache$k <- cache$k + 1
    }
    layers <- forw.prop(E, output$Ac, A)
    grads.g <- back.prop(Y, output$Ac, layers, A, A.prime)
    grads.p <- pen(output$Ac, lambda)
    output <- adadelta(output, grads.g, grads.p)
    j <- j + 1
  }
  cat(j - 1, "\n")
  return(list(Ac = cache$Ac., error = cache$error., j = cache$j.))
}