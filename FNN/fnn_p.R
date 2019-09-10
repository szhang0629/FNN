FNN.p <- function(Y, E, output, Bases, A, A.prime, lambda., P, int){
  I <- length(lambda.)
  if (I > 1) {
    groups <- divide(Y, names = c("subtr", "valid"))
    if (is.list(Y))
      list2env(split(mget(c('Y', 'E', 'Bases')), groups), envir = environment())
    else {
      list2env(split(mget(c('Y', 'E')), groups), envir = environment())
      Bases.subtr <- Bases
      Bases.valid <- Bases
    }
    error. <- rep(0, I)
    for (i in 1:I) {
      Ac <- FNN.p.(Y.subtr, E.subtr, output, Bases.subtr, A, A.prime,
                     lambda.[[i]], P, int)$Ac
      error.[i] <- cost(Y.valid, Ac, E.valid, Bases.valid, A, int = int)
    }
    lambda <- lambda.[which.min(error.)]
  } else lambda <- lambda.
  Ac <- FNN.p.(Y, E, output, Bases, A, A.prime, lambda, P, int)$Ac
  return(list(Ac = Ac, lambda = lambda))
}
FNN.p. <- function(Y, E, output, Bases, A, A.prime, lambda, P, int) {
  j <- 0
  cache <- list(k = 0, error. = Inf)
  while (cache$k < 30) { 
    if (j %% 10 == 0) {
      error <- cost(Y, output$Ac, E, Bases, A, int = int)
      # print(error)
      if (is.nan(error)) break
      if (error < (cache$error. * 0.999))
      # if (error < (cache$error.))
        cache <- list(Ac. = output$Ac, j. = j, k = 0, error. = error)
      else cache$k <- cache$k + 1
    }
    output <- FNN.(Y, E, output, Bases, A, A.prime, lambda, P, int)
    j <- j + 1
  }
  cat(j-1, "\n")
  return(list(Ac = cache$Ac., error = cache$error.))
}
