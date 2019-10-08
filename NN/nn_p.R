NN.p <- function(Y, X, output, A, A.prime, lambda.){
  I <- length(lambda.)
  if (I > 1) {
    groups <- divide(Y$PTID, names = c("subtr", "valid"))
    list2env(split(mget(c('Y', 'X')), groups), envir = environment())
    error. <- rep(0, I)
    for (i in 1:I) {
      Ac <- NN.p.(Y.subtr, X.subtr, output, A, A.prime, lambda.[[i]])$Ac
      error.[i] <- cost(Y.valid, Ac, X.valid, A)
    }
    lambda <- lambda.[which.min(error.)]
  } else lambda <- lambda.
  Ac <- NN.p.(Y, X, output, A, A.prime, lambda)$Ac
  return(list(Ac = Ac, lambda = lambda))
}
NN.p. <- function(Y, E, output, A, A.prime, lambda) {
  Y$PTID <- NULL
  Y <- as.matrix(Y)
  j <- 0
  cache <- list(k = 0, error. = Inf)
  while (cache$k < 30) { 
    if (j %% 10 == 0) {
      error <- cost(Y, output$Ac, E, A)
      # print(error)
      if (error < (cache$error. * 0.999))
      # if (error < (cache$error.))
        cache <- list(Ac. = output$Ac, j. = j, k = 0, error. = error)
      else cache$k <- cache$k + 1
    }
    output <- NN.(Y, E, output, A, A.prime, lambda)
    j <- j + 1
  }
  cat(j-1, "\n")
  return(list(Ac = cache$Ac., error = cache$error.))
}