cost <- function(Y, Ac, E, A){
  Y.hat <- pred(Ac, E, A)
  if (is.data.frame(Y)) Y <- Y$Y
  result <- cost.(Y, Y.hat)
  return(result)
}
Error.nn <- function(Y.train, E.train, Y.test, E.test, Bases, A, lambda.) {
  nn.p <- NN.p(Y.train, E.train, Bases, A = A, lambda. = lambda.)
  Y.train. <- pred(nn.p$Ac, E.train, A)
  Y.test. <- pred(nn.p$Ac, E.test, A)
  error <- Error(Y.train, Y.test, Y.train., Y.test.)
  result <- cbind(error, data.frame(lambda = nn.p$lambda, j = nn.p$j))
  return(result)
}
NN.p <- function(Y, X, Bases, A, lambda.){
  I <- length(lambda.)
  if (I > 1) {
    if (is.list(Y)) groups <- divide(Y$PTID, names = c("subtr", "valid"))
    else groups <- divide(1:nrow(Y), names = c("subtr", "valid"))
    list2env(sep(mget(c('Y', 'X')), groups), envir = environment())
    error. <- rep(0, I)
    for (i in 1:I) {
      Ac <- NN.p.(Y.subtr, X.subtr, Bases, A, lambda.[[i]])$Ac
      error.[i] <- cost(Y.valid, Ac, X.valid, A)
    }
    lambda <- lambda.[which.min(error.)]
  } else lambda <- lambda.
  nn.p <- NN.p.(Y, X, Bases, A, lambda)
  return(list(Ac = nn.p$Ac, lambda = lambda, j = nn.p$j))
}