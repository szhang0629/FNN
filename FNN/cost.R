cost <- function(Y, Ac, X, G, Bases, A, pos = NULL, int = T){
  # Bases <- Bf2m(Bases, pos, Y$loc)
  Y.hat <- pred(Ac, X, G, A, Bases = Bases, pos, int = int, 
                loc = subset(Y, select = -Y))
  return(cost.(Y$Y, Y.hat))
}
Error.fnn <- function(Y.train, X.train, G.train, Y.test, X.test, G.test, 
                      pos, Bases, A, lambda.) {
  fnn.p <- FNN.p(Y.train, X.train, G.train, pos, A = A, Bases, lambda. = lambda.)
  Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, 
                   loc = subset(Y.train, select = -Y))
  Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, 
                  loc = subset(Y.test, select = -Y))
  train <- cost.(Y.train$Y, Y.train.)
  test <- cost.(Y.test$Y, Y.test.)
  cor1 <- cor(Y.train$Y, Y.train.)
  cor2 <- cor(Y.test$Y, Y.test.)
  result <- data.frame(train = train, test = test, cor1 = cor1, cor2 = cor2, 
                       lambda = fnn.p$lambda, j = fnn.p$j)
  return(result)
}
FNN.p <- function(Y, X, G, pos, A, Bases, lambda.){
  P <- pen.Bases(Bases)
  int <- int.fun(Bases)
  # Bases <- Bf2m(Bases., pos, subset(Y, select = -c(PTID, Y)))
  Ac <- init(Bases, seed = 1, X.col = ncol(X))
  e0 <- init(Bases, X.col = ncol(X), para.ub = 0)
  output <- list(Ac = Ac, eg. = e0, ex. = e0)
  A.prime <- lapply(A, Deriv)
  I <- length(lambda.)
  if (I > 1) {
    groups <- divide(Y$PTID, type = "name", names = c("subtr", "valid"))
    list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
    valid. <- rep(0, I)
    subtr. <- rep(0, I)
    j. <- rep(0, I)
    for (i in 1:I) {
      fnn.p <- FNN.p.(Y.subtr, X.subtr, G.subtr, output, Bases, pos, A, A.prime,
                   lambda.[[i]], P, int)
      valid.[i] <- cost(Y.valid, fnn.p$Ac, X.valid, G.valid, Bases, A, pos, int)
      subtr.[i] <- cost(Y.subtr, fnn.p$Ac, X.subtr, G.subtr, Bases, A, pos, int)
      j.[i] <- fnn.p$j
    }
    print(data.frame(j., subtr., valid.))
    lambda <- lambda.[which.min(valid.)]
  } else lambda <- lambda.
  fnn.p <- FNN.p.(Y, X, G, output, Bases, pos, A, A.prime, lambda, P, int)
  return(list(Ac = fnn.p$Ac, lambda = lambda, j = fnn.p$j))
}
