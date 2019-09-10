cost <- function(Y, Ac, E, Bases, A, pos = NULL, loc = NULL, X = NULL, int = T){
  Bases <- Bf2m(Bases, pos, loc)
  Y.hat <- pred(Ac, E, A, Bases = Bases, int = int)
  cost <- cost.(Y, Y.hat)
  return(cost)
}
Error.fnn <- function(Y.train, E.train, loc.train, Y.test, E.test, loc.test,
                      pos, Bases, A, lambda.) {
  fnn.p <- FNN(Y.train, E.train, loc.train, pos, Bases, A = A, 
               lambda. = lambda.)
  Y.train. <- pred(fnn.p$Ac, E.train, A, Bases, pos, loc.train)
  Y.test. <- pred(fnn.p$Ac, E.test, A, Bases, pos, loc.test)
  train <- cost.(Y.train, Y.train.)
  test <- cost.(Y.test, Y.test.)
  cor1 <- cor(c(unlist(Y.train)), c(unlist(Y.train.)))
  cor2 <- cor(c(unlist(Y.test)), c(unlist(Y.test.)))
  result <- data.frame(train = train, test = test, lambda = fnn.p$lambda, 
                       cor1 = cor1, cor2 = cor2)
  return(result)
}