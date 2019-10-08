cost <- function(Y, Ac, E, A){
  Y.hat <- pred(Ac, E, A)
  if ("PTID" %in% colnames(Y)) {
    Y$PTID <- NULL
    Y <- as.matrix(Y)
  }
  result <- cost.(Y, Y.hat)
  return(result)
}
Error.nn <- function(Y.train, E.train, Y.test, E.test, Bases, A, lambda.) {
  fnn.p <- NN(Y.train, E.train, Bases, A = A, lambda. = lambda.)
  Y.train. <- pred(fnn.p$Ac, E.train, A)
  Y.test. <- pred(fnn.p$Ac, E.test, A)
  Y.train$PTID <- NULL
  Y.train <- as.matrix(Y.train)
  Y.test$PTID <- NULL
  Y.test <- as.matrix(Y.test)
  train <- cost.(Y.train, Y.train.)
  test <- cost.(Y.test, Y.test.)
  cor1 <- cor(c(Y.train), c(Y.train.))
  cor2 <- cor(c(unlist(Y.test)), c(unlist(Y.test.)))
  result <- data.frame(train = train, test = test, lambda = fnn.p$lambda, 
                       cor1 = cor1, cor2 = cor2)
  return(result)
}
NN. <- function(Y, E, output, A, A.prime, lambda = 0) {
  lambda <- lambda <- lambda/nrow(Y)
  layers <- forw.prop(E, output$Ac, A)
  grads <- back.prop(Y, output$Ac, layers, A, A.prime)
  output <- adadelta(output, grads, lambda)
  return(output)
}
NN <- function(Y, E, Bases, A, lambda. = 0, Ac.init = NULL){
  Ac <- init(Bases, Ac.init, seed = 1)
  e0 <- init(Bases, para.ub = 0)
  output <- list(Ac = Ac, eg. = e0, ex. = e0)
  A.prime <- lapply(A, Deriv)
  return(NN.p(Y, E, output, A, A.prime, lambda.))
}