cost <- function(Y, Ac, X, G, Bases, A, pos = NULL){
  # Bases <- Bf2m(Bases, pos, Y$loc)
  Y.hat <- pred(Ac, X, G, A, Bases = Bases, pos, 
                loc = subset(Y, select = -Y))
  return(cost.(Y$Y, Y.hat))
}
Error.fnn <- function(Y.train, X.train, G.train, Y.test, X.test, G.test, 
                      pos, A, lambda.) {
  fnn.p <- FNN.p(Y.train, X.train, G.train, pos, A = A, lambda. = lambda.)
  Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, pos = pos, 
                   loc = subset(Y.train, select = -Y))
  Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, pos = pos, 
                  loc = subset(Y.test, select = -Y))
  train <- cost.(Y.train$Y, Y.train.)
  test <- cost.(Y.test$Y, Y.test.)
  cor1 <- cor(Y.train$Y, Y.train.)
  cor2 <- cor(Y.test$Y, Y.test.)
  result <- data.frame(train = train, test = test, cor1 = cor1, cor2 = cor2, 
                       lambda = do.call(paste, as.list(fnn.p$lambda)), 
                       j = fnn.p$j)
  # result$lambda = list(fnn.p$lambda)
  return(result)
}
FNN. <- function(Y, X, G, output, Bases, A, A.prime, Pc, PA) {
  layers <- forw.prop(output$Ac, X, G, Bases, A, Y[, c("PTID", "loc")])
  grads <- back.prop(Y, output$Ac, layers, Bases, A, A.prime)
  output <- adadelta(output, grads, Pc, PA)
  return(output)
}