flm1 <- function(X, bb0, pos, para, Z = NULL) {
  B0 <- eval.basis(pos, bb0)
  X. <- X %*% B0
  if (is.null(Z)) Y <- rep.row(para$theta, nrow(X))
  else Y <- cbind(1, Z) %*% para$theta
  Y <- Y + X. %*% para$W
  return(Y)
}
solve.flm1 <- function(Y, X, bb0, pos, Z = NULL, lambda = 0) {
  P <- bsplinepen(bb0)
  # coef <- ncol(P) / sum(diag(P))
  # P <- coef * P
  n <- nrow(X)
  B0 <- eval.basis(pos, bb0)
  D <- X %*% B0
  if (is.null(Z)) C <- rep.row(1, n)
  else C <- cbind(1, Z)
  M0 <- ginv(t(C) %*% C) %*% t(C)
  M <- C %*% M0
  L <- t(D) - t(D) %*% M
  W <- ginv(L %*% D + lambda * P) %*% L %*% Y
  theta <- M0 %*% (Y - D %*% W)
  return(list(W = W, theta = theta))
}
Error.flm1 <- function(Y.train, G.train, Y.test, G.test, bb0, pos, 
                       X.train = NULL, X.test = NULL, 
                       lambda. = c(0, 1e-3, 1e-2, 1e-1)){
  nbreaks <- 10
  folds <- cut(seq(1,nrow(Y.train)),breaks = nbreaks, labels = FALSE)
  n <-length(lambda.)
  results <- rep(0, n)
  for (j in 1:n) {
    for (i in 1:nbreaks) {
      validIndexes <- which(folds == i, arr.ind = TRUE)
      G.valid <- G.train[validIndexes, ]
      G.subtr <- G.train[-validIndexes, ]
      X.valid <- X.train[validIndexes, ,drop = FALSE]
      X.subtr <- X.train[-validIndexes, ,drop = FALSE]
      Y.valid <- Y.train[validIndexes, ,drop = FALSE]
      Y.subtr <- Y.train[-validIndexes, ,drop = FALSE]
      para <- solve.flm1(Y.subtr, G.subtr, bb0, pos, X.subtr, lambda.[j])
      output <- cost.(Y.valid, flm1(G.valid, bb0, pos, para, X.valid))
      results[j] <- output/nbreaks + results[j]
    }
  }
  index <- which.min(results)
  lambda <- lambda.[index]
  para <- solve.flm1(Y.train, G.train, bb0, pos, X.train, lambda)
  Y.train. <- flm1(G.train, bb0, pos, para, X.train)
  Y.test. <- flm1(G.test, bb0, pos, para, X.test)
  train <- cost.(Y.train, Y.train.)
  test <- cost.(Y.test, Y.test.)
  cor1 <- cor(c(unlist(Y.train)), c(unlist(Y.train.)))
  cor2 <- cor(c(unlist(Y.test)), c(unlist(Y.test.)))
  result <- data.frame(train = train, test = test, lambda = lambda, 
                       cor1 = cor1, cor2 = cor2)
  return(result)
}