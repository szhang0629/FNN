pred.flm1 <- function(G, bb0, pos, para, Z = NULL) {
  B0 <- eval.basis(pos, bb0)
  D <- G %*% B0
  X <- cbind(rep(1, nrow(G)), Z)
  Y <-  X %*% para$theta + D %*% para$W
  return(Y)
}
solve.flm1 <- function(Y, G, bb0, pos, Z = NULL, lambda = 0) {
  if (is.list(Y)) Y <- Y$Y
  P <- bsplinepen(bb0)
  # coef <- ncol(P) / sum(diag(P))
  # P <- coef * P
  B0 <- eval.basis(pos, bb0)
  D <- G %*% B0
  X <- cbind(rep(1, nrow(G)), Z)
  M <- ginv(t(X) %*% X) %*% t(X)
  L <- t(D) - t(D) %*% X %*% M
  W <- ginv(L %*% D + lambda * P) %*% L %*% Y
  theta <- M %*% (Y - D %*% W)
  return(list(W = W, theta = theta))
}
flm1 <- function(Y.train, G.train, bb0, pos, X.train = NULL, 
                 lambda. = c(1e-3, 1e-2, 1e-1, 1)){
  n <- length(lambda.)
  if (n > 1) {
    nbreaks <- 10
    results <- rep(0, n)
    if (!is.null(rownames(Y.train))) groups <- divide.cv(rownames(Y.train), 
                                                         n = nbreaks)
    else groups <- divide.cv(1:nrow(Y.train), n = nbreaks)
    for (j in 1:n) {
      for (i in 1:nbreaks) {
        G.valid <- G.train[groups[[i]], ]
        G.subtr <- G.train[-groups[[i]], ]
        X.valid <- X.train[groups[[i]], ,drop = FALSE]
        X.subtr <- X.train[-groups[[i]], ,drop = FALSE]
        Y.valid <- Y.train[groups[[i]], ,drop = FALSE]
        Y.subtr <- Y.train[-groups[[i]], ,drop = FALSE]
        para <- solve.flm1(Y.subtr, G.subtr, bb0, pos, X.subtr, lambda.[j])
        output <- cost.(c(Y.valid), c(pred.flm1(G.valid, bb0, pos, para, X.valid)))
        results[j] <- output/nbreaks + results[j]
      }
    }
    print(data.frame(lambda., results))
    index <- which.min(results)
    lambda <- lambda.[index]
  } else lambda <- lambda.
  para <- solve.flm1(Y.train, G.train, bb0, pos, X.train, lambda)
  return(list(para = para, lambda = lambda))
}