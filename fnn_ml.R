FNN.ml <- function(Y1, X1, G1, Y2, X2, G2, Bases, pos, A, lambda) {
  # idy <- idx.names(subset(Y1, select = -c(PTID, Y)))
  # if (is.null(Y1$loc)) {
  #   fBases <- c(rep(list(.bb), length(lambda)), list(matrix(1, 1)))
  #   Bases1 <- Bf2m(fBases, pos, NULL)
  #   id1 <- NULL
  # } else {
  #   fBases <- rep(list(.bb), length(lambda))
  #   Bases1 <- Bf2m(fBases, pos, 
  #                 subs(subset(Y1, select = -c(PTID, Y)), match(1:max(idy), idy)))
  #   idx <- idx.names(Y1$PTID, rownames(G))
  #   id1 <- data.frame(idx = idx, idy = idy)
  # }
  # Y1 <- Y1$Y
  # 
  # idy <- idx.names(subset(Y2, select = -c(PTID, Y)))
  # if (is.null(Y2$loc)) {
  #   fBases <- c(rep(list(.bb), length(lambda)), list(matrix(1, 1)))
  #   Bases2 <- Bf2m(fBases, pos, NULL)
  #   id2 <- NULL
  # } else {
  #   fBases <- rep(list(.bb), length(lambda))
  #   Bases2 <- Bf2m(fBases, pos, 
  #                 subs(subset(Y2, select = -c(PTID, Y)), match(1:max(idy), idy)))
  #   idx <- idx.names(Y2$PTID, rownames(G))
  #   id2 <- data.frame(idx = idx, idy = idy)
  # }
  # Y2 <- Y2$Y
  
  P <- pen.Bases(Bases)
  int <- lapply(Bases, is.basis)
  Bases1 <- Bf2m(Bases, pos, NULL)
  Bases2 <- Bf2m(Bases, pos, NULL)
  e1 <- init(Bases1, X.col = ncol(X1), para.ub = 0)
  output1 <- list(Ac = init(Bases1, seed = 1, X.col = ncol(X1)), 
                  eg = e1, ex = e1)
  e2 <- init(Bases2, X.col = ncol(X2), para.ub = 0)
  output2 <- list(Ac = init(Bases2, seed = 1, X.col = ncol(X2)), 
                  eg = e2, ex = e2)
  
  A.prime <- lapply(A, Deriv)
  j <- 0
  k <- 0
  cache <- list(error. = Inf)
  while (k < 30) { 
    if (j %% 10 == 0) {
      error1.g <- sum((Y1 - pred(output1$Ac, X1, G1, A, Bases1, int = int))^2)
      error1.p <- Error.p(output1$Ac, P, lambda)
      error2.g <- sum((Y2 - pred(output2$Ac, X2, G2, A, Bases2, int = int))^2)
      error2.p <- Error.p(output2$Ac, P, lambda)
      error <- sum(error1.g, error1.p, error2.g, error2.p)
      if (j %% 1000 == 0) {
        cat(error1.g, error1.p, error2.g, error2.p, j)
        cat("\n")
      }
      if (is.nan(error))
        break
      # if (error < (cache$error. * 0.999))
      if (error < (cache$error.))
        cache <- list(Ac1 = output1$Ac, Ac2 = output2$Ac, j. = j, error. = error)
      else k <- k + 1
    }
    layers1 <- forw.prop(output1$Ac, X1, G1, Bases1, A, int = int)
    layers2 <- forw.prop(output2$Ac, X2, G2, Bases2, A, int = int)
    grads1.g <- back.prop(Y1, output1$Ac, layers1, Bases1, A, A.prime, int = int)
    grads1.p <- pen(output1$Ac, P, lambda)
    grads2.g <- back.prop(Y2, output2$Ac, layers2, Bases2, A, A.prime, int = int)
    grads2.p <- pen(output2$Ac, P, lambda)
    grads0.g <- grads1.g[[1]] + grads2.g[[1]]
    grads0.p <- grads1.p[[1]] + grads2.p[[1]]
    grads1.g[[1]] <- grads0.g
    grads2.g[[1]] <- grads0.g
    grads1.p[[1]] <- grads0.p
    grads2.p[[1]] <- grads0.p
    output1 <- adadelta(output1, grads1.g, grads1.p)
    output2 <- adadelta(output2, grads2.g, grads2.p)
    j <- j + 1
  }
  return(cache)
}