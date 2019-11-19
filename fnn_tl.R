FNN.tl <- function(Y, X, G, Bases, A, lambda, pos, Ac.) {
  # idy <- idx.names(subset(Y, select = -c(PTID, Y)))
  # if (is.null(Y$loc)) {
  #   fBases <- c(rep(list(.bb), length(lambda)), list(matrix(1, 1)))
    # Bases <- Bf2m(fBases, pos, NULL)
  #   id <- NULL
  # } else {
  #   fBases <- rep(list(.bb), length(lambda))
  #   Bases <- Bf2m(fBases, pos, 
  #                 subs(subset(Y, select = -c(PTID, Y)), match(1:max(idy), idy)))
  #   idx <- idx.names(Y$PTID, rownames(G))
  #   id <- data.frame(idx = idx, idy = idy)
  # }
  # Y <- Y$Y
  
  P <- pen.Bases(Bases)
  int <- lapply(Bases, is.basis)
  Bases <- Bf2m(Bases, pos, NULL)
  e0 <- init(Bases, X.col = ncol(X), para.ub = 0)
  # e0 <- rep(list(0), length(A))
  output <- list(Ac = Ac., eg = e0, ex = e0)
  A.prime <- lapply(A, Deriv)
  
  j <- 0
  cache <- list(k = 0, error. = Inf)
  while (cache$k < 30) { 
    if (j %% 10 == 0) {
      error.g <- sum((Y - pred(output$Ac, X, G, A, Bases, int = int))^2)
      error.p <- Error.p(output$Ac, P, lambda)
      error <- sum(error.g, error.p)
      if (j %% 200 == 0) {
        cat(error.g, error.p, j)
        cat("\n")
      }
      if (is.nan(error))
        break
      # if (error < (cache$error. * 0.999))
      if (error < (cache$error.))
        cache <- list(Ac. = output$Ac, j. = j, k = 0, error. = error)
      else cache$k <- cache$k + 1
    }
    layers <- forw.prop(output$Ac, X, G, Bases, A, int = int)
    grads.g <- back.prop(Y, output$Ac, layers, Bases, A, A.prime, int = int)
    grads.g[[1]][grads.g[[1]] != 0] <- 0
    grads.p <- pen(output$Ac, P, lambda)
    grads.p[[1]][grads.p[[1]] != 0] <- 0
    output <- adadelta(output, grads.g, grads.p)
    j <- j + 1
  }
  return(list(Ac = cache$Ac., error = cache$error., j = cache$j.))
}