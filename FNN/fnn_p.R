FNN.p. <- function(Y, X, G, output, Bases, pos, A, A.prime, lambda, P, int) {
  idy <- idx.names(subset(Y, select = -c(PTID, Y)))
  Bases <- Bf2m(Bases, pos, subs(subset(Y, select = -c(PTID, Y)), 
                                 match(1:max(idy), idy)))
  idx <- idx.names(Y$PTID, rownames(G))
  # Y <- data.frame(Y = Y$Y, idx = idx, idy = idy)
  id <- data.frame(idx = idx, idy = idy)
  Y <- Y$Y
  
  j <- 0
  lambda <- lambda/nrow(id)
  cache <- list(k = 0, error. = Inf)
  while (cache$k < 30) { 
    if (j %% 10 == 0) {
      error <- cost.(Y, pred.(output$Ac, X, G, A, Bases, int, id))
      if (j %% 1000 == 0) {
        cat(paste(error, j))
        cat("\n")
      }
      if (is.nan(error))
        break
      if (error < (cache$error. * 0.999))
      # if (error < (cache$error.))
        cache <- list(Ac. = output$Ac, j. = j, k = 0, error. = error)
      else cache$k <- cache$k + 1
    }
    layers <- forw.prop(output$Ac, X, G, Bases, A, int = int, id)
    grads <- back.prop(Y, output$Ac, layers, Bases, A, A.prime, int = int, id)
    output <- adadelta(output, grads, lambda, P)
    j <- j + 1
  }
  cat(j - 1, "\n")
  return(list(Ac = cache$Ac., error = cache$error., j = cache$j.))
}