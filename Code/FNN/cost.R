cost <- function(Y, Ac, X, G, Bases, A, pos = NULL, loc = NULL){
  Y.hat <- pred(Ac, X, G, A, Bases = Bases, pos, loc)
  if (is.list(Y)) return(cost.(Y$Y, Y.hat))
  else return(cost.(Y, Y.hat))
}
FNN.p <- function(Y, X, G, Bases, A, lambda., pos, loc = NULL){
  I <- length(lambda.)
  if (I > 1) {
    if (!is.null(rownames(Y))) 
      groups <- divide(rownames(Y), type = "name", names = c("subtr", "valid"))
    else groups <- divide(1:nrow(Y), names = c("subtr", "valid"))
    list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
    valid. <- rep(0, I)
    subtr. <- rep(0, I)
    j. <- rep(0, I)
    for (i in 1:I) {
      fnn.p <- FNN.p.(Y.subtr, X.subtr, G.subtr, Bases, A, lambda.[[i]], pos, 
                      loc)
      if (is.null(loc)) {
        valid.[i] <- cost(Y.valid[, 1, drop = FALSE], fnn.p$Ac, X.valid, 
                          G.valid, Bases, A, pos, Y.valid[, -1, drop = FALSE])
        subtr.[i] <- cost(Y.subtr[, 1, drop = FALSE], fnn.p$Ac, X.subtr, 
                          G.subtr, Bases, A, pos, Y.subtr[, -1, drop = FALSE])
      } else {
        valid.[i] <- cost(Y.valid, fnn.p$Ac, X.valid, 
                          G.valid, Bases, A, pos, loc)
        subtr.[i] <- cost(Y.subtr, fnn.p$Ac, X.subtr, 
                          G.subtr, Bases, A, pos, loc)
      }
      j.[i] <- fnn.p$j
    }
    print(data.frame(j., subtr., valid.))
    lambda <- lambda.[[which.min(valid.)]]
  } else lambda <- lambda.[[1]]
  fnn.p <- FNN.p.(Y, X, G, Bases, A, lambda, pos, loc)
  return(list(Ac = fnn.p$Ac, lambda = lambda, j = fnn.p$j))
}