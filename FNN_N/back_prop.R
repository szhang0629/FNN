back.prop <- function(Y, Ac, layers, Bases, A, A.prime){
  grads <- list()
  n <- length(layers)
  Y.hat <- layers[[n]]$G
  M <- c(Y.hat - Y$Y) * Bases[[length(Bases)]][Y$loc, ]/length(Y.hat)
  # odr <- as.character(rownames(layers[[1]]$G))
  # M. <- outer(odr, Y$PTID, "==") * 1
  M. <- outer(1:max(Y$PTID), Y$PTID, "==") * 1
  M <- M. %*% M
  grads[[n]] <- back.prop.(M, layers[[n]]$D, Ac[[n]])
  for (i in (n - 1):1) {
    M <- integ((M %*% t(Ac[[i + 1]][-1, ]) %*% t(Bases[[i + 1]])) *
                 A.prime[[i]](layers[[i]]$H), Bases[[i + 1]])
    grads[[i]] <- back.prop.(M, layers[[i]]$D, Ac[[i]])
  }
  return(grads)
}
back.prop. <- function(M, D, Ac) {
  N <- dim(M)[1]
  dc. <- colMeans(M)
  dA. <- t(D) %*% M/N
  grads. <- rbind(dc., dA.)
  return(grads.)
}