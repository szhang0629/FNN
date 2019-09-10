back.prop <- function(Y, Ac, layers, Bases, A, A.prime, int = F){
  grads <- list()
  n <- length(layers)
  Y.hat <- layers[[n]]$E
  M <- integ(Y.hat - Y, Bases[[length(Bases)]], T)
  grads[[n]] <- back.prop.(M, layers[[n]]$D, Ac[[n]])
  for (i in (n - 1):1) {
    M <- integ((M %*% t(Ac[[i + 1]][-1, ]) %*% t(Bases[[i + 1]])) *
                 A.prime[[i]](layers[[i]]$H), Bases[[i + 1]], int[[i]])
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