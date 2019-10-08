forw.prop <- function(E, Ac, A){
  layers <- list()
  for (i in 1:length(A)) {
    cache <- forw.prop.(E, Ac[[i]], A[[i]])
    E <- cache$E
    layers[[i]] <- cache
  }
  return(layers)
}
forw.prop. <- function(E, Ac, A){
  D <- E
  H <- D %*% Ac[-1, ] + rep.row(Ac[1, ], nrow(D))
  E <- A(H)
  cache <- mget(c('E', 'D', 'H'))
  return(cache)
}
pred <- function(Ac, E, A){
  layers <- forw.prop(E, Ac, A)
  result <- layers[[length(layers)]]$E
  return(result)
}
back.prop <- function(Y, Ac, layers, A, A.prime){
  grads <- list()
  n <- length(layers)
  Y.hat <- layers[[n]]$E
  M <- Y.hat - Y
  grads[[n]] <- back.prop.(M, layers[[n]]$D, Ac[[n]])
  for (i in (n - 1):1) {
    M <- (M %*% t(Ac[[i + 1]][-1, ]) * A.prime[[i]](layers[[i]]$H))
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