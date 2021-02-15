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
  D <- cbind(matrix(rep(1, nrow(E)), nrow = nrow(E)), E)
  H <- D %*% Ac
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
  M <- (Y.hat - Y) * 2
  grads[[n]] <- t(layers[[n]]$D) %*% M
  if (n > 1) {
    for (i in (n - 1):1) {
      M <- (M %*% t(Ac[[i + 1]][-1, ]) * A.prime[[i]](layers[[i]]$H))
      grads[[i]] <- t(layers[[i]]$D) %*% M
    }
  }
  return(grads)
}
