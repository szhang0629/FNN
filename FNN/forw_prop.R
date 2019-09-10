forw.prop <- function(E, Ac, Bases, A, int = T){
  layers <- list()
  for (i in 1:length(A)) {
    cache <- forw.prop.(E, Ac[[i]], Bases[[i]], Bases[[i + 1]], 
                        A[[i]], int = int[[i]])
    E <- cache$E
    layers[[i]] <- cache
  }
  return(layers)
}
forw.prop. <- function(E, Ac, B0, B, A, int = T){
  if (is.list(E)) {
    G <- integ(E$G, B0, int = int)
    D <- cbind(E$X, G)
  }  else D <- integ(E, B0, int = int)
  C <- D %*% Ac[-1, ] + rep.row(Ac[1, ], nrow(D))
  if (is.list(B)) {
    E <- list()
    for (i in 1:length(B)) 
      E[[i]] <- A(C[i, , drop = FALSE] %*% t(B[[i]]))
    cache <- mget(c('E', 'D'))
  } else {
    H <- C %*% t(B)
    E <- A(H)
    cache <- mget(c('E', 'D', 'H'))
  }
  return(cache)
}