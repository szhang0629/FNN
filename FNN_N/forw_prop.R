forw.prop <- function(Ac, X, G, Bases, A, loc = NULL){
  layers <- list()
  dp <- length(A)
  layers[[1]] <- forw.prop.(Ac[[1]], G, Bases[[1]], Bases[[2]], A[[1]], X)
  if (dp > 2) {
    for (i in 2:(dp - 1)) 
      layers[[i]] <- forw.prop.(Ac[[i]], layers[[i-1]]$G, Bases[[i]], 
                                Bases[[i + 1]], A[[i]])
  }
  layers[[dp]] <- forw.prop.(Ac[[dp]], layers[[dp-1]]$G, Bases[[dp]], 
                             Bases[[dp + 1]], A[[dp]], loc = loc)
  return(layers)
}
forw.prop. <- function(Ac, G, B0, B, A, X = NULL, loc = NULL){
  G <- integ(G, B0)
  D <- cbind(X, G)
  C <- D %*% Ac[-1, ] + rep.row(Ac[1, ], nrow(D))
  if (is.null(loc)) 
    H <- C %*% t(B)
  else 
    H <- rowSums(C[loc$PTID, ] * B[loc$loc, ])
  G <- A(H)
  cache <- mget(c('G', 'D', 'H'))
  return(cache)
}