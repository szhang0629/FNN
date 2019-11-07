forw.prop <- function(Ac, X, G, Bases, A, int, loc = NULL){
  layers <- list()
  dp <- length(A)
  layers[[1]] <- forw.prop.(Ac[[1]], G, Bases[[1]], Bases[[2]], 
                            A[[1]], int = int[[1]], X)
  if (dp > 2) {
    for (i in 2:(dp - 1)) 
      layers[[i]] <- forw.prop.(Ac[[i]], layers[[i - 1]]$G, Bases[[i]], 
                                Bases[[i + 1]], A[[i]], int[[i]])
  }
  layers[[dp]] <- forw.prop.(Ac[[dp]], layers[[dp - 1]]$G, Bases[[dp]], 
                             Bases[[dp + 1]], A[[dp]], int[[dp]], loc = loc)
  return(layers)
}
forw.prop. <- function(Ac, G, B0, B, A, int, X = NULL, loc = NULL){
  D <- cbind(matrix(rep(1, nrow(G)), nrow = nrow(G)), X, 
             integ(G, B0, int = int))
  C <- D %*% Ac
  if (is.null(loc)) H <- C %*% t(B)
  else 
    H <- rowSums(C[loc$idx, , drop = FALSE] * B[loc$idy, , drop = FALSE])
  G <- A(H)
  cache <- mget(c('G', 'D', 'H'))
  return(cache)
}
pred <- function(Ac, X, G, A, Bases = NULL, pos = NULL, loc = NULL, int = NULL){
  if (!is.null(pos)) {
    int <- lapply(Bases, is.basis)
    if (is.list(loc)) {
      idy <- idx.names(subset(loc, select = -PTID))
      idx <- idx.names(loc$PTID, rownames(G))
      Bases <- Bf2m(Bases, pos, subs(subset(loc, select = -PTID), 
                                     match(1:max(idy), idy)))
      loc <- data.frame(idx = idx, idy = idy)
    } else {
      Bases <- Bf2m(Bases, pos, loc)
      loc <- NULL
    }
  }
  layers <- forw.prop(Ac, X, G, Bases, A, int, loc)
  return(layers[[length(layers)]]$G)
}
back.prop <- function(Y, Ac, layers, Bases, A, A.prime, int = F, id){
  grads <- list()
  n <- length(layers)
  Y.hat <- layers[[n]]$G
  if (!is.null(id)) {
    M <- 2 * c(Y.hat - Y) * Bases[[length(Bases)]][id$idy, ]#/length(Y.hat)
    M. <- outer(1:max(id$idx), id$idx, "==") * 1
    M <- M. %*% M
  } else {
    M <- 2 * (Y.hat - Y) %*% Bases[[length(Bases)]]
  }
  grads[[n]] <- t(layers[[n]]$D) %*% M
  for (i in (n - 1):1) {
    M <- integ((M %*% t(Ac[[i + 1]][-1, ]) %*% t(Bases[[i + 1]])) *
                 A.prime[[i]](layers[[i]]$H), Bases[[i + 1]], int[[i]])
    grads[[i]] <- t(layers[[i]]$D) %*% M
  }
  return(grads)
}