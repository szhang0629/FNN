forw.prop <- function(Ac, X, G, Bases, A, int){
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
                             Bases[[dp + 1]], A[[dp]], int[[dp]])
  return(layers)
}
forw.prop. <- function(Ac, G, B0, B, A, int, X = NULL){
  D <- cbind(matrix(rep(1, nrow(G)), nrow = nrow(G)), X, 
             integ(G, B0, int = int))
  C <- D %*% Ac
  if (is.null(rownames(B))) H <- C %*% t(B)
  else H <- rowSums(C[as.numeric(rownames(B)), , drop = FALSE] * B)
  G <- A(H)
  cache <- mget(c('G', 'D', 'H'))
  return(cache)
}
pred <- function(Ac, X, G, A, Bases, pos = NULL, loc = NULL, int = NULL){
  if (!is.null(pos)) {
    int <- lapply(Bases, is.basis)
    Bases <- Bf2m(Bases, pos, loc)
    if (length(loc) > 0) 
      rownames(Bases[[length(Bases)]]) <- idx.names(rownames(loc), rownames(G))
  }
  layers <- forw.prop(Ac, X, G, Bases, A, int)
  return(layers[[length(layers)]]$G)
}
back.prop <- function(Y, Ac, layers, Bases, A, A.prime, int){
  grads <- list()
  n <- length(layers)
  Y.hat <- layers[[n]]$G
  B. <- Bases[[length(Bases)]]
  if (is.null(rownames(B.))) M <- 2 * (Y.hat - Y) %*% B.
  else {
    id <- as.numeric(rownames(B.))
    M <- 2 * c(Y.hat - Y) * B.#/length(Y.hat)
    M. <- outer(1:max(id), id, "==") * 1
    M <- M. %*% M
  }
  grads[[n]] <- t(layers[[n]]$D) %*% M
  for (i in (n - 1):1) {
    M <- integ((M %*% t(Ac[[i + 1]][-1, ]) %*% t(Bases[[i + 1]])) *
                 A.prime[[i]](layers[[i]]$H), Bases[[i + 1]], int[[i]])
    grads[[i]] <- t(layers[[i]]$D) %*% M
  }
  return(grads)
}