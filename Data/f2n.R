f2n <- function(Y, X, G, loc) { # list
  G. <- list()
  X. <- list()
  Y. <- list()
  V. <- list()
  for (i in 1:length(loc)) {
    loc. <- loc[[i]]
    Y.. <- Y[[i]]
    G.. <- G[i, ]
    X.. <- c(X[i, ], loc.[1])
    n.loc <- length(loc.)
    elem <- loc.[2:length(loc.)] - loc.[1]
    l.elem <- length(elem)
    a <- matrix(rep(elem, l.elem), l.elem, l.elem)
    a[lower.tri(a)] <- t(a)[lower.tri(a)]
    # V.[[i]] <- solve
    V.[[i]] <- ginv(a)
    if (!is.null(G)) {
      G.. <- G[i, ]
      G.[[i]] <- rep.row(G.., n.loc - 1)
    }
    X.[[i]] <- cbind(rep.row(X.., n.loc - 1), loc.[2:n.loc])
    Y.[[i]] <- matrix(Y..[2:n.loc] - rep(Y..[1], n.loc - 1), ncol = 1)
  }
  if (!is.null(G))
    G <- do.call(rbind, G.)
  X <- do.call(rbind, X.)
  Y <- do.call(rbind, Y.)
  V <- do.call(bdiag, V.)
  return(list(G = G, X = X, Y = Y, V = V))
}