FNN.p. <- function(Y, X, G, Bases, A, lambda, pos, loc = null) {
  P <- pen.Bases(Bases)
  int <- lapply(Bases, is.list)
  # if (ncol(Y) > 1) {
  #   idy <- idx.names(Y)
  #   idx <- idx.names(rownames(Y), rownames(G))
  #   id <- data.frame(idx = idx, idy = idy)
  #   Bases <- Bf2m(Bases, pos, loc[match(1:max(idy), idy), , drop = FALSE])
  # } else {
  #   if (is.matrix(loc)) id <- NA
  #   else id <- NULL
  #   Bases <- Bf2m(Bases, pos, loc)
  # } 
  if (is.null(loc) && length(Y[, -1, drop = FALSE]) > 0) {
    Bases <- Bf2m(Bases, pos, Y[, -1, drop = FALSE])
    Y <- Y[, 1, drop = FALSE]
    rownames(Bases[[length(Bases)]]) <- idx.names(rownames(Y), rownames(G))
  } else Bases <- Bf2m(Bases, pos, loc)
  e0 <- init(Bases, X.col = ncol(X), para.ub = 0)
  output <- list(Ac = init(Bases, seed = 1, X.col = ncol(X)), eg = e0, ex = e0)
  A.prime <- lapply(A, Deriv)
  
  j <- 0
  k <- 0
  cache <- list(error = Inf)
  while (k < 30) { 
    if (j %% 10 == 0) {
      error.g <- sum((Y - pred(output$Ac, X, G, A, Bases, int = int))^2)
      error.p <- Error.p(output$Ac, P, lambda)
      error <- sum(error.g, error.p)
      if (j %% 1000 == 0) {
        cat(error.g, error.p, j)
        cat("\n")
      }
      if (is.nan(error))
        break
      if (error < (cache$error)) {
        cache <- list(Ac = output$Ac, j = j, error = error)
        k <- 0
      } else k <- k + 1
    }
    layers <- forw.prop(output$Ac, X, G, Bases, A, int)
    grads.g <- back.prop(Y, output$Ac, layers, Bases, A, A.prime, int)
    grads.p <- pen(output$Ac, P, lambda)
    output <- adadelta(output, grads.g, grads.p)
    j <- j + 1
  }
  cat(cache$j, "\n")
  return(cache)
}
idx.names <- function(df, nms = NULL) {
  if (is.matrix(df) && ncol(df) > 2)
    return(1:nrow(df))
  if (is.null(nms)) 
    df <- idx.names(df, unique(df))
  else {
    nms <- unlist(nms)
    df. <- 1:length(nms)
    names(df.) <- nms
    df <- df.[as.character(unlist(df))]
  }
  return(df)
}