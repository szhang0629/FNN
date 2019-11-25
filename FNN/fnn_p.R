FNN.p. <- function(Y, X, G, Bases, A, lambda, pos, loc = NULL) {
  P <- pen.Bases(Bases)
  int <- lapply(Bases, is.list)
  
  if (is.list(Y)) {
    idy <- idx.names(subset(Y, select = -c(PTID, Y)))
    idx <- idx.names(Y$PTID, rownames(G))
    # Y <- data.frame(Y = Y$Y, idx = idx, idy = idy)
    id <- data.frame(idx = idx, idy = idy)
    Bases <- Bf2m(Bases, pos, subs(subset(Y, select = -c(PTID, Y)), 
                                   match(1:max(idy), idy)))
    Y <- Y$Y
  } else {
    Bases <- Bf2m(Bases, pos, loc)
    id <- NULL
  }
  
  e0 <- init(Bases, X.col = ncol(X), para.ub = 0)
  output <- list(Ac = init(Bases, seed = 1, X.col = ncol(X)), eg = e0, ex = e0)
  A.prime <- lapply(A, Deriv)
  
  j <- 0
  cache <- list(k = 0, error. = Inf)
  while (cache$k < 30) { 
    if (j %% 10 == 0) {
      error.g <- sum((Y - pred(output$Ac, X, G, A, Bases, loc = id, int = int))^2)
      error.p <- Error.p(output$Ac, P, lambda)
      error <- sum(error.g, error.p)
      if (j %% 1000 == 0) {
        cat(error.g, error.p, j)
        cat("\n")
      }
      if (is.nan(error))
        break
      # if (error < (cache$error. * 0.999))
      if (error < (cache$error.))
        cache <- list(Ac. = output$Ac, j. = j, k = 0, error. = error)
      else cache$k <- cache$k + 1
    }
    layers <- forw.prop(output$Ac, X, G, Bases, A, int = int, id)
    grads.g <- back.prop(Y, output$Ac, layers, Bases, A, A.prime, int = int, id)
    grads.p <- pen(output$Ac, P, lambda)
    output <- adadelta(output, grads.g, grads.p)
    j <- j + 1
  }
  return(list(Ac = cache$Ac., error = cache$error., j = cache$j.))
}
idx.names <- function(df, nms = NULL) {
  if (length(df) == 2 && is.data.frame(df))
    df <- data.frame(paste(df[, 1], df[, 2]), 
                     stringsAsFactors = FALSE)
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