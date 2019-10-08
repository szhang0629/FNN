FNN.p <- function(Y, X, G, pos, A, lambda.){
  if (is.list(lambda.)) {
    groups <- divide(Y$PTID, type = "name", names = c("subtr", "valid"))
    list2env(split(mget(c('Y', 'X', 'G')), groups), envir = environment())
    I <- length(lambda.)
    error. <- rep(0, I)
    for (i in 1:I) {
      Ac <- FNN.p.(Y.subtr, X.subtr, G.subtr, pos, A, lambda.[[i]])$Ac
      print("done")
      error.[i] <- cost(Y.valid, Ac, X.valid, G.valid, Bases = NULL, A, pos)
    }
    lambda <- lambda.[[which.min(error.)]]
  } else lambda <- lambda.
  fnn.p <- FNN.p.(Y, X, G, pos, A, lambda)
  return(list(Ac = fnn.p$Ac, lambda = lambda, j = fnn.p$j))
}
FNN.p. <- function(Y, X, G, pos, A, lambda) {
  Bases <- Bf2m(rep(list(.bb), length(lambda)), 
                pos, subset(Y, select = -c(PTID, Y)))
  Ac <- init(Bases, seed = 1, X.col = ncol(X))
  e0 <- init(Bases, X.col = ncol(X), para.ub = 0)
  output <- list(Ac = Ac, eg. = e0, ex. = e0)
  # A.prime <- list()
  # for (i in 1:length(A))
  #   A.prime[[i]] <- Deriv(A[[i]])
  A.prime <- lapply(A, Deriv)
  Y$PTID <- idx.names(Y$PTID, rownames(G))
  Y$loc <- idx.names(subset(Y, select = -c(Y, PTID)))
  Y <- subset(Y, select = c(PTID, loc, Y))
  lambda <- lambda/nrow(Y)
  Pc <- list()
  PA <- list()
  for (i in 1:length(lambda)) {
    Pc[[i]] <- lambda[i + 1] * .P2
    PA[[i]] <- lambda[i + 1] * .P02 + lambda[i] * .P20
  }
  
  j <- 0
  cache <- list(k = 0, error. = Inf)
  while (cache$k < 30) { 
    if (j %% 10 == 0) {
      error <- cost(Y, output$Ac, X, G, Bases, A)
      if (j %% 100 == 0) cat(error, j, "\n")
      if (is.nan(error)) break
      if (error < (cache$error. * 0.999)) # if (error < (cache$error.))
        cache <- list(Ac. = output$Ac, j. = j, k = 0, error. = error)
      else cache$k <- cache$k + 1
    }
    output <- FNN.(Y, X, G, output, Bases, A, A.prime, Pc, PA)
    j <- j + 1
  }
  return(list(Ac = cache$Ac., error = cache$error., j = cache$j.))
}
idx.names <- function(df, nms = NULL) {
  if (is.null(nms)) {
    if (is.data.frame(df)) {
      if (length(df) > 1) df <- 1:nrow(df)
      else df <- idx.names(unlist(df), unique(unlist(df)))
    } else df <- idx.names(df, unique(df))
  } else {
    nms <- unlist(nms)
    df. <- 1:length(nms)
    names(df.) <- nms
    df <- df.[as.character(df)]
  }
  return(df)
}